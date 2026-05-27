"""Configuration loading for the UTI89 pipeline.

The project config is YAML, but the first wrapper should not require a package
install just to prepare a run directory. We use PyYAML when available and fall
back to a small parser that supports the simple maps/lists used by the example
config.
"""

from pathlib import Path
from typing import Any, Dict, List, Tuple, Union


def load_config(path: Union[str, Path]) -> Dict[str, Any]:
    """Load a pipeline YAML config file."""
    config_path = Path(path)
    try:
        import yaml  # type: ignore
    except ImportError:
        return _load_simple_yaml(config_path)

    with config_path.open() as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Config file must contain a mapping: {config_path}")
    return data


def write_yaml(path: Union[str, Path], data: Dict[str, Any]) -> None:
    """Write a simple YAML file using stable key order."""
    output_path = Path(path)
    output_path.write_text(_dump_mapping(data), encoding="utf-8")


def get_required(config: Dict[str, Any], dotted_key: str) -> Any:
    """Read a required nested config value using dot notation."""
    current: Any = config
    for part in dotted_key.split("."):
        if not isinstance(current, dict) or part not in current:
            raise KeyError(f"Missing required config key: {dotted_key}")
        current = current[part]
    return current


def _load_simple_yaml(path: Path) -> Dict[str, Any]:
    lines = []  # type: List[Tuple[int, str]]
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            if not raw_line.strip() or raw_line.lstrip().startswith("#"):
                continue
            indent = len(raw_line) - len(raw_line.lstrip(" "))
            lines.append((indent, raw_line.strip()))

    value, index = _parse_block(lines, 0, 0)
    if index != len(lines):
        raise ValueError(f"Could not parse entire config file: {path}")
    if not isinstance(value, dict):
        raise ValueError(f"Config file must contain a mapping: {path}")
    return value


def _parse_block(
    lines: List[Tuple[int, str]], start: int, indent: int
) -> Tuple[Any, int]:
    if start >= len(lines):
        return {}, start

    actual_indent, first_text = lines[start]
    if actual_indent < indent:
        return {}, start
    if actual_indent != indent:
        raise ValueError(f"Unexpected indentation near: {first_text}")

    if first_text.startswith("- "):
        return _parse_list(lines, start, indent)
    return _parse_mapping(lines, start, indent)


def _parse_mapping(
    lines: List[Tuple[int, str]], start: int, indent: int
) -> Tuple[Dict[str, Any], int]:
    result = {}  # type: Dict[str, Any]
    index = start
    while index < len(lines):
        line_indent, text = lines[index]
        if line_indent < indent:
            break
        if line_indent > indent:
            raise ValueError(f"Unexpected nested line: {text}")
        if text.startswith("- "):
            break
        if ":" not in text:
            raise ValueError(f"Expected key/value line: {text}")

        key, raw_value = text.split(":", 1)
        key = key.strip()
        raw_value = raw_value.strip()
        index += 1
        if raw_value:
            result[key] = _parse_scalar(raw_value)
            continue
        if index >= len(lines) or lines[index][0] <= indent:
            result[key] = {}
            continue
        child_indent = lines[index][0]
        result[key], index = _parse_block(lines, index, child_indent)
    return result, index


def _parse_list(
    lines: List[Tuple[int, str]], start: int, indent: int
) -> Tuple[List[Any], int]:
    result = []  # type: List[Any]
    index = start
    while index < len(lines):
        line_indent, text = lines[index]
        if line_indent < indent:
            break
        if line_indent != indent or not text.startswith("- "):
            break
        raw_value = text[2:].strip()
        index += 1
        if raw_value:
            result.append(_parse_scalar(raw_value))
            continue
        if index >= len(lines) or lines[index][0] <= indent:
            result.append(None)
            continue
        child_indent = lines[index][0]
        child, index = _parse_block(lines, index, child_indent)
        result.append(child)
    return result, index


def _parse_scalar(raw_value: str) -> Any:
    if raw_value in {"null", "Null", "NULL", "~"}:
        return None
    if raw_value in {"true", "True", "TRUE"}:
        return True
    if raw_value in {"false", "False", "FALSE"}:
        return False
    if (
        len(raw_value) >= 2
        and raw_value[0] == raw_value[-1]
        and raw_value[0] in {"'", '"'}
    ):
        return raw_value[1:-1]
    try:
        return int(raw_value)
    except ValueError:
        pass
    try:
        return float(raw_value)
    except ValueError:
        return raw_value


def _dump_mapping(data: Dict[str, Any], indent: int = 0) -> str:
    lines = []  # type: List[str]
    prefix = " " * indent
    for key, value in data.items():
        if isinstance(value, dict):
            lines.append(f"{prefix}{key}:")
            lines.append(_dump_mapping(value, indent + 2).rstrip())
        elif isinstance(value, list):
            lines.append(f"{prefix}{key}:")
            for item in value:
                if isinstance(item, dict):
                    lines.append(f"{prefix}  -")
                    lines.append(_dump_mapping(item, indent + 4).rstrip())
                else:
                    lines.append(f"{prefix}  - {_format_scalar(item)}")
        else:
            lines.append(f"{prefix}{key}: {_format_scalar(value)}")
    return "\n".join(lines) + "\n"


def _format_scalar(value: Any) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)
