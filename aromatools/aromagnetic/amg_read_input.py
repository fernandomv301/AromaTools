#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

BoolLike = Union[bool, str]
Scalar = Union[int, float, str, bool]
Value = Union[Scalar, List[Scalar], None]


_TRUE = {"YES", "Y", "TRUE", "T", "ON", "1"}
_FALSE = {"NO", "N", "FALSE", "F", "OFF", "0"}

_KEEP_STRING_KEYS = {
    "KEYWORD_LINE",
    "PROGRAM",
    "MOLECULAR_COORDINATES",
    "SCALE",
}

_MODULE_KEYS = {
    "GRID_MODULE",
    "VERTICAL_PROFILE",
    "HORIZONTAL_PROFILE",
}


def _strip_comment(line: str) -> str:
    if "#" in line:
        line = line.split("#", 1)[0]
    return line.strip()


def _as_bool(token: str) -> Optional[bool]:
    t = token.strip().upper()
    if t in _TRUE:
        return True
    if t in _FALSE:
        return False
    return None


def _parse_value(tokens: List[str], *, keep_as_string: bool = False) -> Value:
    if not tokens:
        return None

    if keep_as_string:
        return " ".join(tokens).strip()

    if len(tokens) == 1:
        tok = tokens[0].strip()
        b = _as_bool(tok)
        if b is not None:
            return b
        try:
            return int(tok)
        except ValueError:
            pass
        try:
            return float(tok)
        except ValueError:
            pass
        return tok

    parsed: List[Scalar] = []
    for tok in tokens:
        tok = tok.strip()
        if not tok:
            continue
        b = _as_bool(tok)
        if b is not None:
            parsed.append(b)
            continue
        try:
            parsed.append(int(tok))
            continue
        except ValueError:
            pass
        try:
            parsed.append(float(tok))
            continue
        except ValueError:
            pass
        parsed.append(tok)
    return parsed


def read_input(filepath: Union[str, Path] = "INPUT.txt"):

    path = Path(filepath).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"INPUT file not found: {path}")

    general: Dict[str, Any] = {}
    modules: Dict[str, bool] = {}
    blocks: Dict[str, Dict[str, Any]] = {}
    post: Dict[str, Any] = {}
    orca_raw_lines: List[str] = []
    current_block: Optional[str] = None

    orca_block_name = "ORCA_BLOCK"

    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for lineno, raw in enumerate(lines, start=1):
        line = _strip_comment(raw)
        if not line:
            continue
        upper = line.strip().upper()

        # ---- block open ----
        if upper.startswith("BEGIN_"):
            blk = upper[len("BEGIN_"):].strip()
            current_block = blk
            if blk == orca_block_name:
                orca_raw_lines = []
            else:
                blocks.setdefault(blk, {})
            continue

        # ---- block close ----
        if upper.startswith("END_"):
            blk = upper[len("END_"):].strip()
            if current_block != blk:
                raise ValueError(
                    f"INPUT.txt block mismatch at line {lineno}: "
                    f"END_{blk} while inside {current_block}"
                )
            if blk == orca_block_name:
                blocks[orca_block_name] = {
                    "raw": "\n".join(orca_raw_lines).rstrip()
                    + ("\n" if orca_raw_lines else "")
                }
            current_block = None
            continue

        # ---- raw ORCA block lines ----
        if current_block == orca_block_name:
            orca_raw_lines.append(raw.rstrip("\n"))
            continue

        # ---- key-value pair ----
        parts = line.split()
        key = parts[0].strip()
        val_tokens = parts[1:]
        key_upper = key.upper()
        keep_as_string = key_upper in _KEEP_STRING_KEYS
        value = _parse_value(val_tokens, keep_as_string=keep_as_string)

        if current_block is None:
            # Module switches: GRID_MODULE, VERTICAL_PROFILE, HORIZONTAL_PROFILE
            if key_upper in _MODULE_KEYS or key_upper.endswith("_MODULE"):
                b: Optional[bool] = None
                if isinstance(value, bool):
                    b = value
                elif isinstance(value, str):
                    b = _as_bool(value)
                if b is None:
                    raise ValueError(
                        f"Invalid module switch value for {key} at line {lineno}: {value!r}"
                    )
                modules[key_upper] = b
            else:
                general[key_upper] = value

        else:
            blk_dict = blocks.setdefault(current_block, {})
            if current_block == "ORBITALS":
                if value is None:
                    blk_dict[key_upper] = []
                elif isinstance(value, list):
                    ints: List[int] = []
                    for v in value:
                        if isinstance(v, int):
                            ints.append(v)
                        elif isinstance(v, str) and v.strip().isdigit():
                            ints.append(int(v))
                        else:
                            raise ValueError(
                                f"Invalid MO list for {key} in ORBITALS block "
                                f"at line {lineno}: {value!r}"
                            )
                    blk_dict[key_upper] = ints
                elif isinstance(value, int):
                    blk_dict[key_upper] = [value]
                else:
                    raise ValueError(
                        f"Invalid ORBITALS entry for {key} at line {lineno}: {value!r}"
                    )
            elif key_upper == "PREDI_XY_PATH":
                blk_dict.setdefault(key_upper, []).append(value)
            else:
                blk_dict[key_upper] = value

    # ---- post-processing keys ----
    post_keys = {
        "RING_CURRENT_STRENGHT",
        "INTEGRAL_NICS",
        "BINDZ_VALUES",
    }
    for k in list(general.keys()):
        if k in post_keys:
            post[k] = general.pop(k)

    config: Dict[str, Any] = {
        "general": general,
        "modules": modules,
        "_meta": {"path": str(path)},
    }

    # Expose known blocks at top level for easy access
    # New names: VP (Vertical Profile) and HP (Horizontal Profile)
    if "GRID" in blocks:
        config["GRID"] = blocks["GRID"]
    if "VP" in blocks:
        config["VP"] = blocks["VP"]
    if "HP" in blocks:
        config["HP"] = blocks["HP"]
    if "ORBITALS" in blocks:
        config["ORBITALS"] = blocks["ORBITALS"]
    if "ORCA_BLOCK" in blocks:
        config["ORCA_BLOCK"] = blocks["ORCA_BLOCK"]
    if post:
        config["post"] = post

    return config


if __name__ == "__main__":
    cfg = read_input("INPUT.txt")
