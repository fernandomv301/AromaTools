#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import re
import sys
from typing import List, Tuple, Optional

SUPPORTED_EXTS = (".inp",)
ELEM_RE = re.compile(r"^\s*([A-Za-z]{1,2}:?)(?=\s)")
NUM_RE = re.compile(r"[-+]?\d*\.\d+|[-+]?\d+\.?|\d+\.?\d*(?:[Ee][-+]?\d+)?")


def is_geom_line(line: str) -> Optional[str]:
    m = ELEM_RE.match(line)
    if not m:
        return None

    elem_raw = m.group(1)
    nums = NUM_RE.findall(line)
    if len(nums) < 3:
        return None

    return elem_raw

def normalize_elem(elem_raw: str) -> str:
    e = elem_raw.strip()
    tagged = e.endswith(":")
    core = e[:-1] if tagged else e

    if len(core) == 0:
        return e

    core_norm = core[0].upper() + core[1:].lower()
    if core_norm.lower() == "bq":
        core_norm = "Bq"

    return core_norm + (":" if tagged else "")

def elem_core(elem_norm: str) -> str:
    return elem_norm[:-1] if elem_norm.endswith(":") else elem_norm

def is_tagged(elem_norm: str) -> bool:
    return elem_norm.endswith(":")

def find_geom_block(lines: List[str]) -> Optional[Tuple[int, int]]:
    start = None
    for i, line in enumerate(lines):
        if start is None:
            if is_geom_line(line) is not None:
                start = i
        else:
            if line.strip() == "":
                return (start, i)
    if start is not None:
        return (start, len(lines))
    return None

def check_only_CH(geom_lines: List[str]) -> Tuple[bool, List[str]]:
    species = []
    bad = set()
    for ln in geom_lines:
        eraw = is_geom_line(ln)
        if eraw is None:
            continue
        en = normalize_elem(eraw)
        if is_tagged(en):
            continue
        if elem_core(en) == "Bq":
            continue
        c = elem_core(en)
        species.append(c)
        if c not in ("C", "H"):
            bad.add(c)

    return (len(bad) == 0), sorted(set(species))

def pseudo_pi(geom_lines: List[str]) -> List[str]:
    out = []
    for ln in geom_lines:
        eraw = is_geom_line(ln)
        if eraw is None:
            out.append(ln)
            continue
        en = normalize_elem(eraw)
        if is_tagged(en):
            out.append(ln)
            continue
        core = elem_core(en)
        if core == "Bq":
            out.append(ln)
            continue
        if core == "H":
            continue
        if core == "C":
            parts = ln.rstrip("\n").split()
            if not parts:
                continue
            parts[0] = "H"
            out.append("  ".join(parts) + "\n")
            continue
        out.append(ln)

    return out

def main():

    files = []
    for ext in SUPPORTED_EXTS:
        files.extend(sorted(glob.glob(f"*{ext}")))
    if not files:
        print("No *.inp files found in the current directory.")
        sys.exit(0)
    changed = 0
    skipped = 0
    for path in files:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
        block = find_geom_block(lines)
        if block is None:
            print(f"[skip] {path}: geometry block not detected.")
            skipped += 1
            continue
        s, e = block
        geom = lines[s:e]
        ok, species = check_only_CH(geom)
        if not ok:
            print(
                f"[NO pseudo-π] {path}: detected REAL-atom species = {species}. "
                f"Pseudo-π requires CH-only inputs (excluding tagged dummies like H:)."
            )
            skipped += 1
            continue

        new_geom = pseudo_pi(geom)

        if new_geom == geom:
            print(f"[skip] {path}: no geometry changes applied.")
            skipped += 1
            continue

        out_lines = lines[:s] + new_geom + lines[e:]
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.writelines(out_lines)

        print(f"[ok] {path}: pseudo-π applied.")
        changed += 1

    print(f"\nDone. Modified={changed}, Skipped={skipped}, Total={len(files)}")

if __name__ == "__main__":
    main()

