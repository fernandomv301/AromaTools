#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import sys

ORCA_OK = "****ORCA TERMINATED NORMALLY****"
G16_OK = "Normal termination of Gaussian 16 at"
G09_OK = "Normal termination of Gaussian 09 at"

def read_tail(path, max_bytes=4096):
    with open(path, "rb") as f:
        f.seek(0, 2) 
        size = f.tell()
        offset = min(size, max_bytes)
        f.seek(-offset, 1)
        data = f.read()
    return data.decode(errors="ignore")

def check_file(path: Path) -> str:
    tail = read_tail(path)
    lines = tail.splitlines()

    if any(ORCA_OK in line for line in lines):
        return "OK ORCA 6.1"

    if any(G09_OK in line for line in lines):
        return "OK Gaussian 09"

    if any(G16_OK in line for line in lines):
        return "OK Gaussian 16"

    return "Fail :c"

def main():
    if len(sys.argv) > 1:
        files = [Path(p) for p in sys.argv[1:]]
    else:
        files = sorted(Path(".").glob("*.out"))

    if not files:
        print("There's no file in this directory!")
        return

    for f in files:
        if not f.is_file():
            continue
        estado = check_file(f)
        print(f"{f.name}: {estado}")

if __name__ == "__main__":
    main()

