#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np


def read_vtk_scalar(path):
    with open(path, "r") as f:
        lines = f.readlines()
    dims = None
    for l in lines:
        if l.startswith("DIMENSIONS"):
            dims = tuple(map(int, l.split()[1:4]))
            break
    if dims is None:
        raise SystemExit("ERROR: There's no DIMENSIONS in vtk file!")
    pd_idx = None
    for i, l in enumerate(lines):
        if l.startswith("POINT_DATA"):
            pd_idx = i
            npts = int(l.split()[1])
            break
    if pd_idx is None:
        raise SystemExit("ERROR: Does not find POINT_DATA in vtk file!")
    data_start = None
    for i in range(pd_idx, len(lines)):
        if "LOOKUP_TABLE" in lines[i]:
            data_start = i + 1
            break
    if data_start is None:
        raise SystemExit("ERROR: Does not find LOOKUP_TABLE in vtk file!")
    data = []
    for l in lines[data_start:]:
        if l.strip():
            data.append(float(l.split()[0]))
    if len(data) != npts:
        raise SystemExit("ERROR: Number of points does not match POINT_DATA")
    # Return only the header lines (up to and including LOOKUP_TABLE)
    # so write_vtk_scalar doesn't duplicate the old data block.
    return lines[:data_start], dims, np.array(data)


def write_vtk_scalar(path, header_lines, data):
    with open(path, "w") as f:
        for l in header_lines:
            f.write(l)
        for v in data:
            f.write(f"{v:.8e}\n")


def read_cube(path):
    with open(path, "r") as f:
        lines = f.readlines()
    header = lines[:6]
    nx = int(header[3].split()[0])
    ny = int(header[4].split()[0])
    nz = int(header[5].split()[0])
    npts = nx * ny * nz
    data = []
    for l in lines[6:]:
        data.extend([float(x) for x in l.split()])
    if len(data) != npts:
        raise SystemExit("ERROR: Number of points does not match in cube file")
    return lines[:6], (nx, ny, nz), np.array(data)


def write_cube(path, header, data):
    with open(path, "w") as f:
        for l in header:
            f.write(l)
        for i, v in enumerate(data):
            f.write(f"{v:13.6e}")
            if (i + 1) % 6 == 0:
                f.write("\n")
        if len(data) % 6 != 0:
            f.write("\n")


def main():
    if len(sys.argv) != 4:
        raise SystemExit("Use: sum_rest.py [sum|rest] file1 file2")
    op = sys.argv[1].lower()
    f1 = sys.argv[2]
    f2 = sys.argv[3]
    if op not in ("sum", "rest"):
        raise SystemExit("ERROR: The operation must be 'sum' or 'rest'")
    ext1 = os.path.splitext(f1)[1]
    ext2 = os.path.splitext(f2)[1]
    if ext1 != ext2:
        raise SystemExit("ERROR: Can't mix VTK with CUBE")
    if ext1 == ".vtk":
        h1, dims1, d1 = read_vtk_scalar(f1)
        h2, dims2, d2 = read_vtk_scalar(f2)
        if dims1 != dims2:
            raise SystemExit("ERROR: The vtk files do not have the same dimensions!")
        if d1.size != d2.size:
            raise SystemExit("ERROR: The vtk files do not have the same point count!")
        result = d1 + d2 if op == "sum" else d1 - d2
        out = f"{op}_{os.path.basename(f1)}_{os.path.basename(f2)}"
        write_vtk_scalar(out, h1, result)
        print(f"File generated!: {out}")
    elif ext1 == ".cube":
        h1, dims1, d1 = read_cube(f1)
        h2, dims2, d2 = read_cube(f2)
        if dims1 != dims2:
            raise SystemExit("ERROR: The cube files do not have the same dimensions!")
        if d1.size != d2.size:
            raise SystemExit("ERROR: The cube files do not have the same point count!")
        result = d1 + d2 if op == "sum" else d1 - d2
        out = f"{op}_{os.path.basename(f1)}_{os.path.basename(f2)}"
        write_cube(out, h1, result)
        print(f"File generated!: {out}")
    else:
        raise SystemExit("ERROR: This script only supports vtk and cube files!")
