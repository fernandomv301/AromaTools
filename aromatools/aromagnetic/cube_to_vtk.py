#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os

BOHR2ANG = 0.52917721092  # Bohr to Angstrom

def read_cube(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        title = f.readline().rstrip("\n")
        comment = f.readline().rstrip("\n")
        parts = f.readline().split()
        if len(parts) < 4:
            raise ValueError(f"Bad atom/origin line in {path}")
        natoms = int(parts[0])
        ox, oy, oz = map(float, parts[1:4])

        parts = f.readline().split()
        if len(parts) < 4:
            raise ValueError(f"Bad X-axis line in {path}")
        nx = int(parts[0])
        dx_x, dx_y, dx_z = map(float, parts[1:4])

        parts = f.readline().split()
        if len(parts) < 4:
            raise ValueError(f"Bad Y-axis line in {path}")
        ny = int(parts[0])
        dy_x, dy_y, dy_z = map(float, parts[1:4])

        parts = f.readline().split()
        if len(parts) < 4:
            raise ValueError(f"Bad Z-axis line in {path}")
        nz = int(parts[0])
        dz_x, dz_y, dz_z = map(float, parts[1:4])

        for _ in range(natoms):
            f.readline()

        values = []
        for line in f:
            for tok in line.split():
                values.append(float(tok))

    expected = nx * ny * nz
    if len(values) != expected:
        raise ValueError(
            f"{path}: expected {expected} grid values, got {len(values)}"
        )

    dx = (dx_x**2 + dx_y**2 + dx_z**2) ** 0.5
    dy = (dy_x**2 + dy_y**2 + dy_z**2) ** 0.5
    dz = (dz_x**2 + dz_y**2 + dz_z**2) ** 0.5

    origin_bohr = (ox, oy, oz)
    steps_bohr = (dx, dy, dz)
    shape = (nx, ny, nz)

    return title, comment, natoms, shape, origin_bohr, steps_bohr, values

def cube_to_vtk(cube_path, vtk_path=None):
    if not os.path.isfile(cube_path):
        raise FileNotFoundError(f"{cube_path} not found")

    title, comment, natoms, (nx, ny, nz), origin_bohr, steps_bohr, values = read_cube(cube_path)

    ox_b, oy_b, oz_b = origin_bohr
    dx_b, dy_b, dz_b = steps_bohr

    ox = ox_b * BOHR2ANG
    oy = oy_b * BOHR2ANG
    oz = oz_b * BOHR2ANG
    sx = dx_b * BOHR2ANG
    sy = dy_b * BOHR2ANG
    sz = dz_b * BOHR2ANG

    ntot = nx * ny * nz

    if vtk_path is None:
        vtk_path = cube_path + "_.vtk"

    scalar_name = os.path.basename(cube_path)

    grid = [[[0.0 for _ in range(nz)] for _ in range(ny)] for _ in range(nx)]
    idx = 0
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                grid[ix][iy][iz] = values[idx]
                idx += 1

    with open(vtk_path, "w", encoding="utf-8") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write(f"{scalar_name}  vtk_file from CUBE conversion\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {nx:d} {ny:d} {nz:d}\n")
        f.write(f"ORIGIN {ox:12.6f} {oy:12.6f} {oz:12.6f}\n")
        f.write(f"SPACING {sx:12.6f} {sy:12.6f} {sz:12.6f}\n")
        f.write(f"POINT_DATA {ntot:d}\n")
        f.write(f"SCALARS {scalar_name} float 1\n")
        f.write("LOOKUP_TABLE default\n")

        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    f.write(f"{grid[ix][iy][iz]: .6e}\n")

    print(f"[ok] {cube_path} -> {vtk_path}")

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    if not argv:
        print("Use: cube_to_vtk.py archivo1.cube [archivo2.cube ...]")
        sys.exit(1)

    for cube in argv:
        cube_to_vtk(cube)

if __name__ == "__main__":
    main()

