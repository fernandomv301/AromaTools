#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import os
import re
import sys
import math
import subprocess
import matplotlib.pyplot as plt
import numpy as np
from ase import io as ase_io
from pathlib import Path

from . import parser_mo
from .amg_read_input import read_input
from .get_cube import *
from .get_nics_scan import (
    _parse_heights,
    best_fit_plane_axes,
    compute_half_range,
    build_scan_positions
)
from .get_integral import romberg

# ======================================== For MO ========================================

MU_0 = 4.0 * math.pi * 1.0e-7
BOHR_M = 5.29177210903e-11
PROP_1_OVER_MU0 = (1.0 / MU_0) * BOHR_M * 1.0e3

def integrate_series(x, y, tol=1e-6):
    x = np.array(x, float)
    y = np.array(y, float)

    def f(t):
        return float(np.interp(t, x, y))

    val, it = romberg(f, float(x.min()), float(x.max()), tol=tol)
    return float(val), int(it)


def ring_current_from_bindz(integral_bindz, symmetry=4.0):
    factor = PROP_1_OVER_MU0 * float(symmetry)
    return float(-factor * integral_bindz)

# ======================================== For MO ========================================

def convert_cube_vtk(cube_files):
    script = Path(__file__).with_name("cube_to_vtk.py")

    if not script.exists():
        print("WARNING: cube_to_vtk.py not found. Skipping VTK conversion.")
        return
    
    for cube in cube_files:
        cube = str(cube)
        try:
            subprocess.run(
                [sys.executable, str(script), cube],
                check=True
            )
            print(f"[ok] VTK generated from {cube}")
        except subprocess.CalledProcessError:
            print(f"WARNING: failed to convert {cube} to VTK")

def find_profile_prefixes(name):
    rx = re.compile(
        rf"^{re.escape(name)}_ring\d{{4}}(?:_\d{{4}})?\.(out|log)$",
        re.IGNORECASE
    )
    prefixes = set()
    for fn in os.listdir(os.getcwd()):
        if rx.match(fn):
            base = re.sub(r"\.(out|log)$", "", fn, flags=re.IGNORECASE)
            base = re.sub(r"_\d{4}$", "", base)
            prefixes.add(base)
    return sorted(prefixes)

def build_profile_distances(profile_cfg):
    n = int(profile_cfg["NUM_GA"])
    scale = str(profile_cfg.get("SCALE", "LINEAR")).upper()

    if scale == "LOGARITHMIC":
        start = float(profile_cfg["START"])
        end = float(profile_cfg["END"])
        return np.logspace(np.log10(start), np.log10(end), n).tolist()

    start = float(profile_cfg.get("START", 0.0))
    delta = float(profile_cfg.get("DELTA_X", profile_cfg.get("DELTA_Z", 0.1)))
    return [start + i * delta for i in range(n)]

def _truthy(value):
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().upper() in {"ON", "TRUE", "YES", "1"}

def _get_general_value(general, *keys, default=None):
    for key in keys:
        if key in general and general[key] not in (None, ""):
            return general[key]
        low = key.lower()
        if low in general and general[low] not in (None, ""):
            return general[low]
        up = key.upper()
        if up in general and general[up] not in (None, ""):
            return general[up]
    return default

def detect_active_module(cfg):
    mods = cfg.get("modules", {})
    active = []

    for key in ("GRID_MODULE", "PROFILE_MODULE", "NICS_SCAN_XY_MODULE"):
        if _truthy(mods.get(key, False)):
            active.append(key)

    if len(active) != 1:
        raise SystemExit(
            f"ERROR: exactly one module must be active in INPUT.txt. "
            f"Detected active modules: {active}"
        )

    return active[0]


def load_common_cfg():
    cfg = read_input("INPUT.txt")
    if not cfg:
        raise SystemExit("ERROR: INPUT.txt does not exist or could not be read.")

    general = cfg.get("general", {})
    program = str(
        _get_general_value(general, "PROGRAM", default="Gaussian")
    ).strip().lower()

    if program != "gaussian":
        raise SystemExit(
            "ERROR: MO decomposition is currently implemented only for Gaussian outputs."
        )

    name = str(
        _get_general_value(
            general,
            "NAME_OUTPUT_FILE",
            "NAME_MOLECULE",
            default="molecule"
        )
    ).strip()

    if not name:
        raise SystemExit(
            "ERROR: NAME_OUTPUT_FILE (or NAME_MOLECULE) is missing in INPUT.txt"
        )

    spec, others_labels = parser_mo.parse_orbitals_from_cfg(cfg)

    if not spec and not others_labels:
        raise SystemExit(
            "ERROR: BEGIN_ORBITALS section not found or empty in INPUT.txt. "
            "Define MO_* groups before running 'aromagnetic mo'."
        )

    labels = list(spec.keys()) + list(others_labels)
    return cfg, name, spec, others_labels, labels

def _validate_vector_lengths(prefix, expected, labels, data_dict, what):
    expected = int(expected)
    bad = []
    for lb in labels:
        arr = data_dict.get(lb)
        if arr is None:
            bad.append(f"{lb}: missing")
            continue
        if len(arr) != expected:
            bad.append(f"{lb}: {len(arr)} != {expected}")

    if bad:
        raise SystemExit(
            f"ERROR: invalid {what} data length for prefix '{prefix}': " + "; ".join(bad)
        )

def _validate_xyz_lengths(prefix, expected, labels, data_dict, what):
    expected = int(expected)
    bad = []
    for lb in labels:
        arrs = data_dict.get(lb)
        if arrs is None or len(arrs) != 3:
            bad.append(f"{lb}: missing xyz tuple")
            continue
        lens = tuple(len(x) for x in arrs)
        if lens != (expected, expected, expected):
            bad.append(f"{lb}: {lens} != ({expected}, {expected}, {expected})")

    if bad:
        raise SystemExit(
            f"ERROR: invalid {what} xyz data length for prefix '{prefix}': " + "; ".join(bad)
        )

def run_mo_profile(cfg, name, spec, others_labels, labels):
    profile = cfg.get("PROFILE", {})
    if "NUM_GA" not in profile:
        raise SystemExit(
            "PROFILE mode: NUM_GA must be defined in INPUT.txt (BEGIN_PROFILE)"
        )

    expected = int(profile["NUM_GA"])
    grid = cfg.get("GRID", {})
    da_max = int(profile.get("GA_MAX", grid.get("GA_MAX", expected)))
    dist = build_profile_distances(profile)
    prefixes = find_profile_prefixes(name)

    if len(dist) != expected:
        raise SystemExit(
            f"ERROR: profile distance array length ({len(dist)}) does not match "
            f"NUM_GA ({expected})."
        )

    if not prefixes:
        raise SystemExit(
            "PROFILE mode: no ring outputs found "
            "(expected name_ringXXXX*.out/.log)"
        )

    block = [
    "==========================================\n",
    "      Integral Analysis Summary (MO)      \n",
    "==========================================\n",
    f"Molecule            : {name}\n",
    "Program             : gaussian\n",
    f"Ring outputs found  : {len(prefixes)}\n",
    "\n"
    ]

    Bext = np.array((0.0, 0.0, 1.0), float)
    for prefix in prefixes:
        print(f"[info] MO-profile: {prefix}")

        nics_mo = parser_mo.parse_gaussian_nics(
            prefix,
            expected,
            spec,
            others_labels=others_labels,
            da_max=da_max
        )
        bindz_mo = parser_mo.parse_gaussian_bind(
            prefix,
            expected,
            spec,
            others_labels=others_labels,
            Bext=Bext,
            da_max=da_max,
            component="z"
        )
        _validate_vector_lengths(prefix, expected, labels, nics_mo, "NICS")
        _validate_vector_lengths(prefix, expected, labels, bindz_mo, "Bind(z)")
        out_csv = f"{prefix}_mo.csv"
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(
                ["distance_A"]
                + [f"NICS_{lb}_ppm" for lb in labels]
                + [f"Bindz_{lb}_ppm" for lb in labels]
            )
            for i, d in enumerate(dist):
                row = [f"{d:.8f}"]
                row += [f"{nics_mo[lb][i]:.8f}" for lb in labels]
                row += [f"{bindz_mo[lb][i]:.8f}" for lb in labels]
                w.writerow(row)
        plt.figure()
        for lb in labels:
            plt.plot(dist, nics_mo[lb], label=lb)
        plt.xlabel("Distance (Å)")
        plt.ylabel("NICS (ppm)")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{prefix}_NICS_mo.png", dpi=200)
        plt.close()
        plt.figure()
        for lb in labels:
            plt.plot(dist, bindz_mo[lb], label=lb)
        plt.xlabel("Distance (Å)")
        plt.ylabel("Bindz (ppm)")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{prefix}_Bindz_mo.png", dpi=200)
        plt.close()
        print(f"[ok] {out_csv}, {prefix}_NICS_mo.png, {prefix}_Bindz_mo.png")
        for lb in labels:
            plt.figure()
            plt.plot(dist, bindz_mo[lb])
            plt.xlabel("Distance (Å)")
            plt.ylabel("Bindz (ppm)")
            plt.tight_layout()
            plt.savefig(f"{prefix}_Bindz_{lb}.png", dpi=200)
            plt.close()
            plt.figure()
            plt.plot(dist, nics_mo[lb])
            plt.xlabel("Distance (Å)")
            plt.ylabel("NICS (ppm)")
            plt.tight_layout()
            plt.savefig(f"{prefix}_NICS_{lb}.png", dpi=200)
            plt.close()
        npts = len(dist)
        if npts < 2:
            block.append(f"{prefix}\n")
            block.append("  Status                 : skipped (not enough points)\n\n")
            continue
        x0 = float(dist[0])
        x1 = float(dist[-1])
        step_A = float(dist[1] - dist[0]) if npts >= 2 else float("nan")
        block.append(f"{prefix}\n")
        block.append("------------------------------------------\n")
        block.append(f"CSV file                : {out_csv}\n")
        block.append(f"Points                  : {npts}\n")
        block.append(f"Distance range [Å]      : {x0:.6f} -> {x1:.6f}\n")
        block.append(f"Step [Å]                : {step_A:.6f}\n")
        for lb in labels:
            block.append(f"\n[MO group: {lb}]\n")
            I_bindz, it = integrate_series(dist, bindz_mo[lb], tol=1e-6)
            current = ring_current_from_bindz(I_bindz, symmetry=4.0)
            block.append(f"Integral Bindz [ppm·Å]  : {I_bindz:.10f}\n")
            block.append(f"Romberg iterations      : {it}\n")
            block.append(f"1/mu0 factor            : {PROP_1_OVER_MU0:.10f}\n")
            block.append(f"Ring current [nA/T]     : {current:.10f}\n")
            nics_vals = np.array(nics_mo[lb], float)
            if not np.all(np.isnan(nics_vals)):
                valid = [(d, n) for d, n in zip(dist, nics_vals) if not np.isnan(n)]
                if len(valid) >= 2:
                    d2 = [p[0] for p in valid]
                    n2 = [p[1] for p in valid]
                    I_nics, itn = integrate_series(d2, n2, tol=1e-6)
                    block.append(f"Integral NICS [ppm·Å]   : {I_nics:.10f}\n")
                    block.append(f"NICS Romberg iterations : {itn}\n")
                else:
                    block.append("Integral NICS [ppm·Å]   : not available\n")
            else:
                block.append("Integral NICS [ppm·Å]   : not available\n")

        with open("int_out.txt", "w", encoding="utf-8") as f:
            f.writelines(block)
        print("[done] Wrote int_out.txt")

    print("Done! MO profile CSV/plots generated!")

def run_mo_nics_scan(cfg, name, spec, others_labels, labels):
    scan = cfg.get("NICS_SCAN_XY", {})
    general = cfg.get("general", {})
    xyz_file = _get_general_value(general, "MOLECULAR_COORDINATES")
    if not xyz_file or not os.path.isfile(xyz_file):
        raise SystemExit("ERROR: MOLECULAR_COORDINATES must exist for NICS_SCAN_XY MO analysis")
    npts = int(scan.get("POINTS", 81))
    margin = float(scan.get("MARGIN", 2.0))
    delta_x = float(scan.get("DELTA_X", 0.25))
    make_y = _truthy(scan.get("MAKE_AXIS_Y", "NO"))
    heights = _parse_heights(scan.get("HEIGHT", 1.0))
    da_max = int(scan.get("GA_MAX", npts))
    ats = ase_io.read(xyz_file, index=0)
    coords = ats.get_positions()
    center, xhat, yhat, zhat = best_fit_plane_axes(coords)
    xs, _ = build_scan_positions(coords, center, xhat, margin, npts, delta_x)
    ys = None
    if make_y:
        ys, _ = build_scan_positions(coords, center, yhat, margin, npts, delta_x)
    print("\n=== NICS-SCAN (MO analyze) ===")
    print(f"Molecule    : {name}")
    print(f"POINTS      : {npts}")
    print(f"HEIGHTS     : {', '.join(str(h) for h in heights)}")
    print(f"MAKE_AXIS_Y : {'YES' if make_y else 'NO'}\n")
    for h in heights:
        tag = f"h{h:.2f}".replace(".", "p")
        prefix_x = f"{name}_nicsX_{tag}"
        print(f"[info] MO NICS-scan X: {prefix_x}")
        nics_x = parser_mo.parse_gaussian_nics(
            prefix_x,
            npts,
            spec,
            others_labels=others_labels,
            da_max=da_max
        )
        _validate_vector_lengths(prefix_x, npts, labels, nics_x, "NICS-scan X")
        csv_x = f"{prefix_x}_mo.csv"
        with open(csv_x, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["scan_coordinate_A"] + [f"NICS_{lb}_ppm" for lb in labels])
            for i, x in enumerate(xs):
                row = [f"{x:.8f}"] + [f"{nics_x[lb][i]:.8f}" for lb in labels]
                w.writerow(row)
        plt.figure()
        for lb in labels:
            plt.plot(xs, nics_x[lb], label=lb)
        plt.xlabel("Scan coordinate (Å)")
        plt.ylabel("NICS (ppm)")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{prefix_x}_mo.png", dpi=200)
        plt.close()
        for lb in labels:
            plt.figure()
            plt.plot(xs, nics_x[lb])
            plt.xlabel("Scan coordinate (Å)")
            plt.ylabel("NICS (ppm)")
            plt.tight_layout()
            plt.savefig(f"{prefix_x}_{lb}.png", dpi=200)
            plt.close()
        print(f"[ok] {csv_x}, {prefix_x}_mo.png")
        if make_y and ys is not None:
            prefix_y = f"{name}_nicsY_{tag}"
            print(f"[info] MO NICS-scan Y: {prefix_y}")
            nics_y = parser_mo.parse_gaussian_nics(
                prefix_y,
                npts,
                spec,
                others_labels=others_labels,
                da_max=da_max
            )
            _validate_vector_lengths(prefix_y, npts, labels, nics_y, "NICS-scan Y")
            csv_y = f"{prefix_y}_mo.csv"
            with open(csv_y, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(["scan_coordinate_A"] + [f"NICS_{lb}_ppm" for lb in labels])
                for i, y in enumerate(ys):
                    row = [f"{y:.8f}"] + [f"{nics_y[lb][i]:.8f}" for lb in labels]
                    w.writerow(row)
            plt.figure()
            for lb in labels:
                plt.plot(ys, nics_y[lb], label=lb)
            plt.xlabel("Scan coordinate (Å)")
            plt.ylabel("NICS (ppm)")
            plt.legend()
            plt.tight_layout()
            plt.savefig(f"{prefix_y}_mo.png", dpi=200)
            plt.close()
            for lb in labels:
                plt.figure()
                plt.plot(ys, nics_y[lb])
                plt.xlabel("Scan coordinate (Å)")
                plt.ylabel("NICS (ppm)")
                plt.tight_layout()
                plt.savefig(f"{prefix_y}_{lb}.png", dpi=200)
                plt.close()
            print(f"[ok] {csv_y}, {prefix_y}_mo.png")
    print("\nDone! MO NICS-scan CSV/plots generated!")

def run_mo_grid(cfg, name, spec, others_labels, labels):
    grid = cfg.get("GRID", {})
    if "NUM_POINTS" not in grid:
        raise SystemExit(
            "ERROR: NUM_POINTS Nx Ny Nz must be defined in INPUT.txt (BEGIN_GRID)"
        )

    Bext = np.array(grid.get("EXTERNAL_FIELD", (0.0, 0.0, 1.0)), float)
    da_max = int(grid.get("GA_MAX", 30))
    Nx, Ny, Nz = map(int, grid["NUM_POINTS"])
    expected = Nx * Ny * Nz
    gi = read_grid_info("grid_info.txt")
    calc_xyz = gi.get("calc_xyz") or f"{name}in_calcgrid.xyz"
    full_xyz = gi.get("full_xyz") or f"{name}in_symgrid.xyz"
    mol_file = f"{name}_oriented.xyz"

    if not Path(mol_file).exists():
        raise SystemExit(f"ERROR: expected oriented XYZ not found: {mol_file}")
    if not Path(calc_xyz).exists():
        raise SystemExit(f"ERROR: reduced grid xyz file not found: {calc_xyz}")
    if not Path(full_xyz).exists():
        raise SystemExit(f"ERROR: full grid xyz file not found: {full_xyz}")

    atoms_rows = load_atoms_bohr(mol_file)
    _, calc_pts = read_xpts(calc_xyz)
    _, full_pts = read_xpts(full_xyz)
    idx_full, shape_full, origin_ang, steps_ang = build_full_mapping(calc_pts, full_pts)
    ox, oy, oz = np.array(origin_ang) * ANG2BOHR
    dx, dy, dz = np.array(steps_ang) * ANG2BOHR
    bind_data = parser_mo.parse_gaussian_bind(
        name,
        expected,
        spec,
        others_labels=others_labels,
        Bext=Bext,
        da_max=da_max,
        component="xyz"
    )

    nics_data = parser_mo.parse_gaussian_nics(
        name,
        expected,
        spec,
        others_labels=others_labels,
        da_max=da_max
    )

    _validate_xyz_lengths(name, expected, labels, bind_data, "Bind")
    _validate_vector_lengths(name, expected, labels, nics_data, "NICS")

    generated_cubes = []
    for label in labels:
        Bx_red, By_red, Bz_red = bind_data[label]
        Bx = [Bx_red[i] for i in idx_full]
        By = [By_red[i] for i in idx_full]
        Bz = [Bz_red[i] for i in idx_full]
        NICS = [nics_data[label][i] for i in idx_full]
        shape_out, origin_out, steps_out, Bx_out = pad_vis(
            shape_full, (ox, oy, oz), (dx, dy, dz), Bx, pad_to=3
        )
        _, _, _, By_out = pad_vis(
            shape_full, (ox, oy, oz), (dx, dy, dz), By, pad_to=3
        )
        _, _, _, Bz_out = pad_vis(
            shape_full, (ox, oy, oz), (dx, dy, dz), Bz, pad_to=3
        )
        _, _, _, NICS_out = pad_vis(
            shape_full, (ox, oy, oz), (dx, dy, dz), NICS, pad_to=3
        )

        cube_bindx = f"{name}_Bindx_{label}.cube"
        cube_bindy = f"{name}_Bindy_{label}.cube"
        cube_bindz = f"{name}_Bindz_{label}.cube"
        cube_nics  = f"{name}_NICS_{label}.cube"

        write_cube(
            cube_bindx,
            atoms_rows,
            origin_out,
            steps_out,
            shape_out,
            Bx_out,
            per_line=1
        )
        write_cube(
            cube_bindy,
            atoms_rows,
            origin_out,
            steps_out,
            shape_out,
            By_out,
            per_line=1
        )
        write_cube(
            cube_bindz,
            atoms_rows,
            origin_out,
            steps_out,
            shape_out,
            Bz_out,
            per_line=1
        )
        write_cube(
            cube_nics,
            atoms_rows,
            origin_out,
            steps_out,
            shape_out,
            NICS_out,
            per_line=1
        )

        generated_cubes.extend([cube_bindx, cube_bindy, cube_bindz, cube_nics])
    
    vtk_files = grid.get("VTK_FILES", False)

    if _truthy(vtk_files):
        print("[info] Converting MO cube files to VTK...")
        convert_cube_vtk(generated_cubes)
        
    print("Done! Cube files generated!")

def main():
    cfg, name, spec, others_labels, labels = load_common_cfg()
    active_module = detect_active_module(cfg)

    if active_module == "PROFILE_MODULE":
        run_mo_profile(cfg, name, spec, others_labels, labels)
        return

    if active_module == "GRID_MODULE":
        run_mo_grid(cfg, name, spec, others_labels, labels)
        return

    if active_module == "NICS_SCAN_XY_MODULE":
        run_mo_nics_scan(cfg, name, spec, others_labels, labels)
        return

    raise SystemExit("ERROR: unknown active module")

if __name__ == "__main__":
    main()

