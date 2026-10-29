#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# get_results.py
# Unified analysis script for all profile types.
#
#   VERTICAL_PROFILE   → parses ring outputs, writes CSVs + plots,
#                         computes Romberg integrals and ring currents,
#                         writes int_out.txt summary.
#
#   HORIZONTAL_PROFILE → parses hpX / hpY outputs, writes CSVs + plots.

import os
import re
import math
import numpy as np
import matplotlib.pyplot as plt
from ase import io as ase_io

from pathlib import Path as _Path

from .amg_read_input import read_input
from .tools import _truthy, _parse_heights
from .check_out import check_file as _check_file


def _warn_bad_terminations(out_files):
    for f in out_files:
        if not _check_file(_Path(f)).startswith("OK"):
            print(f"[WARN] {f}: calculation did not terminate normally")

# ============================================================
#  Physical constants  (used for ring-current integration)
# ============================================================

MU_0              = 4.0 * math.pi * 1.0e-7
BOHR_M            = 5.29177210903e-11
PROP_1_OVER_MU0   = (1.0 / MU_0) * BOHR_M * 1.0e3   # nA/T per (ppm·Å)


# ============================================================
#  Output-file finders
# ============================================================

def find_ring_outputs(name):
    """Return {prefix: [sorted output files]} for VP ring outputs."""
    rx = re.compile(
        rf"^({re.escape(name)}_ring\d{{4}})(?:_(\d{{4}}))?\.(out|log)$",
        re.IGNORECASE,
    )
    groups = {}
    for fn in os.listdir(os.getcwd()):
        m = rx.match(fn)
        if not m:
            continue
        prefix    = m.group(1)
        chunk     = m.group(2)
        chunk_idx = int(chunk) if chunk is not None else -1
        groups.setdefault(prefix, []).append((chunk_idx, fn))
    return {
        prefix: [fn for _, fn in sorted(items, key=lambda x: x[0])]
        for prefix, items in groups.items()
    }


def find_hp_outputs(name, axis, tag):
    """
    Return sorted output files for one HP prefix.
    axis : 'X' or 'Y'
    tag  : e.g. 'h0p00'
    """
    prefix = f"{name}_hp{axis}_{tag}"
    rx = re.compile(
        rf"^{re.escape(prefix)}(?:_(\d{{4}}))?\.(out|log)$",
        re.IGNORECASE,
    )
    items = []
    for fn in os.listdir(os.getcwd()):
        m = rx.match(fn)
        if not m:
            continue
        chunk = m.group(1)
        idx   = int(chunk) if chunk is not None else -1
        items.append((idx, fn))
    items.sort(key=lambda x: x[0])
    return [fn for _, fn in items], prefix


def find_ring_csvs(name):
    """Return sorted list of VP ring CSV files."""
    rx = re.compile(rf"^({re.escape(name)}_ring\d{{4}})\.csv$", re.IGNORECASE)
    csvs = sorted(fn for fn in os.listdir(os.getcwd()) if rx.match(fn))
    return csvs


# ============================================================
#  QC output parsers  (ORCA and Gaussian, shared by VP and HP)
# ============================================================

_float_rx = re.compile(r"[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?")


def _parse_orca_tensors(filename):
    """Return list of (xx, yy, zz) diagonal shielding for every atom found."""
    diags = []
    with open(filename, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if "Total shielding tensor (ppm):" in lines[i]:
            tensor_lines = []
            for j in range(i + 1, i + 4):
                if j < len(lines):
                    t = lines[j].strip()
                    if t:
                        tensor_lines.append(t)
            if len(tensor_lines) == 3:
                r1 = _float_rx.findall(tensor_lines[0])
                r2 = _float_rx.findall(tensor_lines[1])
                r3 = _float_rx.findall(tensor_lines[2])
                if len(r1) >= 3 and len(r2) >= 3 and len(r3) >= 3:
                    diags.append((float(r1[0]), float(r2[1]), float(r3[2])))
            i += 4
        else:
            i += 1
    return diags


def read_orca_dummy_diags(out_files, n_real_atoms):
    """Return (xx, yy, zz) tuples for dummy atoms only (VP use)."""
    diags = []
    for fn in out_files:
        all_diags = _parse_orca_tensors(fn)
        diags.extend(all_diags[n_real_atoms:])
    return diags


def read_orca_dummy_isotropic(out_files, n_real_atoms):
    """Return isotropic shielding values for dummy atoms only (HP use)."""
    iso = []
    for fn in out_files:
        all_diags = _parse_orca_tensors(fn)
        for xx, yy, zz in all_diags[n_real_atoms:]:
            iso.append((xx + yy + zz) / 3.0)
    return iso


def _parse_gaussian_one_file(filename):
    """Return list of {'zz': float, 'iso': float|None} for every Bq atom."""
    num = r"[-+]?\d+(?:\.\d+)?(?:[EeDd][-+]?\d+)?"
    rx_zz      = re.compile(rf"ZZ=\s*({num})")
    rx_bq_iso  = re.compile(
        rf"^\s*(?:\d+\s+)?Bq\b.*?Isotropic\s*=\s*({num})", re.IGNORECASE
    )
    with open(filename, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    iso_list = []
    for line in lines:
        m = rx_bq_iso.search(line)
        if m:
            iso_list.append(float(m.group(1).replace("D", "E").replace("d", "e")))

    bq_data = []
    for i, line in enumerate(lines):
        if "Bq" not in line:
            continue
        zz = None
        for j in range(i + 1, min(i + 4, len(lines))):
            if "ZZ=" in lines[j]:
                mzz = rx_zz.search(lines[j])
                if mzz:
                    zz = float(mzz.group(1).replace("D", "E").replace("d", "e"))
                break
        if zz is not None:
            bq_data.append({"zz": zz, "iso": None})

    for k in range(min(len(bq_data), len(iso_list))):
        bq_data[k]["iso"] = iso_list[k]

    return bq_data


def read_gaussian_dummy_diags(out_files):
    """Return list of {'zz', 'iso'} dicts for all Bq atoms (VP use)."""
    data = []
    for fn in out_files:
        data.extend(_parse_gaussian_one_file(fn))
    return data


def read_gaussian_bq_isotropic(out_files):
    """Return isotropic shielding values for all Bq atoms (HP use)."""
    num = r"[-+]?\d+(?:\.\d+)?(?:[EeDd][-+]?\d+)?"
    rx  = re.compile(
        rf"^\s*(?:\d+\s+)?Bq\b.*?Isotropic\s*=\s*({num})", re.IGNORECASE
    )
    iso = []
    for fn in out_files:
        with open(fn, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                m = rx.search(line)
                if m:
                    iso.append(float(m.group(1).replace("D", "E").replace("d", "e")))
    return iso


# ============================================================
#  Distance / scan-coordinate builders
# ============================================================

def build_vp_distances(vp_block):
    """Reconstruct the distance array used when generating VP inputs."""
    scale = str(vp_block.get("SCALE", "Linear")).strip().upper()
    if scale == "LINEAR":
        start = float(vp_block.get("START", 0.0))
        dz    = float(vp_block.get("DELTA_Z", 0.1))
        npts  = int(vp_block.get("NUM_GA", 100))
        return [start + i * dz for i in range(npts)]
    if scale == "LOGARITHMIC":
        start = float(vp_block.get("START", 0.1))
        end   = float(vp_block.get("END", 5.0))
        npts  = int(vp_block.get("NUM_GA", 100))
        return list(np.logspace(np.log10(start), np.log10(end), npts))
    raise SystemExit("ERROR: VP SCALE must be Linear or Logarithmic")


def build_hp_xcoords(coords, hp_block):
    """
    Reconstruct X (and Y) scan coordinate arrays used when generating HP inputs.
    Returns (xs, ys_or_None, heights, make_y).
    """
    npts         = int(hp_block.get("POINTS", 81))
    margin       = float(hp_block.get("MARGIN", 2.0))
    delta_x_user = float(hp_block.get("DELTA_X", 0.25))
    make_y       = _truthy(hp_block.get("MAKE_AXIS_Y", "NO"))
    heights      = _parse_heights(hp_block.get("HEIGHT", 1.0))

    center = coords.mean(axis=0)
    xhat   = np.array([1.0, 0.0, 0.0])
    yhat   = np.array([0.0, 1.0, 0.0])

    def _half_range(axis_hat):
        rel  = coords - center
        proj = rel @ axis_hat
        geom = float(np.max(np.abs(proj)) + margin)
        step = 0.5 * (npts - 1) * delta_x_user
        return max(geom, step)

    hr_x = _half_range(xhat)
    xs   = np.linspace(-hr_x, hr_x, npts).tolist()

    ys = None
    if make_y:
        hr_y = _half_range(yhat)
        ys   = np.linspace(-hr_y, hr_y, npts).tolist()

    return xs, ys, heights, make_y


# ============================================================
#  CSV writers
# ============================================================

def save_vp_csv(prefix, distances, bindz, nics):
    """Write VP result CSV.  All magnetic quantities in ppm."""
    csv_name = f"{prefix}.csv"
    with open(csv_name, "w", encoding="utf-8") as f:
        f.write("distance_A,Bindz_ppm,NICS_ppm\n")
        for d, b, n in zip(distances, bindz, nics):
            if np.isnan(n):
                f.write(f"{d:.8f},{b:.8f},\n")
            else:
                f.write(f"{d:.8f},{b:.8f},{n:.8f}\n")
    return csv_name


def save_hp_csv(prefix, xvals, nics_vals):
    """Write HP result CSV.  NICS in ppm."""
    csv_name = f"{prefix}.csv"
    with open(csv_name, "w", encoding="utf-8") as f:
        f.write("x_A,NICS_ppm\n")
        for x, n in zip(xvals, nics_vals):
            f.write(f"{x:.8f},{n:.8f}\n")
    return csv_name


def read_vp_csv(csv_file):
    """Read a VP CSV back for integration."""
    dist, bindz, nics = [], [], []
    with open(csv_file, "r", encoding="utf-8", errors="ignore") as f:
        f.readline()                            # skip header
        for line in f:
            if not line.strip():
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 2:
                continue
            dist.append(float(parts[0]))
            bindz.append(float(parts[1]))
            nics.append(float(parts[2]) if len(parts) >= 3 and parts[2] else float("nan"))
    return dist, bindz, nics


# ============================================================
#  Plotting
# ============================================================

def plot_vp(prefix, x, y, ylabel):
    """Line plot for a vertical-profile quantity vs distance."""
    png = f"{prefix}_{ylabel.replace(' ', '_')}.png"
    plt.figure()
    plt.plot(x, y)
    plt.xlabel("Distance (Å)")
    plt.ylabel(ylabel)
    plt.title(f"{prefix}: {ylabel} vs Distance")
    plt.tight_layout()
    plt.savefig(png, dpi=200)
    plt.close()
    return png


def plot_hp(prefix, xvals, nics_vals):
    """Line plot for a horizontal-profile NICS scan."""
    png = f"{prefix}.png"
    plt.figure()
    plt.plot(xvals, nics_vals)
    plt.xlabel("Scan coordinate (Å)")
    plt.ylabel("NICS (ppm)")
    plt.title(prefix)
    plt.tight_layout()
    plt.savefig(png, dpi=200)
    plt.close()
    return png


# ============================================================
#  Romberg integration + ring-current
# ============================================================

def romberg(f, a, b, tol=1e-6, max_iter=20):
    R = [[0.5 * (b - a) * (f(a) + f(b))]]
    n = 1
    for iteration in range(1, max_iter + 1):
        h     = (b - a) / (2 ** n)
        sum_f = sum(f(a + (2 * k - 1) * h) for k in range(1, 2 ** (n - 1) + 1))
        R.append([0.5 * R[n - 1][0] + h * sum_f])
        for m in range(1, n + 1):
            R[n].append(R[n][m - 1] + (R[n][m - 1] - R[n - 1][m - 1]) / (4 ** m - 1))
        if abs(R[n][n] - R[n - 1][n - 1]) < tol:
            return R[n][n], iteration
        n += 1
    return R[-1][-1], max_iter


def integrate_series(x, y, tol=1e-6):
    x = np.array(x, float)
    y = np.array(y, float)
    val, it = romberg(lambda t: float(np.interp(t, x, y)),
                      float(x.min()), float(x.max()), tol=tol)
    return float(val), int(it)


def ring_current_from_bindz(integral_bindz, symmetry=4.0):
    return float(-PROP_1_OVER_MU0 * float(symmetry) * integral_bindz)


# ============================================================
#  VP analysis
# ============================================================

def analyze_vertical_profile(cfg, name, program, xyz_file, n_real):
    vp_block  = cfg.get("VP", {})
    distances = build_vp_distances(vp_block)
    groups    = find_ring_outputs(name)

    if not groups:
        raise SystemExit(
            f"ERROR: No ring output files found. Expected {name}_ring####.out/.log"
        )

    print(f"\n=== VERTICAL PROFILE (analyze) ===")
    print(f"Molecule : {name}")
    print(f"Program  : {program}")
    print(f"XYZ      : {xyz_file}  (real atoms: {n_real})")
    print("Convention: Bext=(0,0,1)  →  Bindz = -sigma_zz\n")

    summary = [
        "==========================================\n",
        "        Integral Analysis Summary         \n",
        "==========================================\n",
        f"Molecule            : {name}\n",
        f"Program             : {program}\n",
        f"XYZ                 : {xyz_file}\n",
        f"CSV files found     : {len(groups)}\n",
        "\n",
    ]

    for prefix, out_files in sorted(groups.items()):
        print(f"[info] {prefix}: reading {len(out_files)} file(s)")
        _warn_bad_terminations(out_files)

        # --- extract tensor diagonals ---
        if "orca" in program:
            diags = read_orca_dummy_diags(out_files, n_real)
            if not diags:
                print(f"[warn] {prefix}: no dummy tensors found. Skipping.")
                summary += [f"{prefix}\n", "  Status : skipped (no tensors)\n\n"]
                continue
            n_use = min(len(diags), len(distances))
            if len(diags) != len(distances):
                print(f"[warn] {prefix}: distances={len(distances)} tensors={len(diags)} using {n_use}")
            bindz = [-float(zz)                            for xx, yy, zz in diags[:n_use]]
            nics  = [-(float(xx)+float(yy)+float(zz))/3.0 for xx, yy, zz in diags[:n_use]]

        elif "gaussian" in program:
            data = read_gaussian_dummy_diags(out_files)
            if not data:
                print(f"[warn] {prefix}: no Bq ZZ found. Skipping.")
                summary += [f"{prefix}\n", "  Status : skipped (no Bq data)\n\n"]
                continue
            n_use = min(len(data), len(distances))
            if len(data) != len(distances):
                print(f"[warn] {prefix}: distances={len(distances)} tensors={len(data)} using {n_use}")
            bindz, nics = [], []
            for item in data[:n_use]:
                bindz.append(-float(item["zz"]))
                nics.append(-float(item["iso"]) if item["iso"] is not None else float("nan"))
        else:
            raise SystemExit("ERROR: PROGRAM must be ORCA or GAUSSIAN.")

        # --- write CSV ---
        save_vp_csv(prefix, distances[:n_use], bindz, nics)
        print(f"[ok] wrote {prefix}.csv")

    # --- read CSVs back and integrate ---
    csv_files = find_ring_csvs(name)
    if not csv_files:
        raise SystemExit(
            f"ERROR: No ring CSV files found for integration. Expected {name}_ring####.csv"
        )

    for csv_file in csv_files:
        prefix = os.path.splitext(csv_file)[0]
        print(f"[info] integrating {prefix}")
        dist, bindz, nics = read_vp_csv(csv_file)
        npts = len(dist)

        summary.append(f"{prefix}\n")
        summary.append("------------------------------------------\n")
        summary.append(f"CSV file                : {csv_file}\n")
        summary.append(f"Points                  : {npts}\n")

        if npts < 2:
            print(f"[warn] {prefix}: not enough points for integration. Skipping.")
            summary.append("  Status                 : skipped (not enough points)\n\n")
            continue

        x0     = float(dist[0])
        x1     = float(dist[-1])
        step_A = float(dist[1] - dist[0])
        summary.append(f"Distance range [Å]      : {x0:.6f} -> {x1:.6f}\n")
        summary.append(f"Step [Å]                : {step_A:.6f}\n")

        # plots
        png_bindz = plot_vp(prefix, dist, bindz, "Bindz (ppm)")
        print(f"[ok] plot: {png_bindz}")
        summary.append(f"Bindz plot              : {png_bindz}\n")

        if not all(np.isnan(v) for v in nics):
            png_nics = plot_vp(prefix, dist, nics, "NICS (ppm)")
            print(f"[ok] plot: {png_nics}")
            summary.append(f"NICS plot               : {png_nics}\n")
        else:
            print(f"[info] {prefix}: NICS all NaN, skipping NICS plot.")
            summary.append("NICS plot               : not available\n")

        # Romberg integration of Bindz
        I_bindz, it      = integrate_series(dist, bindz)
        factor_rc        = PROP_1_OVER_MU0 * 2.0
        current          = ring_current_from_bindz(I_bindz, symmetry=4.0)
        summary.append(f"Integral Bindz [ppm·Å]  : {I_bindz:.10f}\n")
        summary.append(f"Romberg iterations      : {it}\n")
        summary.append(f"1/mu0 factor            : {PROP_1_OVER_MU0:.10f}\n")
        summary.append(f"Total factor            : {factor_rc:.10f}\n")
        summary.append(f"=========== Ring current [nA/T]   : {current:.10f} ===========\n")

        # Romberg integration of NICS
        nics_valid = [(d, n) for d, n in zip(dist, nics) if not np.isnan(n)]
        if len(nics_valid) >= 2:
            d2, n2    = zip(*nics_valid)
            I_nics, itn = integrate_series(list(d2), list(n2))
            summary.append(f"=========== Integral NICS [ppm·Å] : {I_nics:.10f} ===========\n")
            summary.append(f"NICS Romberg iterations : {itn}\n")
        else:
            summary.append("Integral NICS [ppm·Å]   : not available (missing NICS values)\n")

        summary.append("\n")

    with open("int_out.txt", "w", encoding="utf-8") as f:
        f.writelines(summary)
    print("[done] wrote int_out.txt")


# ============================================================
#  HP analysis
# ============================================================

def analyze_horizontal_profile(cfg, name, program, xyz_file, n_real, interactive=False):
    hp_block   = cfg.get("HP", {})
    predixy_on = _truthy(hp_block.get("PREDI_XY", "NO"))
    nics_xy_on = _truthy(hp_block.get("NICS_SCAN_XY", "NO"))

    if predixy_on:
        from .predixy import run_predixy_analyze
        run_predixy_analyze(cfg, name, program, xyz_file, n_real, interactive=interactive)

    if not nics_xy_on:
        return

    ats                 = ase_io.read(xyz_file, index=0)
    coords              = ats.get_positions()
    xs, ys, heights, make_y = build_hp_xcoords(coords, hp_block)

    print(f"\n=== HORIZONTAL PROFILE (analyze) ===")
    print(f"Molecule    : {name}")
    print(f"Program     : {program}")
    print(f"Heights (Å) : {', '.join(str(h) for h in heights)}")
    print(f"MAKE_AXIS_Y : {'YES' if make_y else 'NO'}\n")

    for h in heights:
        tag = f"h{h:.2f}".replace(".", "p")

        axes_to_scan = [("X", xs)] + ([("Y", ys)] if make_y else [])
        for axis, scan_coords in axes_to_scan:
            out_files, prefix = find_hp_outputs(name, axis, tag)
            if not out_files:
                print(f"[warn] No outputs found for {prefix}.* (skipping)")
                continue
            _warn_bad_terminations(out_files)

            # --- extract isotropic shielding ---
            if "gaussian" in program:
                iso = read_gaussian_bq_isotropic(out_files)
            elif "orca" in program:
                iso = read_orca_dummy_isotropic(out_files, n_real)
            else:
                raise SystemExit("ERROR: PROGRAM must be ORCA or GAUSSIAN.")

            n_use = min(len(scan_coords), len(iso))
            if len(iso) != len(scan_coords):
                print(
                    f"[warn] {prefix}: expected {len(scan_coords)} points, "
                    f"got {len(iso)} isotropic values. Using {n_use}."
                )

            nics = [-float(v) for v in iso[:n_use]]
            csv_name = save_hp_csv(prefix, scan_coords[:n_use], nics)
            png_name = plot_hp(prefix, scan_coords[:n_use], nics)
            print(f"[ok] {prefix}: {csv_name}, {png_name}")

    print("\nDone.")


# ============================================================
#  main()
# ============================================================

def main(interactive=False):
    cfg = read_input("INPUT.txt")
    if not cfg:
        raise SystemExit("ERROR: INPUT.txt not found")

    general  = cfg.get("general", {})

    name     = str(general.get("NAME_MOLECULE", "molecule")).strip()
    program  = str(general.get("PROGRAM", "orca")).strip().lower()
    xyz_file = general.get("MOLECULAR_COORDINATES")

    if not xyz_file or not os.path.isfile(xyz_file):
        raise SystemExit("ERROR: MOLECULAR_COORDINATES must exist for profile analysis.")

    n_real = len(ase_io.read(xyz_file, index=0))

    vp_on = "VP" in cfg and bool(cfg["VP"])
    hp_on = "HP" in cfg and bool(cfg["HP"])

    if not vp_on and not hp_on:
        raise SystemExit(
            "ERROR: INPUT.txt must contain BEGIN_VP...END_VP "
            "or BEGIN_HP...END_HP (or both) for profile analysis."
        )

    if vp_on:
        analyze_vertical_profile(cfg, name, program, xyz_file, n_real)

    if hp_on:
        analyze_horizontal_profile(cfg, name, program, xyz_file, n_real, interactive=interactive)


if __name__ == "__main__":
    main()
