#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import subprocess
import numpy as np
import networkx as nx
from ase import io as ase_io
from ase.data import atomic_numbers, covalent_radii

# AromaTools modules
from .tools import orient, _truthy, _parse_heights
from .amg_read_input import read_input
from .amg_write_inputs import make_inputs, flatten_cfg
from . import pseudo_pi
from .predixy import run_predixy_prepare

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def _only_ch(ase_atoms):
    syms = [s.upper() for s in ase_atoms.get_chemical_symbols()]
    return all(s in ("C", "H") for s in syms)


def _gen_used(keyword_line: str) -> bool:
    if not keyword_line:
        return False
    return re.search(r"\bgen(ecp)?\b", keyword_line, flags=re.IGNORECASE) is not None


def _apply_pseudo_pi(ats, label=""):
    if _only_ch(ats):
        try:
            pseudo_pi.main()
            print(f"[ok] pseudo_pi applied to generated inputs{' (' + label + ')' if label else ''}.")
        except Exception as e:
            raise SystemExit(f"\n[ERROR] pseudo_pi failed: {e}\n")
    else:
        print(
            f"[info] PSEUDO_PI=YES but molecule is not CH-only. "
            f"Skipping pseudo_pi{' (' + label + ')' if label else ''}."
        )


def _insert_basis(script_dir):
    try:
        insert_script = os.path.join(script_dir, "insert_basis.py")
        subprocess.check_call([sys.executable, insert_script])
        print("[ok] basis.txt appended automatically to generated inputs (GEN/GENECP).")
    except Exception as e:
        raise SystemExit(
            f"\n[ERROR] GEN basis insertion failed: {e}\n"
            "You can run it manually: python insert_basis.py\n"
        )


# ============================================================
#  VERTICAL PROFILE helpers
# ============================================================

def _setup_mol_graph(symbols, coords, bond_tol=0.4):
    G = nx.Graph()
    n = len(symbols)
    G.add_nodes_from(range(n))
    Z = [atomic_numbers[s] for s in symbols]
    Rcov = [covalent_radii[z] for z in Z]
    for i in range(n):
        for j in range(i):
            d_act = float(np.linalg.norm(coords[i] - coords[j]))
            d_cov = float(Rcov[i] + Rcov[j])
            if (d_act - d_cov) <= bond_tol:
                G.add_edge(i, j)
    return G


def get_cycles(symbols, coords, bond_tol=0.4):
    G = _setup_mol_graph(symbols, coords, bond_tol=bond_tol)
    return nx.minimum_cycle_basis(G)


def ring_center(coords, atom_indices):
    ring_coords = coords[np.array(atom_indices, dtype=int)]
    return ring_coords.mean(axis=0)


def ring_normal(atom_coords):
    if len(atom_coords) < 3:
        raise ValueError("At least 3 ring atoms are required to compute a normal vector.")
    v1 = atom_coords[1] - atom_coords[0]
    v2 = atom_coords[2] - atom_coords[0]
    nvec = np.cross(v1, v2)
    norm = float(np.linalg.norm(nvec))
    if norm == 0.0:
        raise ValueError("Normal vector is zero (degenerate ring geometry).")
    return nvec / norm


def build_vp_points(center, normal_vec, vp_block):
    """Build ghost-atom positions along the ring normal (vertical profile)."""
    pts = []
    zhat = np.array([0.0, 0.0, 1.0], float)
    if float(np.dot(normal_vec, zhat)) < 0.0:
        normal_vec = -normal_vec

    scale = str(vp_block.get("SCALE", "Linear")).strip().upper()
    if scale == "LINEAR":
        start = float(vp_block.get("START", 0.0))
        dz = float(vp_block.get("DELTA_Z", 0.1))
        npts = int(vp_block.get("NUM_GA", 100))
        if npts < 1:
            raise ValueError("NUM_GA must be >= 1 for Linear scale.")
        for i in range(npts):
            pts.append(center + (start + i * dz) * normal_vec)
    elif scale == "LOGARITHMIC":
        start = float(vp_block.get("START", 0.1))
        end = float(vp_block.get("END", 5.0))
        npts = int(vp_block.get("NUM_GA", 100))
        if start <= 0.0 or end <= 0.0:
            raise ValueError("START and END must be > 0 for Logarithmic scale.")
        if end <= start:
            raise ValueError("END must be > START for Logarithmic scale.")
        if npts < 2:
            raise ValueError("NUM_GA must be >= 2 for Logarithmic scale.")
        for z in np.logspace(np.log10(start), np.log10(end), npts):
            pts.append(center + float(z) * normal_vec)
    else:
        raise ValueError("SCALE must be Linear or Logarithmic (case-insensitive).")

    return [np.array(p, float) for p in pts]


def write_vp_xyz(path, ase_atoms, all_pts):
    """Write oriented molecule + all vertical profile ghost atoms."""
    real_syms = ase_atoms.get_chemical_symbols()
    real_xyz = ase_atoms.get_positions()
    n_total = len(real_syms) + len(all_pts)
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{n_total}\n")
        f.write("Real atoms + vertical profile points (X)\n")
        for sym, (x, y, z) in zip(real_syms, real_xyz):
            f.write(f"{sym:<2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")
        for (x, y, z) in all_pts:
            f.write(f"X  {x:14.8f} {y:14.8f} {z:14.8f}\n")


# ============================================================
#  HORIZONTAL PROFILE helpers
# ============================================================

def _compute_half_range(coords, center, axis_hat, margin):
    rel = coords - center
    proj = rel @ axis_hat
    return float(np.max(np.abs(proj)) + float(margin))


def _generate_hp_points(center, axis_hat, zhat, half_range, npts, height):
    """Return (scan_coords_1d, list_of_3d_points) for one height."""
    xs = np.linspace(-half_range, half_range, int(npts))
    pts = [center + x * axis_hat + float(height) * zhat for x in xs]
    return xs.tolist(), [np.array(p, float) for p in pts]


def write_hp_xyz(path, ase_atoms, all_pts):
    """Write oriented molecule + all horizontal scan ghost atoms."""
    real_syms = ase_atoms.get_chemical_symbols()
    real_xyz = ase_atoms.get_positions()
    n_total = len(real_syms) + len(all_pts)
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{n_total}\n")
        f.write("Real atoms + horizontal profile points (X)\n")
        for sym, (x, y, z) in zip(real_syms, real_xyz):
            f.write(f"{sym:<2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")
        for p in all_pts:
            x, y, z = p
            f.write(f"X  {x:14.8f} {y:14.8f} {z:14.8f}\n")


# ============================================================
#  Entry-point: VERTICAL PROFILE
# ============================================================

def run_vertical_profile(cfg, ats, coords, name, gen_flag, pseudo_pi_flag):
    vp_block = cfg.get("VP", {})
    bond_tol = float(vp_block.get("BOND_TOL", 0.4))
    symbols = ats.get_chemical_symbols()

    cycles = get_cycles(symbols, coords, bond_tol=bond_tol)
    if not cycles:
        raise SystemExit(
            "ERROR: No rings detected. VERTICAL_PROFILE requires at least one ring."
        )

    print(f"\n=== VERTICAL PROFILE (prepare) ===")
    print(f"Rings detected : {len(cycles)}")
    print(f"Bond tol (Å)   : {bond_tol}")
    print("Profile normal is oriented with fixed reference field (0,0,1).")

    cfg_flat = flatten_cfg(cfg, section="VP")
    all_profile_pts = []

    for rid, atom_indices in enumerate(cycles):
        center = ring_center(coords, atom_indices)
        ring_coords = coords[np.array(atom_indices, dtype=int)]
        nvec = ring_normal(ring_coords)
        pts = build_vp_points(center, nvec, vp_block)
        all_profile_pts.extend([p.tolist() for p in pts])

        cfg_ring = dict(cfg_flat)
        cfg_ring["da_max"] = max(1, int(len(pts)))
        cfg_ring["name_molecule"] = f"{name}_ring{rid:04d}"
        n_written = make_inputs(cfg_ring["name_molecule"], ats, pts, cfg_ring)
        print(
            f"[ok] ring {rid:04d}: atoms={len(atom_indices)} "
            f"profile_pts={len(pts)} inputs={n_written}"
        )

    profiles_xyz = f"{name}_vprofiles.xyz"
    write_vp_xyz(profiles_xyz, ats, all_profile_pts)
    print(f"[ok] Wrote vertical profiles XYZ: {profiles_xyz}")


# ============================================================
#  Entry-point: HORIZONTAL PROFILE
# ============================================================

def run_horizontal_profile(cfg, ats, coords, name, gen_flag, pseudo_pi_flag):
    hp_block    = cfg.get("HP", {})
    nics_xy_on  = _truthy(hp_block.get("NICS_SCAN_XY", "NO"))
    predixy_on  = _truthy(hp_block.get("PREDI_XY",     "NO"))

    if not nics_xy_on and not predixy_on:
        raise SystemExit(
            "ERROR: BEGIN_HP block found but neither NICS_SCAN_XY nor PREDI_XY is YES."
        )

    # ---- NICS_SCAN_XY ----
    if nics_xy_on:
        npts         = int(hp_block.get("POINTS", 81))
        margin       = float(hp_block.get("MARGIN", 2.0))
        delta_x_user = float(hp_block.get("DELTA_X", 0.25))
        make_y       = _truthy(hp_block.get("MAKE_AXIS_Y", "NO"))

        height_raw = hp_block.get("HEIGHT", 1.0)
        heights = _parse_heights(height_raw)

        if npts < 2:
            raise SystemExit("ERROR: POINTS must be >= 2 for NICS_SCAN_XY")

        center = coords.mean(axis=0)
        xhat   = np.array([1.0, 0.0, 0.0])
        yhat   = np.array([0.0, 1.0, 0.0])
        zhat   = np.array([0.0, 0.0, 1.0])

        half_range_geom_x  = _compute_half_range(coords, center, xhat, margin)
        half_range_from_step = 0.5 * (npts - 1) * delta_x_user
        half_range_x       = max(half_range_geom_x, half_range_from_step)
        delta_x_eff        = (2.0 * half_range_x) / (npts - 1)

        half_range_y = None
        if make_y:
            half_range_geom_y = _compute_half_range(coords, center, yhat, margin)
            half_range_y      = max(half_range_geom_y, half_range_from_step)

        print(f"\n=== NICS SCAN XY (prepare) ===")
        print(f"Plane          : inertia-oriented Cartesian XY")
        print(f"HEIGHT         : {height_raw} Å (along +Z)")
        print(f"MARGIN         : {margin} Å")
        print(f"POINTS         : {npts}")
        print(f"DELTA_X (input): {delta_x_user} Å")
        print(f"DELTA_X (used) : {delta_x_eff:.6f} Å")
        if half_range_x > half_range_from_step + 1e-9:
            print(
                "[warn] POINTS and DELTA_X did not cover the molecule+MARGIN; "
                "range was expanded, DELTA_X adjusted."
            )

        cfg_flat          = flatten_cfg(cfg, section="GRID")
        cfg_flat["da_max"] = npts
        all_pts_for_xyz   = []

        for h in heights:
            xs, ptsX = _generate_hp_points(center, xhat, zhat, half_range_x, npts, h)
            all_pts_for_xyz.extend([p.tolist() for p in ptsX])
            cfg_x = dict(cfg_flat)
            cfg_x["name_molecule"] = f"{name}_hpX_h{h:.2f}".replace(".", "p")
            cfg_x["da_max"]        = npts
            make_inputs(cfg_x["name_molecule"], ats, ptsX, cfg_x)
            print(f"[ok] Inputs for X-scan generated (height={h} Å) → {cfg_x['name_molecule']}")

            if make_y and half_range_y is not None:
                ys, ptsY = _generate_hp_points(center, yhat, zhat, half_range_y, npts, h)
                all_pts_for_xyz.extend([p.tolist() for p in ptsY])
                cfg_y = dict(cfg_flat)
                cfg_y["name_molecule"] = f"{name}_hpY_h{h:.2f}".replace(".", "p")
                cfg_y["da_max"]        = npts
                make_inputs(cfg_y["name_molecule"], ats, ptsY, cfg_y)
                print(f"[ok] Inputs for Y-scan generated (height={h} Å) → {cfg_y['name_molecule']}")

        xyz_out = f"{name}_hprofiles.xyz"
        write_hp_xyz(xyz_out, ats, all_pts_for_xyz)
        print(f"[ok] Wrote horizontal profiles XYZ: {xyz_out}")

    # ---- PREDI_XY ----
    if predixy_on:
        run_predixy_prepare(cfg, ats, coords, name,
                            gen_flag=gen_flag,
                            pseudo_pi_flag=pseudo_pi_flag)


# ============================================================
#  main()  – called by aromagnetic.py dispatcher
# ============================================================

def main():
    cfg = read_input("INPUT.txt")
    if not cfg:
        raise SystemExit("ERROR: INPUT.txt not found in current directory")

    general = cfg.get("general", {})

    name = str(general.get("NAME_MOLECULE", "molecule")).strip()
    xyz_file = general.get("MOLECULAR_COORDINATES", None)

    if not xyz_file:
        raise SystemExit("ERROR: MOLECULAR_COORDINATES must be defined in INPUT.txt")
    if not os.path.isfile(xyz_file):
        raise SystemExit(f"ERROR: XYZ file not found: {xyz_file}")

    kline = str(general.get("KEYWORD_LINE", "")).strip()
    gen_flag = _gen_used(kline)
    if gen_flag and not os.path.isfile("basis.txt"):
        raise SystemExit(
            "\n[ERROR] KEYWORD_LINE contains 'GEN'/'GENECP' but 'basis.txt' was not found.\n"
            "Create 'basis.txt' in the current directory with your generic basis block.\n"
        )

    pseudo_pi_flag = _truthy(general.get("PSEUDO_PI", "NO"))

    # Orient molecule once; both profile types share the same oriented structure
    ats, _, _, _, _ = orient(xyz_file, margin=0.0)
    oriented_xyz = f"{name}_oriented.xyz"
    ase_io.write(oriented_xyz, ats, format="xyz")
    coords = ats.get_positions()

    print(f"\n=== GET PROFILE ===")
    print(f"Molecule     : {name}")
    print(f"XYZ original : {xyz_file}")
    print(f"XYZ oriented : {oriented_xyz}")

    vp_on = "VP" in cfg and bool(cfg["VP"])
    hp_on = "HP" in cfg and bool(cfg["HP"])

    if not vp_on and not hp_on:
        raise SystemExit(
            "ERROR: INPUT.txt must contain BEGIN_VP...END_VP "
            "or BEGIN_HP...END_HP (or both)."
        )

    if vp_on:
        run_vertical_profile(cfg, ats, coords, name, gen_flag, pseudo_pi_flag)

    if hp_on:
        run_horizontal_profile(cfg, ats, coords, name, gen_flag, pseudo_pi_flag)

    # Post-processing — run once after all inputs are generated to avoid
    # inserting the basis or applying pseudo-π multiple times when both
    # VP and HP are active.
    if gen_flag:
        _insert_basis(SCRIPT_DIR)
    if pseudo_pi_flag:
        _apply_pseudo_pi(ats)

    print("\nDone. Profile inputs are ready.")


if __name__ == "__main__":
    main()
