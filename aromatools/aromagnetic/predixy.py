#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import math
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ase import io as ase_io
from ase.data import atomic_numbers, covalent_radii

from pathlib import Path as _Path

from .amg_write_inputs import make_inputs, flatten_cfg
from .amg_read_input import read_input
from .tools import orient, _parse_heights
from .check_out import check_file as _check_file

# ============================================================
#  Constants
# ============================================================

_STEP    = 0.1   # Å — fixed, matches Gershoni-Poranne 2020
_COV_TOL = 1.3   # covalency factor (same as Predi-XY default)

# ============================================================
#  Small helpers
# ============================================================

def _prefix(name, height, path_idx=None):
    tag = f"h{height:.2f}".replace(".", "p")
    if path_idx is not None:
        return f"{name}_prediXY_p{path_idx:02d}_{tag}"
    return f"{name}_prediXY_{tag}"


def _parse_path_input(raw):
    """
    Parse a single PREDI_XY_PATH value.
    Returns "AUTO" or list[int].
    """
    if raw is None or isinstance(raw, bool):
        return "AUTO"
    if isinstance(raw, str) and raw.strip().upper() == "AUTO":
        return "AUTO"
    if isinstance(raw, (list, tuple)):
        try:
            return [int(x) for x in raw]
        except (TypeError, ValueError):
            return "AUTO"
    if isinstance(raw, int):
        return [raw]
    try:
        return [int(x) for x in str(raw).split()]
    except ValueError:
        return "AUTO"


def _parse_all_paths(raw):
    """
    Parse PREDI_XY_PATH from the config, handling single or multiple entries.

    The parser accumulates repeated PREDI_XY_PATH lines into a list, so raw is:
      - "AUTO"              → single AUTO path
      - [[0,1,2]]           → single manual path (new parser, one line)
      - [[0,1,2],[3,2,1]]   → two manual paths
      - [0,1,2]             → single manual path (old format, backward compat)

    Returns list of path specs, each "AUTO" or list[int].
    """
    if raw is None or (isinstance(raw, str) and raw.strip().upper() == "AUTO"):
        return ["AUTO"]
    if isinstance(raw, int):
        return [[raw]]
    if isinstance(raw, list) and raw:
        first = raw[0]
        if isinstance(first, list):
            return [_parse_path_input(p) for p in raw]
        if isinstance(first, (int, float)):
            return [[int(x) for x in raw]]
        if isinstance(first, str) and first.strip().upper() == "AUTO":
            return ["AUTO"]
    return ["AUTO"]


# ============================================================
#  Molecular graph + ring detection
# ============================================================

def _build_atom_graph(symbols, coords):
    """Atom-level NetworkX graph using covalent radii * _COV_TOL."""
    G    = nx.Graph()
    n    = len(symbols)
    G.add_nodes_from(range(n))
    Z    = [atomic_numbers[s] for s in symbols]
    Rcov = [covalent_radii[z] for z in Z]
    for i in range(n):
        for j in range(i + 1, n):
            d     = float(np.linalg.norm(coords[i] - coords[j]))
            d_cov = (Rcov[i] + Rcov[j]) * _COV_TOL
            if d <= d_cov:
                G.add_edge(i, j)
    return G


def _get_rings_and_edges(symbols, coords):
    """
    Returns
    -------
    cycles     : list[list[int]]
    centers    : (N,3) ndarray
    atom_G     : nx.Graph
    bond_edges : list[tuple]
    """
    atom_G     = _build_atom_graph(symbols, coords)
    bond_edges = list(atom_G.edges())
    # minimum_cycle_basis guarantees the smallest rings (actual ring faces),
    # unlike cycle_basis which may return large non-minimal cycles whose
    # centroids fall outside any individual ring.
    cycles     = nx.minimum_cycle_basis(atom_G)
    centers    = np.array(
        [coords[np.array(c, dtype=int)].mean(axis=0) for c in cycles],
        dtype=float,
    )
    # Sort rings left→right (X), then bottom→top (Y) for predictable numbering
    order   = sorted(range(len(cycles)), key=lambda i: (centers[i, 0], centers[i, 1]))
    cycles  = [cycles[i] for i in order]
    centers = centers[np.array(order)]
    return cycles, centers, atom_G, bond_edges


# ============================================================
#  Knot graph + AUTO path (Renana)
# ============================================================

def _build_knot_graph(cycles):
    n = len(cycles)
    G = nx.Graph()
    G.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            if set(cycles[i]) & set(cycles[j]):
                G.add_edge(i, j)
    return G


def find_knot_path_auto(cycles, centers):
    """
    AUTO path: Gershoni-Poranne 2020 strategy.
    Only works for linear catacondensed systems (exactly 2 terminal rings).
    Raises SystemExit with a helpful message if the molecule is branched.
    """
    knot_G    = _build_knot_graph(cycles)
    terminals = [v for v, d in knot_G.degree() if d == 1]

    if len(terminals) != 2:
        raise SystemExit(
            f"\nERROR: PREDI_XY_PATH AUTO only works for linear catacondensed systems.\n"
            f"Found {len(terminals)} terminal ring(s) — expected exactly 2.\n\n"
            f"For branched or pericondensed systems (e.g. pyrene, coronene),\n"
            f"run  'amg rings'  first to see the ring numbering, then set:\n"
            f"  PREDI_XY_PATH   0 1 3 2   (your desired order)\n"
            f"in the BEGIN_HP ... END_HP block of INPUT.txt."
        )

    start, end = terminals[0], terminals[1]
    path = nx.shortest_path(knot_G, source=start, target=end)

    # Orient left → right by X coordinate (Renana's core.py convention)
    if centers[path[0]][0] > centers[path[-1]][0]:
        path = list(reversed(path))

    return path


def _validate_manual_path(path_indices, n_rings, cycles):
    """
    Validate a user-supplied path:
    - All indices must be in range [0, n_rings)
    - No duplicates
    - Every consecutive pair must be physically adjacent (share atoms)
    """
    if not path_indices:
        raise SystemExit("ERROR: PREDI_XY_PATH is empty.")

    seen = set()
    for idx in path_indices:
        if idx < 0 or idx >= n_rings:
            raise SystemExit(
                f"ERROR: PREDI_XY_PATH contains index {idx} but there are "
                f"only {n_rings} rings (0-{n_rings - 1}).\n"
                f"Run 'amg rings' to see the correct numbering."
            )
        if idx in seen:
            raise SystemExit(
                f"ERROR: PREDI_XY_PATH contains duplicate index {idx}."
            )
        seen.add(idx)

    knot_G = _build_knot_graph(cycles)
    for k in range(len(path_indices) - 1):
        a, b = path_indices[k], path_indices[k + 1]
        if not knot_G.has_edge(a, b):
            raise SystemExit(
                f"ERROR: rings {a} and {b} in PREDI_XY_PATH are not adjacent "
                f"(they share no atoms).\n"
                f"Run 'amg rings' to see which rings are neighbours."
            )


# ============================================================
#  Physical path of key points  (Renana get_aroma_path)
# ============================================================

def _shared_atoms(cycles, i, j):
    return sorted(set(cycles[i]) & set(cycles[j]))


def _opposing_edge(cycle_atoms, shared, coords):
    n      = len(cycle_atoms)
    n_half = n // 2

    if n == 5:
        mid        = coords[np.array(shared)].mean(axis=0)
        candidates = [a for a in cycle_atoms if a not in shared]
        best       = max(candidates,
                         key=lambda a: float(np.linalg.norm(coords[a] - mid)))
        return [best]

    a_idx = cycle_atoms.index(shared[0])
    b_idx = cycle_atoms.index(shared[1])
    if max(a_idx, b_idx) == n - 1 and min(a_idx, b_idx) == 0:
        start_index = 0
    else:
        start_index = max(a_idx, b_idx)

    opp = [(start_index + n_half) % n, (start_index + n_half - 1) % n]
    return [cycle_atoms[o] for o in opp]


def build_aroma_path(cycles, centers, path_indices, coords):
    """
    Build the sequence of key 3-D points:
        opp_edge_0 → center_0 → shared_01 → center_1 → ... → center_N → opp_edge_N
    """
    key_pts = []

    for step, ring_idx in enumerate(path_indices):
        cyc = cycles[ring_idx]

        if step == 0:
            next_shared = _shared_atoms(cycles, ring_idx, path_indices[1])
            opp         = _opposing_edge(cyc, next_shared, coords)
            key_pts.append(coords[np.array(opp)].mean(axis=0).copy())

        key_pts.append(centers[ring_idx].copy())

        if step < len(path_indices) - 1:
            nxt    = path_indices[step + 1]
            shared = _shared_atoms(cycles, ring_idx, nxt)
            key_pts.append(coords[np.array(shared)].mean(axis=0).copy())
        else:
            prev_shared = _shared_atoms(cycles, path_indices[-2], ring_idx)
            opp         = _opposing_edge(cyc, prev_shared, coords)
            key_pts.append(coords[np.array(opp)].mean(axis=0).copy())

    return key_pts


# ============================================================
#  Interpolate at _STEP spacing
# ============================================================

def interpolate_path(key_pts, step=_STEP):
    pts      = []
    cum_dist = []
    d_acc    = 0.0

    for seg in range(len(key_pts) - 1):
        A       = np.asarray(key_pts[seg],     dtype=float)
        B       = np.asarray(key_pts[seg + 1], dtype=float)
        seg_vec = B - A
        seg_len = float(np.linalg.norm(seg_vec))

        if seg_len < 1e-10:
            continue

        n_steps = max(1, int(round(seg_len / step)))
        ts      = np.linspace(0.0, 1.0, n_steps + 1)

        for k, t in enumerate(ts):
            if k == len(ts) - 1 and seg < len(key_pts) - 2:
                continue
            p = A + t * seg_vec
            pts.append(np.array([p[0], p[1], p[2]], dtype=float))
            cum_dist.append(d_acc + t * seg_len)

        d_acc += seg_len

    return np.array(pts, dtype=float), np.array(cum_dist, dtype=float)


def lift_points(pts_plane, height):
    lifted       = pts_plane.copy()
    lifted[:, 2] += float(height)
    return lifted


# ============================================================
#  Visualization XYZ
# ============================================================

def write_predixy_xyz(path, ase_atoms, all_pts):
    real_syms = ase_atoms.get_chemical_symbols()
    real_xyz  = ase_atoms.get_positions()
    n_total   = len(real_syms) + len(all_pts)
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{n_total}\n")
        f.write("Real atoms + Predi-XY ghost atoms (X)\n")
        for sym, (x, y, z) in zip(real_syms, real_xyz):
            f.write(f"{sym:<2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")
        for (x, y, z) in all_pts:
            f.write(f"X  {x:14.8f} {y:14.8f} {z:14.8f}\n")


# ============================================================
#  Rings info PNG  (amg rings)
# ============================================================

def plot_rings_info(name, ase_atoms, coords, bond_edges, centers):
    """
    Generate a PNG showing the molecule (bonds in XY projection) with each
    ring center numbered.  Saved as {name}_rings.png.
    """
    png_out = f"{name}_rings.png"

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_aspect("equal")
    ax.axis("off")

    # Draw bonds
    for i, j in bond_edges:
        xi, yi = coords[i, 0], coords[i, 1]
        xj, yj = coords[j, 0], coords[j, 1]
        ax.plot([xi, xj], [yi, yj], color="0.35", linewidth=1.0, zorder=1)

    # Draw ring centers with index labels
    for idx, c in enumerate(centers):
        cx, cy = c[0], c[1]
        ax.plot(cx, cy, "o", color="crimson", markersize=6, zorder=3)
        ax.text(
            cx, cy + 0.25, str(idx),
            ha="center", va="bottom",
            fontsize=9, fontweight="bold",
            color="crimson", zorder=4,
        )

    ax.set_title(f"{name} — ring indices for PREDI_XY_PATH", fontsize=10, pad=8)
    plt.tight_layout()
    plt.savefig(png_out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return png_out


def rings_info():
    """
    Entry point for  amg rings.
    Reads INPUT.txt, orients the molecule, detects rings, and writes
    {name}_rings.png with each ring center numbered.
    """
    cfg = read_input("INPUT.txt")
    if not cfg:
        raise SystemExit("ERROR: INPUT.txt not found.")

    general  = cfg.get("general", {})
    name     = str(general.get("NAME_MOLECULE", "molecule")).strip()
    xyz_file = general.get("MOLECULAR_COORDINATES")

    if not xyz_file or not os.path.isfile(xyz_file):
        raise SystemExit("ERROR: MOLECULAR_COORDINATES not found.")

    ats, _, _, _, _ = orient(xyz_file, margin=0.0)
    coords   = ats.get_positions()
    symbols  = ats.get_chemical_symbols()

    cycles, centers, atom_G, bond_edges = _get_rings_and_edges(symbols, coords)

    if not cycles:
        raise SystemExit("ERROR: No rings detected in the molecule.")

    png = plot_rings_info(name, ats, coords, bond_edges, centers)

    print(f"\n=== RINGS INFO ===")
    print(f"Molecule      : {name}")
    print(f"Rings found   : {len(cycles)}")
    print(f"PNG           : {png}")
    print()
    print("Ring centers (oriented XYZ):")
    for idx, c in enumerate(centers):
        print(f"  Ring {idx:3d}:  x={c[0]:8.4f}  y={c[1]:8.4f}  z={c[2]:8.4f}"
              f"  atoms: {sorted(cycles[idx])}")
    print()
    print("To set a manual path, add in BEGIN_HP ... END_HP:")
    print(f"  PREDI_XY_PATH   {' '.join(str(i) for i in range(len(cycles)))}"
          f"   # change this order as needed")


# ============================================================
#  Metadata
# ============================================================

def _save_meta(name, path_indices, key_pts, heights, batch_size, n_pts_total, batch_records, path_idx=None):
    suffix    = f"_p{path_idx:02d}" if path_idx is not None else ""
    meta_file = f"{name}_prediXY{suffix}_path.txt"
    with open(meta_file, "w", encoding="utf-8") as f:
        f.write("# Predi-XY path metadata — DO NOT EDIT\n")
        f.write(f"step {_STEP}\n")
        f.write("path_indices " + " ".join(str(i) for i in path_indices) + "\n")
        f.write("heights "      + " ".join(f"{h:.6f}" for h in heights)  + "\n")
        f.write(f"n_keypts {len(key_pts)}\n")
        f.write(f"batch_size {batch_size}\n")
        f.write(f"n_pts_total {n_pts_total}\n")
        for k, p in enumerate(key_pts):
            f.write(f"keypt {k} {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}\n")
        for h, b, start, end in batch_records:
            h_tag = f"h{h:.2f}".replace(".", "p")
            f.write(f"batch {h_tag} {b} {start} {end}\n")
    return meta_file


def _load_meta(name, path_idx=None):
    suffix    = f"_p{path_idx:02d}" if path_idx is not None else ""
    meta_file = f"{name}_prediXY{suffix}_path.txt"
    if not os.path.isfile(meta_file):
        raise SystemExit(
            f"ERROR: {meta_file} not found. "
            "Run 'amg run' first to generate Predi-XY inputs."
        )
    key_pts       = []
    path_indices  = []
    heights       = []
    step          = _STEP
    batch_size    = 80
    n_pts_total   = 0
    batch_records = []

    with open(meta_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            key   = parts[0]
            if key == "step":
                step = float(parts[1])
            elif key == "path_indices":
                path_indices = [int(x) for x in parts[1:]]
            elif key == "heights":
                heights = [float(x) for x in parts[1:]]
            elif key == "batch_size":
                batch_size = int(parts[1])
            elif key == "n_pts_total":
                n_pts_total = int(parts[1])
            elif key == "keypt":
                key_pts.append(
                    np.array([float(parts[2]), float(parts[3]), float(parts[4])])
                )
            elif key == "batch":
                batch_records.append((parts[1], int(parts[2]), int(parts[3]), int(parts[4])))

    return key_pts, path_indices, heights, step, batch_size, n_pts_total, batch_records


# ============================================================
#  PREPARE  (amg run)
# ============================================================

def run_predixy_prepare(cfg, ats, coords, name,
                        gen_flag=False, pseudo_pi_flag=False):
    """
    Called from get_profile.py when PREDI_XY = YES inside BEGIN_HP.
    """
    hp_block   = cfg.get("HP", {})
    heights    = _parse_heights(hp_block.get("HEIGHT", 1.0))
    path_raw   = hp_block.get("PREDI_XY_PATH", "AUTO")
    all_paths  = _parse_all_paths(path_raw)
    batch_size = int(hp_block.get("BATCH_SIZE", 80))
    multi_path = len(all_paths) > 1
    symbols    = ats.get_chemical_symbols()

    cycles, centers, atom_G, bond_edges = _get_rings_and_edges(symbols, coords)

    if len(cycles) < 2:
        raise SystemExit(
            "ERROR: PREDI_XY requires at least 2 rings. "
            "Check MOLECULAR_COORDINATES."
        )

    print(f"\n=== PREDI-XY (prepare) ===")
    print(f"Rings detected : {len(cycles)}")
    print(f"Paths          : {len(all_paths)}")
    print(f"Heights (Å)    : {heights}")
    print(f"Batch size     : {batch_size} DA/input")

    cfg_flat = flatten_cfg(cfg, section="GRID")

    for p_idx, path_spec in enumerate(all_paths):
        p_label = p_idx if multi_path else None

        if path_spec == "AUTO":
            print(f"\n[path {p_idx}] AUTO → Gershoni-Poranne 2020")
            path_indices = find_knot_path_auto(cycles, centers)
        else:
            path_indices = path_spec
            print(f"\n[path {p_idx}] manual → {path_indices}")
            _validate_manual_path(path_indices, len(cycles), cycles)

        key_pts = build_aroma_path(cycles, centers, path_indices, coords)
        pts_plane, cum_dist = interpolate_path(key_pts, step=_STEP)

        n_total   = len(pts_plane)
        n_batches = math.ceil(n_total / batch_size)

        print(f"  Ring path  : {' → '.join(str(i) for i in path_indices)}")
        print(f"  DA/scan    : {n_total}  ({n_batches} batch(es) per height)")
        print(f"  Length     : {cum_dist[-1]:.4f} Å")

        all_pts_for_xyz = []
        batch_records   = []

        for h in heights:
            pts_3d = lift_points(pts_plane, h)
            all_pts_for_xyz.extend(pts_3d.tolist())
            prefix = _prefix(name, h, path_idx=p_label)

            for b in range(n_batches):
                start     = b * batch_size
                end       = min(start + batch_size, n_total)
                batch_pts = pts_3d[start:end]

                batch_prefix           = f"{prefix}_b{b:04d}"
                cfg_h                  = dict(cfg_flat)
                cfg_h["name_molecule"] = batch_prefix
                cfg_h["da_max"]        = len(batch_pts)
                make_inputs(batch_prefix, ats, batch_pts, cfg_h)
                batch_records.append((h, b, start, end - 1))
                print(
                    f"  [ok] {batch_prefix}"
                    f"  (h={h:.2f} Å, DA {start}–{end - 1})"
                )

        xyz_label = f"_p{p_idx:02d}" if multi_path else ""
        xyz_out   = f"{name}_prediXY{xyz_label}.xyz"
        write_predixy_xyz(xyz_out, ats, all_pts_for_xyz)
        print(f"  [ok] reference XYZ → {xyz_out}")

        meta = _save_meta(
            name, path_indices, key_pts, heights,
            batch_size, n_total, batch_records, path_idx=p_label,
        )
        print(f"  [ok] metadata → {meta}")


# ============================================================
#  Overlay plot  (multiple paths on the same axes)
# ============================================================

_OVERLAY_COLORS = [
    "steelblue", "crimson", "forestgreen", "darkorange",
    "purple", "saddlebrown", "teal", "deeppink",
]


def _plot_predixy_overlay(name, all_collected, interactive=False):
    """
    all_collected : list of (path_idx, collected)
    collected     : list of (h, dist, nics)
    One subplot per unique height; all paths drawn on each subplot.
    """
    all_heights = sorted({h for _, col in all_collected for h, _, _ in col})

    for h in all_heights:
        h_tag  = f"h{h:.2f}".replace(".", "p")
        prefix = f"{name}_prediXY_overlay_{h_tag}"

        fig, ax = plt.subplots(figsize=(9, 4))
        ax.axhline(0, color="black", linewidth=0.6, linestyle="--")

        for color, (path_idx, collected) in zip(_OVERLAY_COLORS, all_collected):
            p_label = f"p{path_idx:02d}" if path_idx is not None else "path"
            for ch, dist, nics, *_ in collected:
                if abs(ch - h) < 1e-6:
                    ax.plot(dist, nics, color=color, linewidth=1.5, label=p_label)
                    break

        ax.set_xlabel("Path distance (Å)")
        ax.set_ylabel("NICS (ppm)")
        ax.set_title(f"{name} — PREDI-XY overlay  (h = {h:.2f} Å)")
        ax.legend(loc="best", fontsize=8)
        fig.tight_layout()
        png_name = f"{prefix}.png"
        fig.savefig(png_name, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"[ok] overlay → {png_name}")

        if interactive:
            import plotly.graph_objects as go
            fig_px = go.Figure()
            fig_px.add_hline(y=0, line=dict(color="black", width=1, dash="dash"))
            for color, (path_idx, collected) in zip(_OVERLAY_COLORS, all_collected):
                p_label = f"p{path_idx:02d}" if path_idx is not None else "path"
                for ch, dist, nics, *_ in collected:
                    if abs(ch - h) < 1e-6:
                        fig_px.add_trace(go.Scatter(
                            x=dist, y=nics, mode="lines", name=p_label,
                            line=dict(width=2),
                            hovertemplate=(
                                f"[{p_label}] dist: %{{x:.3f}} Å"
                                "<br>NICS: %{y:.3f} ppm<extra></extra>"
                            ),
                        ))
                        break
            fig_px.update_layout(
                xaxis_title="Path distance (Å)",
                yaxis_title="NICS (ppm)",
                title=f"{name} — PREDI-XY overlay  (h = {h:.2f} Å)",
            )
            html_name = f"{prefix}.html"
            fig_px.write_html(html_name)
            print(f"[ok] overlay (interactive) → {html_name}")


# ============================================================
#  3D plot helper
# ============================================================

_3D_COLORS = [
    "steelblue", "crimson", "forestgreen", "darkorange",
    "purple", "saddlebrown", "teal", "deeppink",
]


def _plot_predixy_3d(name, collected, interactive=False, mol_data=None):
    """
    collected : list of (height, dist_list, nics_list, xy_arr).
    X,Y = actual molecular XY coordinates of each scan point.
    Z   = NICS value.
    mol_data  : optional dict with keys 'coords', 'bonds' for the molecular skeleton.
    """
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    all_nics_flat = [n for _, _, nics, _ in collected for n in nics]
    z_mol = min(all_nics_flat) - 1.0

    png_name = f"{name}_prediXY_3D.png"
    fig_mpl  = plt.figure(figsize=(9, 7))
    ax       = fig_mpl.add_subplot(111, projection="3d")

    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.fill = False
        pane.set_edgecolor("none")
    ax.grid(False)

    if mol_data is not None:
        coords_m = np.array(mol_data["coords"])
        for i, j in mol_data["bonds"]:
            ax.plot(
                [coords_m[i, 0], coords_m[j, 0]],
                [coords_m[i, 1], coords_m[j, 1]],
                [z_mol, z_mol],
                color="0.45", linewidth=0.9, alpha=0.7,
            )

    for i, (h, dist, nics, xy) in enumerate(collected):
        color = _3D_COLORS[i % len(_3D_COLORS)]
        ax.plot(xy[:, 0], xy[:, 1], np.array(nics),
                color=color, linewidth=1.5, label=f"h = {h:.2f} Å")

    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_zlabel("NICS (ppm)")
    ax.set_title(f"{name} — PREDI-XY 3D")
    if len(collected) > 1:
        ax.legend(loc="best", fontsize=8)
    fig_mpl.tight_layout()
    fig_mpl.savefig(png_name, dpi=200, bbox_inches="tight")
    plt.close(fig_mpl)

    if interactive:
        import plotly.graph_objects as go

        fig_px = go.Figure()

        if mol_data is not None:
            coords_m = np.array(mol_data["coords"])
            bx, by, bz = [], [], []
            for i, j in mol_data["bonds"]:
                bx += [float(coords_m[i, 0]), float(coords_m[j, 0]), None]
                by += [float(coords_m[i, 1]), float(coords_m[j, 1]), None]
                bz += [z_mol, z_mol, None]
            fig_px.add_trace(go.Scatter3d(
                x=bx, y=by, z=bz,
                mode="lines",
                line=dict(color="gray", width=3),
                name="molecule",
                showlegend=False,
                hoverinfo="skip",
            ))

        for i, (h, dist, nics, xy) in enumerate(collected):
            nics_arr = np.array(nics)
            if len(collected) == 1:
                line_kw = dict(
                    color=nics_arr.tolist(),
                    colorscale="RdBu_r",
                    reversescale=True,
                    width=6,
                    colorbar=dict(title="NICS (ppm)"),
                )
            else:
                line_kw = dict(color=_3D_COLORS[i % len(_3D_COLORS)], width=5)
            fig_px.add_trace(go.Scatter3d(
                x=xy[:, 0].tolist(),
                y=xy[:, 1].tolist(),
                z=nics_arr.tolist(),
                mode="lines",
                name=f"h = {h:.2f} Å",
                line=line_kw,
                hovertemplate=(
                    f"h = {h:.2f} Å<br>"
                    "X: %{{x:.3f}} Å<br>Y: %{{y:.3f}} Å<br>"
                    "NICS: %{{z:.3f}} ppm<extra></extra>"
                ),
            ))

        _ax_clean = dict(
            showgrid=False,
            showbackground=False,
            showspikes=False,
            showline=True,
            linecolor="black",
            linewidth=2,
            showticklabels=True,
        )
        fig_px.update_layout(
            scene=dict(
                xaxis=dict(**_ax_clean, title="X (Å)"),
                yaxis=dict(**_ax_clean, title="Y (Å)"),
                zaxis=dict(**_ax_clean, title="NICS (ppm)"),
            ),
            title=f"{name} — PREDI-XY 3D",
        )
        html_3d = f"{name}_prediXY_3D.html"
        fig_px.write_html(html_3d)
        return png_name, html_3d

    return png_name, None


# ============================================================
#  ANALYZE  (amg get)
# ============================================================

def _analyze_one_path(name, path_idx, prog, n_real, interactive,
                      read_gaussian_fn, read_orca_fn):
    """Analyze a single PREDI-XY path. Returns list of (h, dist, nics)."""
    key_pts, path_indices, heights, step, batch_size, n_pts_total, batch_records = \
        _load_meta(name, path_idx=path_idx)

    if not batch_records:
        print(f"[warn] metadata for path {path_idx} has no batch records — skipping")
        return []

    pts_plane, cum_dist = interpolate_path(key_pts, step=step)
    p_label = f"p{path_idx:02d}" if path_idx is not None else None

    print(f"\n{'--- ' if p_label else ''}PREDI-XY path: "
          f"{' → '.join(str(i) for i in path_indices)}")

    collected = []

    for h in heights:
        h_tag  = f"h{h:.2f}".replace(".", "p")
        prefix = _prefix(name, h, path_idx=path_idx)

        h_batches = sorted(
            [(b, start, end) for (ht, b, start, end) in batch_records if ht == h_tag],
            key=lambda x: x[0],
        )
        if not h_batches:
            print(f"  [warn] No batches for height {h:.2f} Å — skipping")
            continue

        iso_all = []
        for b, start, end in h_batches:
            batch_prefix = f"{prefix}_b{b:04d}"
            out_file = None
            for ext in (".log", ".out"):
                candidate = batch_prefix + ext
                if os.path.isfile(candidate):
                    out_file = candidate
                    break
            if out_file is None:
                raise SystemExit(
                    f"ERROR: output file not found for {batch_prefix} (tried .log and .out)"
                )
            if not _check_file(_Path(out_file)).startswith("OK"):
                print(f"[WARN] {out_file}: calculation did not terminate normally")
            if "gaussian" in prog:
                iso_all.extend(read_gaussian_fn([out_file]))
            elif "orca" in prog:
                iso_all.extend(read_orca_fn([out_file], n_real))
            else:
                raise SystemExit("ERROR: PROGRAM must be GAUSSIAN or ORCA.")

        n_use = min(len(cum_dist), len(iso_all))
        nics  = [-float(v) for v in iso_all[:n_use]]
        dist  = list(cum_dist[:n_use])
        xy    = pts_plane[:n_use, :2].copy()
        collected.append((h, dist, nics, xy))

        csv_name = f"{prefix}.csv"
        with open(csv_name, "w", encoding="utf-8") as f:
            f.write("dist_A,NICS_ppm\n")
            for d, n in zip(dist, nics):
                f.write(f"{d:.8f},{n:.8f}\n")

        png_name = f"{prefix}.png"
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(dist, nics, color="steelblue", linewidth=1.5)
        ax.axhline(0, color="black", linewidth=0.6, linestyle="--")
        ax.set_xlabel("Path distance (Å)")
        ax.set_ylabel("NICS (ppm)")
        ax.set_title(f"{name} — PREDI-XY  (h = {h:.2f} Å)"
                     + (f"  [{p_label}]" if p_label else ""))
        fig.tight_layout()
        fig.savefig(png_name, dpi=200, bbox_inches="tight")
        plt.close(fig)

        out_files_str = f"{csv_name}, {png_name}"
        if interactive:
            import plotly.graph_objects as go
            fig_px = go.Figure()
            fig_px.add_trace(go.Scatter(
                x=dist, y=nics, mode="lines",
                line=dict(color="steelblue", width=2),
                hovertemplate="dist: %{x:.3f} Å<br>NICS: %{y:.3f} ppm<extra></extra>",
            ))
            fig_px.add_hline(y=0, line=dict(color="black", width=1, dash="dash"))
            fig_px.update_layout(
                xaxis_title="Path distance (Å)",
                yaxis_title="NICS (ppm)",
                title=f"{name} — PREDI-XY  (h = {h:.2f} Å)"
                      + (f"  [{p_label}]" if p_label else ""),
            )
            html_name = f"{prefix}.html"
            fig_px.write_html(html_name)
            out_files_str += f", {html_name}"

        print(f"  [ok] h={h:.2f} Å → {out_files_str}  ({n_use} pts)")

    return collected


def run_predixy_analyze(cfg, name, program, xyz_file, n_real, interactive=False):
    from .get_results import read_gaussian_bq_isotropic, read_orca_dummy_isotropic

    prog = program.lower()

    print(f"\n=== PREDI-XY (analyze) ===")
    print(f"Molecule : {name}  |  Program : {program}")

    mol_data = None
    if xyz_file and os.path.isfile(str(xyz_file)):
        ats_viz, _, _, _, _ = orient(str(xyz_file), margin=0.0)
        coords_viz  = ats_viz.get_positions()
        symbols_viz = ats_viz.get_chemical_symbols()
        _, _, _, bond_edges_viz = _get_rings_and_edges(symbols_viz, coords_viz)
        mol_data = {"coords": coords_viz, "bonds": bond_edges_viz}

    # Discover metadata files: single-path (old naming) or multi-path (p00, p01, …)
    all_path_indices = []
    if os.path.isfile(f"{name}_prediXY_path.txt"):
        all_path_indices.append(None)
    p = 0
    while os.path.isfile(f"{name}_prediXY_p{p:02d}_path.txt"):
        all_path_indices.append(p)
        p += 1

    if not all_path_indices:
        raise SystemExit(
            f"ERROR: No metadata file found for '{name}'.\n"
            "Run 'amg run' first."
        )

    print(f"Paths found : {len(all_path_indices)}")

    all_collected = []
    for path_idx in all_path_indices:
        collected = _analyze_one_path(
            name, path_idx, prog, n_real, interactive,
            read_gaussian_bq_isotropic, read_orca_dummy_isotropic,
        )
        if collected:
            all_collected.append((path_idx, collected))

    for path_idx, collected in all_collected:
        p_label = f"p{path_idx:02d}" if path_idx is not None else None
        suffix  = f"_{p_label}" if p_label else ""
        png_3d, html_3d = _plot_predixy_3d(
            f"{name}{suffix}", collected, interactive=interactive, mol_data=mol_data
        )
        msg = f"[ok] 3D{suffix} → {png_3d}"
        if html_3d:
            msg += f", {html_3d}"
        print(msg)

    if len(all_collected) > 1:
        _plot_predixy_overlay(name, all_collected, interactive=interactive)

    print("\nDone.")
