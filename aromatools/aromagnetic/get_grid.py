#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Python modules
import os
import re
import sys
import subprocess
import importlib.util
import numpy as np

from ase import io
from pymatgen.core import Molecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer

# AromaTools modules
from . import pseudo_pi
from .amg_read_input import read_input
from .amg_write_inputs import make_inputs

script_dir = os.path.dirname(os.path.abspath(__file__))

if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

def only_ch(mol_file):
    try:
        with open(mol_file, "r", encoding="utf-8", errors="ignore") as fh:
            lines = fh.readlines()
    except FileNotFoundError:
        return False
    if len(lines) < 3:
        return False
    for line in lines[2:]:
        parts = line.strip().split()
        if not parts:
            continue
        elem = parts[0].upper()
        if elem in ("X", "XX", "BQ", "DUM", "DUMMY"):
            continue
        if elem not in ("C", "H"):
            return False
    return True

# ---------------- IMF helpers (raw scan) ----------------
def _imf_find_tokens(key, imf_path="INPUT.txt"):
    key_u = key.strip().upper()
    try:
        with open(imf_path, "r", encoding="utf-8", errors="ignore") as fh:
            for raw in fh:
                s = raw.strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split()
                if not parts:
                    continue
                if parts[0].upper() == key_u:
                    return parts[1:]
    except FileNotFoundError:
        return None
    return None

def _imf_truthy(tok):
    if tok is None:
        return False
    t = str(tok).strip().upper()
    return t in ("YES", "Y", "TRUE", "T", "1", "ON")

def _imf_get_bool(key, default=False, imf_path="INPUT.txt"):
    toks = _imf_find_tokens(key, imf_path=imf_path)
    if not toks:
        return default
    return _imf_truthy(toks[0])

def _imf_get_float(key, default=None, imf_path="INPUT.txt"):
    toks = _imf_find_tokens(key, imf_path=imf_path)
    if not toks:
        return default
    try:
        return float(toks[0])
    except Exception:
        return default

def orient_molecule(xyz_path, margin=4.0):
    ats = io.read(xyz_path, index=-1)
    m = ats.get_masses()
    R = ats.get_positions()
    cm = np.average(R, axis=0, weights=m)
    Rc = R - cm
    I = np.zeros((3, 3))
    for mi, ri in zip(m, Rc):
        x, y, z = ri
        r2 = x*x + y*y + z*z
        I += mi * (r2 * np.eye(3) - np.outer(ri, ri))
    w, V = np.linalg.eigh(I)
    idx = np.argsort(w)
    V = V[:, idx]
    Ror = Rc @ V
    ats.set_positions(Ror)
    xmin, ymin, zmin = Ror.min(axis=0)
    xmax, ymax, zmax = Ror.max(axis=0)
    bbox_nomargin = ((xmin, xmax), (ymin, ymax), (zmin, zmax))
    xmin -= margin
    ymin -= margin
    zmin -= margin
    xmax += margin
    ymax += margin
    zmax += margin
    bbox = ((xmin, xmax), (ymin, ymax), (zmin, zmax))
    return ats, Ror, bbox, bbox_nomargin

def _snap_vertical_mirror(coords, pga, atol=1e-3):
    """
    For molecules with degenerate in-plane moments of inertia (D3h, C3v,
    C6v, …) the principal-axis orientation leaves the XY rotation
    undetermined.  As a result, σv planes land at arbitrary angles that
    don't coincide with X or Y, and axis-aligned mirror detection misses
    them.

    This function finds all σv plane normals (those lying in the XY plane)
    from the pymatgen PGA, then applies the smallest rotation around Z that
    snaps the closest one onto the X or Y axis. After rotation the
    axis-aligned detector can find it.

    Only acts on planar molecules (all |z| < 0.1 Å). Returns
    (new_coords, angle_rad); if no snap is needed angle_rad == 0.
    """
    # Guard: only snap for planar molecules
    if np.max(np.abs(coords[:, 2])) > 0.1:
        return coords, 0.0

    v_normals = []
    for op in pga.get_symmetry_operations():
        Rm = op.rotation_matrix
        if abs(np.linalg.det(Rm) + 1.0) > 0.15:
            continue                         # not improper
        if not np.allclose(Rm @ Rm, np.eye(3), atol=atol):
            continue                         # not a pure reflection (e.g. Sn)
        # Normal of this σ = eigenvector of R with eigenvalue -1
        S = 0.5 * (np.eye(3) - Rm)
        w, V = np.linalg.eigh(S)
        n = V[:, np.argmax(w)]
        n /= np.linalg.norm(n)
        if abs(n[2]) < 0.15:               # σv: normal is in the XY plane
            v_normals.append(n)

    if not v_normals:
        return coords, 0.0                   # no σv → nothing to do

    # Find the smallest rotation around Z that maps any σv normal onto
    # one of the four Cartesian half-axes in XY (0°, 90°, 180°, 270°).
    best_angle, best_score = 0.0, float('inf')
    for n in v_normals:
        theta = np.arctan2(float(n[1]), float(n[0]))
        for target in [0.0, np.pi / 2, np.pi, -np.pi / 2]:
            angle = target - theta
            angle = (angle + np.pi) % (2 * np.pi) - np.pi   # → [-π, π]
            if abs(angle) < best_score:
                best_score = abs(angle)
                best_angle = angle

    if best_score < 1e-8:
        return coords, 0.0                   # already aligned

    c, s = np.cos(best_angle), np.sin(best_angle)
    Rz = np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
    return coords @ Rz.T, best_angle


def detect_mirror_axes(coords, species, sym_tol=1e-3):
    mol = Molecule(species, coords)
    pga = PointGroupAnalyzer(mol, tolerance=sym_tol)
    if hasattr(pga, "sch_symbol"):
        pg_symbol = pga.sch_symbol                 
    elif hasattr(pga, "get_pointgroup"):
        pg_symbol = pga.get_pointgroup()          
    else:
        pg_symbol = "UNKNOWN"
    axes = set()
    for op in pga.get_symmetry_operations():
        R = op.rotation_matrix
        detR = np.linalg.det(R)
        if detR < -0.9:
            if not np.allclose(R @ R, np.eye(3), atol=1e-3):
                continue
            if np.allclose(R, np.diag([-1, 1, 1]), atol=1e-3):
                axes.add('x')
            elif np.allclose(R, np.diag([1, -1, 1]), atol=1e-3):
                axes.add('y')
            elif np.allclose(R, np.diag([1, 1, -1]), atol=1e-3):
                axes.add('z')

    return pg_symbol, axes

def choose_grid_symmetry(axes):
    n = len(axes)
    if n == 0:
        return "FULL"
    elif n == 1:
        return "HALF"
    elif n == 2:
        return "QUARTER"
    else:
        return "OCT"

def make_full_box(bbox, cube_min=4.0):
    (xmin, xmax), (ymin, ymax), (zmin, zmax) = bbox
    Lx = xmax - xmin
    Ly = ymax - ymin
    Lz = zmax - zmin
    Lmax = max(Lx, Ly, Lz, cube_min)
    half = 0.5 * Lmax
    return (-half, half), (-half, half), (-half, half)

def fundamental_ranges(full_box, axes):
    (xmin, xmax), (ymin, ymax), (zmin, zmax) = full_box
    def half_pos(lo, hi):
        return (0.0, hi)
    rx = (xmin, xmax)
    ry = (ymin, ymax)
    rz = (zmin, zmax)
    if 'x' in axes:
        rx = half_pos(*rx)
    if 'y' in axes:
        ry = half_pos(*ry)
    if 'z' in axes:
        rz = half_pos(*rz)

    return rx, ry, rz

def build_regular_grid(xrng, yrng, zrng, Nx, Ny, Nz):
    (xmin, xmax) = xrng
    (ymin, ymax) = yrng
    (zmin, zmax) = zrng
    def step(lo, hi, n):
        if n <= 1:
            return 0.0
        return (hi - lo) / (n - 1)

    dx = step(xmin, xmax, Nx)
    dy = step(ymin, ymax, Ny)
    dz = step(zmin, zmax, Nz)

    pts = []
    x = xmin
    for ix in range(Nx):
        y = ymin
        for iy in range(Ny):
            z = zmin
            for iz in range(Nz):
                pts.append((x, y, z))
                z += dz
            y += dy
        x += dx

    return np.array(pts, float), dx, dy, dz

def replicate_grid(P, axes):
    
    if not axes:
        return P
    axis_to_idx = {'x': 0, 'y': 1, 'z': 2}
    mirrored_axes = [axis_to_idx[a] for a in sorted(axes)]
    combos = []
    m = len(mirrored_axes)
    for mask in range(2 ** m):
        s = [1, 1, 1]
        for bit in range(m):
            idx = mirrored_axes[bit]
            if (mask >> bit) & 1:
                s[idx] = -1
        combos.append(tuple(s))

    pts = []
    for x, y, z in P:
        for sx, sy, sz in combos:
            xx = sx * x
            yy = sy * y
            zz = sz * z
            pts.append((xx, yy, zz))

    full = np.array(pts, float)
    full = np.unique(full, axis=0)
    idx = np.lexsort((full[:, 2], full[:, 1], full[:, 0]))
    return full[idx]

def write_xyz(path, base_atoms, pts):
   
    symbols = base_atoms.get_chemical_symbols()
    coords = base_atoms.get_positions()
    n_atoms = len(symbols)
    n_pts = len(pts)
    total = n_atoms + n_pts
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{total}\n")
        f.write("Grid generated by amg_new (symmetry-aware)\n")
        for sym, (x, y, z) in zip(symbols, coords):
            f.write(f"{sym:<2s} {x:15.8f} {y:15.8f} {z:15.8f}\n")
        for x, y, z in pts:
            f.write(f"X  {x:15.8f} {y:15.8f} {z:15.8f}\n")

# ---------------- main ----------------

def main():
    cfg = read_input("INPUT.txt")
    if not cfg:
        raise SystemExit("error: INPUT.txt not found in current directory")

    general = cfg.get("general", {})
    grid_blk = cfg.get("GRID", {})
    name = general.get("NAME_MOLECULE", "molecule")
    mol_file = general.get("MOLECULAR_COORDINATES") or f"{name}.xyz"
    margin = float(grid_blk.get("GRID_MARGIN", 4.0))
    cube_min = float(grid_blk.get("CUBE_MIN", 4.0))
    pseudo_pi_used = False # si se uso pseudo_pi
    kline = (general.get("KEYWORD_LINE") or "").strip()
    gen_used = bool(re.search(r"\bgen(ecp)?\b", kline, flags=re.IGNORECASE))
    basis_file = "basis.txt"

    if gen_used:
        if not os.path.isfile(basis_file):
            raise SystemExit(
                "\n[ERROR] KEYWORD_LINE contains 'GEN'/'GENECP' but 'basis.txt' was not found.\n"
                "Create 'basis.txt' in the current directory with your generic basis block, then rerun.\n"
                "Tip: basis.txt will be appended to every generated *.inp input.\n"
            )

    num_points = grid_blk.get("NUM_POINTS", None)
    if not num_points or not isinstance(num_points, (list, tuple)) or len(num_points) != 3:
        raise SystemExit("ERROR: NUM_POINTS Nx Ny Nz must be defined in INPUT.txt within BEGIN_GRID ... END_GRID")

    Nx, Ny, Nz = map(int, num_points)

    imf = {
        "name_molecule": name,
        "mol_file": mol_file,
        "grid_margin": margin,
        "cube_min": cube_min,
        "keyword_line": general.get("KEYWORD_LINE", ""),
        "program": general.get("PROGRAM", "ORCA"),
        "procs": int(general.get("PROCS", 1) or 1),
        "memory": int(general.get("MEMORY_GB", 16) or 16),
        "da_max": int(grid_blk.get("GA_MAX", 20) or 20),
        "charge": int(general.get("CHARGE", 0) or 0),
        "multiplicity": int(general.get("MULTIPLICITY", 1) or 1),
        "orca_block": (cfg.get("ORCA_BLOCK", {}).get("raw", "") if isinstance(cfg.get("ORCA_BLOCK", {}), dict) else ""),
        "num_points": (Nx, Ny, Nz),
    }

    if Nx == 1 and Ny > 1 and Nz > 1:
        grid_type = "YZ plane"
        print("[info] A YZ plane will be constructed (Nx = 1)!")
    elif Ny == 1 and Nx > 1 and Nz > 1:
        grid_type = "XZ plane"
        print("[info] A XZ plane will be constructed (Ny = 1)!")
    elif Nz == 1 and Nx > 1 and Ny > 1:
        grid_type = "XY plane"
        print("[info] A XY plane will be constructed (Nz = 1)!")
    elif Nx == 1 and Ny == 1 and Nz > 1:
        grid_type = "Z line"
        print("[warning] Grid degenerates into a Z line (Nx = Ny = 1).")
    elif Nx == 1 and Nz == 1 and Ny > 1:
        grid_type = "Y line"
        print("[warning] Grid degenerates into a Y line (Nx = Nz = 1).")
    elif Ny == 1 and Nz == 1 and Nx > 1:
        grid_type = "X line"
        print("[warning] Grid degenerates into an X line (Ny = Nz = 1).")
    else:
        grid_type = "3D grid"
        print("[info] A full 3D grid will be constructed!")

    Nx = int(Nx); Ny = int(Ny); Nz = int(Nz)
    if Nx < 1 or Ny < 1 or Nz < 1:
        raise SystemExit("ERROR: NUM_POINTS must be positive integers")

    ats, coords, bbox0, bbox_nomargin = orient_molecule(mol_file, margin=margin)
    species = ats.get_chemical_symbols()

    # 2a) Snap σv planes onto Cartesian axes when in-plane moments are
    #     degenerate (D3h, C3v, C6v, …). Inertia orientation leaves the
    #     XY rotation undetermined in those cases, so σv planes land at
    #     arbitrary angles. A minimal rotation around Z fixes this.
    _mol_pre  = Molecule(species, coords)
    _pga_pre  = PointGroupAnalyzer(_mol_pre, tolerance=1e-3)
    coords, snap_angle = _snap_vertical_mirror(coords, _pga_pre)
    if abs(snap_angle) > 1e-8:
        ats.set_positions(coords)
        # Recompute molecule bbox without margin after rotation
        _cmin = coords.min(axis=0)
        _cmax = coords.max(axis=0)
        bbox_nomargin = (
            (float(_cmin[0]), float(_cmax[0])),
            (float(_cmin[1]), float(_cmax[1])),
            (float(_cmin[2]), float(_cmax[2])),
        )
        print(f"[info] σv alignment: applied {np.degrees(snap_angle):.3f}° rotation around Z")

    oriented_xyz = f"{name}_oriented.xyz"
    io.write(oriented_xyz, ats, format="xyz")

    # 2b) Detect point group and axis-aligned reflection planes
    pg, axes = detect_mirror_axes(coords, species)

    # 3) Choose a grid
    gsym = choose_grid_symmetry(axes)

    # 4) All box
    full_box = make_full_box(bbox0, cube_min=cube_min)

    # 5) Fundamental interval based on symmetry
    xrng, yrng, zrng = fundamental_ranges(full_box, axes)

    # 6) Build a grid
    P, dx, dy, dz = build_regular_grid(xrng, yrng, zrng, Nx, Ny, Nz)

    # 7) Complete grid!
    Pf = replicate_grid(P, axes)

    # 8) Write XYZ
    calc_xyz = f"{name}in_calcgrid.xyz"
    full_xyz = f"{name}in_symgrid.xyz"

    write_xyz(calc_xyz, ats, P)
    write_xyz(full_xyz, ats, Pf)

    # 9) Input's count
    n_inp = 0
    try:
        n_inp = make_inputs(name, ats, P, imf)
        if gen_used:
            try:
                subprocess.check_call([sys.executable, os.path.join(script_dir, "insert_basis.py")])
                print("[ok] basis.txt was appended automatically to all generated *.inp files.")
            except Exception as e:
                raise SystemExit(
                    f"\n[ERROR] GEN was requested and basis.txt exists, but automatic insertion failed: {e}\n"
                    "You can run it manually: python insert_basis.py\n"
                )
    except Exception as e:
        print(f"[WARNING] could not generate inputs: {e}")

    # 9b) Optional pseudo-π model (edits generated *.inp files)
    if _imf_get_bool("PSEUDO_PI", default=False):
        if not only_ch(mol_file):
            print("[WARNING] PSEUDO_PI requested, but molecule contains atoms other than C/H. Skipping pseudo-π (only valid for C/H systems).")
        else:
            try:
                pseudo_pi.main()
                pseudo_pi_used = True
                # print(f"Pseudo-π model was applied!")
            except Exception as e:
                print(f"[WARNING] PSEUDO_PI YES but could not apply pseudo-π: {e}")

    # 10) grid_info.txt
    origin = tuple(P[0].tolist())
    Lx = xrng[1] - xrng[0]
    Ly = yrng[1] - yrng[0]
    Lz = zrng[1] - zrng[0]

    info = []
    info.append("=== GRID INFO ===")
    info.append(f"NAME_MOLECULE            : {name}")
    info.append(f"POINT_GROUP              : {pg}")
    info.append(f"GRID_SYMMETRY            : {gsym}")
    info.append(f"MIRROR PLANES            : {sorted(axes) if axes else 'none'}")
    info.append(f"CALC_XYZ                 : {calc_xyz}")
    info.append(f"FULL_XYZ                 : {full_xyz}")
    info.append(f"NEW MOLECULAR FILE       : {oriented_xyz}")
    info.append(f"SHAPE(N)                 : {Nx} {Ny} {Nz}")
    info.append(f"GRID_TYPE                : {grid_type}")
    info.append(f"STEP(angstroms)          : {dx:.6f} {dy:.6f} {dz:.6f}")
    info.append(f"AXIS_LENGTHS(angstroms)  : {Lx:.6f} {Ly:.6f} {Lz:.6f}")
    info.append(f"GRID ORIGIN (angstroms)  : {origin[0]:.6f} {origin[1]:.6f} {origin[2]:.6f}")
    info.append(f"GRID_MARGIN (angstroms)  : {margin:.2f}")
    # info.append(f"CUBE_MIN (angstroms)     : {cube_min:.2f}")
    info.append("\n")
    (x0, x1), (y0, y1), (z0, z1) = bbox_nomargin
    info.append("BBOX_MOLECULE (angstroms)")
    info.append(f"  x : {x0:.6f}  {x1:.6f}   Lx = {x1 - x0:.6f}")
    info.append(f"  y : {y0:.6f}  {y1:.6f}   Ly = {y1 - y0:.6f}")
    info.append(f"  z : {z0:.6f}  {z1:.6f}   Lz = {z1 - z0:.6f}")
    info.append("\n")
    (x0, x1), (y0, y1), (z0, z1) = bbox0
    info.append("BBOX_WITH_MARGIN (angstroms)")
    info.append(f"  x : {x0:.6f}  {x1:.6f}   Lx = {x1 - x0:.6f}")
    info.append(f"  y : {y0:.6f}  {y1:.6f}   Ly = {y1 - y0:.6f}")
    info.append(f"  z : {z0:.6f}  {z1:.6f}   Lz = {z1 - z0:.6f}")
    info.append("\n")
    if pseudo_pi_used:
        info.append("PSEUDO-π model was applied!")
        
    with open("grid_info.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(info) + "\n")

    # 11) Log
    print("\n=== GRID SUMMARY (symmetry-aware simple) ===")
    print(f"Molecule          : {name}")
    print(f"Point group       : {pg}")
    print(f"Mirror planes     : {sorted(axes) if axes else 'none'}")
    print(f"Grid symmetry     : {gsym}")
    print(f"Calc grid         : N=({Nx},{Ny},{Nz})  step(Å)=({dx:.6f},{dy:.6f},{dz:.6f})  points={P.shape[0]}")
    print(f"Full grid         : points={Pf.shape[0]}")
    print(f"Files             : {calc_xyz}, {full_xyz}, grid_info.txt, {oriented_xyz}")
    print(f"Inputs            : {n_inp} .inp files")
    if pseudo_pi_used:
        print("PSEUDO-π model was applied!")

if __name__ == "__main__":
    main()

