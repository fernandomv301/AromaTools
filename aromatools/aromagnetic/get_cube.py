#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Python modules
import os
import re
import glob
import numpy as np
import subprocess
import sys
from ase.data import atomic_numbers

# AromaTools modules
from .amg_read_input import read_input
from .check_out import check_file as _check_file

script_dir = os.path.dirname(os.path.abspath(__file__))

def _truthy(val):
    return str(val).strip().upper() in ("YES", "Y", "TRUE", "T", "1", "ON")

ANG2BOHR = 1.8897259886 # Value getting form NIST!!

def znum(sym: str):
    s = sym.strip()
    m = re.match(r"[A-Za-z]+", s)
    if not m:
        raise ValueError("Unknow elemnt symbol for ASE")
    elem = m.group(0).title()
    try:
        return atomic_numbers[elem]
    except KeyError:
        raise ValueError(f"Unknow elemnt symbol for ASE: {elem!r} (from {sym!r})")

def load_atoms_bohr(xyz_path):
    if not os.path.exists(xyz_path):
        raise FileNotFoundError(f"Missing XYZ: {xyz_path}")
    rows = []
    with open(xyz_path, "r", encoding="utf-8", errors="ignore") as f:
        try:
            n = int(f.readline().split()[0])
        except Exception:
            raise ValueError(f"Bad first line in {xyz_path!r}")
        _ = f.readline()  # comentario
        for _ in range(n):
            ln = f.readline()
            if not ln:
                break
            parts = ln.split()
            if len(parts) < 4:
                continue
            symU = parts[0].upper()
            if symU in {"X", "BQ"}:
                continue
            x, y, z = map(float, parts[1:4])
            core_match = re.match(r"[A-Za-z]+", parts[0])
            core = core_match.group(0) if core_match else parts[0]
            Z = znum(core)
            rows.append((Z, 0.0, x*ANG2BOHR, y*ANG2BOHR, z*ANG2BOHR))
    if not rows:
        raise ValueError(f"No real atoms found in {xyz_path!r}")
    return rows

def read_xpts(xyz_path):
    if not os.path.exists(xyz_path):
        raise FileNotFoundError(f"Missing {xyz_path}")
    pts = []
    first = None
    with open(xyz_path, "r", encoding="utf-8", errors="ignore") as f:
        try:
            n = int(f.readline().split()[0])
        except Exception:
            raise ValueError(f"Bad header in {xyz_path}")
        _ = f.readline()
        for _ in range(n):
            ln = f.readline()
            if not ln:
                break
            parts = ln.split()
            if len(parts) < 4:
                continue
            if parts[0].upper() == 'X':
                x, y, z = map(float, parts[1:4])
                if first is None:
                    first = (x, y, z)
                pts.append((x, y, z))
    if first is None or not pts:
        raise ValueError(f"No 'X' points in {xyz_path}")
    return first, pts

_floatlabs = ("XX","XY","XZ","YX","YY","YZ","ZX","ZY","ZZ") #Gaussian shielding tensors!

def _labvals(line):
    out = {}
    for lab in _floatlabs:
        m = re.search(rf"{lab}\s*=\s*([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)", line)
        if m:
            out[lab] = float(m.group(1))
    return out

def find_outs(prefix):
    files = glob.glob(f"{prefix}_*.out") + glob.glob(f"{prefix}_*.log")
    def idx(p):
        m = re.search(r"_(\d{4})\.(?:out|log)$", os.path.basename(p), re.I)
        return int(m.group(1)) if m else -1
    files = [p for p in files if idx(p) >= 0]
    files.sort(key=idx)
    return files

# ================ GAUSSIAN 16 Parser ================

def read_tens_g16(prefix, expected=None): 
    outs = find_outs(prefix)
    if not outs:
        raise FileNotFoundError(f"No Gaussian outputs found for '{prefix}_*.out/.log'")
    head = re.compile(r"SCF GIAO Magnetic shielding tensor \(ppm\):", re.I)
    bq = re.compile(r"^\s*(?:\d+|[A-Za-z]{1,2})?\s*Bq\b.*Isotropic", re.I)
    T = []
    for path in outs:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            lines = fh.readlines()
        pos = [i for i, s in enumerate(lines) if head.search(s)]
        if not pos:
            continue
        i = pos[0] + 1
        n = len(lines)
        while i < n:
            if bq.search(lines[i]):
                if i + 3 >= n:
                    break
                r1, r2, r3 = lines[i+1], lines[i+2], lines[i+3]
                v1, v2, v3 = _labvals(r1), _labvals(r2), _labvals(r3)
                try:
                    XX, YX, ZX = v1["XX"], v1["YX"], v1["ZX"]
                    XY, YY, ZY = v2["XY"], v2["YY"], v2["ZY"]
                    XZ, YZ, ZZ = v3["XZ"], v3["YZ"], v3["ZZ"]
                except KeyError:
                    nums = []
                    for r in (r1, r2, r3):
                        nums += [float(x) for x in re.findall(
                            r"[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?", r)]
                    if len(nums) < 9:
                        i += 1
                        continue
                    XX, YX, ZX, XY, YY, ZY, XZ, YZ, ZZ = nums[:9]
                sigma = np.array([
                    [XX, XY, XZ],
                    [YX, YY, YZ],
                    [ZX, ZY, ZZ]
                ], float)
                T.append(sigma)
                i += 4
                continue
            i += 1
    if expected is not None and len(T) != expected:
        raise RuntimeError(f"Found {len(T)} Bq blocks, expected {expected}")
    if not T:
        raise RuntimeError("No Bq blocks found in outputs.")
    return T

# ================ ORCA 6.1 Parser ================

def real_atoms(xyz_path):
    if not os.path.exists(xyz_path):
        raise FileNotFoundError(f"Missing XYZ: {xyz_path}")
    with open(xyz_path,"r", encoding="utf-8", errors="ignore") as f:
        line = f.readline()
    try:
        return int(line.split()[0])
    except Exception:
        raise ValueError("There's no number of atoms in the XZY file!!")
    
def read_tens_orca(prefix, expected=None, n_real_atoms=0):
    outs = find_outs(prefix)
    if not outs:
        raise FileNotFoundError("No ORCA outputs found :c")
    T_all = []
    for path in outs:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            lines = fh.readlines()
        n = len(lines)
        in_section = False
        sigmas_this_file = []
        i = 0
        while i < n:
            line = lines[i]
            if not in_section:
                if "CHEMICAL SHIELDINGS (ppm)" in line:
                    in_section = True
                i += 1
                continue
            if line.strip().startswith("Nucleus"):
                j = i + 1
                while j < n and "Total shielding tensor (ppm)" not in lines[j]:
                    j += 1
                if j >= n or "Total shielding tensor (ppm)" not in lines[j]:
                    i += 1
                    continue
                if j + 3 >= n:
                    break
                try:
                    row1 = [float(x) for x in lines[j+1].split()[:3]]
                    row2 = [float(x) for x in lines[j+2].split()[:3]]
                    row3 = [float(x) for x in lines[j+3].split()[:3]]
                except ValueError:
                    i = j + 4
                    continue
                sigma = np.array([row1, row2, row3], float)
                sigmas_this_file.append(sigma)
                i = j + 4
                continue
            i += 1
        if len(sigmas_this_file) <= n_real_atoms:
            continue
        dummies = sigmas_this_file[n_real_atoms:]
        T_all.extend(dummies)
    if expected is not None and len(T_all) != expected:
        raise RuntimeError(
            f"ORCA: Found {len(T_all)} dummy shielding tensors, expected {expected}"
        )
    if not T_all:
        raise RuntimeError("ORCA: No shielding tensors found for dummy atoms.")
    return T_all

# ---------------- NICS - IMF ----------------

def nics(tensors):
    # NICS = -(1/3) * (σxx + σyy + σzz) = - σ_iso
    return [-(np.trace(sigma) / 3.0) for sigma in tensors]

# ---------------- B_ind = -σ·B_ext ----------------
def bind_from_T(tensors, B):
    Bvec = np.array(B, float)
    N = len(tensors)
    bx = np.empty(N, float)
    by = np.empty(N, float)
    bz = np.empty(N, float)

    for i, sigma in enumerate(tensors):
        Bind = - sigma @ Bvec
        bx[i], by[i], bz[i] = Bind
    bx[np.abs(bx) < 5e-13] = 0.0
    by[np.abs(by) < 5e-13] = 0.0
    bz[np.abs(bz) < 5e-13] = 0.0
    return bx.tolist(), by.tolist(), bz.tolist()

# ---------------- make CUBE ----------------
def write_cube(path, atoms_rows, origin_bohr, steps_bohr, shape, values, per_line=1):
    Nx, Ny, Nz = shape
    if len(values) != Nx*Ny*Nz:
        raise ValueError(f"values has {len(values)} elements, expected {Nx*Ny*Nz}")
    ox, oy, oz = origin_bohr
    dx, dy, dz = steps_bohr
    with open(path, "w", encoding="utf-8") as f:
        f.write("CUBE FILE from CINVESTAV MDA\n")
        f.write("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n")
        f.write(f"{len(atoms_rows):5d} {ox:12.6f} {oy:12.6f} {oz:12.6f}\n")
        f.write(f"{Nx:5d} {dx:12.6f} {0.0:12.6f} {0.0:12.6f}\n")
        f.write(f"{Ny:5d} {0.0:12.6f} {dy:12.6f} {0.0:12.6f}\n")
        f.write(f"{Nz:5d} {0.0:12.6f} {0.0:12.6f} {dz:12.6f}\n")
        for (Z, q, xb, yb, zb) in atoms_rows:
            f.write(f"{Z:5d} {q:12.6f} {xb:12.6f} {yb:12.6f} {zb:12.6f}\n")
        line = []
        for ix in range(Nx):
            for iy in range(Ny):
                base = ix*(Ny*Nz) + iy*Nz
                for iz in range(Nz):
                    line.append(f"{values[base+iz]:13.6e}")
                    if len(line) == per_line:
                        f.write(" ".join(line) + "\n")
                        line = []
        if line:
            f.write(" ".join(line) + "\n")

# ---------------- grid_info.txt ----------------
def read_grid_info(path="grid_info.txt"):
    info = {}
    if not os.path.exists(path):
        return info
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            s = raw.strip()
            if not s or s.startswith("===") or ":" not in s:
                continue
            k, v = s.split(":", 1)
            k = k.strip().lower()
            v = v.strip()
            if k == "calc_xyz":
                info["calc_xyz"] = v
            elif k == "full_xyz":
                info["full_xyz"] = v
            elif k == "mol_file":
                info["mol_file"] = v
    return info

# -------- mapear grid reducida a grid completa --------
def build_full_mapping(calc_pts, full_pts, tol=1e-6):
    calc = np.array(calc_pts, float)
    full = np.array(full_pts, float)

    Xc = np.unique(calc[:, 0]); Yc = np.unique(calc[:, 1]); Zc = np.unique(calc[:, 2])
    Xf = np.unique(full[:, 0]); Yf = np.unique(full[:, 1]); Zf = np.unique(full[:, 2])

    def is_mirrored(c, f):
        return (c.min() >= -tol) and (f.min() < -tol) and np.isclose(c.max(), f.max(), atol=tol)

    mx = is_mirrored(Xc, Xf)
    my = is_mirrored(Yc, Yf)
    mz = is_mirrored(Zc, Zf)

    def key(v):
        return round(v/tol)*tol

    def k3(x, y, z):
        if mx: x = abs(x)
        if my: y = abs(y)
        if mz: z = abs(z)
        return (key(x), key(y), key(z))

    idx_map = {}
    for i, (x, y, z) in enumerate(calc_pts):
        idx_map[k3(x, y, z)] = i

    idx_full = []
    for (x, y, z) in full_pts:
        k = k3(x, y, z)
        if k not in idx_map:
            raise RuntimeError(
                f"full-grid point ({x:.6f},{y:.6f},{z:.6f}) not found in calc grid (after symmetry)."
            )
        idx_full.append(idx_map[k])

    Nx_full, Ny_full, Nz_full = len(Xf), len(Yf), len(Zf)
    Ox, Oy, Oz = Xf[0], Yf[0], Zf[0]
    dx = Xf[1] - Xf[0] if Nx_full > 1 else 0.0
    dy = Yf[1] - Yf[0] if Ny_full > 1 else 0.0
    dz = Zf[1] - Zf[0] if Nz_full > 1 else 0.0

    return idx_full, (Nx_full, Ny_full, Nz_full), (Ox, Oy, Oz), (dx, dy, dz)

def pad_vis(shape, origin_bohr, steps_bohr, values_list, pad_to=3):
    Nx, Ny, Nz = map(int, shape)
    ox, oy, oz = map(float, origin_bohr)
    dx, dy, dz = map(float, steps_bohr)

    if (Nx != 1) and (Ny != 1) and (Nz != 1):
        return (Nx, Ny, Nz), (ox, oy, oz), (dx, dy, dz), values_list

    nonzero = [s for s in (dx, dy, dz) if abs(s) > 0.0]
    step_pad = min(nonzero) if nonzero else 1e-3
    arr = np.asarray(values_list, dtype=float).reshape((Nx, Ny, Nz), order="C")

    if Nx == 1:
        dx = step_pad
        ox = ox - dx
        arr = np.repeat(arr, pad_to, axis=0)
        Nx = pad_to

    if Ny == 1:
        dy = step_pad
        oy = oy - dy
        arr = np.repeat(arr, pad_to, axis=1)
        Ny = pad_to

    if Nz == 1:
        dz = step_pad
        oz = oz - dz
        arr = np.repeat(arr, pad_to, axis=2)
        Nz = pad_to

    values2 = arr.ravel(order="C").tolist()
    return (Nx, Ny, Nz), (ox, oy, oz), (dx, dy, dz), values2

# ---------------- main ----------------
def main():
    cfg = read_input("INPUT.txt")
    if not cfg:
        raise SystemExit("ERROR: INPUT.txt not found in current directory")

    general = cfg.get("general", {})
    grid_blk = cfg.get("GRID", {})
    name = general.get("NAME_MOLECULE", "molecule")
    prog = (general.get("PROGRAM") or general.get("program") or "Gaussian").lower()
    cube_to_vtk_flag = grid_blk.get("VTK_FILES", "NO")
    B = grid_blk.get("EXTERNAL_FIELD", (0.0, 0.0, 0.0))
    num_points = grid_blk.get("NUM_POINTS", None)
    if not num_points or not isinstance(num_points, (list, tuple)) or len(num_points) != 3:
        raise SystemExit(
            "ERROR: NUM_POINTS Nx Ny Nz must be defined in INPUT.txt within BEGIN_GRID ... END_GRID"
        )

    Nx, Ny, Nz = map(int, num_points)
    Nx = int(Nx); Ny = int(Ny); Nz = int(Nz)
    if Nx < 1 or Ny < 1 or Nz < 1:
        raise SystemExit('error: NUM_POINTS must be positive integers')
    gi = read_grid_info("grid_info.txt")
    calc_xyz = gi.get("calc_xyz") or f"{name}in_calcgrid.xyz"
    full_xyz = gi.get("full_xyz") or f"{name}in_symgrid.xyz"
    mol_file = f"{name}_oriented.xyz"
    n_real = real_atoms(mol_file)
    _, calc_pts = read_xpts(calc_xyz)
    expected = Nx * Ny * Nz
    if len(calc_pts) != expected:
        raise SystemExit(
            f"error: calc grid points ({len(calc_pts)}) do not match NUM_POINTS "
            f"({Nx}*{Ny}*{Nz}={expected}). Regenerate grid with python amg_new.py"
        )

    print("\n=== GET CUBE (symmetry-aware) ===")
    print(f"Molecule        : {name}")
    print(f"Program         : {prog}")
    print(f"Calc grid XYZ   : {calc_xyz}")
    print(f"Full grid XYZ   : {full_xyz}")
    print(f"Number of atoms : {n_real}")
    print(f"Shape (fundamental) : (Nx,Ny,Nz)=({Nx},{Ny},{Nz})  points={expected}")
    print(f"External field: B=({B[0]}, {B[1]}, {B[2]})")

    # check output files before parsing tensors
    out_files = find_outs(name)
    bad = [f for f in out_files if not _check_file(f).startswith("OK")]
    for f in bad:
        print(f"[WARN] {f}: calculation did not terminate normally")
    if bad:
        print(
            f"[WARN] {len(bad)} output file(s) terminated abnormally — "
            "tensor count may be wrong. Run 'amg check' for details."
        )

    # tensores
    try:
        if "gaussian" in prog:
            T = read_tens_g16(name, expected=expected)
        elif "orca" in prog:
            T = read_tens_orca(name, expected=expected, n_real_atoms=n_real)
        else:
            raise SystemExit("This software only has support with Gaussian 16 and Orca version 6.1")
    except RuntimeError as e:
        raise SystemExit(
            f"ERROR: {e}\n"
            "One or more calculations terminated abnormally (see [WARN] above).\n"
            "Fix the failing jobs and rerun 'amg get'."
        )
    
    print(f"[ok] Found {len(T)} shielding tensors!")

    Bx_red, By_red, Bz_red = bind_from_T(T, B)
    _, full_pts = read_xpts(full_xyz)
    idx_full, shape_full, origin_ang, steps_ang = build_full_mapping(calc_pts, full_pts)
    NICS_read = nics(T)
    NICS_FULL = [NICS_read[i] for i in idx_full]
    Nx_full, Ny_full, Nz_full = shape_full
    Ox_ang, Oy_ang, Oz_ang = origin_ang
    dx_ang, dy_ang, dz_ang = steps_ang

    print(f"Full grid shape : (Nx,Ny,Nz)=({Nx_full},{Ny_full},{Nz_full})  points={len(full_pts)}")

    Bx_full = [Bx_red[i] for i in idx_full]
    By_full = [By_red[i] for i in idx_full]
    Bz_full = [Bz_red[i] for i in idx_full]

    total_full = Nx_full * Ny_full * Nz_full
    if len(Bx_full) != total_full:
        raise RuntimeError(
            f"Inconsistent full grid: {total_full} expected values vs {len(Bx_full)} mapped values"
        )

    ox = Ox_ang * ANG2BOHR
    oy = Oy_ang * ANG2BOHR
    oz = Oz_ang * ANG2BOHR
    dx = dx_ang * ANG2BOHR
    dy = dy_ang * ANG2BOHR
    dz = dz_ang * ANG2BOHR
    shape0 = (Nx_full, Ny_full, Nz_full)
    origin0 = (ox, oy, oz)
    steps0 = (dx, dy, dz)
    shape_out, origin_out, steps_out, NICS_out = pad_vis(shape0, origin0, steps0, NICS_FULL, pad_to=3)
    _, _, _, Bx_out = pad_vis(shape0, origin0, steps0, Bx_full,  pad_to=3)
    _, _, _, By_out = pad_vis(shape0, origin0, steps0, By_full,  pad_to=3)
    _, _, _, Bz_out = pad_vis(shape0, origin0, steps0, Bz_full,  pad_to=3)
    Nx_out, Ny_out, Nz_out = shape_out
    ox_out, oy_out, oz_out = origin_out
    dx_out, dy_out, dz_out = steps_out
    atoms_rows = load_atoms_bohr(mol_file)
    per_line = 1
    outn = f"{name}NICS.cube"
    outx = f"{name}B_ind_x.cube"
    outy = f"{name}B_ind_y.cube"
    outz = f"{name}B_ind_z.cube"
    write_cube(outn, atoms_rows, (ox_out, oy_out, oz_out), (dx_out, dy_out, dz_out),
               (Nx_out, Ny_out, Nz_out), NICS_out, per_line=per_line)
    write_cube(outx, atoms_rows, (ox_out, oy_out, oz_out), (dx_out, dy_out, dz_out),
               (Nx_out, Ny_out, Nz_out), Bx_out, per_line=per_line)
    write_cube(outy, atoms_rows, (ox_out, oy_out, oz_out), (dx_out, dy_out, dz_out),
               (Nx_out, Ny_out, Nz_out), By_out, per_line=per_line)
    write_cube(outz, atoms_rows, (ox_out, oy_out, oz_out), (dx_out, dy_out, dz_out),
               (Nx_out, Ny_out, Nz_out), Bz_out, per_line=per_line)
    if _truthy(cube_to_vtk_flag):
        vtk_script = os.path.join(script_dir, "cube_to_vtk.py")
        if not os.path.isfile(vtk_script):
            raise SystemExit(
                "\n[ERROR] VTK_FILES = YES but cube_to_vtk.py was not found in this directory.\n"
            )
        try:
            cube_files = [
                f"{name}B_ind_x.cube",
                f"{name}B_ind_y.cube",
                f"{name}B_ind_z.cube",
                f"{name}NICS.cube",
            ]
            subprocess.check_call(
                [sys.executable, vtk_script] + cube_files
            )
            print("[ok] cube_to_vtk.py executed automatically.")
        except Exception as e:
            raise SystemExit(
                f"\n[ERROR] Automatic CUBE -> VTK conversion failed: {e}\n"
                "You can run it manually: python cube_to_vtk.py\n"
            )
    print(f"Done!!! {outx}, {outy}, {outz}, {outn}")

if __name__ == "__main__":
    main()

