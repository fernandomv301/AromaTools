#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
from ase import io
from ase.atom import Atom
from pymatgen.symmetry.analyzer import PointGroupAnalyzer

# ---------------- orientation / bbox ----------------

def orient(fname, margin=0.0, fmt=None):
    ats = io.read(fname, index=-1, format=fmt)
    m = ats.get_masses()
    R = ats.get_positions()
    cm = np.average(R, axis=0, weights=m)
    R0 = R - cm
    ats.set_positions(R0)
    moments, axes = ats.get_moments_of_inertia(vectors=True)
    E = axes.T
    R1 = R0.dot(E)
    ats.set_positions(R1)
    mi, ma = R1.min(axis=0), R1.max(axis=0)
    if isinstance(margin, (int, float)):
        mi -= margin; ma += margin
    elif hasattr(margin, '__len__') and len(margin) == 3:
        mi -= np.array(margin); ma += np.array(margin)
    else:
        raise ValueError("margin must be scalar or len-3 sequence")
    L = ma - mi
    bbox = {
        'xmin': float(mi[0]), 'xmax': float(ma[0]),
        'ymin': float(mi[1]), 'ymax': float(ma[1]),
        'zmin': float(mi[2]), 'zmax': float(ma[2]),
        'lengths': (float(L[0]), float(L[1]), float(L[2]))
    }
    base = os.path.basename(fname)
    name, _ = os.path.splitext(base)
    return ats, name, moments, axes, bbox

# ---------------- symmetry utils ----------------

def is_reflection(R, atol=1e-6):
    return (np.isclose(np.linalg.det(R), -1.0, atol=atol) and
            np.allclose(R, R.T, atol=atol) and
            np.allclose(R @ R, np.eye(3), atol=atol) and
            np.isclose(np.trace(R), 1.0, atol=atol))

def get_reflections(pga, atol=1e-6):
    Rs, normals = [], []
    for op in pga.get_symmetry_operations():
        R = op.rotation_matrix
        if is_reflection(R, atol=atol):
            S = 0.5*(np.eye(3) - R)
            w, V = np.linalg.eigh(S)
            n = V[:, int(np.argmax(w))]
            n = n / np.linalg.norm(n)
            Rs.append(R); normals.append(n)
    uniq_n, uniq_R = [], []
    for i, n in enumerate(normals):
        if all(abs(np.dot(n, m)) <= 1 - 1e-6 for m in uniq_n):
            uniq_n.append(n); uniq_R.append(Rs[i])
    return uniq_R, uniq_n

def axis_of(v):
    eX = np.array([1.,0.,0.]); eY = np.array([0.,1.,0.]); eZ = np.array([0.,0.,1.])
    dots = {'X': abs(np.dot(v, eX)), 'Y': abs(np.dot(v, eY)), 'Z': abs(np.dot(v, eZ))}
    return max(dots, key=dots.get)

def pick_half_axis(normals):
    if not normals:
        return 'X'
    bases = {'X': np.array([1,0,0.]), 'Y': np.array([0,1,0.]), 'Z': np.array([0,0,1.])}
    best, score = 'X', -1.0
    for A, u in bases.items():
        s = max(abs(np.dot(n, u)) for n in normals)
        if s > score:
            best, score = A, s
    return best

def pick_quarter_axis(normals):
    if len(normals) < 2:
        return 'X', None
    bases = {'X': np.array([1,0,0.]), 'Y': np.array([0,1,0.]), 'Z': np.array([0,0,1.])}
    bestA, bestS, best_pair = 'X', 0.0, (normals[0], normals[1])
    N = len(normals)
    for i in range(N):
        for j in range(i+1, N):
            a = np.cross(normals[i], normals[j])
            na = np.linalg.norm(a)
            if na < 1e-6: continue
            a /= na
            for A, u in bases.items():
                s = abs(np.dot(a, u))
                if s > bestS:
                    bestA, bestS, best_pair = A, s, (normals[i], normals[j])
    return bestA, best_pair

def axis_aligned(R, tol=1e-6):
    return (
        np.allclose(R.dot(R.T), np.eye(3), atol=tol) and
        (np.count_nonzero(np.abs(R) > tol) <= 3)
    )

# ---------------- grid building ----------------

def cut_ranges(bbox, grid_sym, half_axis=None, quarter_axis=None, cube_min=6.0, isotropic=True):
    def in_zero(lo, hi): return (min(lo, 0.0), max(hi, 0.0))
    def pos_half(lo, hi): lo2, hi2 = in_zero(lo, hi); return (0.0, hi2)
    def span(lo, hi): return hi - lo
    def center(lo, hi): return 0.5*(lo + hi)
    def inflate(lo, hi, L):
        c = center(lo, hi); h = 0.5*float(L); return (c - h, c + h)
    xmin, xmax = in_zero(bbox['xmin'], bbox['xmax'])
    ymin, ymax = in_zero(bbox['ymin'], bbox['ymax'])
    zmin, zmax = in_zero(bbox['zmin'], bbox['zmax'])
    if isotropic:
        sx, sy, sz = span(xmin, xmax), span(ymin, ymax), span(zmin, zmax)
        target = max(cube_min, sx, sy, sz)
        xmin, xmax = inflate(xmin, xmax, target)
        ymin, ymax = inflate(ymin, ymax, target)
        zmin, zmax = inflate(zmin, zmax, target)
        xmin, xmax = in_zero(xmin, xmax)
        ymin, ymax = in_zero(ymin, ymax)
        zmin, zmax = in_zero(zmin, zmax)
    mode = (grid_sym or '').upper()
    if mode.startswith('FULL'):
        return (xmin, xmax), (ymin, ymax), (zmin, zmax)
    if mode.startswith('HALF') and half_axis:
        ax = half_axis.upper()
        if ax == 'X': return pos_half(xmin, xmax), (ymin, ymax), (zmin, zmax)
        if ax == 'Y': return (xmin, xmax), pos_half(ymin, ymax), (zmin, zmax)
        if ax == 'Z': return (xmin, xmax), (ymin, ymax), pos_half(zmin, zmax)
    if mode.startswith('QUARTER') and quarter_axis:
        keep = quarter_axis.upper()
        if keep == 'X': return (xmin, xmax), pos_half(ymin, ymax), pos_half(zmin, zmax)
        if keep == 'Y': return pos_half(xmin, xmax), (ymin, ymax), pos_half(zmin, zmax)
        if keep == 'Z': return pos_half(xmin, xmax), pos_half(ymin, ymax), (zmin, zmax)
    if mode.startswith('OCT'):
        return pos_half(xmin, xmax), pos_half(ymin, ymax), pos_half(zmin, zmax)
    return (xmin, xmax), (ymin, ymax), (zmin, zmax)

def linsp(bmin, bmax, N):
    if N <= 1: return np.array([0.5*(bmin+bmax)])
    return np.linspace(bmin, bmax, N)

def gen_grid(xr, yr, zr, Nx, Ny, Nz):
    X = linsp(xr[0], xr[1], Nx)
    Y = linsp(yr[0], yr[1], Ny)
    Z = linsp(zr[0], zr[1], Nz)
    dx = (X[-1]-X[0])/(Nx-1) if Nx>1 else 0.0
    dy = (Y[-1]-Y[0])/(Ny-1) if Ny>1 else 0.0
    dz = (Z[-1]-Z[0])/(Nz-1) if Nz>1 else 0.0
    pts = [(x,y,z) for x in X for y in Y for z in Z]
    return np.array(pts, float), dx, dy, dz

# ---------------- replication ----------------

def _reflect(P, fx=False, fy=False, fz=False, tol=1e-9):
    out = []
    for x,y,z in P:
        if fx: x = -x
        if fy: y = -y
        if fz: z = -z
        out.append((x,y,z))
    key = lambda q: (round(q[0]/tol)*tol, round(q[1]/tol)*tol, round(q[2]/tol)*tol)
    uniq = {}
    for q in out:
        uniq[key(q)] = q
    return np.array(list(uniq.values()), float)

def replicate(P, grid_sym, axis):
    pts = [P]
    m = (grid_sym or '').upper()
    if m.startswith('FULL'):
        full = P
    elif m.startswith('HALF'):
        if axis == 'X': pts.append(_reflect(P, fx=True))
        elif axis == 'Y': pts.append(_reflect(P, fy=True))
        else: pts.append(_reflect(P, fz=True))
        full = np.vstack(pts)
    elif m.startswith('QUARTER'):
        if axis == 'X':
            pts += [_reflect(P, fy=True), _reflect(P, fz=True), _reflect(P, fy=True, fz=True)]
        elif axis == 'Y':
            pts += [_reflect(P, fx=True), _reflect(P, fz=True), _reflect(P, fx=True, fz=True)]
        else:
            pts += [_reflect(P, fx=True), _reflect(P, fy=True), _reflect(P, fx=True, fy=True)]
        full = np.vstack(pts)
    elif m.startswith('OCT'):
        pts += [
            _reflect(P, fx=True), _reflect(P, fy=True), _reflect(P, fz=True),
            _reflect(P, fx=True, fy=True), _reflect(P, fx=True, fz=True),
            _reflect(P, fy=True, fz=True), _reflect(P, fx=True, fy=True, fz=True)
        ]
        full = np.vstack(pts)
    else:
        full = P
    idx = np.lexsort((full[:,2], full[:,1], full[:,0]))
    return full[idx]

# ---------------- xyz writer ----------------

def write_xyz(path, base_atoms, pts):
    symbols = base_atoms.get_chemical_symbols()
    coords = base_atoms.get_positions()
    n_atoms = len(symbols)
    n_pts = len(pts)
    total = n_atoms + n_pts

    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{total}\n")
        f.write("Grid generated by AroMagnetic\n")
        for sym, (x, y, z) in zip(symbols, coords):
            f.write(f"{sym:<2s} {x:15.8f} {y:15.8f} {z:15.8f}\n")
        for x, y, z in pts:
            f.write(f"X  {x:15.8f} {y:15.8f} {z:15.8f}\n")


# ----------------------------------------------------------------
#  Shared small helpers used across multiple modules
# ----------------------------------------------------------------

def _truthy(value):
    """Return True when *value* represents a yes/on/true flag."""
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().upper() in {"YES", "Y", "TRUE", "T", "ON", "1"}


def _parse_heights(height_raw):
    """Normalise a HEIGHT value from INPUT.txt to list[float]."""
    if isinstance(height_raw, (list, tuple)):
        return [float(h) for h in height_raw]
    return [float(height_raw)]

