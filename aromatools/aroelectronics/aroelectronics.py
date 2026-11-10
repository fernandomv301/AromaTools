#!/usr/bin/env python3

import os
import sys
import math
import warnings
import numpy as np
import networkx as nx

warnings.filterwarnings("ignore")

from ase.data import covalent_radii, atomic_numbers

_BOHR_TO_ANG = 0.529177210903
_COV_TOL     = 1.3

SUPPORTED_PARTITIONS = {"meta-lowdin", "nao", "iao", "mulliken", "lowdin", "qtaim"}

BANNER = r"""
======================================================================
                         AROELECTRONICS
          Electronic Criteria for Aromaticity Evaluation
======================================================================

 Computes electron-delocalization aromaticity indices via ESIpy:

   Iring     Multicenter bond index (Giambiagi)
   MCI       Multicenter index (Bultinck)
   PDI       Para-delocalization index (6-membered rings only)
   FLU       Aromatic fluctuation index
   AV1245    4-center index (>= 6-membered rings)
   BOA       Bond-order alternation

 Three calculation modes:
   • PySCF (from scratch): provide xyz + functional + basis set
   • Wavefunction file  : provide a .fchk or .molden file
   • QTAIM (AIMAll)     : PARTITION qtaim + folder with the atomic .int files

----------------------------------------------------------------------
 Developed by:
   Fernando Martínez-Villarino · Joan Grebol-Tomas · Pedro Salvador ·
   Mesías Orozco-Ic · Eduard Matito · Gabriel Merino

 Citation:
   Martínez-Villarino, F.; Grebol-Tomas, J.; Salvador, P.; Orozco-Ic, M.; Matito, E.; Merino, G.
   AromaTools 1.0, Cinvestav Mérida, Mérida, Mexico, 2025.

 Cinvestav – Mérida, Universitat de Girona and Donostia International Physics Center, 2026
======================================================================
"""


def _read_input(path):
    data = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            i = line.find("#")
            if i != -1:
                line = line[:i]
            line = line.strip()
            if not line:
                continue
            parts = line.split(None, 1)
            data[parts[0].upper()] = parts[1].strip() if len(parts) > 1 else ""
    return data


def _truthy(val):
    return str(val).strip().upper() in ("YES", "Y", "TRUE", "1", "ON")


def _read_xyz(path):
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    n = int(lines[0].strip())
    syms, coords = [], []
    for line in lines[2:2 + n]:
        p = line.split()
        syms.append(p[0])
        coords.append([float(p[1]), float(p[2]), float(p[3])])
    return syms, coords


def _detect_rings(symbols, coords):
    G = nx.Graph()
    n = len(symbols)
    G.add_nodes_from(range(n))
    for i in range(n):
        zi = atomic_numbers.get(symbols[i], 0)
        ri = covalent_radii[zi] if zi else 0.0
        for j in range(i + 1, n):
            zj = atomic_numbers.get(symbols[j], 0)
            rj = covalent_radii[zj] if zj else 0.0
            d = math.sqrt(sum((coords[i][k] - coords[j][k]) ** 2 for k in range(3)))
            if d <= (ri + rj) * _COV_TOL:
                G.add_edge(i, j)
    cycles = nx.minimum_cycle_basis(G)
    return [[a + 1 for a in c] for c in cycles]  # keep the cyclic order of the ring


def _run_pyscf(symbols, coords, functional, basis, charge, mult):
    from pyscf import gto, dft
    atom_str = "\n".join(f"{s}  {c[0]}  {c[1]}  {c[2]}" for s, c in zip(symbols, coords))
    mol = gto.Mole()
    mol.atom     = atom_str
    mol.basis    = basis.lower()
    mol.charge   = charge
    mol.spin     = mult - 1
    mol.verbose  = 3
    mol.build()
    is_open = mult != 1
    mf = dft.UKS(mol) if is_open else dft.RKS(mol)
    mf.xc = functional.lower()
    mf.kernel()
    return mol, mf


def _parse_wavefunction(wf_file):
    from esipy import parse
    return parse(wf_file)


def _run_esi(mol, mf, rings, partition):
    from esipy import ESI
    return ESI(mol=mol, mf=mf, rings=rings, partition=partition)


def _read_qtaim(path):
    # The AOMs of every atom are read from the AIMAll .int files with qtaim_parser.
    try:
        from .qtaim_parser import read_qtaim
    except ImportError:
        from qtaim_parser import read_qtaim
    from esipy.atomicfiles import read_molinfo

    data = read_qtaim(path)
    molinfo = read_molinfo(path)  # symbols and geometry (from the .wfx) used by FLU and the ring search
    if data["is_open_shell"]:
        aom = data["aoms_by_spin"]       # [alpha, beta]
    else:
        aom = data["aoms_by_spin"][0]
    return aom, molinfo


def _run_esi_qtaim(aom, rings, molinfo):
    from esipy import ESI
    return ESI(aom=aom, rings=rings, partition="qtaim", molinfo=molinfo)


def _fmt(val, decimals=6):
    if val is None:
        return "N/A"
    try:
        return f"{val:.{decimals}f}"
    except (TypeError, ValueError):
        return str(val)


def _print_results(esi, rings, want_iring, want_mci, want_pdi, want_flu, want_av1245, want_boa=False):
    for i, (ring, ind) in enumerate(zip(rings, esi.indicators)):
        ring_str = "-".join(str(a) for a in ring)
        print(f"\n  Ring {i + 1}: [{ring_str}]")
        if want_iring:
            print(f"    Iring    = {_fmt(ind.iring)}   (^1/n = {_fmt(ind.iring_n)})")
        if want_mci:
            if ind.mci is not None:
                print(f"    MCI      = {_fmt(ind.mci)}   (^1/n = {_fmt(ind.mci_n)})")
            else:
                print(f"    MCI      = N/A  (ring > 12 atoms)")
        if want_pdi:
            print(f"    PDI      = {_fmt(ind.pdi) if ind.pdi is not None else 'N/A (not 6-membered)'}")
        if want_flu:
            print(f"    FLU      = {_fmt(ind.flu)}")
        if want_av1245:
            av  = getattr(ind, "av1245", None)
            avm = getattr(ind, "avmin",  None)
            print(f"    AV1245   = {_fmt(av)}   AVmin = {_fmt(avm)}")
        if want_boa:
            print(f"    BOA      = {_fmt(getattr(ind, 'boa', None))}   BOAc = {_fmt(getattr(ind, 'boa_c', None))}")


def _write_csv(fname, header, rows):
    with open(fname, "w", encoding="utf-8") as f:
        f.write(header + "\n")
        for row in rows:
            f.write(",".join(str(x) for x in row) + "\n")


def _write_txt(name, partition, data, wf_file, xyz_file, rings, esi,
               want_iring, want_mci, want_pdi, want_flu, want_av1245, csv_files, want_boa=False):
    inds = esi.indicators
    lines = []
    lines.append("=" * 60)
    lines.append(f"  AroElectronics — {name}")
    lines.append("=" * 60)
    lines.append("")
    lines.append("Calculation setup")
    lines.append("-" * 30)
    if partition == "qtaim":
        lines.append(f"  Mode              : QTAIM (AIMAll atomic files)")
        lines.append(f"  Atomic files      : {data.get('ATOMICFILES_PATH', '')}")
    elif wf_file:
        lines.append(f"  Mode              : wavefunction file")
        lines.append(f"  Wavefunction file : {wf_file}")
    else:
        lines.append(f"  Mode              : PySCF (from scratch)")
        lines.append(f"  Coordinates       : {xyz_file}")
        lines.append(f"  Functional        : {data.get('FUNCTIONAL', '')}")
        lines.append(f"  Basis set         : {data.get('BASIS_SET', '')}")
        lines.append(f"  Charge            : {data.get('CHARGE', 0)}")
        lines.append(f"  Multiplicity      : {data.get('MULTIPLICITY', 1)}")
    lines.append(f"  Partition         : {partition}")
    lines.append("")
    lines.append("Rings")
    lines.append("-" * 30)
    for i, ring in enumerate(rings):
        lines.append(f"  Ring {i + 1}: {'-'.join(str(a) for a in ring)}")
    lines.append("")
    lines.append("Aromaticity indices")
    lines.append("-" * 30)
    for i, (ring, ind) in enumerate(zip(rings, inds)):
        ring_str = "-".join(str(a) for a in ring)
        lines.append(f"\n  Ring {i + 1}: [{ring_str}]")
        if want_iring:
            lines.append(f"    Iring    = {_fmt(ind.iring)}   (^1/n = {_fmt(ind.iring_n)})")
        if want_mci:
            if ind.mci is not None:
                lines.append(f"    MCI      = {_fmt(ind.mci)}   (^1/n = {_fmt(ind.mci_n)})")
            else:
                lines.append(f"    MCI      = N/A  (ring > 12 atoms)")
        if want_pdi:
            lines.append(f"    PDI      = {_fmt(ind.pdi) if ind.pdi is not None else 'N/A (not 6-membered)'}")
        if want_flu:
            lines.append(f"    FLU      = {_fmt(ind.flu)}")
        if want_av1245:
            av  = getattr(ind, "av1245", None)
            avm = getattr(ind, "avmin",  None)
            lines.append(f"    AV1245   = {_fmt(av)}   AVmin = {_fmt(avm)}")
        if want_boa:
            lines.append(f"    BOA      = {_fmt(getattr(ind, 'boa', None))}   BOAc = {_fmt(getattr(ind, 'boa_c', None))}")
    lines.append("")
    lines.append("Output files")
    lines.append("-" * 30)
    for f in csv_files:
        lines.append(f"  {f}")
    lines.append("")
    lines.append("=" * 60)

    fname = f"{name}_{partition}_results.txt"
    with open(fname, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return fname


def _save_results(name, partition, rings, esi, want_iring, want_mci, want_pdi, want_flu, want_av1245, want_boa=False):
    inds  = esi.indicators
    files = []

    def ring_label(r):
        return "-".join(str(a) for a in r)

    if want_iring:
        rows  = [(i + 1, ring_label(r), _fmt(ind.iring), _fmt(ind.iring_n))
                 for i, (r, ind) in enumerate(zip(rings, inds))]
        fname = f"{name}_{partition}_iring.csv"
        _write_csv(fname, "Ring,Atoms,Iring,Iring**(1/n)", rows)
        files.append(fname)

    if want_mci:
        rows  = [(i + 1, ring_label(r), _fmt(ind.mci), _fmt(ind.mci_n))
                 for i, (r, ind) in enumerate(zip(rings, inds)) if ind.mci is not None]
        fname = f"{name}_{partition}_mci.csv"
        _write_csv(fname, "Ring,Atoms,MCI,MCI**(1/n)", rows)
        files.append(fname)

    if want_pdi:
        rows  = [(i + 1, _fmt(ind.pdi) if ind.pdi is not None else "N/A")
                 for i, (_, ind) in enumerate(zip(rings, inds))]
        fname = f"{name}_{partition}_pdi.csv"
        _write_csv(fname, "Ring,PDI", rows)
        files.append(fname)

    if want_flu:
        rows  = [(i + 1, ring_label(r), _fmt(ind.flu))
                 for i, (r, ind) in enumerate(zip(rings, inds))]
        fname = f"{name}_{partition}_flu.csv"
        _write_csv(fname, "Ring,Atoms,FLU", rows)
        files.append(fname)

    if want_av1245:
        rows  = [(i + 1, _fmt(getattr(ind, "av1245", None)), _fmt(getattr(ind, "avmin", None)))
                 for i, (_, ind) in enumerate(zip(rings, inds))]
        fname = f"{name}_{partition}_av1245.csv"
        _write_csv(fname, "Ring,AV1245,AVmin", rows)
        files.append(fname)

    if want_boa:
        rows  = [(i + 1, ring_label(r), _fmt(getattr(ind, "boa", None)), _fmt(getattr(ind, "boa_c", None)))
                 for i, (r, ind) in enumerate(zip(rings, inds))]
        fname = f"{name}_{partition}_boa.csv"
        _write_csv(fname, "Ring,Atoms,BOA,BOAc", rows)
        files.append(fname)

    return files


def main():
    print(BANNER)

    input_path = sys.argv[1] if len(sys.argv) > 1 else "INPUT.txt"
    if not os.path.isfile(input_path):
        raise SystemExit(f"ERROR: {input_path} not found.")

    data = _read_input(input_path)

    name      = data.get("NAME_MOLECULE", "molecule")
    partition = data.get("PARTITION", "meta-lowdin").strip().lower()
    wf_file   = data.get("WAVEFUNCTION_FILE", "").strip()
    xyz_file  = data.get("MOLECULAR_COORDINATES", "").strip()
    qtaim_path = data.get("ATOMICFILES_PATH", "").strip()

    if partition not in SUPPORTED_PARTITIONS:
        raise SystemExit(
            f"ERROR: PARTITION '{partition}' not supported.\n"
            f"Options: {', '.join(sorted(SUPPORTED_PARTITIONS))}"
        )

    want_iring  = _truthy(data.get("IRING",  "NO"))
    want_mci    = _truthy(data.get("MCI",    "NO"))
    want_pdi    = _truthy(data.get("PDI",    "NO"))
    want_flu    = _truthy(data.get("FLU",    "NO"))
    want_av1245 = _truthy(data.get("AV1245", "NO"))
    want_boa    = _truthy(data.get("BOA",    "NO"))

    if not any([want_iring, want_mci, want_pdi, want_flu, want_av1245, want_boa]):
        raise SystemExit(
            "ERROR: No index requested. Set at least one of IRING, MCI, PDI, FLU, AV1245, BOA to YES."
        )

    if partition == "qtaim":
        if not os.path.isdir(qtaim_path):
            raise SystemExit("ERROR: PARTITION qtaim needs ATOMICFILES_PATH, the folder with the AIMAll .int files.")
        print(f"[info] Reading QTAIM atomic files: {qtaim_path}")
        aom, molinfo = _read_qtaim(qtaim_path)
        symbols = list(molinfo["symbols"])
        coords  = None if molinfo.get("geom") is None else (molinfo["geom"] * _BOHR_TO_ANG).tolist()
    elif wf_file:
        if not os.path.isfile(wf_file):
            raise SystemExit(f"ERROR: WAVEFUNCTION_FILE not found: {wf_file}")
        print(f"[info] Reading wavefunction: {wf_file}")
        mol, mf = _parse_wavefunction(wf_file)
        symbols = list(mol.atomic_symbols)
        coords  = (mol.atom_coords() * _BOHR_TO_ANG).tolist()
    else:
        for key in ("MOLECULAR_COORDINATES", "FUNCTIONAL", "BASIS_SET"):
            if not data.get(key, "").strip():
                raise SystemExit(f"ERROR: {key} is required when WAVEFUNCTION_FILE is not provided.")
        if not os.path.isfile(xyz_file):
            raise SystemExit(f"ERROR: MOLECULAR_COORDINATES not found: {xyz_file}")
        charge = int(data.get("CHARGE", 0))
        mult   = int(data.get("MULTIPLICITY", 1))
        print(f"[info] Running PySCF: {data['FUNCTIONAL']}/{data['BASIS_SET']}")
        symbols, coords = _read_xyz(xyz_file)
        mol, mf = _run_pyscf(symbols, coords, data["FUNCTIONAL"], data["BASIS_SET"], charge, mult)

    ring_raw = data.get("RING", "").strip()
    if ring_raw:
        rings = [[int(x) for x in ring_raw.split()]]
        print(f"[info] Manual ring: {rings[0]}")
    else:
        if coords is None:
            raise SystemExit(
                "ERROR: The .wfx file was not found, so the rings cannot be detected.\n"
                "Add  RING  <atom1> <atom2> ...  to INPUT.txt to specify it manually."
            )
        rings = _detect_rings(symbols, coords)
        if not rings:
            raise SystemExit(
                "ERROR: No rings detected automatically.\n"
                "Add  RING  <atom1> <atom2> ...  to INPUT.txt to specify it manually."
            )
        print(f"[info] Auto-detected {len(rings)} ring(s)")

    print(f"[info] Partition: {partition}")
    print(f"[info] Indices  : "
          + "  ".join(k for k, v in [("Iring", want_iring), ("MCI", want_mci),
                                      ("PDI", want_pdi), ("FLU", want_flu),
                                      ("AV1245", want_av1245), ("BOA", want_boa)] if v))

    if partition == "qtaim":
        esi = _run_esi_qtaim(aom, rings, molinfo)
    else:
        esi = _run_esi(mol, mf, rings, partition)

    print(f"\n{'='*50}")
    print(f"  {name}  |  partition: {partition}")
    print("="*50)
    _print_results(esi, rings, want_iring, want_mci, want_pdi, want_flu, want_av1245, want_boa)

    files = _save_results(name, partition, rings, esi,
                          want_iring, want_mci, want_pdi, want_flu, want_av1245, want_boa)

    txt = _write_txt(name, partition, data, wf_file, xyz_file, rings, esi,
                     want_iring, want_mci, want_pdi, want_flu, want_av1245, files, want_boa=want_boa)
    files.append(txt)

    print(f"\n[ok] Output files:")
    for f in files:
        print(f"  {f}")


if __name__ == "__main__":
    main()
