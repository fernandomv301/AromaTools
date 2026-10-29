#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys
from pathlib import Path

from .amg_read_input import read_input
from .reset_calc import main as reset_calc_main
from .check_out import main as check_out_main
from .sum_rest import main as sum_rest_main

BANNER = r"""
======================================================================
                           AROMAGNETIC
              Magnetic Criteria for Aromaticity Evaluation
======================================================================

AroMagnetic implements widely used magnetic descriptors of aromaticity,
allowing detailed analysis of magnetic response properties derived from
the nuclear magnetic shielding tensor.

Main capabilities:

 • NICS-based analysis:
      NICS-scan-XY and Predi-XY
      NICS-2D
      NICS-3D
      NICS-SP
      INICS (Integral NICS)

 • Evaluation of ring current strength via numerical integration

 • Reconstruction of the induced magnetic field (Bind) from the
   nuclear magnetic shielding tensor under an external magnetic field (Bext)

 • Molecular orbital decomposition of the shielding tensor

 • Pseudo-π model implementation for π-electron magnetic response analysis

 • Generation of scalar fields (.cube) and visualization grids (.vtk)

These tools enable quantitative characterization of aromaticity,
antiaromaticity and magnetic response in molecular systems.

----------------------------------------------------------------------
Developed by:

   Fernando Martínez-Villarino
   Mesías Orozco-Ic
   Gabriel Merino

   Cinvestav - Merida
   Centro de Investigacion y de Estudios Avanzados del IPN
   Merida, Yucatan, Mexico

----------------------------------------------------------------------
Citation:

   Martinez-Villarino, F.; Orozco-Ic, M.; Merino, G.
   AromaTools 1.0, Cinvestav Merida, Merida, Mexico, 2025.

----------------------------------------------------------------------
Available workflows:

   amg run            -> generate inputs from INPUT.txt
   amg get            -> compute NICS / Bind / integrals
   amg mo             -> MO decomposition analysis
   amg rings          -> show ring numbering for PREDI_XY_PATH
   amg reset          -> remove generated files
   amg check          -> verify Gaussian/ORCA calculations
   amg add / amg sub  -> algebra with .cube or .vtk fields

======================================================================
"""

def print_banner():
    print(BANNER)


def _read_cfg():
    cfg = read_input("INPUT.txt")
    if not cfg:
        raise SystemExit("ERROR: INPUT.txt was not found or could not be read.")
    return cfg


def _detect_build_mode(cfg):
    has_grid = "GRID" in cfg and bool(cfg["GRID"])
    has_vp   = "VP"   in cfg and bool(cfg["VP"])
    has_hp   = "HP"   in cfg and bool(cfg["HP"])

    if not has_grid and not has_vp and not has_hp:
        raise SystemExit(
            "ERROR: INPUT.txt must contain at least one of:\n"
            "  BEGIN_GRID … END_GRID\n"
            "  BEGIN_VP   … END_VP\n"
            "  BEGIN_HP   … END_HP"
        )

    if has_grid and (has_vp or has_hp):
        raise SystemExit(
            "ERROR: BEGIN_GRID cannot be used together with BEGIN_VP or BEGIN_HP.\n"
            "Use separate INPUT.txt files for each calculation type."
        )

    if has_grid:
        return "GRID"

    return "PROFILE"   # covers VP, HP, or both


def keyword_needs_generic_basis(cfg):
    general = cfg.get("general", {})
    kline = str(general.get("KEYWORD_LINE", "") or "").strip()
    return bool(re.search(r"\bGEN(ECP)?\b", kline, flags=re.IGNORECASE))


def ensure_basis_if_needed(cfg):
    if not keyword_needs_generic_basis(cfg):
        return
    try:
        with open("basis.txt", "r", encoding="utf-8") as f:
            content = f.read().strip()
    except FileNotFoundError:
        raise SystemExit(
            "ERROR: KEYWORD_LINE contains GEN or GENECP "
            "but basis.txt was not found."
        )
    if not content:
        raise SystemExit("ERROR: basis.txt exists but is empty.")


def _run_check_out_main(args=None):
    args = [] if args is None else list(args)
    old_argv = sys.argv[:]
    try:
        sys.argv = ["check_outs"] + args
        check_out_main()
    finally:
        sys.argv = old_argv


def _run_sum_rest_main(op, file1, file2):
    mapped = {"add": "sum", "sub": "rest"}
    old_argv = sys.argv[:]
    try:
        sys.argv = ["sum_rest", mapped[op], file1, file2]
        sum_rest_main()
    finally:
        sys.argv = old_argv


# ──────────────────────────────────────────────
#  Command implementations
# ──────────────────────────────────────────────

def run_work():
    """amg run — generate QC inputs from INPUT.txt."""
    cfg    = _read_cfg()
    mode   = _detect_build_mode(cfg)
    ensure_basis_if_needed(cfg)

    if mode == "GRID":
        from .get_grid import main as _main
        _main()

    elif mode == "PROFILE":
        # get_profile handles VP, HP, or both simultaneously
        from .get_profile import main as _main
        _main()


def run_analysis(interactive=False):
    """amg get — parse QC outputs, generate CSVs, plots, and integrals."""
    cfg  = _read_cfg()
    mode = _detect_build_mode(cfg)

    if mode == "GRID":
        from .get_cube import main as _main
        _main()

    elif mode == "PROFILE":
        from .get_results import main as _main
        _main(interactive=interactive)


def run_rings():
    """amg rings — show ring numbering PNG for PREDI_XY_PATH."""
    from .predixy import rings_info
    rings_info()


def run_mo():
    """amg mo — MO decomposition analysis (Gaussian only)."""
    cfg     = _read_cfg()
    general = cfg.get("general", {})
    program = str(general.get("PROGRAM", "GAUSSIAN")).strip().upper()

    if program != "GAUSSIAN":
        raise SystemExit(
            "ERROR: 'amg mo' is currently implemented only for Gaussian.\n"
            "Set PROGRAM GAUSSIAN in INPUT.txt."
        )

    from .get_mo import main as _main
    _main()


# ──────────────────────────────────────────────
#  Dispatcher
# ──────────────────────────────────────────────

_USAGE = """\
Usage: amg <command> [options]

Commands:
  run            Generate QC inputs from INPUT.txt
  get [--interactive]  Parse outputs, generate CSVs / plots / integrals
  mo             MO decomposition analysis (Gaussian only)
  rings          Show ring numbering PNG (for PREDI_XY_PATH)
  reset          Remove files generated by a previous run
  check          Verify that Gaussian/ORCA calculations finished correctly
  add <f1> <f2>  Add two .cube or .vtk scalar fields
  sub <f1> <f2>  Subtract two .cube or .vtk scalar fields
"""


def dispatch(args):
    """
    args = sys.argv[1:]   (everything after the 'amg' executable name)
    """
    if not args:
        print_banner()
        return

    cmd  = args[0].lower()
    rest = args[1:]

    if cmd == "run":
        run_work()

    elif cmd == "get":
        run_analysis(interactive="--interactive" in rest)

    elif cmd == "rings":
        run_rings()

    elif cmd == "mo":
        run_mo()

    elif cmd == "reset":
        reset_calc_main()

    elif cmd == "check":
        _run_check_out_main(rest)

    elif cmd in {"add", "sub"}:
        if len(rest) != 2:
            raise SystemExit(
                f"ERROR: 'amg {cmd}' requires exactly 2 files.\n"
                f"Usage: amg {cmd} file1.vtk file2.vtk\n"
                f"   or: amg {cmd} file1.cube file2.cube"
            )
        _run_sum_rest_main(cmd, rest[0], rest[1])

    else:
        raise SystemExit(f"Unknown command: '{cmd}'\n\n{_USAGE}")


def main():
    dispatch(sys.argv[1:])


if __name__ == "__main__":
    main()
