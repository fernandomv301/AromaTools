# AromaTools

AromaTools is a Python package designed for the **analysis of aromaticity** using three main approaches.

---

## Aromaticity Criteria

### Magnetic Criteria (AroMagnetic)

- NICS-based analysis:
  - NICS-SP
  - NICS-2D
  - NICS-3D
  - NICS-scan
  - INICS (Integral NICS)

- Evaluation of ring current strength via numerical integration

- Reconstruction of the induced magnetic field (Bind) from the 
  nuclear magnetic shielding tensor under an external magnetic field (Bext)

- Molecular orbital decomposition of the shielding tensor

- Pseudo-pi model implementation for pi-electron magnetic response analysis

- Generation of scalar fields (.cube) and visualization grids (.vtk)

These tools enable quantitative and qualitative characterization of aromaticity, antiaromaticity and magnetic response in molecular systems.

### Geometric Criteria (AroGeometric)

- Bond Length Alternation (BLA)
- HOMA index
- HOMER index (excited state aromaticity)
- HOMAc index

Interpretation:

- HOMA = 1 → aromatic
- HOMA = 0 → non-aromatic
- HOMAc = −1 → antiaromatic

---

### Electronic Criteria (AroElectronics)

- Iring (multicenter bond index)
- MCI (multicenter index)
- PDI (para-delocalization index)
- FLU (aromatic fluctuation index)
- AV1245 and AVmin
- BOA and BOAc (bond-order alternation)
- DI (delocalization indices of the ring bonds) and atomic populations

The indices are computed with ESIpy from a PySCF calculation, from a
wavefunction file (.fchk, .molden) or from AIMAll atomic files (QTAIM), using the
meta-Löwdin, NAO, IAO, Löwdin or Mulliken atomic partitions (or QTAIM):

```bash
ael run
```

(`ael` and `aroelectronics` are the same command; `ael run` reads `INPUT.txt` from the current folder.)

---

## Installation

AromaTools requires Python 3.10 or newer.

Install using pip:

```bash
pip install aromatools
```

Or install it from the Git repository:

```bash
pip install git+https://github.com/fernandomv301/AromaTools.git
```

You can also clone the repository and install it from the local folder:

```bash
git clone https://github.com/fernandomv301/AromaTools.git
cd AromaTools
pip install .
```

See the manual (Installation page) for updating, editable installs and the list of requirements.
