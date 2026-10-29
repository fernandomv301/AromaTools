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

Coming soon...

---

## Installation

Install using pip:

```bash
pip install aromatools