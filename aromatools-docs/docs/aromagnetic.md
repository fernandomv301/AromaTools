# AroMagnetic

AroMagnetic is a Python toolkit for evaluating **magnetic criteria of aromaticity** from quantum chemistry calculations.

The program automates the generation of ghost atom grids, reconstruction of the induced magnetic field (Bind), calculation of NICS-based descriptors, ring-current strength evaluation, and molecular orbital decomposition of magnetic response properties.

AroMagnetic integrates workflows typically performed manually when using Gaussian or ORCA.

---

## Capabilities

AroMagnetic implements the following magnetic descriptors:

- $B^{ind}_{z}$ (induced magnetic field)
- NICS-SP
- NICS-scan
- NICS-2D
- NICS-3D
- INICS (integral NICS)
- Ring current strength (Ampere–Maxwell law)
- MO decomposition of magnetic response (only with Gaussian)
- pseudo-π model

### Automation features

- Automatic grid construction
- Automatic molecular orientation
- Automatic ring detection
- Symmetry-adapted grid reduction
- Automatic cube generation
- cube → vtk conversion
- Numerical integration of profiles
- Algebra with scalar fields (cube/vtk)

---

## General workflow

AroMagnetic operates using an input file:

```text
INPUT.txt
```

Typical workflow:  

- Prepare calcuation inputs
- Run Gaussian or Orca
- Analyze magnetic response

## Grid module
### Input format

```text
# ================================
#        AroMagnetic INPUT
# ================================

# --------------------------------
# General Information
# --------------------------------

NAME_MOLECULE             mol
PROGRAM                   Gaussian            # ORCA or GAUSSIAN
KEYWORD_LINE              Functional/Basis NMR INTEGRAL(SUPERFINEGRID) SCF=XQC

PROCS                     16 #If you use Orca this value doesn't apear in the input files
MEMORY_GB                 32

CHARGE                    0
MULTIPLICITY              1

MOLECULAR_COORDINATES     mol.xyz
PSEUDO_PI                 NO              # YES or NO


# --------------------------------
# Module Switchboard
# --------------------------------
# Turn ON only what you need.

GRID_MODULE                 ON
PROFILE_MODULE             OFF
NICS_SCAN_XY_MODULE        OFF


# --------------------------------
# Grid (Grid generation)
# --------------------------------
# Required if GRID_MODULE = ON

BEGIN_GRID

GA_MAX            80               # Maximum number of ghost atoms per input
EXTERNAL_FIELD    0.0 0.0 1.0      # Bx By Bz
NUM_POINTS        19 19 19         # Nx Ny Nz
GRID_MARGIN       2.0              # Å
VTK_FILES         YES              # YES or NO

END_GRID
```

## Profile module
### Input format

```text
# ================================
#        AroMagnetic INPUT
# ================================

# --------------------------------
# General Information
# --------------------------------

NAME_MOLECULE             mol
PROGRAM                   Gaussian            # ORCA or GAUSSIAN
KEYWORD_LINE              Functional/Basis NMR INTEGRAL(SUPERFINEGRID) SCF=XQC

PROCS                     16 #If you use Orca this value doesn't apear in the input files
MEMORY_GB                 32

CHARGE                    0
MULTIPLICITY              1

MOLECULAR_COORDINATES     mol.xyz
PSEUDO_PI                 NO              # YES or NO


# --------------------------------
# Module Switchboard
# --------------------------------
# Turn ON only what you need.

GRID_MODULE                OFF
PROFILE_MODULE             ON
NICS_SCAN_XY_MODULE        OFF

# --------------------------------
# Profile (1D normal to ring)
# --------------------------------
# Required if PROFILE_MODULE = ON
#
# Linear scale requires:
#   START, DELTA_Z, NUM_GA
#
# Logarithmic scale requires:
#   START, END, NUM_GA

BEGIN_PROFILE

SCALE       Linear               # Linear or Logarithmic

START       0.0
DELTA_Z     0.1
NUM_GA      100

# Only for Logarithmic scale:
# END      5.0

END_PROFILE
```

## NICS-scan module
### Input format

```text
# ================================
#        AroMagnetic INPUT
# ================================

# --------------------------------
# General Information
# --------------------------------

NAME_MOLECULE             mol
PROGRAM                   Gaussian            # ORCA or GAUSSIAN
KEYWORD_LINE              Functional/Basis NMR INTEGRAL(SUPERFINEGRID) SCF=XQC

PROCS                     16 #If you use Orca this value doesn't apear in the input files
MEMORY_GB                 32

CHARGE                    0
MULTIPLICITY              1

MOLECULAR_COORDINATES     mol.xyz
PSEUDO_PI                 NO              # YES or NO


# --------------------------------
# Module Switchboard
# --------------------------------
# Turn ON only what you need.

GRID_MODULE                OFF
PROFILE_MODULE             OFF
NICS_SCAN_XY_MODULE        ON

# --------------------------------
# NICS Scan XY (Transversal scan)
# --------------------------------
# Required if NICS_SCAN_XY_MODULE = ON
#
# If both POINTS and DELTA_X are defined:
#   Code checks consistency. If inconsistent → error.

BEGIN_NICS_SCAN_XY

POINTS        81
MARGIN        2.0
HEIGHT        0.0 1.0 1.5
DELTA_X       0.1

MAKE_AXIS_Y   YES            # YES or NO

END_NICS_SCAN_XY
```

## Interfaz with Orca
### Orca section

```text
# --------------------------------
# ORCA Block (Optional)
# --------------------------------
# Used only if PROGRAM = ORCA.
# If present, this block overrides PROCS/MEMORY_GB
# and adds advanced ORCA-specific settings.

BEGIN_ORCA_BLOCK

%pal
  nprocs 16
end

%maxcore 32000

END_ORCA_BLOCK
```

## How to run AroMagnetic?

1. Prepare your `INPUT.txt` file and your molecular structure file (XYZ format).
    - If you are using generic databases in your calculations, you must create the *basis.txt* file and enter your database, leaving two spaces at the end of the file.
2. Generate the quantum chemistry input files:
```text
work_amg
```
3. Run the generated input files using Gaussian or ORCA.
4. Once the calculations are finished, analyze the results:
```text
analyze_amg
```
This command generates .cube files (and .vtk files if selected in INPUT.txt), .csv files and/or graph, depending on the module you've chosen.

**Molecular orbital decomposition of the shielding tensor (Gaussian only)**

To compute molecular orbital contributions to the shielding tensor, follow these steps (This module works with NBO6).

1. Include the following keywords in the KEYWORD_LINE section of your INPUT.txt:
```text
functional/basis NMR INTEGRAL(SUPERFINEGRID) SCF=XQC POP=NBO6read
```
2. Generate the Gaussian input files:
```text
work_amg
```
3. Run the Gaussian calculations.
4. After the calculations are completed, perform the molecular orbital analysis:
```text
mo_amg
```

