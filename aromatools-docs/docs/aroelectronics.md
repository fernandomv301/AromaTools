# AroElectronics

AroElectronics is the AromaTools module for evaluating aromaticity using **electronic criteria**, that is, indices that quantify how many electrons are shared (delocalized) among the atoms of a ring.

AroElectronics **does not reimplement the indices**: it imports from [ESIpy](https://github.com/jgrebol/ESIpy) the functions that build the atomic overlap matrices and compute the indices (`ESI`, the wavefunction reader and the reader of AIMAll files). AroElectronics prepares the wavefunction, defines the ring and the atomic partition, runs ESIpy and writes the results, so the values reported are those of ESIpy. If you use this module, please cite ESIpy and the works of the indices and partitions you used (see [Methods and references](methods.md)). The theory of the indices and the comparison between partitions are described in the ESIpy paper (Grèbol-Tomàs *et al.*, *Chem. Eur. J.* **2024**, *30*, e202401282) and are not repeated here.

---

## Overview

Electronic indices are built from the **delocalization index** (DI) between pairs of atoms,

$$
\delta(A,B)=4\, \sum_{i,j}{S}_{ij}\left({A}\right){S}_{ij}\left({B}\right)
$$

where ${S}^{A}$ is the atomic overlap matrix (AOM) of atom $A$ in the basis of the occupied molecular orbitals. The AOMs depend on how an atom in the molecule is defined, so **the value of every index depends on the chosen partition** (see [Partitions](#partitions)).

AroElectronics automatically:

- Runs a Kohn–Sham calculation with PySCF, or reads a wavefunction file (see [Wavefunction files](#wavefunction-files))
- Detects the rings of the molecule (or uses the ring given by the user)
- Builds the AOMs with the requested partition
- Computes the requested indices for every ring
- Saves one `.csv` file per index and a summary `.txt` file

---

## Implemented indices

| Index | Description | Aromatic character |
|-------|-------------|--------------------|
| $I_{ring}$ | Multicenter bond index (Giambiagi) | Higher is more aromatic |
| MCI | Multicenter index (Bultinck) | Higher is more aromatic |
| PDI | Para-delocalization index | Higher is more aromatic |
| FLU | Aromatic fluctuation index | **Lower** is more aromatic (0 = ideal) |
| AV1245 | Average of the 4-center indices of 1-2-4-5 fragments | Higher is more aromatic |
| AVmin | Minimum of the 4-center indices | Higher is more aromatic |
| BOA | Bond-order alternation: difference between the mean DI of the alternate bonds of the ring | **Lower** is more aromatic (0 = full equalization) |
| BOAc | Variant of BOA built from the difference between **consecutive** bonds (see below) | **Lower** is more aromatic (0 = full equalization) |
| DI | Delocalization index $\delta(A,B)$ of each bond of the ring | Not an index by itself: the DI of the C–C bonds of an aromatic ring are close to each other |
| Populations | Atomic population $N(A)$, the number of electrons assigned to each atom of the ring | The atomic charge is $Z_A - N(A)$ |

**DI and populations.** They are printed for the atoms of each ring (the DI for every bond between consecutive atoms, including the one that closes the ring) and are the quantities from which the indices above are built. As with the indices, they depend on the chosen partition.

**BOA and BOAc.** Both are the electronic counterpart of the bond-length alternation (BLA), using the DI of the bonds of the ring, $\delta_i = \delta(A_i, A_{i+1})$, instead of the bond lengths:

$$
BOA = \left| \frac{1}{n_1}\sum_{i\ \mathrm{odd}} \delta_i - \frac{1}{n_2}\sum_{i\ \mathrm{even}} \delta_i \right|
\qquad
BOA_c = \frac{1}{n}\sum_{i=1}^{n} \left| \delta_i - \delta_{i+1} \right|
$$

BOA compares the two sets of alternate bonds, so it can vanish when the bonds alternate in a pattern that cancels out in the averages, and it is not well defined for rings with an odd number of atoms. BOAc compares each bond with the next one around the ring, so it detects any irregularity and is valid for rings of any size. In a perfectly regular ring both are zero.

$I_{ring}$, MCI and AV1245 are reported also in the normalized form $x^{1/n}$, where $n$ is the number of atoms of the ring, so that rings of different size can be compared.

**Limitations by ring size**

- PDI is defined only for six-membered rings.
- AV1245 and AVmin require rings of at least six atoms.
- The cost of MCI grows factorially with the number of ring atoms; it is reported as `N/A` for rings of more than 12 atoms. AV1245 and AVmin are the alternatives for large rings.
- FLU needs a reference DI for every type of bond in the ring (CC, CN, BN, NN, CS, CO).

---

## Quick start

AroElectronics runs from the terminal with a single input file. Go to the folder that contains `INPUT.txt` and run:

```bash
ael run
```

The program can be called with `ael` or with `aroelectronics`; the two names are equivalent. With no input file given, `ael run` reads the file `INPUT.txt` of the current folder, whatever its capitalization (`input.txt`, `Input.txt`, ...). The input file does not have to be called `INPUT.txt`: to use a file with another name, write it after the command:

```bash
ael run benzene.txt
```

**Available commands**

| Command | Description |
|---------|-------------|
| `ael run [INPUT]` | Computes the indices requested in `INPUT.txt` (or in the file `INPUT`, if one is given) |
| `ael help` | Shows the list of commands (the same happens when `ael` is run with no arguments) |

The form `aroelectronics INPUT.txt`, without the word `run`, also works.

There are three ways of obtaining the atomic overlap matrices: from a PySCF calculation, from a wavefunction file, or from the atomic files of a QTAIM calculation with AIMAll (see [QTAIM](#qtaim)).

**1. PySCF (from scratch).** Provide the coordinates, the functional, the basis set, the charge and the multiplicity:

```text
NAME_MOLECULE         benzene
MOLECULAR_COORDINATES benzene.xyz
FUNCTIONAL            B3LYP
BASIS_SET             aug-cc-pVDZ
CHARGE                0
MULTIPLICITY          1

PARTITION             meta-lowdin
RING                  1 2 6 10 7 3

IRING                 YES
MCI                   YES
PDI                   YES
FLU                   YES
AV1245                YES
BOA                   YES
DI                    YES
POPULATIONS           YES
```

**2. Wavefunction file.** Use a file generated by another program (for example `.fchk` or `.molden`). The keyword `WAVEFUNCTION_FILE` replaces the coordinates, functional, basis set, charge and multiplicity:

```text
NAME_MOLECULE         benzene
WAVEFUNCTION_FILE     benzene.fchk

PARTITION             nao
IRING                 YES
MCI                   YES
```

---

## Input keywords

Everything after `#` in a line is ignored. Keywords are not case sensitive.

| Keyword | Description | Default |
|---------|-------------|---------|
| `NAME_MOLECULE` | Prefix of the output files | `molecule` |
| `MOLECULAR_COORDINATES` | `.xyz` file (PySCF mode) | required in PySCF mode |
| `FUNCTIONAL` | Exchange–correlation functional (PySCF mode) | required in PySCF mode |
| `BASIS_SET` | Basis set (PySCF mode) | required in PySCF mode |
| `CHARGE` | Total charge | `0` |
| `MULTIPLICITY` | Spin multiplicity. `1` runs a restricted calculation (RKS); any other value runs an unrestricted one (UKS) | `1` |
| `WAVEFUNCTION_FILE` | Wavefunction file (file mode), see [Wavefunction files](#wavefunction-files) | none |
| `PARTITION` | Atomic partition (see below) | `meta-lowdin` |
| `ATOMICFILES_PATH` | Folder with the AIMAll `.int` files (only with `PARTITION qtaim`), see [QTAIM](#qtaim) | none |
| `RING` | Atom numbers of **one** ring, **in cyclic order**, starting at 1 | all rings are detected |
| `IRING`, `MCI`, `PDI`, `FLU`, `AV1245`, `BOA`, `DI`, `POPULATIONS` | `YES` to compute the index (`AV1245` also gives AVmin, and `BOA` also gives BOAc). `DI` gives the delocalization indices of the ring bonds and `POPULATIONS` the atomic populations of the ring atoms | `NO` |

At least one index must be requested.

---

## Wavefunction files

Wavefunction files are read through [IOData](https://github.com/theochem/iodata), so any format that IOData can load is accepted (for example `.fchk` and `.molden`, written by Gaussian, ORCA, PySCF and other programs). The file is converted to an `.fchk` file named `<file>.esipy.fchk`, which is created in the same folder as the original file; a Gaussian `.fchk` file is used directly. The molecule and the wavefunction are then read from it, so no additional calculation is performed.

The file must contain the molecular orbitals and the basis set. A closed-shell wavefunction is treated as restricted and an open-shell one as unrestricted.

---

## Partitions

| `PARTITION` | Description |
|-------------|-------------|
| `meta-lowdin` | Löwdin orthogonalization of the atomic orbitals using the occupied-space weights (meta-Löwdin) |
| `nao` | Natural atomic orbitals |
| `iao` | Intrinsic atomic orbitals |
| `lowdin` | Löwdin symmetric orthogonalization |
| `mulliken` | Mulliken partition |
| `qtaim` | Real-space partition of the Quantum Theory of Atoms in Molecules, read from AIMAll files (see [QTAIM](#qtaim)) |

The three first partitions (`meta-lowdin`, `nao`, `iao`) are recommended: their indices are largely independent of the basis set. Mulliken and Löwdin depend strongly on it and are included mainly for comparison.

### QTAIM

With `PARTITION qtaim` the AOMs are not computed by AromaTools: they are read from the atomic files of an AIMAll calculation, so no PySCF calculation or wavefunction file is needed.

```text
NAME_MOLECULE         benzene
PARTITION             qtaim
ATOMICFILES_PATH      benzene_atomicfiles

IRING                 YES
MCI                   YES
PDI                   YES
FLU                   YES
AV1245                YES
BOA                   YES
```

- `ATOMICFILES_PATH` is the folder generated by AIMAll with one `.int` file per atom (for example `benzene_atomicfiles`). The files of **all** the atoms of the molecule must be in the folder, including the hydrogens, because the atom numbers are taken from them.
- Each `.int` file must contain the atomic overlap matrix. Only the occupied orbitals are used, so the `.wfx` file that the `.int` files refer to is also needed; it is searched in the same folder and up to three folders above.
- The wavefunction type is detected from the files: restricted (closed-shell) or unrestricted (open-shell). Wavefunctions with fractional occupations (natural orbitals) are not supported.
- The geometry for the automatic detection of rings is read from the `.wfx` file. If it is not found, give the ring with `RING`.
- `MOLECULAR_COORDINATES`, `FUNCTIONAL`, `BASIS_SET`, `CHARGE`, `MULTIPLICITY` and `WAVEFUNCTION_FILE` are not used.

Real-space partitions are much more expensive than the Hilbert-space ones, since the AOMs come from a numerical integration over the atomic basins, and the quality of the integration affects the indices.

!!! note
    FLU uses reference DI values that are specific to each partition and bond type. They are taken from ESIpy.

---

## Rings

If `RING` is not given, all the rings of the molecule are detected from the covalent radii (a pair of atoms is bonded if their distance is below 1.3 times the sum of their covalent radii) and the indices are computed for each of them. The atoms of each ring are listed in cyclic order.

To analyze a specific ring, give its atoms with `RING`. The atoms must be written **in the order in which they are connected around the ring**, because the multicenter indices depend on that order. For benzene, for example, `1 2 6 10 7 3` is a valid order only if atom 1 is bonded to 2, 2 to 6, and so on, and 3 is bonded back to 1.

---

## Open-shell systems

For `MULTIPLICITY` different from 1 the calculation is unrestricted and the indices, the DI and the populations are the sum of the α and β contributions. For closed-shell systems the single-spin result is doubled by ESIpy for Iring, MCI, PDI and AV1245, whereas FLU is not doubled; the same convention is used in the unrestricted case, where FLU is obtained as the sum of the α and β values.

---

## Results

For a molecule named `benzene` with the `meta-lowdin` partition the program writes:

```
benzene_meta-lowdin_iring.csv
benzene_meta-lowdin_mci.csv
benzene_meta-lowdin_pdi.csv
benzene_meta-lowdin_flu.csv
benzene_meta-lowdin_av1245.csv
benzene_meta-lowdin_boa.csv
benzene_meta-lowdin_di.csv
benzene_meta-lowdin_populations.csv
benzene_meta-lowdin_results.txt
```

Only the files of the requested indices are created. The `.txt` file summarizes the setup (mode, partition, rings) and all the indices of every ring. The screen output has the following format (one block per ring):

```
Ring 1: [1-2-6-10-7-3]
  Atomic populations (electrons):
    C1 = <value>
    ...
  Delocalization indices (ring bonds):
    C1-C2 = <value>
    ...
  Iring    = <value>   (^1/n = <value>)
  MCI      = <value>   (^1/n = <value>)
  PDI      = <value>
  FLU      = <value>
  AV1245   = <value>   AVmin = <value>
  BOA      = <value>   BOAc = <value>
```

---

## Requirements

AroElectronics needs [ESIpy](https://github.com/jgrebol/ESIpy) 1.1.0 or newer (and therefore Python 3.10+) and, to run calculations from scratch, [PySCF](https://pyscf.org/). Both are installed with AromaTools (see [Installation](installation.md)). The QTAIM mode needs the files of an [AIMAll](https://aim.tkgristmill.com/) calculation.

If you install them by hand, note that the package is called `esipython` on PyPI, but it is imported as `esipy`:

```bash
pip install "esipython>=1.1.0" pyscf
```

Do not install `esipy`: that is a different project (an EVE Online API client).

---

