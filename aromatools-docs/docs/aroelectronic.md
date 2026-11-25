# AroElectronic

AroElectronic is the AromaTools module for evaluating aromaticity using **electronic criteria**, that is, indices that quantify how many electrons are shared (delocalized) among the atoms of a ring.

AroElectronic **does not reimplement the indices**: it imports the functions that build the atomic overlap matrices and compute the indices from [ESIpy](https://github.com/jgrebol/ESIpy) (`ESI`, the wavefunction reader, and the AIMAll file reader). AroElectronic prepares the wavefunction, defines the ring and the atomic partition, runs ESIpy and writes the results, so the reported values are those computed by ESIpy. If you use this module, please cite ESIpy and the original works describing the indices and partitions you used (see [Methods and references](methods.md)). The theory of the indices and the comparison between partitions are described in the ESIpy paper (Grèbol-Tomàs *et al.*, *Chem. Eur. J.* **2024**, *30*, e202401282) and are not repeated here.

---

## Overview

Electronic indices are built from the **delocalization index** (DI) between pairs of atoms,

$$
\delta(A,B)=4\, \sum_{i,j}{S}_{ij}\left({A}\right){S}_{ij}\left({B}\right)
$$

where $S_{ij}(A)$ are the elements of the atomic overlap matrix (AOM) of atom $A$ in the basis of the occupied molecular orbitals. The AOMs depend on how an atom in the molecule is defined, so **the value of every index depends on the chosen partition** (see [Partitions](#partitions)).

AroElectronic automatically:

- Runs a Kohn–Sham calculation with PySCF, or reads a wavefunction file (see [Wavefunction files](#wavefunction-files))
- Detects the rings of the molecule (or uses the ring specified by the user)
- Builds the AOMs with the requested partition
- Computes the requested indices for every ring
- Saves all the results in a `.txt` file and, optionally (`CSV YES`), one `.csv` file per index

---

## Implemented indices

| Index | Description | Aromatic character |
|-------|-------------|--------------------|
| $I_{ring}$ | Multicenter bond index (Giambiagi) | Higher is more aromatic |
| MCI | Multicenter index (Bultinck) | Higher is more aromatic |
| PDI | Para-delocalization index | Higher is more aromatic |
| FLU | Aromatic fluctuation index | **Closer to 0** is more aromatic (0 = ideal) |
| AV1245 | Average of the 4-center indices of 1-2-4-5 fragments | Higher is more aromatic |
| AV<sub>min</sub> | Minimum of the 4-center indices | Higher is more aromatic |
| BOA | Bond-order alternation: difference between the mean DIs of the two sets of alternating bonds in the ring | **Lower** is more aromatic (0 = full equalization) |
| BOAc | Variant of BOA built from the difference between **consecutive** bonds (see below) | **Lower** is more aromatic (0 = full equalization) |
| DI | Delocalization index $\delta(A,B)$ of each bond of the ring | Not an index by itself: the DI of an aromatic C–C bond is close to 1.4, and the C–C bonds of an aromatic ring have similar DI values |
| Populations | Atomic population $N(A)$, the number of electrons assigned to each atom of the ring | The atomic charge is $Z_A - N(A)$ |

**DI and populations.** Both are printed for each ring (the DI for every bond between consecutive atoms, including the bond that closes the ring). They are the quantities from which the indices above are built. As with the indices, they depend on the chosen partition.

**BOA and BOAc.** Both are electronic counterparts of the bond-length alternation (BLA), using the DI of the bonds of the ring, $\delta_i = \delta(A_i, A_{i+1})$, instead of the bond lengths:

$$
BOA = \left| \frac{1}{n_1}\sum_{i\ \mathrm{odd}} \delta_i - \frac{1}{n_2}\sum_{i\ \mathrm{even}} \delta_i \right|
\qquad
BOA_c = \frac{1}{n}\sum_{i=1}^{n} \left| \delta_i - \delta_{i+1} \right|
$$

BOA compares the two sets of alternate bonds, so it can vanish when the bonds alternate in a pattern that cancels out in the averages, and it is not well defined for rings with an odd number of atoms. BOAc compares each bond with the next one around the ring, so it detects any irregularity and is valid for rings of any size. In a perfectly regular ring both are zero.

$I_{ring}$ and MCI are also reported in the normalized form $x^{1/n}$, where $n$ is the number of atoms in the ring, so that rings of different sizes can be compared. AV1245 and AV<sub>min</sub> are not normalized.

**Limitations by ring size**

- PDI is defined only for six-membered rings.
- AV1245 and AV<sub>min</sub> require rings of at least six atoms; for smaller rings the program reports `N/A`.
- The cost of MCI grows factorially with the number of ring atoms; it is reported as `N/A` for rings with more than 12 atoms. AV1245 and AV<sub>min</sub> are the recommended alternatives for large rings.
- FLU requires a reference DI for every type of bond in the ring (CC, CN, BN, NN, CS, CO).

---

## Quick start

AroElectronic runs from the terminal with a single input file. Go to the folder that contains `INPUT.txt` and run:

```bash
ael run
```

The program can be invoked as either `ael` or `aroelectronic`. If no input file is specified, `ael run` reads the `INPUT.txt` file in the current folder, regardless of capitalization (`input.txt`, `Input.txt`, etc.). The input file does not have to be called `INPUT.txt`: to use a file with a different name, specify it after the command:

```bash
ael run benzene.txt
```

### Available commands

| Command | Description |
|---------|-------------|
| `ael run [INPUT]` | Computes the indices requested in `INPUT.txt` (or in the specified file, if one is given) |
| `ael help` | Displays the list of available commands (equivalent to running `ael` with no arguments) |

The shorter form `aroelectronic INPUT.txt` (without `run`) also works.

### Input modes

There are three ways of obtaining the atomic overlap matrices: from a PySCF calculation, from a wavefunction file, or from the atomic files of a QTAIM calculation with AIMAll (see [QTAIM](#qtaim)).

**1. PySCF (from scratch).** Provide the coordinates, the functional, the basis set, the charge and the multiplicity:

```text
NAME_MOLECULE         benzene
MOLECULAR_COORDINATES benzene.xyz
FUNCTIONAL            B3LYP
BASIS_SET             aug-cc-pVDZ
CHARGE                0
MULTIPLICITY          1
CSV                   YES

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

**2. Wavefunction file.** Use a file generated by another program (for example, `.fchk` or `.molden`). The keyword `WAVEFUNCTION_FILE` replaces the coordinates, functional, basis set, charge and multiplicity:

```text
NAME_MOLECULE         benzene
WAVEFUNCTION_FILE     benzene.fchk

PARTITION             nao
IRING                 YES
MCI                   YES
```

---

## Input keywords

Everything after `#` on a line is ignored. Keywords are case-insensitive.

| Keyword | Description | Default |
|---------|-------------|---------|
| `NAME_MOLECULE` | Prefix of the output files | `molecule` |
| `MOLECULAR_COORDINATES` | `.xyz` file (PySCF mode) | required in PySCF mode |
| `FUNCTIONAL` | Exchange–correlation functional (PySCF mode) | required in PySCF mode |
| `BASIS_SET` | Basis set (PySCF mode) | required in PySCF mode |
| `CHARGE` | Total charge | `0` |
| `MULTIPLICITY` | Spin multiplicity. `1` runs a restricted calculation (RKS); any other value runs an unrestricted one (UKS) | `1` |
| `WAVEFUNCTION_FILE` | Wavefunction file (file mode), see [Wavefunction files](#wavefunction-files) | none |
| `PARTITION` | Atomic partition (see below). Several partitions can be written in the same line, or `ALL`, see [Several partitions at once](#several-partitions-at-once) | `meta-lowdin` |
| `ATOMICFILES_PATH` | Folder with the AIMAll `.int` files (only with `PARTITION qtaim`), see [QTAIM](#qtaim) | none |
| `RING` | Atom numbers of **one** ring (or of any path of atoms that you want to study), **in the order to follow**, starting at 1 | all rings are detected |
| `IRING`, `MCI`, `PDI`, `FLU`, `AV1245`, `BOA`, `DI`, `POPULATIONS` | `YES` to compute the index (`AV1245` also gives AV<sub>min</sub>, and `BOA` also gives BOAc). `DI` gives the delocalization indices of the ring bonds and `POPULATIONS` the atomic populations of the ring atoms | `NO` |
| `CSV` | `YES` to also write the results as `.csv` files, one per requested index and partition (see [Results](#results)). The `.txt` file with all the results is always written | `NO` |

At least one index must be requested.

---

## Wavefunction files

Wavefunction files are read through [IOData](https://github.com/theochem/iodata), so any format that IOData can load is accepted (for example `.fchk` and `.molden`, written by Gaussian, ORCA, PySCF and other programs). Files in other formats are converted to an `.fchk` file named `<file>.esipy.fchk` in the same folder as the original; Gaussian `.fchk` files are used directly. The molecule and the wavefunction are then read from it, so no additional calculation is performed.

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

The first three partitions (`meta-lowdin`, `nao`, `iao`) are recommended because the resulting indices are largely independent of the basis set. The Mulliken and Löwdin partitions are strongly basis-set dependent and are included mainly for comparison.

### Several partitions at once

Because the value of every index depends on the partition, it is good practice to check that your conclusions do not change from one partition to another. You do not need one `INPUT.txt` for each partition: write all of them in the `PARTITION` line, separated by spaces or commas.

```text
NAME_MOLECULE         benzene
MOLECULAR_COORDINATES benzene.xyz
FUNCTIONAL            B3LYP
BASIS_SET             aug-cc-pVDZ
CHARGE                0
MULTIPLICITY          1

PARTITION             meta-lowdin nao iao      # or: meta-lowdin, nao, iao
RING                  1 2 6 10 7 3

IRING                 YES
MCI                   YES
AV1245                YES
```

The keyword `ALL` is a shortcut for `meta-lowdin nao iao lowdin mulliken`; `qtaim` is not included in it, because it needs the AIMAll files, so write it explicitly if you want it (it can be combined with the others, see below).

- The PySCF calculation (or the reading of the wavefunction file) is done **only once** and shared by all the partitions, so each additional partition only costs the construction of its atomic overlap matrices and the indices.
- The results of all the partitions are collected in a single text file, `benzene_results.txt`, and the screen shows one block per partition. With `CSV YES`, each partition also has its own `.csv` files (`benzene_meta-lowdin_iring.csv`, `benzene_nao_iring.csv`, ...), exactly as if you had run each partition separately.
- When more than one partition is computed, the program prints on the screen, and writes in the `.txt` file, a table with the requested Iring, MCI, PDI, FLU, AV1245/AVmin and BOA/BOAc of every ring for every partition, so they can be compared at a glance. With `CSV YES` the same table is also saved as `<NAME_MOLECULE>_partitions_summary.csv`.
- If a partition fails, the program prints a warning, skips it and continues with the others. If a single partition was requested, the error stops the program as usual.
- To combine `qtaim` with other partitions, give both the data of the other partitions (`MOLECULAR_COORDINATES`, `FUNCTIONAL` and `BASIS_SET`, or `WAVEFUNCTION_FILE`) and `ATOMICFILES_PATH`. The QTAIM values come from the AIMAll calculation and the rest from the wavefunction, so make sure that both were obtained for the same geometry and level of theory.

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

- `ATOMICFILES_PATH` is the folder generated by AIMAll with one `.int` file per atom (for example `benzene_atomicfiles`). The folder must contain the files for **all** atoms in the molecule, including hydrogens, because the atom numbering is taken from them.
- Each `.int` file must contain the atomic overlap matrix. Only the occupied orbitals are used, so the `.wfx` file that the `.int` files refer to is also needed; the program searches for it in the same folder and up to three parent folders.
- The wavefunction type is detected from the files: restricted (closed-shell) or unrestricted (open-shell). Wavefunctions with fractional occupations (natural orbitals) are not supported.
- The geometry used for automatic ring detection is read from the `.wfx` file. If this file is not found, specify the ring with `RING`.
- `MOLECULAR_COORDINATES`, `FUNCTIONAL`, `BASIS_SET`, `CHARGE`, `MULTIPLICITY` and `WAVEFUNCTION_FILE` are not used.

Real-space partitions are much more computationally expensive than Hilbert-space partitions, since the AOMs come from a numerical integration over the atomic basins, and the quality of the integration affects the indices.

!!! note
    FLU uses reference DI values that are specific to each partition and bond type. They are taken from ESIpy.

---

## Rings

If `RING` is not specified, all the rings of the molecule are detected from the covalent radii (a pair of atoms is bonded if their distance is below 1.3 times the sum of their covalent radii) and the indices are computed for each of them. The atoms of each ring are listed in cyclic order.

To analyze a specific ring, list its atoms with `RING`. The atoms must be written **in the order in which they are connected around the ring**, because the multicenter indices depend on that order. For benzene, for example, `1 2 6 10 7 3` is a valid order only if atom 1 is bonded to atom 2, atom 2 to atom 6, and so on, and atom 3 is bonded back to atom 1.

`RING` is not limited to rings of the molecule: you can also list any path of atoms that you want to follow (for example, the atoms of a conjugated chain), and the indices are computed for the atoms in the order written.

---

## Open-shell systems

When `MULTIPLICITY` is not 1, the calculation is unrestricted, and the indices, DIs, and populations are obtained as the sum of the α and β contributions. For closed-shell systems, ESIpy doubles the single-spin result for Iring, MCI, PDI, and AV1245, but not for FLU. The same convention is applied in the unrestricted case, where FLU is obtained as the sum of the α and β values.

---

## Results

The results of **all** the partitions are always written in a **single** text file, `<NAME_MOLECULE>_results.txt`. It summarizes the setup (mode, partitions, rings), then gives a block with all the indices of every ring for each partition and, when there are several partitions, the comparison table between them.

With `CSV YES`, the program also writes one `.csv` file per requested index and partition. For a molecule named `benzene` with the `meta-lowdin` partition and all the indices requested, it writes:

```
benzene_meta-lowdin_iring.csv
benzene_meta-lowdin_mci.csv
benzene_meta-lowdin_pdi.csv
benzene_meta-lowdin_flu.csv
benzene_meta-lowdin_av1245.csv
benzene_meta-lowdin_boa.csv
benzene_meta-lowdin_di.csv
benzene_meta-lowdin_populations.csv
benzene_results.txt
```

The `.csv` files carry the name of the partition, and they are created only for the requested indices. If several partitions are requested, the same set of `.csv` files is written for each of them, plus `benzene_partitions_summary.csv` with the comparison table. Without `CSV YES` (the default), only the `.txt` file is written. The screen output has the following format (one block per ring):

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

AroElectronic requires [ESIpy](https://github.com/jgrebol/ESIpy) 1.1.0 or newer (and therefore Python 3.10+) and, for calculations from scratch, [PySCF](https://pyscf.org/). Both are installed with AromaTools (see [Installation](installation.md)). The QTAIM mode requires the output files of an [AIMAll](https://aim.tkgristmill.com/) calculation.

If you install them manually, note that the package is named `esipython` on PyPI but is imported as `esipy`:

```bash
pip install "esipython>=1.1.0" pyscf
```

---

