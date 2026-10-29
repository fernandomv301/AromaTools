# AroGeometric

AroGeometric is the AromaTools module for evaluating aromaticity using **geometric descriptors** derived from molecular structure.

The module provides an automated workflow for detecting molecular rings, evaluating bond-length equalization, and computing aromaticity indices commonly used in computational chemistry.

---

## Overview

AroGeometric automatically:

- Detects molecular rings
- Constructs the molecular graph
- Computes bond distances
- Evaluates geometric aromaticity indices
- Generates tables with bond-length information
- Saves results as CSV files
- Produces graphical representations of cyclic systems

Only molecular geometry is required as input.

Supported file formats:

- Gaussian output files (`.log`)
- ORCA output files (`.out`)
- XYZ files (`.xyz`)

---

## Implemented indices

The following geometric descriptors are implemented:

- HOMA93
- HOMAc
- HOMER
- BLA (Bond Length Alternation)

These indices quantify aromaticity through the degree of bond-length equalization in cyclic systems.

---

## Quick start

AroGeometric runs directly from the terminal. No scripting required.

### Run the program

```bash
arogeometric
```

Or, if you are working directly with the script:

```bash
python3 arogeometric.py
```

### Provide the required inputs

The program will ask for two things:

```
Output file (Gaussian or ORCA): benzene.log
Index type [HOMA/HOMER]: HOMA
```

| Input | Description |
|-------|-------------|
| Output file | Path to your Gaussian (`.log`), ORCA (`.out`), or XYZ (`.xyz`) file |
| Index type | `HOMA` for ground state (S₀) or `HOMER` for excited state (T₁) |

### Results

After the calculation, the results are printed to the screen and saved automatically as CSV files in the same folder as the input file:

```
benzene_indices.csv
benzene_bonds.csv
```

The program then asks if you want to open the interactive 3D visualization:

```
Show 3D visualization? [y/n]:
```

---

## Workflow

The module follows these steps:

1. Read molecular geometry from the input file.
2. Build the molecular graph based on covalent radii.
3. Detect rings (up to 8-membered rings).
4. Compute bond distances.
5. Evaluate aromaticity indices.
6. Print and save results.
7. Generate visualization (optional).

---

## Theoretical background

### HOMA index

The harmonic oscillator model of aromaticity is defined as:

$$
HOMA = 1 - \frac{1}{n}\sum_{i=1}^{n} \alpha (R_i - R_{opt})^2
$$

where:

- $n$ is the number of bonds
- $R_i$ is the bond length
- $R_{opt}$ is the optimal aromatic bond length
- $\alpha$ is a normalization constant

HOMA values close to 1 indicate aromatic character, values near 0 correspond to non-aromatic systems, and negative values indicate antiaromaticity.

---

### Parameter sets

**HOMA93 parameters**

| Bond | Ropt (Å) | α (Å⁻²) |
|------|----------|---------|
| CC   | 1.388    | 257.70  |
| CN   | 1.334    | 93.52   |
| NN   | 1.309    | 130.33  |
| CO   | 1.265    | 157.38  |

Reference: [HOMA article](https://www.sciencedirect.com/science/article/pii/S0040403901941759)

**HOMAc parameters (ground state S₀)**

| Bond | Ropt (Å) | α (Å⁻²) |
|------|----------|---------|
| CC   | 1.392    | 153.37  |
| CN   | 1.333    | 111.83  |
| NN   | 1.318    | 98.99   |
| CO   | 1.315    | 335.16  |

Reference: [HOMAc article](https://pubs.acs.org/doi/full/10.1021/acs.joc.4c02475)

**HOMER parameters (excited state T₁)**

| Bond | Ropt (Å) | α (Å⁻²) |
|------|----------|---------|
| CC   | 1.437    | 950.74  |
| CN   | 1.390    | 506.43  |
| NN   | 1.375    | 187.36  |
| CO   | 1.379    | 164.96  |

Reference: [HOMER article](https://pubs.rsc.org/en/content/articlehtml/2023/cp/d3cp00842h)

---

## Output

### Aromaticity indices

Printed to the screen and saved as `<filename>_indices.csv`.

| Column | Description |
|--------|-------------|
| Ring | Ring identifier |
| HOMA93 | HOMA index with experimental parameters |
| HOMAc | HOMA index with computational parameters |
| HOMER | HOMER index (only when `index_type="HOMER"`) |
| BLA | Bond Length Alternation |

HOMA93 and HOMAc are both reported when `index_type="HOMA"`. Only HOMER is reported when `index_type="HOMER"`.


### Bond lengths

Printed to the screen and saved as `<filename>_bonds.csv`.

| Column | Description |
|--------|-------------|
| Ring | Ring identifier |
| Type_Bond | Bond type (e.g. CC, CN, CO) |
| Length (Å) | Bond length in Ångströms |


### Graphical output

If the 3D visualization is requested, an interactive plot opens in the browser showing:

- Molecular connectivity
- Detected rings with their index values
- Ring labels

