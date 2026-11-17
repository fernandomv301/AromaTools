# Installation

AromaTools requires **Python 3.10 or newer**. It can be installed from PyPI or directly from its Git repository. In both cases, working in a virtual environment is recommended. A virtual environment isolates AromaTools and its dependencies from the rest of your Python installation, which avoids version conflicts with other programs.

```bash
python3 -m venv aromatools-env
source aromatools-env/bin/activate
```

Alternatively, with conda: `conda create -n aromatools python=3.10` followed by `conda activate aromatools`.

## Install via pip

```bash
pip install aromatools
```

All the required libraries (including ESIpy, PySCF, and IOData) are installed automatically; there is no need to install them manually.

## Install from Git

The source code is available at [github.com/fernandomv301/AromaTools](https://github.com/fernandomv301/AromaTools). Git must be installed on your computer.

**Directly with pip** (no manual download required):

```bash
pip install git+https://github.com/fernandomv301/AromaTools.git
```

**By cloning the repository** (download the code and install it from the local folder):

```bash
git clone https://github.com/fernandomv301/AromaTools.git
cd AromaTools
pip install .
```

To modify the code and apply the changes without reinstalling, use an editable installation (`pip install -e .`) instead of `pip install .`.

To update an installation made from Git:

```bash
pip install --upgrade git+https://github.com/fernandomv301/AromaTools.git
```

If you cloned the repository, run `git pull` in the repository folder and then `pip install --upgrade .`.

## Check the installation

The package can be imported from Python:

```bash
python -c "import aromatools"
```

and the three modules are installed as commands: `aromagnetic` (also available as `amg`), `arogeometric`, and `aroelectronic` (also available as `ael`). Running `amg` or `ael` without arguments lists their available commands, `arogeometric` prompts for the calculation output file, and `ael run` runs AroElectronic using the `INPUT.txt` file in the current folder. See each module's page for details.

### Quick reference

| Module | Command | Description |
|--------|---------|-------------|
| AroGeometric | `arogeometric` | Geometric indices (HOMA, HOMAc, HOMER, BLA); asks for the output file of the calculation |
| AroMagnetic | `aromagnetic` or `amg` | Magnetic criteria; `amg run`, `amg get`, ... (see the AroMagnetic page) |
| AroElectronic | `aroelectronic` or `ael` | Electronic criteria; `ael run` (see the AroElectronic page) |

## Requirements

Python 3.10+

Required libraries (installed automatically with AromaTools):

- numpy
- pandas
- networkx
- matplotlib
- plotly
- ase
- pymatgen
- esipython 1.1.0 or newer (required by AroElectronic). It is imported as `esipy` and installs [PySCF](https://pyscf.org/) and [IOData](https://github.com/theochem/iodata) as dependencies.

---
