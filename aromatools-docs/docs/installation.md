# Installation

AromaTools requires **Python 3.10 or newer**. It can be installed from PyPI or directly from its Git repository. In both cases it is recommended to work in a virtual environment.

```bash
python3 -m venv aromatools-env
source aromatools-env/bin/activate
```

With conda: `conda create -n aromatools python=3.10` and `conda activate aromatools`.

## Install via pip

```bash
pip install aromatools
```

## Install from Git

The source code is available at [github.com/fernandomv301/AromaTools](https://github.com/fernandomv301/AromaTools). Git must be installed on your computer.

**Directly with pip** (no need to download anything by hand):

```bash
pip install git+https://github.com/fernandomv301/AromaTools.git
```

**Cloning the repository.** Download the code and install it from the local folder:

```bash
git clone https://github.com/fernandomv301/AromaTools.git
cd AromaTools
pip install .
```

To modify the code and use the changes without reinstalling, use an editable installation (`pip install -e .`) instead of `pip install .`.

To update an installation made from Git:

```bash
pip install --upgrade git+https://github.com/fernandomv301/AromaTools.git
```

If you cloned the repository, run `git pull` in the folder and then `pip install --upgrade .`.

## Check the installation

The package can be imported from Python:

```bash
python -c "import aromatools"
```

and the three modules are installed as commands: `aromagnetic` (also available as `amg`), `arogeometric`, and `aroelectronics` (also available as `ael`). `amg` and `ael` with no arguments list their commands; `arogeometric` asks for the output file of the calculation; and `ael run` runs AroElectronics with the input file `INPUT.txt` of the current folder (see the pages of each module).

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
- esipython 1.1.0 or newer (AroElectronics). It is imported as `esipy`, and installs [PySCF](https://pyscf.org/) and [IOData](https://github.com/theochem/iodata) as dependencies.

---
