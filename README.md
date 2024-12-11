# About

SA like score for TMCs.

A virtual library of reference correct TMCs is used to build a dictionary of allowed chemical features. The chemical features of input molecules are compared against this dictionary to yield a familiarity score of a given TMC.

Based on work by:
[Kerstjens, A., De Winter, H. Molecule auto-correction to facilitate molecular design. J Comput Aided Mol Des 38, 10 (2024).](https://doi.org/10.1007/s10822-024-00549-1)

# Installation

### Prerequisites

Ensure the following dependencies are installed:

- [RDKit](https://rdkit.org/)
- [Molpert](https://github.com/AlanKerstjens/Molpert)
- [Boost](https://www.boost.org/). You already have this if you installed the RDKit. If you'd like to build the Python bindings make sure Boost.Python is installed.
- [CMake](https://cmake.org/)

### Instructions

The following instructions are for GNU+Linux. For alternative operating systems you'll have to adapt these commands slightly.

First install the [env.yml](env.yml) file into a conda env.
Then install MoleculeAutoCorrect and Molpert.

Then run the install script:

```shell
./install.sh
```

To be able to import the library from Python add `${MOLECULE_AUTO_CORRECT}/lib` and `${MOLPERT}/lib` to your `${PYTHONPATH}`. Consider doing so in your `bash_profile` file. Otherwise you'll have to manually extend `${PYTHONPATH}` everytime you open a new shell.

```shell
export PYTHONPATH="${PYTHONPATH}:${MOLECULE_AUTO_CORRECT}/lib"
export PYTHONPATH="${PYTHONPATH}:${MOLPERT}/lib"
```

### Troubleshooting

CMake will try to find the rest of the dependencies for you. To avoid problems ensure you build the software with the same Boost and Python versions that you used to build Molpert and the RDKit. If CMake finds a different Boost or Python installation you'll need to point it to the correct one, as described [here](https://cmake.org/cmake/help/latest/module/FindBoost.html) and [here](https://cmake.org/cmake/help/latest/module/FindPython.html).

CMake will search for the RDKit in the active Anaconda environment (if you have one) and at `${RDBASE}` if set. If neither of these are the case you need to specify the path to the RDKit yourself. Replace the above CMake command with the one below, substituting the `<placeholder/path>` with your paths.

```shell
cmake -DRDKit_ROOT=<path/to/rdkit> -DMolpert_INCLUDE_DIRS=<path/to/molpert> ..
```

# Quick start

Get your hands on a virtual library of molecules you would like to use as reference of correct chemistry (here `tmc.smi`). Then use this library to create a dictionary of chemical features (here `tmc.dict`). You can specify the radius of circular atomic environments as the last argument (here `1`).

```shell
bin/MakeChemicalDictionary tmc.smi tmc.dict 1
```

Given the SMILES string of a TMC `CCCN(C)[Mo](<-[C]1N(CC)C=CN1CC)(N(C)CCC)N(C)CCC` you can inspect if it has any issues:

```shell
bin/HighlightMoleculeErrors tmc.dict "CCCN(C)[Mo](<-[C]1N(CC)C=CN1CC)(N(C)CCC)N(C)CCC" molecule_errors.svg
```

If it has issues you can proceed to try correcting them:

```shell
python AutoCorrectMolecule.py tmc.dict "CCCN(C)[Mo](<-[C]1N(CC)C=CN1CC)(N(C)CCC)N(C)CCC"
```

You can experiment with different settings, including tree policies. Access the `--help` for more information.

```shell
python AutoCorrectMolecule.py --help
```
