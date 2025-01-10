# TMC_SA

SA like score for TMCs.

A virtual library of reference correct TMCs is used to build a dictionary of allowed chemical features. The chemical features of input molecules are compared against this dictionary to yield a familiarity score of a given TMC.

Based on work by:
[Kerstjens, A., De Winter, H. Molecule auto-correction to facilitate molecular design. J Comput Aided Mol Des 38, 10 (2024).](https://doi.org/10.1007/s10822-024-00549-1)
[MoleculeAutoCorrect](https://github.com/AlanKerstjens/MoleculeAutoCorrect)

# Installation

The following instructions are for GNU+Linux. For alternative operating systems you'll have to adapt these commands slightly.

*IMPORTANT:*
The MoleculeAutoCorrect and Molpert repos depend on the source C++ code from RDKit. In newer versions of RDKit, certain header files have been removed from the conda install.
Therefore, these repos need to be compiled in a conda environment with RDKit version 2022.09.5.

Once the binaries from the MoleculeAutoCorrect repo has been compiled, these binaries can be called from any conda env. E.g a conda env that uses an updated version of RDKit.

To get started, first install the conda env: [env.yml](env.yml)

```shell
conda env create --file ./env.yml
```

Then install MoleculeAutoCorrect and Molpert by running the following:

```shell
./install.sh
```

To be able to import the library from Python add `${MOLECULE_AUTO_CORRECT}/lib` and `${MOLPERT}/lib` to your `${PYTHONPATH}`. Consider doing so in your `.bash_profile` file. Otherwise you'll have to manually extend `${PYTHONPATH}` everytime you open a new shell.

```shell
export PYTHONPATH="${PYTHONPATH}:${MOLECULE_AUTO_CORRECT}/lib"
export PYTHONPATH="${PYTHONPATH}:${MOLPERT}/lib"
```

Now you should be able to run the commands given in [Quick start](#quick-start).

## Quick start

We provide python wrapper functions that call the compiled binaries from MoleculeAutoCorrect.

Get your hands on a virtual library of molecules you would like to use as reference of correct chemistry (here `tmc.smi`). Then use this library to create a dictionary of chemical features (here `tmc.dict`). You can specify the radius of circular atomic environments using the --environment_radius argument (here `1`).

Creating the tmc.dict by calling the MoleculeAutoCorrect binaries with Python using the scripts highlighted below:
NB! When creating the dict, it is important that you use the conda env installed above ([env.yml](env.yml)). Otherwise you will get an error. After the .dict has been created you can switch to an environment with a newer version of RDKit.

```shell
python ./create_chemical_dictionary_from_smiles.py --smiles_data ./data/tmc.smi --dict_name ./dicts/tmc.dict --environment_radius 1
```

Use the tmc.dict to get familiarity scores for a given TMC SMILES:

```shell
python ./get_sa_from_tmc_smiles.py --smiles "CN(C)C=O->[Ni+2]123<-[N-](C(=O)CN->1(CC(=O)[N-]->2c1ccccc1)CC(=O)[N-]->3c1ccccc1)c1ccccc1" --reference_dict ./dicts/tmc.dict
```

You can also inspect the output of HighlightMoleculeErrors directly running the binary:

```shell
bin/HighlightMoleculeErrors ./dicts/tmc.dict "CCCN(C)[Mo](<-[C]1N(CC)C=CN1CC)(N(C)CCC)N(C)CCC" molecule_errors.svg
```

`get_sa_from_tmc_smiles` contains the `get_familiarity` function which returns the calculated familiarities. This function can be imported in other scripts and then used
as an SA score calculator.

<!-- If it has issues you can proceed to try correcting them: -->

<!---->

<!-- ```shell -->

<!-- python AutoCorrectMolecule.py tmc.dict "CCCN(C)[Mo](<-[C]1N(CC)C=CN1CC)(N(C)CCC)N(C)CCC" -->

<!-- ``` -->

<!---->

<!-- You can experiment with different settings, including tree policies. Access the `--help` for more information. -->

<!---->

<!-- ```shell -->

<!-- python AutoCorrectMolecule.py --help -->

<!-- ``` -->
