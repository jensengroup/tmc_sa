import argparse
import logging
import os
import time

import pandas as pd
from rdkit import Chem

import MoleculeAutoCorrect as mac
from get_sa_from_tmc_smiles import run_shell

logger = logging.getLogger(__name__)


def ParseArgs(arg_list=None):
    parser = argparse.ArgumentParser(
        description="Create chemical dictionary", fromfile_prefix_chars="+"
    )
    parser.add_argument(
        "--smiles_data",
        type=str,
        default=None,
        required=True,
        help="Path to .smi file",
    )
    parser.add_argument(
        "--dict_name",
        type=str,
        default=None,
        required=True,
        help="Output dict name",
    )
    parser.add_argument(
        "--environment_radius",
        type=int,
        help="Radius of the environments used in the fragmentation of the TMC atomic environments",
        default=1,
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Set debug mode",
    )
    return parser.parse_args(arg_list)


def get_dict():
    # Get the current environment variables
    current_env = os.environ.copy()

    # # Add/Modify the environment variable
    # current_env["MOLECULE_AUTO_CORRECT"] = "/home/magstr/git/MoleculeAutoCorrect"
    # current_env["MOLPERT"] = "/home/magstr/git/Molpert"

    # Define the command
    command = [
        f"{current_env['MOLECULE_AUTO_CORRECT']}/bin/MakeChemicalDictionary",
        f"{args.smiles_data}",
        f"{args.dict_name}",
        f"{args.environment_radius}",  # Atomic enironment radius
    ]

    # Run the subprocess
    output = run_shell(command, current_env)

    if "Error" in output:
        logger.error("Could not create dictionary")
        return
    logger.info(
        "Succesfully created single entry dict. Adding remaning entries in Python."
    )

    logger.info(f"Loading the single entry dict: {args.dict_name}")
    dictionary = mac.ChemicalDictionary(args.dict_name)

    # Get mols from smiles data
    smiles = pd.read_csv(args.smiles_data, header=None)[0]
    logger.info(f"Loaded data: \n {smiles.head(5).to_string()}")

    if args.debug:
        logger.info("Running in debug mode. Loading up to the first 200 SMILES")
        smiles = smiles[0:200]
        logger.info(f"{len(smiles)} SMILES loaded")

    mols = [Chem.MolFromSmiles(smi) for smi in smiles if smi]

    logger.info(f"Adding {len(mols)} mols to the dictionary")

    start = time.time()
    for i, mol in enumerate(mols[1::]):  # We dont have to add the first entry again
        if i % 1000 == 1:
            logger.info(f"Processed {i} of {len(mols)} mols")
        if not mol:
            continue
        dictionary.AddMolecule(mol)
    logger.info(f"Processed all {len(mols)}")
    end = time.time()

    dictionary.number_of_molecules = len(mols)
    dictionary.dataset = args.smiles_data
    dictionary.Save(args.dict_name)

    logger.info(f"Created dictionary in {end-start} seconds")

    return


if __name__ == "__main__":
    args = ParseArgs()
    level = "DEBUG" if args.debug else "INFO"
    # Initialize logging
    logging.basicConfig(format="%(asctime)s - %(levelname)s - %(message)s", level=level)
    get_dict()
