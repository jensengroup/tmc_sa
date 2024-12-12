import os
import argparse
import rdkit
import json
from rdkit import Chem
import subprocess
import re
import logging
from database import DatabaseManager

logger = logging.getLogger(__name__)


def Familiarity2(n_foreign_atoms, n_foreign_bonds, n_foreign_environments):

    n_foreign_keys = n_foreign_atoms + n_foreign_bonds + n_foreign_environments

    familiarity = 1 / (n_foreign_keys + 1)

    return familiarity


def Familiarity1(m, n_foreign_atoms, n_foreign_bonds, n_foreign_environments):

    n_foreign_keys = n_foreign_atoms + n_foreign_bonds + n_foreign_environments

    n_keys = m.GetNumAtoms() * 2 + m.GetNumBonds()

    familiarity = (n_keys - n_foreign_keys) / n_keys

    return familiarity


def run_shell(command, current_env):
    logger.debug(f"Running subprocess with command: {command}")
    result = subprocess.run(
        command,
        env=current_env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    output = result.stdout
    logger.debug(f"Raw subprocess output: \n{output}")

    if result.stderr:
        logger.error(f"Error output: {result.stderr}")

    return output


def parse_subprocess_output(command, current_env):

    output = run_shell(command, current_env)

    parsed_data = {
        "input_molecule": command[2],
        "foreign_atom_keys": 0,
        "atoms": [],
        "foreign_bond_keys": 0,
        "bonds": [],
        "foreign_atomic_environments": [],
    }

    # Parse the input molecule
    match = re.search(r"Input molecule:\s*(.+)", output)
    if match:
        parsed_data["input_molecule"] = match.group(1)

    # Parse foreign atom keys
    match = re.search(r"Molecule has (\d+) foreign atom keys", output)
    if match:
        parsed_data["foreign_atom_keys"] = int(match.group(1))

    # Parse atom details
    atom_matches = re.finditer(r"Atom (\d+):\s*(\S+)", output)
    for atom_match in atom_matches:
        parsed_data["atoms"].append(
            {
                "atom_id": int(atom_match.group(1)),
                "descriptor": atom_match.group(2),
            }
        )
    # Parse foreign bond keys
    match = re.search(r"Molecule has (\d+) foreign bond keys", output)
    if match:
        parsed_data["foreign_bond_keys"] = int(match.group(1))

    bond_matches = re.finditer(r"(Bond \d+):\s*(.*)", output)
    for bond_match in bond_matches:
        parsed_data["bonds"].append(
            {
                "bond_id": bond_match.group(
                    1
                ),  # Keep the bond ID as a string, e.g., "Bond 12"
                "details": bond_match.group(
                    2
                ).strip(),  # Store everything after the colon as a string
            }
        )

    # Parse foreign atomic environments
    env_matches = re.finditer(r"Environment of atom (\d+):\s*(.+)", output)
    for env_match in env_matches:
        parsed_data["foreign_atomic_environments"].append(
            {
                "atom_id": int(env_match.group(1)),
                "environment": env_match.group(2),
            }
        )

    logger.debug(f"Parsed subprocess output: {parsed_data}")
    return parsed_data


def ParseArgs(arg_list=None):
    parser = argparse.ArgumentParser(
        description="Run GA algorithm", fromfile_prefix_chars="+"
    )
    parser.add_argument("--smiles", type=str, default=None)
    parser.add_argument(
        "--log_level",
        type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Set the logging level",
    )
    parser.add_argument(
        "--db-path",
        type=str,
        default="familiarity_scores.db",
        help="Path to the SQLite database file",
    )
    return parser.parse_args(arg_list)


def get_familiarity(smiles):

    # Get the current environment variables
    current_env = os.environ.copy()

    # Add/Modify the environment variable
    # current_env["MOLECULE_AUTO_CORRECT"] = "/home/magstr/git/MoleculeAutoCorrect"
    # current_env["MOLPERT"] = "/home/magstr/git/Molpert"

    # Define the command
    command = [
        f"{current_env['MOLECULE_AUTO_CORRECT']}/bin/HighlightMoleculeErrors",
        "/home/magstr/git/MoleculeAutoCorrect/dicts/csd_smiles_both_agree.dict",
        smiles,
        "/tmp/image.svg",
    ]
    logger.info(f"Executing command for SMILES: {smiles}")

    # Run the subprocess
    output = parse_subprocess_output(command, current_env)
    # logger.debug(pformat(output))
    logger.debug(json.dumps(output, indent=4))

    logger.info("Getting familiarity")

    n_foreign_atoms = int(output["foreign_atom_keys"])
    n_foreign_bonds = int(output["foreign_bond_keys"])
    n_foreign_environments = int(len(output["foreign_atomic_environments"]))

    m = Chem.MolFromSmiles(smiles)

    fam1 = Familiarity1(m, n_foreign_atoms, n_foreign_bonds, n_foreign_environments)
    fam2 = Familiarity2(n_foreign_atoms, n_foreign_bonds, n_foreign_environments)

    logger.info(f"Familiarity1: {fam1} | Familiarity2: {fam2}")

    return fam1, fam2


if __name__ == "__main__":
    args = ParseArgs()
    # Initialize logging
    logging.basicConfig(format="%(asctime)s - %(levelname)s - %(message)s")
    logger.setLevel(getattr(logging, args.log_level))

    # Initialize database manager
    db_manager = DatabaseManager(args.db_path)

    # Check if the familiarity score already exists
    existing_score = db_manager.get_existing_familiarity_scores(args.smiles)
    if existing_score is not None:
        logger.info(f"Familiarity score for SMILES already exists: {existing_score}")

    else:
        # Compute or retrieve familiarity score
        fam1, fam2 = get_familiarity(args.smiles)
        logger.info(f"Familiarity1: {fam1} | Familiarity2: {fam2}")

        # Store the familiarity score in the database
        db_manager.store_familiarity_scores(args.smiles, fam1, fam2)
