import os
import argparse
import rdkit
from pathlib import Path
import json
from rdkit import Chem
import subprocess
import re
import logging
from database import DatabaseManager

logger = logging.getLogger(__name__)

TRANSITION_METALS_NUM = [
    21,
    22,
    23,
    24,
    25,
    26,
    27,
    57,
    28,
    29,
    30,
    39,
    40,
    41,
    42,
    43,
    44,
    45,
    46,
    47,
    48,
    71,
    72,
    73,
    74,
    75,
    76,
    77,
    78,
    79,
    80,
]

TRANSITION_METALS = [
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "La",
    "Ni",
    "Cu",
    "Zn",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
]
tm = f"[{','.join(TRANSITION_METALS)}]"


def Familiarity2(n_foreign_atoms, n_foreign_bonds, n_foreign_environments):

    n_foreign_keys = n_foreign_atoms + n_foreign_bonds + n_foreign_environments

    familiarity = 1 / (n_foreign_keys + 1)

    return familiarity


def Familiarity1(m, n_foreign_atoms, n_foreign_bonds, n_foreign_environments):

    n_foreign_keys = n_foreign_atoms + n_foreign_bonds + n_foreign_environments

    n_keys = m.GetNumAtoms() * 2 + m.GetNumBonds()

    familiarity = (n_keys - n_foreign_keys) / n_keys

    return familiarity


def Familiarity1_bonds(m, n_foreign_tm_bonds):

    neighbors = m.GetAtomWithIdx(
        m.GetSubstructMatch(Chem.MolFromSmarts(tm))[0]
    ).GetNeighbors()

    n_keys = len(neighbors)

    familiarity = (n_keys - n_foreign_tm_bonds) / n_keys

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
        "foreign_tm_bonds": 0,
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

    pattern = re.compile(
        r"\((-?\d+,-?\d+,-?\d+,-?\d+,-?\d+)\)->\((-?\d+,-?\d+,-?\d+,-?\d+,-?\d+)\)"
    )
    count = 0
    for match in re.finditer(pattern, output):
        left_tuple_str = match.group(1)
        right_tuple_str = match.group(2)

        # Convert to integers
        left_tuple = list(map(int, left_tuple_str.split(",")))
        right_tuple = list(map(int, right_tuple_str.split(",")))

        # The "third position" means index 2 (0-based indexing)
        left_third = left_tuple[2]
        right_third = right_tuple[2]

        left_is_transition = left_third in TRANSITION_METALS_NUM
        right_is_transition = right_third in TRANSITION_METALS_NUM
        if left_is_transition or right_is_transition:
            count += 1
    parsed_data["foreign_tm_bonds"] = count

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
    parser.add_argument(
        "--score_type",
        type=str,
        default="Familiarity1",
        choices=["Familiarity1", "Familiarity1_bonds", "Familiarity2"],
        help="Path to the SQLite database file",
    )
    return parser.parse_args(arg_list)


def get_familiarity(smiles, score_type="Familiarity1", reference_dict=None):

    # Get the current environment variables
    current_env = os.environ.copy()

    # Add/Modify the environment variable
    # current_env["MOLECULE_AUTO_CORRECT"] = "/home/magstr/git/MoleculeAutoCorrect"
    # current_env["MOLPERT"] = "/home/magstr/git/Molpert"

    # Define the command
    command = [
        f"{current_env['MOLECULE_AUTO_CORRECT']}/bin/HighlightMoleculeErrors",
        reference_dict,
        smiles,
        "/tmp/image.svg",
    ]
    logger.debug(f"Executing command for SMILES: {smiles}")

    # Run the subprocess
    output = parse_subprocess_output(command, current_env)
    # logger.debug(pformat(output))
    logger.debug(json.dumps(output, indent=4))

    logger.debug("Getting familiarity")

    n_foreign_atoms = int(output["foreign_atom_keys"])
    n_foreign_bonds = int(output["foreign_bond_keys"])
    n_foreign_environments = int(len(output["foreign_atomic_environments"]))
    n_foreign_tm_bonds = int(output["foreign_tm_bonds"])

    m = Chem.MolFromSmiles(smiles, sanitize=False)

    fam1 = Familiarity1(m, n_foreign_atoms, n_foreign_bonds, n_foreign_environments)
    fam2 = Familiarity2(n_foreign_atoms, n_foreign_bonds, n_foreign_environments)
    fam3 = Familiarity1_bonds(m, n_foreign_tm_bonds)

    logger.debug(
        f"Familiarity1: {fam1} | Familiarity2: {fam2} | Familiarity_bonds: {fam3}"
    )

    return fam1, fam2, fam3


if __name__ == "__main__":
    args = ParseArgs()
    # Initialize logging
    logging.basicConfig(format="%(asctime)s - %(levelname)s - %(message)s")
    logger.setLevel(getattr(logging, args.log_level))

    # Initialize database manager
    db_manager = DatabaseManager(args.db_path)

    # Reference dictionary
    reference_dict = Path("./dicts/csd_smiles_both_agree_train.dict").resolve()
    # Check if the familiarity score already exists
    existing_score = db_manager.get_existing_familiarity_scores(args.smiles)
    if existing_score is not None:
        logger.info(f"Familiarity score for SMILES already exists: {existing_score}")

    else:
        # Compute or retrieve familiarity score
        fam1, fam2, fam3 = get_familiarity(
            args.smiles, reference_dict=str(reference_dict)
        )
        logger.debug(
            f"Familiarity1: {fam1} | Familiarity2: {fam2} | Familiarity_bonds: {fam3}"
        )

        # Store the familiarity score in the database
        db_manager.store_familiarity_scores(args.smiles, fam1, fam2)
