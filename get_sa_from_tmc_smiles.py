import argparse
import json
import logging
import os
import re
import subprocess
from pathlib import Path

from rdkit import Chem

logger = logging.getLogger(__name__)

# fmt: off
TRANSITION_METALS = ["Sc","Ti","V","Cr","Mn","Fe","Co","La","Ni","Cu","Zn",
                     "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","Lu",
                     "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
]

TRANSITION_METALS_NUM = [21,22,23,24,25,26,27,57,28,29,30,39,40,41,
                         42,43,44,45,46,47,48,71,72,73,74,75,76,77,78,79,80,
]
# fmt: on

# SMARTS string for TRANSITION_METALS
tm = f"[{','.join(TRANSITION_METALS)}]"


def familiarity2(n_foreign_atoms, n_foreign_bonds, n_foreign_environments):
    n_foreign_keys = n_foreign_atoms + n_foreign_bonds + n_foreign_environments

    familiarity = 1 / (n_foreign_keys + 1)

    return familiarity


def familiarity1(m, n_foreign_atoms, n_foreign_bonds, n_foreign_environments):
    n_foreign_keys = n_foreign_atoms + n_foreign_bonds + n_foreign_environments

    n_keys = m.GetNumAtoms() * 2 + m.GetNumBonds()

    familiarity = (n_keys - n_foreign_keys) / n_keys

    return familiarity


def familiarity_tm_bonds(m, n_foreign_tm_bonds):
    match = m.GetSubstructMatch(Chem.MolFromSmarts(tm))
    if not match:
        logger.warning(f"No TM metal for this SMILES: {Chem.MolToSmiles(m)}")
        return None
    neighbors = m.GetAtomWithIdx(
        m.GetSubstructMatch(Chem.MolFromSmarts(tm))[0]
    ).GetNeighbors()

    n_keys = len(neighbors)

    familiarity = (n_keys - n_foreign_tm_bonds) / n_keys

    return familiarity


def familiarity4(n_foreign_atoms, n_foreign_bonds, n_foreign_environments):
    familiarity = (
        1
        - 0.2 * n_foreign_atoms
        - 0.10 * n_foreign_bonds
        - 0.05 * n_foreign_environments
    )
    if familiarity < 0:
        familiarity = 0
    return familiarity


def familiarity5(
    n_foreign_atoms, n_foreign_bonds, n_foreign_environments, n_foreign_pd_envs
):
    familiarity = (
        1
        - 0.2 * n_foreign_atoms
        - 0.10 * n_foreign_bonds
        - 0.05 * (n_foreign_environments - n_foreign_pd_envs)
    )

    if familiarity < 0:
        familiarity = 0
    return familiarity


def familiarity6(
    n_foreign_atoms, n_foreign_bonds, n_foreign_environments, n_foreign_pd_envs
):
    familiarity = (
        1
        - 0.1 * n_foreign_atoms
        - 0.05 * n_foreign_bonds
        - 0.01 * (n_foreign_environments - n_foreign_pd_envs)
    )

    if familiarity < 0:
        familiarity = 0

    return familiarity


def run_shell(command, current_env):
    logger.debug(f"Running subprocess with command: {command}")
    result = subprocess.run(
        command,
        env=current_env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=10,
    )

    output = result.stdout
    logger.debug(f"Raw subprocess output: \n{output}")

    if result.stderr:
        logger.error(f"Error output: {result.stderr}")

    return output


def parse_subprocess_output(command, current_env):
    "Function parsing the output returned by the compiled binary"

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

    # Parse the input molecule string
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
                "bond_id": bond_match.group(1),
                "details": bond_match.group(2).strip(),
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

    # Match the bond keys and use to detect TM bonds.
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

    return parsed_data


def ParseArgs(arg_list=None):
    parser = argparse.ArgumentParser(
        description="Get familiarity scores from TMC SMILES", fromfile_prefix_chars="+"
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
        help="Path to the SQLite database file used to store familiarity scores",
    )
    parser.add_argument(
        "--reference_dict",
        type=Path,
        default="./dicts/csd_smiles_both_agree_train.dict",
        help="Path to the TMC reference dict",
    )

    return parser.parse_args(arg_list)


def get_scores_from_output(smiles, output):
    n_foreign_atoms = int(output["foreign_atom_keys"])
    n_foreign_bonds = int(output["foreign_bond_keys"])
    n_foreign_environments = int(len(output["foreign_atomic_environments"]))
    n_foreign_tm_bonds = int(output["foreign_tm_bonds"])
    foreign_envs = output["foreign_atomic_environments"]

    # Get measure for envs with Pd.
    n_foreign_pd_envs = 0
    for d in foreign_envs:
        env = d["environment"]
        if "Pd" in env:
            n_foreign_pd_envs += 1

    m = Chem.MolFromSmiles(smiles, sanitize=False)

    fam1 = familiarity1(m, n_foreign_atoms, n_foreign_bonds, n_foreign_environments)
    fam2 = familiarity2(n_foreign_atoms, n_foreign_bonds, n_foreign_environments)
    fam3 = familiarity_tm_bonds(m, n_foreign_tm_bonds)
    fam4 = familiarity4(n_foreign_atoms, n_foreign_bonds, n_foreign_environments)
    fam5 = familiarity5(
        n_foreign_atoms, n_foreign_bonds, n_foreign_environments, n_foreign_pd_envs
    )
    fam6 = familiarity5(
        n_foreign_atoms, n_foreign_bonds, n_foreign_environments, n_foreign_pd_envs
    )

    logger.debug(
        f"familiarity1: {fam1} | familiarity2: {fam2} | familiarity_tm_bonds: {fam3} | familiarity 4: {fam4} | familiarity5: {fam5}\
        | familiarity6: {fam6}"
    )

    output["familiarity1"] = fam1
    output["familiarity2"] = fam2
    output["familiarity3"] = fam3
    output["familiarity4"] = fam4
    output["familiarity5"] = fam5
    output["familiarity6"] = fam6

    return output


def get_familiarity(smiles, reference_dict=None):
    """Calculate familiarity scores for a given TMC SMILES.

    Args:
        smiles (str): Input TMC SMILES
        reference_dict (str): Path to the TMC reference dict

    Returns:
        output (dict): Dictionary with the calulated familiarities and the parsed output from the MOLECULE_AUTO_CORRECT binary
    """
    current_env = os.environ.copy()

    command = [
        f"{current_env['MOLECULE_AUTO_CORRECT']}/bin/HighlightMoleculeErrors",
        str(reference_dict),
        smiles,
        "/tmp/image.svg",  # This file will contain an image of the molecule with foreign parts highlighted
    ]
    logger.debug(f"Executing command for SMILES: {smiles}")

    # Run the subprocess
    output = parse_subprocess_output(command, current_env)
    logger.debug(json.dumps(output, indent=4))

    # Parse to output and get get familiarity scores.
    output = get_scores_from_output(smiles, output)

    return output


if __name__ == "__main__":
    args = ParseArgs()
    if args.log_level == "DISABLE":
        logging.disable(logging.CRITICAL)
    else:
        logging.basicConfig(format="%(asctime)s - %(levelname)s - %(message)s")
        logger.setLevel(getattr(logging, args.log_level))

    # Initialize database manager
    # db_manager = DatabaseManager(args.db_path)

    # Reference dictionary
    reference_dict = args.reference_dict.resolve()

    # Check if the familiarity score already exists
    # existing_score = db_manager.get_existing_familiarity_scores(args.smiles)
    # if existing_score is not None:
    #     logger.info(f"Familiarity score for SMILES already exists: {existing_score}")

    # Compute or retrieve familiarity score
    output = get_familiarity(args.smiles, reference_dict=str(reference_dict))
    logger.info("Printing the familiarity output:")
    logger.info(json.dumps(output, indent=4))

    # Store the familiarity score in the database
    # db_manager.store_familiarity_scores(args.smiles, fam1, fam2)
