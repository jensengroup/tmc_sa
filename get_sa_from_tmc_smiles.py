import os
import argparse
import rdkit
import json
from rdkit import Chem
import subprocess
import re


def Familiarity2(n_foreign_atoms, n_foreign_bonds, n_foreign_environments):

    n_foreign_keys = n_foreign_atoms + n_foreign_bonds + n_foreign_environments

    familiarity = 1 / (n_foreign_keys + 1)

    return familiarity


def Familiarity1(m, n_foreign_atoms, n_foreign_bonds, n_foreign_environments):

    n_foreign_keys = n_foreign_atoms + n_foreign_bonds + n_foreign_environments

    n_keys = m.GetNumAtoms() * 2 + m.GetNumBonds()

    familiarity = (n_keys - n_foreign_keys) / n_keys

    return familiarity


def parse_subprocess_output(command, current_env):

    result = subprocess.run(
        command,
        env=current_env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    output = result.stdout

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

    return parsed_data


# Example usage
# Replace ['your-command'] with the actual command to run the subprocess
# parsed_output = parse_subprocess_output(['your-command'])
# print(parsed_output)


def ParseArgs(arg_list=None):
    parser = argparse.ArgumentParser(
        description="Run GA algorithm", fromfile_prefix_chars="+"
    )
    parser.add_argument("--smiles", type=str, default=None)
    return parser.parse_args(arg_list)


def get_arguments(arg_list=None) -> argparse.Namespace:
    """

    Args:
        arg_list: Automatically obtained from the commandline if provided.
        Otherwise default arguments are used

    Returns:
        parser.parse_args(arg_list)(Namespace): Dictionary like class that contain the driver arguments.

    """
    parser = argparse.ArgumentParser(
        description="Run GA algorithm", fromfile_prefix_chars="+"
    )
    parser.add_argument(name="smiles", type=str, default=None)
    return parser.parse_args(arg_list)


def get_familiarity(smiles):

    # Get the current environment variables
    current_env = os.environ.copy()

    # Add/Modify the environment variable
    current_env["MOLECULE_AUTO_CORRECT"] = "/home/magstr/git/MoleculeAutoCorrect"
    current_env["MOLPERT"] = "/home/magstr/git/Molpert"
    print(current_env)

    # Define the command
    command = [
        f"{current_env['MOLECULE_AUTO_CORRECT']}/bin/HighlightMoleculeErrors",
        "/home/magstr/git/MoleculeAutoCorrect/dicts/csd_smiles_both_agree.dict",
        smiles,
        "/tmp/image.svg",
    ]
    print(command)

    # Run the subprocess
    output = parse_subprocess_output(command, current_env)

    print("Getting familiarity")

    # m = Chem.MolFromSmiles(smiles)

    n_foreign_atoms = int(output["foreign_atom_keys"])
    n_foreign_bonds = int(output["foreign_bond_keys"])
    n_foreign_environments = int(len(output["foreign_atomic_environments"]))

    # fam1 = Familiarity1(m, n_foreign_atoms, n_foreign_bonds, n_foreign_environments)
    fam2 = Familiarity2(n_foreign_atoms, n_foreign_bonds, n_foreign_environments)

    print(f"Familiarity1: {fam1}\nFamiliarity2: {fam2}")

    return fam2


if __name__ == "__main__":
    args = ParseArgs()
    # Get the current environment variables
    current_env = os.environ.copy()

    # Define the command
    command = [
        f"{current_env['MOLECULE_AUTO_CORRECT']}/bin/HighlightMoleculeErrors",
        "/home/magstr/git/MoleculeAutoCorrect/dicts/csd_smiles_both_agree.dict",
        args.smiles,
        "molecule_errors.svg",
    ]

    # Run the subprocess
    output = parse_subprocess_output(command, current_env)

    print(json.dumps(output, indent=4))

    print("Getting familiarity")

    m = Chem.MolFromSmiles(
        "Cc1cc2ccccc2c2-c3c([O-]->[Ir+2]4(<-[O-]c12)<-N(C)(CCN->4(C)Cc1ccccc1)Cc1ccccc1)c(C)cc1ccccc31"
    )

    n_foreign_atoms = int(output["foreign_atom_keys"])
    n_foreign_bonds = int(output["foreign_bond_keys"])
    n_foreign_environments = int(len(output["foreign_atomic_environments"]))

    fam1 = Familiarity1(m, n_foreign_atoms, n_foreign_bonds, n_foreign_environments)
    fam2 = Familiarity2(n_foreign_atoms, n_foreign_bonds, n_foreign_environments)

    print(f"Familiarity1: {fam1}\nFamiliarity2: {fam2}")
