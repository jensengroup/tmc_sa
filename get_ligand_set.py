import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import GetPeriodicTable, rdchem, rdEHTTools, rdmolops

# import MoleculeAutoCorrect as mac


params = Chem.MolStandardize.rdMolStandardize.MetalDisconnectorOptions()
params.splitAromaticC = True
params.splitGrignards = True
params.adjustCharges = False
MetalNon_Hg = "[#3,#11,#12,#19,#13,#21,#22,#23,#24,#25,#26,#27,#28,#29,#30,#39,#40,#41,#42,#43,#44,#45,#46,#47,#48,#57,#72,#73,#74,#75,#76,#77,#78,#79,#80]~[#1,B,#6,#14,#15,#33,#51,#16,#34,#52,Cl,Br,I,#85]"
mdis = rdMolStandardize.MetalDisconnector(params)
mdis.SetMetalNon(Chem.MolFromSmarts(MetalNon_Hg))


# fmt: off
TRANSITION_METALS_NUM = [21,22,23,24,25,26,27,57,28,29,30,39,40,41,
                         42,43,44,45,46,47,48,71,72,73,74,75,76,77,78,79,80,
]
# fmt: on
#
def disconnect_ligands(smiles):

    mol = Chem.MolFromSmiles(smiles)

    frags = mdis.Disconnect(mol)
    frag_mols = list(rdmolops.GetMolFrags(frags, asMols=True))

    for i, f in enumerate(frag_mols):
        atoms = f.GetAtoms()
        for atom in atoms:
            if atom.GetAtomicNum() in TRANSITION_METALS_NUM:
                frag_mols.pop(i)
                break

    to_write = []
    for f in frag_mols:
        smi = Chem.MolToSmiles(f)
        if smi not in seen:
            seen.add(smi)
            to_write.append(smi)

    with open("ligand_set_csd_both_agree.smi", "a") as f:  # Open the file in write mode
        for item in to_write:
            f.write(f"{item}\n")  # Write each item followed by a newline


if __name__ == "__main__":
    # disconnect_ligands("CCCN(C)[Mo](<-[C]1N(CC)C=CN1CC)(N(C)CCC)N(C)CCC")

    global seen
    seen = set()

    smiles = pd.read_csv(
        "/home/magstr/git/xyz2mol_tm_jensengroup/SMILES_csvs/csd_smiles_both_agree.smi",
        header=None,
    )[0]

    for i, smi in enumerate(smiles):

        if i % 10000 == 0:
            print(f"Processed {i}")
        disconnect_ligands(smi)
