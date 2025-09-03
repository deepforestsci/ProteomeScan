import json
from rdkit import Chem
from rdkit.Chem import AllChem
import os

with open("ligand_SMILES.json", "r") as f:
    data = json.load(f)

# make ligands SDF
os.makedirs('processed', exist_ok=True)
for ligand, smiles in data['ligand_smiles'].items():
    if not smiles:
        continue
    writer = Chem.SDWriter(f'processed/processed_{ligand}.sdf')
    mol = Chem.MolFromSmiles(SMILES=smiles)
    # Add hydrogens to the molecule
    mol_with_H = Chem.AddHs(mol, addCoords=True)
    
    # Optionally, generate 3D coordinates
    AllChem.EmbedMolecule(mol_with_H, AllChem.ETKDG())
    
    # Write the protonated molecule to the output file
    writer.write(mol_with_H)

    writer.close()

# make methylated ligands SDF
os.makedirs('methylated', exist_ok=True)
for ligand, smiles in data['methylated_ligand_smiles'].items():
    if not smiles:
        continue
    writer = Chem.SDWriter(f'methylated/methylated_{ligand}.sdf')
    mol = Chem.MolFromSmiles(SMILES=smiles)
    # Add hydrogens to the molecule
    mol_with_H = Chem.AddHs(mol, addCoords=True)
    
    # Optionally, generate 3D coordinates
    AllChem.EmbedMolecule(mol_with_H, AllChem.ETKDG())
    
    # Write the protonated molecule to the output file
    writer.write(mol_with_H)

    writer.close()


