#!/usr/bin/env python3
import sys

# Cheminformatics stack
from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import AllChem

smiles = sys.argv[1]
name = sys.argv[2]


def smi23dxyz(smiles, name):
    ps = rdDistGeom.srETKDGv3()
    ps.useRandomCoords = False
    try:
        mol = Chem.MolFromSmiles(smiles,False)
        mol.UpdatePropertyCache(strict=False)
        mol_2 = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_2, ps)
        Chem.MolToXYZFile(mol_2, f"./XYZ_Structures_Generated/{name}.xyz")
    except Exception as m:
        raise Exception(m)


smi23dxyz(smiles, name)
