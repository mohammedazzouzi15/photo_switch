#!/usr/env/ python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
import copy
import os

# Cheminformatics stack
import rdkit
from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import AllChem
from rdkit.Chem import FragmentMatcher
from rdkit.Chem import rdFMCS
from glob import glob

# For segfault
import faulthandler
import signal
from subprocess import Popen, PIPE

faulthandler.enable()


def savexyz(mol, name):
    ps = rdDistGeom.ETKDG()
    ps.useRandomCoords = False
    try:
        mol2 = Chem.AddHs(mol)
        # AllChem.EmbedMolecule(mol2, ps)
        # Chem.MolToXYZFile(mol2, f"./XYZ_Structures_Generated/{name}.xyz")
        return mol2
    except Exception as m:
        return None


of = open("Generated_Smiles.csv", "a")
if os.path.exists("finalsmiles.npy") and os.path.exists("finalmols.npy"):
    finalsmiles = np.load("finalsmiles.npy")
    finalmols = np.load("finalmols.npy", allow_pickle=True)
else:
    df = pd.read_csv("Base_Smiles.csv")

    names = df.Structure.values
    smiles = df.SMILES.values
    charges = df.Charge.values

    # Filter non-neutral molecules
    neutral_ids = np.where(charges == 0)

    names = names[neutral_ids]
    smiles = smiles[neutral_ids]
    charges = charges[neutral_ids]

    print(f"A total of {names.shape[0]} neutral smiles were retained!")

    # Filter molecules with weird SMILES
    mols = np.array([Chem.MolFromSmiles(smi) for smi in smiles], dtype=object)
    normal_ids = np.where(mols != None)

    names = names[normal_ids]
    smiles = smiles[normal_ids]
    charges = charges[normal_ids]
    mols = mols[normal_ids]

    print(f"A total of {names.shape[0]} sanitized smiles were retained!")

    # Canonicalization of SMILES to remove duplicates

    csmiles = np.array([Chem.MolToSmiles(mol) for mol in mols])
    _, idx = np.unique(csmiles, return_index=True)
    names = names[idx]
    smiles = smiles[idx]
    charges = charges[idx]
    mols = mols[idx]
    csmiles = csmiles[idx]

    print(f"A total of {names.shape[0]} unique smiles were retained!")

    # Classification rules

    patterns = [
        "NC(N)=O",
        "NC(N)=S",
        "O=C1C(N)=C(N)C1=O",
        "S=C1C(N)=C(N)C1=S",
        "O=C1C(N)=C1N",
        "O=C(C(N)=C(N)C1=O)C1=O",
        "O=S(N)(N)=O",
    ]
    moltype = np.zeros_like(names)

    # Method 1 using SubstructMatch
    mol_patterns = np.array([Chem.MolFromSmiles(smi) for smi in patterns], dtype=object)
    for i, mol in enumerate(mols):
        for j, mp in enumerate(mol_patterns):
            if mol.HasSubstructMatch(mp):
                # print(f"{csmiles[i]} seems to contain {patterns[j]}")
                moltype[i] = j + 1

    # Print sublists from Method 1
    type1 = np.where(moltype == 1)
    type2 = np.where(moltype == 2)
    type3 = np.where(moltype == 3)
    type4 = np.where(moltype == 4)
    type5 = np.where(moltype == 5)
    type6 = np.where(moltype == 6)
    type7 = np.where(moltype == 7)
    types = [type1, type2, type3, type4, type5, type6, type7]
    for i, type in enumerate(types):
        # print(f"\n\n {patterns[i]} list, with {np.count_nonzero(type)} elements:")
        for name, smi, charge in zip(names[type], csmiles[type], charges[type]):
            # print(f"{name},{smi},{charge}")
            pass

    # Method 2 using FragmentMatchers
    p1 = FragmentMatcher.FragmentMatcher()
    p2 = FragmentMatcher.FragmentMatcher()
    p3 = FragmentMatcher.FragmentMatcher()
    p4 = FragmentMatcher.FragmentMatcher()
    p5 = FragmentMatcher.FragmentMatcher()
    p6 = FragmentMatcher.FragmentMatcher()
    p7 = FragmentMatcher.FragmentMatcher()

    pmatchers = [p1, p2, p3, p4, p5, p6, p7]

    for p, pm in zip(patterns, pmatchers):
        pm.Init(p)
        for p2 in patterns:
            if p2 is not p:
                pm.AddExclusion(p2)

    # Generate unique fragment list
    fragsmiles = []
    fragmols = []
    fragaps = []
    fraparent = []
    fraparentmol = []
    for i, mol in enumerate(mols):
        for j, pm in enumerate(pmatchers):
            if pm.HasMatch(mol):
                # print(f"{csmiles[i]} seems to contain {patterns[j]}")
                moltype[i] = j + 1
                matches = pm.GetMatches(mol)
                bonds_to_cut = []
                for idx in matches[0]:
                    for n in mol.GetAtomWithIdx(idx).GetNeighbors():
                        nidx = n.GetIdx()
                        if nidx not in matches[0]:
                            bonds_to_cut.append(mol.GetBondBetweenAtoms(idx, nidx))
                if bonds_to_cut:
                    posi_ids = [b.GetBeginAtomIdx() for b in bonds_to_cut]
                    pose_ids = [b.GetEndAtomIdx() for b in bonds_to_cut]
                    bond_ids = [b.GetIdx() for b in bonds_to_cut]
                    temp_mol_f = Chem.FragmentOnBonds(mol, bond_ids, addDummies=False)
                    fmols = Chem.GetMolFrags(temp_mol_f, asMols=True)
                    fatoms = Chem.GetMolFrags(temp_mol_f)
                    for fmol, fatom in zip(fmols, fatoms):
                        fatom = list(fatom)
                        if not pm.HasMatch(fmol):
                            fragsmiles.append(Chem.MolToSmiles(fmol))
                            fragmols.append(fmol)
                            # print(posi_ids, pose_ids, fatom)
                            if any(a in posi_ids for a in fatom):
                                for pos in posi_ids:
                                    if pos in fatom:
                                        fatom_ap_idx = list(fatom).index(pos)
                            elif any(a in pose_ids for a in fatom):
                                for pos in pose_ids:
                                    if pos in fatom:
                                        fatom_ap_idx = list(fatom).index(pos)
                            fragaps.append(fatom_ap_idx)
                            fraparent.append(names[i])
                            fraparentmol.append(mol)
    fragsmiles, idx = np.unique(np.array(fragsmiles), return_index=True)
    fragmols = np.array(fragmols, dtype=object)[idx]
    fragaps = np.array(fragaps, dtype=int)[idx]
    fraparent = np.array(fraparent)[idx]
    fraparentmol = np.array(fraparentmol)[idx]
    fragidx = np.arange(fragsmiles.shape[0])

    # for i, fragmol in enumerate(fragmols):
    #    mcs = rdFMCS.FindMCS(
    #        [fraparentmol[i], fragmol],
    #        matchValences=True,
    #        ringMatchesRingOnly=True,
    #        completeRingsOnly=True,
    #    )
    #    hit_atoms_frag = list(
    #        fragmol.GetSubstructMatch(Chem.MolFromSmarts(mcs.smartsString))
    #    )
    #    hit_atoms_par = list(
    #        fraparentmol[i].GetSubstructMatch(Chem.MolFromSmarts(mcs.smartsString))
    #    )

    for i, j in zip(fraparent, fragsmiles):
        # print(f"{i},{j}")
        pass

    print(f"Unique substituents detected: {fragsmiles.shape[0]}")
    # Unique smiles, mol objects and attachment points are now stored in fragsmiles, fragmols and fragaps
    # We will use fluorinated cores and replace the fluors by the stored substituents

    patterns_f = [
        "FNC(N(F))=O",
        "FNC(N(F))=S",
        "O=C1C(N(F))=C(N(F))C1=O",
        "S=C1C(N(F))=C(N(F))C1=S",
        "O=C1C(N(F))=C1N(F)",
        "O=C(C(N(F))=C(N(F))C1=O)C1=O",
        "O=S(N(F))(N(F))=O",
    ]
    mol_patterns_f = np.array(
        [Chem.MolFromSmiles(smi) for smi in patterns_f], dtype=object
    )
    finalmols = []
    finalsmiles = []
    count = 0

    for mol_pattern_f in mol_patterns_f:
        print(f"Starting from template {Chem.MolToSmiles(mol_pattern_f)},")
        for i1, i2 in itertools.combinations_with_replacement(fragidx, 2):
            count += 1
            template = copy.deepcopy(mol_pattern_f)
            # print(f"Adding {fragsmiles[i1]}")
            r1 = copy.deepcopy(fragmols[i1])
            r2 = copy.deepcopy(fragmols[i2])
            mod_mol1 = Chem.ReplaceSubstructs(
                template,
                Chem.MolFromSmiles("F"),
                r1,
                False,
                int(fragaps[i1]),
            )[0]
            # print(
            #    f"Substitution 1 done, current molecule SMILES is {Chem.MolToSmiles(mod_mol1)}"
            # )
            # print(f"Adding {fragsmiles[i1]}")
            mod_mol12 = Chem.ReplaceSubstructs(
                mod_mol1,
                Chem.MolFromSmiles("F"),
                r2,
                False,
                int(fragaps[i2]),
            )[0]
            # print(
            #    f"Substitution 2 done, current molecule SMILES is {Chem.MolToSmiles(mod_mol12)}"
            # )
            Chem.SanitizeMol(
                mod_mol12,
                sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
                ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
                ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
            )
            currsmiles = Chem.MolToSmiles(mod_mol12)
            finalsmiles.append(currsmiles)
            finalmols.append(mod_mol12)

    finalsmiles = np.array(finalsmiles)
    finalmols = np.array(finalmols)
    print(
        f"{finalsmiles.shape[0]} molecules generated from substitutions using fragments."
    )
    finalsmiles, idx = np.unique(finalsmiles, return_index=True)
    finalmols = finalmols[idx]
    print(
        f"{finalsmiles.shape[0]} unique molecules generated from substitutions using fragments."
    )
    np.save("finalsmiles.npy", finalsmiles)
    np.save("finalmols.npy", finalmols)

indices = np.ones((finalsmiles.shape[0]), dtype=bool)

for i, smiles in enumerate(finalsmiles):
    if indices[i]:
        name = str(i).rjust(8, "0")
        p = Popen(
            ["python", "gen3d.py", f"{smiles}", f"{name}"], stdout=PIPE, stderr=PIPE
        )
        output, error = p.communicate()
        if p.returncode == 0:
            print(f"{name},{smiles}", file=of, flush=True)
        else:
            print("3D Generation failed %d %s %s" % (p.returncode, output, error))

exit()

## Print sublists from Method 2
#type1 = np.where(moltype == 1)
#type2 = np.where(moltype == 2)
#type3 = np.where(moltype == 3)
#type4 = np.where(moltype == 4)
#type5 = np.where(moltype == 5)
#type6 = np.where(moltype == 6)
#type7 = np.where(moltype == 7)
#types = [type1, type2, type3, type4, type5, type6, type7]
#for i, type in enumerate(types):
#    # print(f"\n\n {patterns[i]} list, with {np.count_nonzero(type)} elements:")
#    for name, smi, charge in zip(names[type], csmiles[type], charges[type]):
#        # print(f"{name},{smi},{charge}")
#        pass
