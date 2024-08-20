import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from photocat_database import building_utils
from photocat_database.calculators import Optimise_ON
import pymongo
import stk
import os


def main():
    mol_smiles_list = pd.read_csv("/media/mohammed/Work/Work/photo_switch/notebooks/data/PSW_database_fragmented.csv")
    photoswitch_E, photoswitch_Z = building_utils.get_photo_switch_EZ_bb(
        mol_smiles_list['linkers_smiles'][0]
    )
    catalytic = building_utils.get_catalytic_bb(mol_smiles_list['smiles_catalytic'][0])
    inhibitor_bb = building_utils.get_inhibitor_bb(mol_smiles_list['smiles_inhibitor'][0])
    constructed_molecule_E, constructed_molecule_Z = (
        building_utils.build_EZ_photocat_3blocks(
            catalytic,
            photoswitch_Z,
            photoswitch_E,
            inhibitor_bb,
        )
    )
    constructed_molecule_E, constructed_molecule_Z = (
        building_utils.optimse_EZ_photocat(
            constructed_molecule_Z, constructed_molecule_E
        )
    )
    print(stk.InchiKey().get_key(constructed_molecule_E))
    print(stk.InchiKey().get_key(constructed_molecule_Z))

    constructed_molecule_E_ON, constructed_molecule_Z_ON = (
        building_utils.run_opt_constrained(
            constructed_molecule_E, constructed_molecule_Z
        )
    )
    print(stk.InchiKey().get_key(constructed_molecule_E_ON))
    print(stk.InchiKey().get_key(constructed_molecule_Z_ON))


if __name__ == "__main__":
    main()
