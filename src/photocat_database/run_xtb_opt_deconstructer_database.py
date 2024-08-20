"""This script constructs and optimizes photocatalyst molecules using specified building blocks.
to run this script, you need to check that the paths to the building blocks are correct.
and then run the script with the following command:
python run_xtb_opt_deconstructer_database.py -p 0 -r 0 -l 0 -i 0
where the arguments are the IDs of the building blocks in the database.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from photocat_database import building_utils
import stk


def main(id_photocat, id_right_block, id_linker, id_inhibitor):
    df_organic_photocat = pd.read_csv(
        "/media/mohammed/Work/Work/photo_switch/notebooks/data/2_murcko_azo_photo_cores_Br_filtered.csv"
    )
    mol_smiles = df_organic_photocat["smiles"][id_photocat]
    photoswitch_E, photoswitch_Z = building_utils.get_photo_switch_EZ_bb(
        mol_smiles
    )

    df_right_block = pd.read_csv(
        "/media/mohammed/Work/Work/photo_switch/notebooks/data/catalyst_right_blocks_unique.csv"
    )
    mol_smiles = df_right_block["smiles"][id_right_block]
    right_block = building_utils.get_right_bb(mol_smiles)

    df_linker = pd.read_csv(
        "/media/mohammed/Work/Work/photo_switch/notebooks/data/catalyst_linkers_unique.csv"
    )
    mol_smiles = df_linker["smiles"][id_linker]
    catalytic = building_utils.get_catalytic_bb(mol_smiles)

    df_inhibitor = pd.read_csv(
        "/media/mohammed/Work/Work/photo_switch/notebooks/data/inhibitor_right_blocks.csv"
    )
    mol_smiles = df_inhibitor["smiles"][id_inhibitor]
    inhibitor_bb = building_utils.get_inhibitor_bb(mol_smiles)

    constructed_molecule_E, constructed_molecule_Z = (
        building_utils.build_EZ_photocat_4blocks(
            inhibitor_bb, photoswitch_Z, photoswitch_E, catalytic, right_block
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
    import argparse as ap

    ap = ap.ArgumentParser(
        description="This script constructs and optimizes photocatalyst molecules using specified building blocks."
    )
    ap.add_argument(
        "-p",
        "--id_photocat",
        type=int,
        help="ID of the photocatalyst in the database",
        default=0,
    )
    ap.add_argument(
        "-r",
        "--id_right_block",
        type=int,
        help="ID of the right block in the database",
        default=0,
    )
    ap.add_argument(
        "-l",
        "--id_linker",
        type=int,
        help="ID of the linker in the database",
        default=0,
    )
    ap.add_argument(
        "-i",
        "--id_inhibitor",
        type=int,
        help="ID of the inhibitor in the database",
        default=0,
    )
    # add help on how to use the script
    # Parse arguments once and store the result
    args = ap.parse_args()

    main(
        args.id_photocat,
        args.id_right_block,
        args.id_linker,
        args.id_inhibitor,
    )
