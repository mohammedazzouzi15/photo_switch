"""list of function to build building objects from smiles
and run xtb calculations on them."""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymongo
import rdkit
import seaborn as sns
import stk
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

from photocat_database.calculators.Optimise_EZ import Optimise_EZ
from photocat_database.calculators.Optimise import Optimise
from photocat_database.calculators import Optimise_Constructed


def get_photo_switch_EZ_bb(mol_smiles):
    """Function to get the photo switch molecule from the smiles.

    Parameters
    ----------
    mol_smiles : :class:`str`
        The smiles of the molecule.

    Returns
    -------
    :class:`stk.BuildingBlock`
        The photo switch molecule.
    """

    calculator = Optimise_EZ(
        client="mongodb://localhost:27017/",
        db_mol="stk_photo_cat_EZ",
    )

    calculator.db_mol_E = "stk_photo_cat_E"  # name of the molecules database
    calculator.db_mol_Z = "stk_photo_cat_Z"  # name of the molecules database
    calculator.load_databaseEZ()
    # Change the xtb path to add the calculation to add to the database
    calculator.xtb_path = "xtb"
    calculator.Db_folder = "data/xtb_calculation/"
    os.makedirs(calculator.Db_folder, exist_ok=True)
    for i in range(1, 10):
        mol_smiles = mol_smiles.replace(f"*:{i}", "Br")
    energy_E, Inchikey_photoswitch_Z, energy_E, Inchikey_photoswitch_E = (
        calculator.evaluate_element(mol_smiles)
    )
    client = pymongo.MongoClient(calculator.client)
    db_polymer = stk.MoleculeMongoDb(
        client,
        database=calculator.db_mol_Z,
    )
    photoswitch_Z = db_polymer.get({"InChIKey": Inchikey_photoswitch_Z})
    db_polymer = stk.MoleculeMongoDb(
        client,
        database=calculator.db_mol_E,
    )
    photoswitch_E = db_polymer.get({"InChIKey": Inchikey_photoswitch_E})
    # display(photoswitch_E.to_rdkit_mol())
    # display(photoswitch_Z.to_rdkit_mol())
    return photoswitch_E, photoswitch_Z


def get_right_bb(mol_smiles):
    """Function to get the right block molecule from the smiles."""
    calculator = Optimise()
    calculator.client = "mongodb://localhost:27017/"
    calculator.db_mol = "right_block"  # name of the molecules database
    calculator.load_database()
    # Change the xtb path to add the calculation to add to the database
    calculator.xtb_path = "xtb"
    calculator.Db_folder = "data/xtb_calculation/"
    os.makedirs(calculator.Db_folder, exist_ok=True)

    mol_smiles = mol_smiles.replace("*:2", "Br")
    energy, inchikey_rb = calculator.evaluate_element(mol_smiles)
    client = pymongo.MongoClient(calculator.client)
    db_polymer = stk.MoleculeMongoDb(
        client,
        database=calculator.db_mol,
    )
    right_block = db_polymer.get({"InChIKey": inchikey_rb})
    # display(right_block.to_rdkit_mol())
    return right_block


def get_catalytic_bb(mol_smiles):
    calculator = Optimise()
    calculator.client = "mongodb://localhost:27017/"
    calculator.db_mol = "catalytic_unit"  # name of the molecules database
    # Change the xtb path to add the calculation to add to the database
    calculator.xtb_path = "xtb"
    calculator.Db_folder = "data/xtb_calculation/"
    calculator.load_database()

    # change the stda bin path if you want to do the calculation of the STDA to add it to the database
    calculator.STDA_bin_path = None  # "/home/mohammed/Work/xtb4stda/"

    os.makedirs(calculator.Db_folder, exist_ok=True)
    calculator.host_IP = "cx1"
    calculator.collection_name = "Precursors"

    mol_smiles = mol_smiles.replace("*:2", "Br")
    mol_smiles = mol_smiles.replace("*:1", "Br")
    energy, inchikey_cat = calculator.evaluate_element(mol_smiles)
    client = pymongo.MongoClient(calculator.client)
    db_polymer = stk.MoleculeMongoDb(
        client,
        database=calculator.db_mol,
    )
    catalytic = db_polymer.get({"InChIKey": inchikey_cat})
    # display(catalytic.to_rdkit_mol())
    return catalytic


def get_inhibitor_bb(mol_smiles):
    calculator = Optimise()
    calculator.client = "mongodb://localhost:27017/"
    calculator.db_mol = "inhibitor_unit"  # name of the molecules database
    # Change the xtb path to add the calculation to add to the database
    calculator.xtb_path = "xtb"
    calculator.Db_folder = "data/xtb_calculation/"
    calculator.load_database()

    # change the stda bin path if you want to do the calculation of the STDA to add it to the database
    calculator.STDA_bin_path = None  # "/home/mohammed/Work/xtb4stda/"

    os.makedirs(calculator.Db_folder, exist_ok=True)
    calculator.host_IP = "cx1"
    calculator.collection_name = "Precursors"

    mol_smiles = mol_smiles.replace("*:2", "Br")
    mol_smiles = mol_smiles.replace("*:1", "Br")
    energy, inchikey_inhibitor = calculator.evaluate_element(mol_smiles)
    client = pymongo.MongoClient(calculator.client)
    db_polymer = stk.MoleculeMongoDb(
        client,
        database=calculator.db_mol,
    )
    inhibitor_bb = db_polymer.get({"InChIKey": inchikey_inhibitor})
    # display(inhibitor_bb.to_rdkit_mol())
    return inhibitor_bb


def build_EZ_photocat_4blocks(
    inhibitor_bb, photoswitch_Z, photoswitch_E, catalytic, right_block
):
    inhibitor_bb = stk.BuildingBlock.init_from_molecule(
        inhibitor_bb, functional_groups=[stk.BromoFactory()]
    )
    photoswitch_Z = stk.BuildingBlock.init_from_molecule(
        photoswitch_Z, functional_groups=[stk.BromoFactory()]
    )
    photoswitch_E = stk.BuildingBlock.init_from_molecule(
        photoswitch_E, functional_groups=[stk.BromoFactory()]
    )

    catalytic = stk.BuildingBlock.init_from_molecule(
        catalytic, functional_groups=[stk.BromoFactory()]
    )
    right_block = stk.BuildingBlock.init_from_molecule(
        right_block, functional_groups=[stk.BromoFactory()]
    )
    constructed_molecule_Z = stk.ConstructedMolecule(
        stk.polymer.Linear(
            building_blocks=[
                inhibitor_bb,
                photoswitch_Z,
                catalytic,
                right_block,
            ],
            repeating_unit="ABCD",
            num_repeating_units=1,
            optimizer=stk.MCHammer(),
        )
    )
    constructed_molecule_E = stk.ConstructedMolecule(
        stk.polymer.Linear(
            building_blocks=[
                inhibitor_bb,
                photoswitch_E,
                catalytic,
                right_block,
            ],
            repeating_unit="ABCD",
            num_repeating_units=1,
            optimizer=stk.MCHammer(),
        )
    )
    # display(constructed_molecule_E.to_rdkit_mol())
    # display(constructed_molecule_Z.to_rdkit_mol())
    return constructed_molecule_E, constructed_molecule_Z


def build_EZ_photocat_3blocks(
    inhibitor_bb, photoswitch_Z, photoswitch_E, catalytic
):
    inhibitor_bb = stk.BuildingBlock.init_from_molecule(
        inhibitor_bb, functional_groups=[stk.BromoFactory()]
    )
    photoswitch_Z = stk.BuildingBlock.init_from_molecule(
        photoswitch_Z, functional_groups=[stk.BromoFactory()]
    )
    photoswitch_E = stk.BuildingBlock.init_from_molecule(
        photoswitch_E, functional_groups=[stk.BromoFactory()]
    )

    catalytic = stk.BuildingBlock.init_from_molecule(
        catalytic, functional_groups=[stk.BromoFactory()]
    )

    constructed_molecule_Z = stk.ConstructedMolecule(
        stk.polymer.Linear(
            building_blocks=[inhibitor_bb, photoswitch_Z, catalytic],
            repeating_unit="ABC",
            num_repeating_units=1,
            optimizer=stk.MCHammer(),
        )
    )
    constructed_molecule_E = stk.ConstructedMolecule(
        stk.polymer.Linear(
            building_blocks=[inhibitor_bb, photoswitch_E, catalytic],
            repeating_unit="ABC",
            num_repeating_units=1,
            optimizer=stk.MCHammer(),
        )
    )
    # display(constructed_molecule_E.to_rdkit_mol())
    # display(constructed_molecule_Z.to_rdkit_mol())
    return constructed_molecule_E, constructed_molecule_Z


def optimse_EZ_photocat(constructed_molecule_Z, constructed_molecule_E):
    calculator = Optimise_Constructed.Optimise_Constructed(isomers="Z")
    calculator.client = "mongodb://localhost:27017/"
    # Change the xtb path to add the calculation to add to the database
    calculator.xtb_path = "xtb"
    # change the stda bin path if you want to do the calculation of the STDA to add it to the database
    calculator.STDA_bin_path = None  # "/home/mohammed/Work/xtb4stda/"

    calculator.host_IP = "cx1"
    calculator.collection_name = "XTB"
    client = pymongo.MongoClient(calculator.client)
    calculator.db_mol = "constructed_photocat_Z"
    calculator.Db_folder = "data/xtb_calculation_Z/"
    calculator.load_database()
    energy_Z, inchikey_photo_cat_Z = calculator.evaluate_constructed(
        constructed_molecule_Z
    )
    db_polymer = stk.ConstructedMoleculeMongoDb(
        client,
        database=calculator.db_mol,
    )
    constructed_molecule_Z = db_polymer.get(
        {"InChIKey": inchikey_photo_cat_Z}
    )
    calculator.db_mol = "constructed_photocat_E"
    calculator.Db_folder = "data/xtb_calculation_E/"
    calculator.isomers = "E"
    calculator.stk_mol_database = stk.ConstructedMoleculeMongoDb(
        client,
        database=calculator.db_mol,
    )
    energy_E, inchikey_photo_cat_E = calculator.evaluate_constructed(
        constructed_molecule_E
    )
    db_polymer = stk.ConstructedMoleculeMongoDb(
        client,
        database=calculator.db_mol,
    )
    constructed_molecule_E = db_polymer.get(
        {"InChIKey": inchikey_photo_cat_E}
    )
    print("energy_Z", energy_Z)
    print("energy_E", energy_E)
    return constructed_molecule_E, constructed_molecule_Z


from photocat_database.calculators import Optimise_ON

def run_opt_constrained(constructed_molecule_E, constructed_molecule_Z):
    calculator = Optimise_ON.Optimise_ON(isomers="Z")
    calculator.client = "mongodb://localhost:27017/"
    # Change the xtb path to add the calculation to add to the database
    calculator.xtb_path = "xtb"
    # change the stda bin path if you want to do the calculation of the STDA to add it to the database
    calculator.STDA_bin_path = None  # "/home/mohammed/Work/xtb4stda/"

    calculator.host_IP = "cx1"
    calculator.collection_name = "XTB"
    client = pymongo.MongoClient(calculator.client)
    calculator.db_mol = "constructed_photocat_Z_ON"
    calculator.Db_folder = "data/xtb_calculation_Z_ON/"
    calculator.load_database()
    energy_Z, inchikey_photo_cat_Z = calculator.evaluate_constructed(
        constructed_molecule_Z
    )
    db_polymer = stk.ConstructedMoleculeMongoDb(
        client,
        database=calculator.db_mol,
    )
    constructed_molecule_Z_ON = db_polymer.get(
        {"InChIKey": inchikey_photo_cat_Z}
    )
    calculator.db_mol = "constructed_photocat_E_ON"
    calculator.Db_folder = "data/xtb_calculation_E_ON/"
    calculator.isomers = "E"
    calculator.stk_mol_database = stk.ConstructedMoleculeMongoDb(
        client,
        database=calculator.db_mol,
    )
    energy_E, inchikey_photo_cat_E = calculator.evaluate_constructed(
        constructed_molecule_E
    )
    db_polymer = stk.ConstructedMoleculeMongoDb(
        client,
        database=calculator.db_mol,
    )
    constructed_molecule_E_ON = db_polymer.get(
        {"InChIKey": inchikey_photo_cat_E}
    )
    print("energy_Z", energy_Z)
    print("energy_E", energy_E)
    return constructed_molecule_E_ON, constructed_molecule_Z_ON
