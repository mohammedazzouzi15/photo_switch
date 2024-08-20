""" This module is used to optimise the E and Z isomers of the molecule
using xtb. The module will take the smiles of the molecule and will
return the energy of the E and Z isomers of the molecule. The module
will also save the results in the database. """


import os
import pymongo
import stk
import stko
from photocat_database.calculators.Optimise import Optimise
from rdkit import Chem




class Optimise_EZ(Optimise):
    """Class to optimise the E and Z isomers of the molecule using xtb.

    The class will take the smiles of the molecule and will return the
    energy of the E and Z isomers of the molecule. The class will also
    save the results in the database.

    Attributes
    ----------
    client : :class:`str`

    db_mol_Z : :class:`str`

    db_mol_E : :class:`str`

    xtb_path : :class:`str`

    Db_folder : :class:`str`

    host_IP : :class:`str`

    collection_name : :class:`str`

    """

    def __init__(self,
                 client = "mongodb://localhost:27017/", 
                 db_mol= "test"):
        """Initialize the class.
        
        """
        super().__init__(client = client, db_mol= db_mol)
        self.db_mol_Z = "Z"
        self.db_mol_E = "E"
        self.load_databaseEZ()

    def load_databaseEZ(self):
        """Function to load the database and collection name.
        """
        client = pymongo.MongoClient(self.client)
        try:
            self.stk_mol_database_Z =  stk.MoleculeMongoDb(
                client,
                database=self.db_mol_Z,
            )
            self.stk_mol_database_E =  stk.MoleculeMongoDb(
                client,
                database=self.db_mol_E,
            )
        except Exception as e:
            print("Database not loaded")
            print(e)    




    def evaluate_element(self, smile,run_opt=True):
        """Function to evaluate the element
        depending on the paths provided (xtb or stda )
        the function will add those calculations to the model.

        """  
        # initialise the database
        client = pymongo.MongoClient(self.client)

        # define the path to xtb and stda
        xtb_path = self.xtb_path

        # define the output directories
        Db_folder = self.Db_folder
        xtb_opt_output_dir = os.path.join(
            Db_folder, "Database_Z_bb", "xtb_opt_output_dir"
        )
        os.makedirs(xtb_opt_output_dir, exist_ok=True)
        # define the database and collection name
        collection_name = self.collection_name
        precursor_Z = self.load_precursors(smile)
        self.stk_mol_database = self.stk_mol_database_Z

        precursor_E = Chem.MolFromSmiles(smile)
        e_db = [b.SetStereo(Chem.BondStereo.STEREOE) for b in precursor_E.GetBonds() if b.GetStereo() == Chem.BondStereo.STEREOZ]
        xtb_opt_output_dir = os.path.join(
            Db_folder, "Database_E_bb", "xtb_opt_output_dir"
        )
        os.makedirs(xtb_opt_output_dir, exist_ok=True)
        precursor_E = self.load_precursors(Chem.MolToSmiles(precursor_E))
        if run_opt:
            self.stk_mol_database = self.stk_mol_database_E
            precursor_E, energy_E = self.run_xtb_opt(
                precursor_E,
                xtb_path,
                xtb_opt_output_dir,
                database=self.db_mol_E,
                collection=collection_name + "_opt_E",
                client=client,
                InchiKey_initial=self.get_inchi_key(precursor_Z),
            )
            self.stk_mol_database = self.stk_mol_database_Z

            precursor_Z, energy = self.run_xtb_opt(
                precursor_Z,
                xtb_path,
                xtb_opt_output_dir,
                database=self.db_mol_Z,
                collection=collection_name + "_opt_Z",
                client=client,
                InchiKey_initial=self.get_inchi_key(precursor_Z),
            )
            Inchikey_Z = stk.InchiKey().get_key(precursor_Z)
            Inchikey_E = stk.InchiKey().get_key(precursor_E)
            return energy, Inchikey_Z, energy_E, Inchikey_E
        Inchikey_Z = stk.InchiKey().get_key(precursor_Z)
        Inchikey_E = stk.InchiKey().get_key(precursor_E)
        

        return 0, Inchikey_Z, 0, Inchikey_E







