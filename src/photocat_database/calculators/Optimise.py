""" This module is used to optimise the molecule
using xtb. The module will take the smiles of the molecule and will
return the lowest energy of the molecule. The module
will also save the results in the database. """

import os
import re
import pymongo
import stk



import stko
import rdkit

class ETKDG(stko.Optimizer):
    """
    Uses ETKDG [3]_ v2 algorithm in :mod:`rdkit` [4]_ to optimize a structure.

    Examples:

        .. code-block:: python

            import stk
            import stko

            mol = stk.BuildingBlock('NCCNCCN')
            etkdg = stko.ETKDG()
            mol = etkdg.optimize(mol)

    References:

        .. [3] http://pubs.acs.org/doi/pdf/10.1021/acs.jcim.5b00654
        .. [4] https://www.rdkit.org/

    """

    def __init__(self, random_seed: int = 12):
        """
        Parameters:

            random_seed:
                The random seed to use.

        """

        self._random_seed = random_seed


    def optimize(self, mol: stk.Molecule) -> stk.Molecule:
        params = rdkit.Chem.rdDistGeom.ETKDGv3()
        params.clearConfs = True
        params.random_seed = self._random_seed
        params.useChirality = True
        rdkit_mol = mol.to_rdkit_mol()

        for i in range(10):
            rdkit.Chem.rdDistGeom.EmbedMolecule(rdkit_mol, params)

            conformer = rdkit_mol.GetConformer()

            mol_opt = mol.with_position_matrix(
                    position_matrix=rdkit_mol.GetConformer().GetPositions()
                )
            if stk.Smiles().get_key(mol_opt) ==stk.Smiles().get_key(mol):
                return mol_opt
        
        print('issue with the ETKDG optimisation of the molecules, could not keep the mol in the same conformation')
        return mol


class Optimise:
    """Class to optimise the molecule using xtb.
    
    The class will take the smiles of the molecule and will return the
    energy of the molecule. The class will also save the results in the
    database.
    
    Attributes
    ----------
    client : :class:`str`
    
    db_mol : :class:`str`
    
    xtb_path : :class:`str`
    
    Db_folder : :class:`str`
    
    host_IP : :class:`str`
    
    collection_name : :class:`str`
    """
    def __init__(self, client = "mongodb://localhost:27017/", 
                 db_mol= "test"):
        """ Initialize the class.
        """
        self.client = client
        self.db_mol= db_mol
        self.xtb_path = None
        self.Db_folder = None
        if self.Db_folder is not None:
            os.makedirs(self.Db_folder, exist_ok=True)
        self.host_IP = "cx1"
        self.collection_name = "Precursors"
        self.load_database()

    def load_database(self):
        """Function to load the database and collection name.
        """
        client = pymongo.MongoClient(self.client)
        try:
            self.stk_mol_database =  stk.MoleculeMongoDb(
                client,
                database=self.db_mol,
            )
        except Exception as e:
            print("Database not loaded")
            print(e)    


    def load_precursors(self, smile):
        """Function to generate stk building block from smiles."""
        return stk.BuildingBlock(
                smile, functional_groups=[
                    stk.BromoFactory()])


    def get_inchi_key(self, molecule):
        """Function to get the InchiKey of the molecule.
        
        Parameters
        ----------
        molecule : stk.BuildingBlock
            The molecule for which the InchiKey is required.
        
        Returns
        -------
        :class:`str`
            The InchiKey of the molecule.

        """
        return stk.InchiKey().get_key(molecule)
    

    def evaluate_element(self, smile):
        """Function to evaluate the element
        depending on the paths provided (xtb )
        the function will add those calculations to the model.

        Parameters
        ----------
        smile : :class:`str`
            The smile of the molecule.

        Returns
        -------
        :class:`float`
            The energy of the molecule.
        :class:`str`
            The InchiKey of the molecule.

        """
        # initialise the database
        client = pymongo.MongoClient(self.client)

        # define the path to xtb and stda
        xtb_path = self.xtb_path

        # define the output directories
        Db_folder = self.Db_folder
        xtb_opt_output_dir = os.path.join(
            Db_folder, "Database", "xtb_opt_output_dir"
        )
        os.makedirs(xtb_opt_output_dir, exist_ok=True)
        # define the database and collection name
        collection_name = self.collection_name
        # print(collection_name)
        precursor = self.load_precursors(smile)
        precursor, energy = self.run_xtb_opt(
                precursor,
                xtb_path,
                xtb_opt_output_dir,
                database=self.db_mol,
                collection=collection_name + "_opt",
                client=client,
                InchiKey_initial=self.get_inchi_key(precursor),
            )
        Inchikey = stk.InchiKey().get_key(precursor)
        return energy, Inchikey
    
    def evaluate_constructed(self, precursor):
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
            Db_folder, "Database", "xtb_opt_output_dir"
        )
        os.makedirs(xtb_opt_output_dir, exist_ok=True)
        # define the database and collection name
        collection_name = self.collection_name
        # print(collection_name)
        precursor, energy = self.run_xtb_opt(
                precursor,
                xtb_path,
                xtb_opt_output_dir,
                database=self.db_mol,
                collection=collection_name + "_opt",
                client=client,
                InchiKey_initial=self.get_inchi_key(precursor),
            )

        Inchikey = stk.InchiKey().get_key(precursor)
        return energy, Inchikey


    def run_ETKDG_opt(
        self,
        polymer,
    ):
        etkdg = stko.OptimizerSequence(
            ETKDG(),
        )
        polymer = etkdg.optimize(polymer)
        self.stk_mol_database.put(polymer)
        return polymer

    def run_xtb_opt(
        self,
        polymer,
        xtb_path,
        xtb_opt_output_dir,
        database="stk_mohammed_BO",
        collection="test",
        client=None,
        InchiKey_initial="",
    ):
        """Function to run xtb optimization on the molecule.
        
        Parameters
        ----------
        polymer : stk.BuildingBlock
            The molecule to be optimized.
            
        xtb_path : :class:`str`
            The path to the xtb executable.
            
            xtb_opt_output_dir : :class:`str`   
            The path to the output directory.
        
        database : :class:`str`, optional
            The name of the database. The default is "stk_mohammed_BO".
            
        collection : :class:`str`, optional
            The name of the collection. The default is "test".

        client : :class:`str`, optional
            The client to the database. The default is None.

        InchiKey_initial : :class:`str`, optional
            The initial InchiKey of the molecule. The default is "".

        Returns
        -------
        
        stk.BuildingBlock
            The optimized molecule.
            
        :class:`float`
            The energy of the molecule.
            
        """
        def save_xtb_opt_calculation(
            polymer, output_dir, collection=None, InchiKey_initial=None
        ) -> None:
            def get_property_value(data, property_name):
                for line in data:
                    if property_name in line:
                        if property_name == "cpu-time":
                            return (
                                re.findall(r"[-+]?(?:\d*\.*\d+)", line)[-3]
                                + " h "
                                + re.findall(r"[-+]?(?:\d*\.*\d+)", line)[-2]
                                + " min "
                                + re.findall(r"[-+]?(?:\d*\.*\d+)", line)[-1]
                                + " s "
                            )
                        return float(
                            re.findall(r"[-+]?(?:\d*\.*\d+)", line)[-1]
                        )  # float(words[3]) #
                return None
            
            polymer_xtb_opt_calc = {
                "InChIKey": stk.InchiKey().get_key(polymer),
                "cal_folder": output_dir,
                "Host IP": self.host_IP,
                "InChIKey_Initial": InchiKey_initial,
            }
            outfile = open(
                os.path.join(
                    polymer_xtb_opt_calc["cal_folder"], "optimization_1.output"
                ),
                encoding="utf8",
            )
            data = outfile.readlines()
            outfile.close()
            polymer_xtb_opt_calc["cpu time"] = get_property_value(
                data, "cpu-time"
            )
            polymer_xtb_opt_calc["total energy (au)"] = get_property_value(
                data, "TOTAL ENERGY"
            )
            polymer_xtb_opt_calc["HOMO-LUMO GAP (eV)"] = get_property_value(
                data, "HOMO-LUMO GAP"
            )
            collection.update_many(
                filter={"InChIKey": self.get_inchi_key(polymer)},
                update={"$set": polymer_xtb_opt_calc},
                upsert=True,
            )
            return polymer_xtb_opt_calc["total energy (au)"] 

        collection = client[database][collection]
        if (
            collection.find_one({"InChIKey_Initial": self.get_inchi_key(polymer)})
            is not None
        ):
            XTB_results = collection.find_one({"InChIKey_Initial": self.get_inchi_key(polymer)})
            InchiKey = XTB_results['InChIKey']
            try:
                molecule = self.stk_mol_database.get({"InChIKey": InchiKey})
                return molecule, XTB_results["total energy (au)"] 
            except:
                print('issue with the dataset')
                

        output_dir = os.path.join(xtb_opt_output_dir, self.get_inchi_key(polymer))
        os.makedirs(output_dir, exist_ok=True)
        xtb = self.optimisation_sequence(xtb_path, output_dir)

        polymer = xtb.optimize(polymer)
        energy = save_xtb_opt_calculation(
            polymer,
            output_dir,
            collection=collection,
            InchiKey_initial=InchiKey_initial,
        )

        self.stk_mol_database.put(polymer)
        
        return polymer, energy
    
    def optimisation_sequence(self, xtb_path, output_dir):
        xtb = stko.OptimizerSequence(
            ETKDG(),
            stko.XTB(
                xtb_path=xtb_path,
                output_dir=output_dir,
                unlimited_memory=True,
                num_cores=24,
            ),
        )
        return xtb 


