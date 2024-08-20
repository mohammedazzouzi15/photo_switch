"""this module is used to optimize the constructer of the database
using xtb. The module will take the smiles of the molecule and will
return the energy of the molecule. The module
will also save the results in the database."""

import pymongo
import stk
import stko
from photocat_database.calculators.Optimise import ETKDG, Optimise
from rdkit import Chem
from rdkit.Chem import FragmentMatcher


class Optimise_Constructed(Optimise):
    """Class to optimize the constructed molecule using xtb.

    The class will take the smiles of the molecule and will return the
    energy of the molecule. The class will also save the results in the
    database.

    Attributes
    ----------
    client : :class:`str`
    """

    def __init__(self,isomers):
        """Initialize the class."""
        super().__init__()
        self.isomers = isomers # can be "E" or "Z"

    def load_database(self):
        client = pymongo.MongoClient(self.client)

        self.stk_mol_database = stk.ConstructedMoleculeMongoDb(
            client,
            database=self.db_mol,
        )

    def optimisation_sequence(self, xtb_path, output_dir):
        xtb = stko.OptimizerSequence(
            ETKDG_isomer(isomers= self.isomers),
            stko.XTB(
                xtb_path=xtb_path,
                output_dir=output_dir,
                unlimited_memory=True,
                num_cores=24,
            ),
        )
        return xtb 
    
    


class ETKDG_isomer(stko.Optimizer):
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

    def __init__(self,isomers, random_seed: int = 12):
        """
        Parameters:

            random_seed:
                The random seed to use.

        """
        self.isomers = "E"  # can be "E" or "Z"
        self._random_seed = random_seed


    def optimize(self, mol: stk.Molecule) -> stk.Molecule:
        params = Chem.rdDistGeom.ETKDGv3()
        params.clearConfs = True
        params.random_seed = self._random_seed
        params.useChirality = True
        rdkit_mol = self.get_constructed_molecule_isomer(mol)
        Chem.rdDistGeom.EmbedMolecule(rdkit_mol, params)
        mol_opt = mol.with_position_matrix(
                position_matrix=rdkit_mol.GetConformer().GetPositions()
            )
        return mol_opt
        
        print('issue with the ETKDG optimisation of the molecules, could not keep the mol in the same conformation')
        return mol
    
    def get_constructed_molecule_isomer(self, stk_mol):
        """Function to get the constructed molecule isomer.

        Parameters
        ----------
        smiles : :class:`str`
            The smiles of the molecule.

        Returns
        -------
        stk.ConstructedMolecule
            The constructed molecule.
        """
        N_id, mol = self.get_NNdoublebond_atoms(stk_mol)
        if N_id is None:
            print("No NN double bond found")
            return None
        Chem.SanitizeMol(mol)
        if self.isomers == "E":
            Chem.rdMolTransforms.SetDihedralDeg(
                mol.GetConformer(),
                N_id[0],
                N_id[1],
                N_id[2],
                N_id[3],
                0,
            )
        if self.isomers == "Z":
            Chem.rdMolTransforms.SetDihedralDeg(
                mol.GetConformer(),
                N_id[0],
                N_id[1],
                N_id[2],
                N_id[3],
                180,
            )
        return mol

    def initialise_matchers(self, patterns):
        pmatchers = [FragmentMatcher.FragmentMatcher() for p in patterns]
        for p, pm in zip(patterns, pmatchers):
            pm.Init(p)

        return pmatchers

    def get_NNdoublebond_atoms(self, stk_mol):
        patterns = [
            "N=N",
            "[N][N]"
        ]
        pmatchers = self.initialise_matchers(patterns)
        mol = stk_mol.to_rdkit_mol()
        for pm in pmatchers:
            if pm.HasMatch(mol):
                match_id = pm.GetMatch(mol)
                atom_type_id = [
                    mol.GetAtomWithIdx(x).GetAtomicNum() for x in match_id
                ]
                NNatom_id = []
                NNatom_id.extend([
                            x.GetIdx()
                            for x in mol.GetAtomWithIdx(match_id[0]).GetNeighbors()
                            if x.GetIdx() not in match_id
                        ])
                NNatom_id.extend(match_id)
                NNatom_id.extend([
                            x.GetIdx()
                            for x in mol.GetAtomWithIdx(match_id[1]).GetNeighbors()
                            if x.GetIdx() not in match_id
                        ])

                # get only unique atoms ids
                #NNatom_id = list(set(NNatom_id))
                #sort the atoms ids
                #NNatom_id.sort()
                print(NNatom_id)
                return NNatom_id, mol
        return None, mol
