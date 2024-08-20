"""this module is used to optimize the constructer of the database
using xtb in the On structure. The module will take the smiles of the molecule and will
return the energy of the molecule. The module
will also save the results in the database."""

from photocat_database.calculators.Optimise import Optimise
import pymongo
import stk
from rdkit.Chem import rdDistGeom
import rdkit.DistanceGeometry as DG
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms
import numpy as np
import stko
from rdkit.Chem import rdMMPA
import os
from rdkit.Chem import FragmentMatcher
from rdkit import Chem
import re
from photocat_database.calculators.Optimise_Constructed import ETKDG_isomer


def initialise_matchers(patterns):
    pmatchers = [FragmentMatcher.FragmentMatcher() for p in patterns]
    for p, pm in zip(patterns, pmatchers):
        pm.Init(p)

    return pmatchers


class ETKDG_constaint(ETKDG_isomer):
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

    def __init__(self, random_seed: int = 12,isomers = "E"):
        """
        Parameters:

            random_seed:
                The random seed to use.

        """

        self._random_seed = random_seed
        self.isomers = isomers

    def optimize(self, mol: stk.Molecule) -> stk.Molecule:
        return self.generate_constained_initial_geom(mol)

    def get_atom_position_to_induce_constain_for_3bb(self, stk_mol):
        patterns_catalytic = [
            "NC(N)=O",
            "NC(N)=S",
            "O=C1C(N)=C(N)C1=O",
            "S=C1C(N)=C(N)C1=S",
            "O=C1C(N)=C1N",
            "O=C(C(N)=C(N)C1=O)C1=O",
            "O=S(N)(N)=O",
        ]
        patterns_inhibitor = [
            "*[N+]([O-])=O",
            "*[N+]=O([O-])",
        ]
        pmatchers = initialise_matchers(patterns_catalytic)

        def get_catalytic_bb_hydrogen_atoms(bb, patterns):
            pmatchers = initialise_matchers(patterns)
            mol = bb.to_rdkit_mol()
            for pm in pmatchers:
                if pm.HasMatch(mol):
                    match_id = pm.GetMatch(mol)
                    atom_type_id = [
                        mol.GetAtomWithIdx(x).GetAtomicNum() for x in match_id
                    ]
                    hydrogen_atom_id = []
                    for id in match_id:
                        hydrogen_atom_id.extend(
                            [
                                x.GetIdx()
                                for x in mol.GetAtomWithIdx(id).GetNeighbors()
                                if x.GetAtomicNum() == 1
                            ]
                        )
                    return hydrogen_atom_id
            print(" no substructure found")
            return []

        def get_inhibitor_bb_oxygen_atoms(bb, patterns):
            pmatchers = initialise_matchers(patterns)
            mol = bb.to_rdkit_mol()
            for pm in pmatchers:
                if pm.HasMatch(mol):
                    match_id = pm.GetMatch(mol)
                    atom_type_id = [
                        mol.GetAtomWithIdx(x).GetAtomicNum() for x in match_id
                    ]
                    oxygen_atom_id = []
                    for id in match_id:
                        oxygen_atom_id.extend(
                            [
                                x.GetIdx()
                                for x in mol.GetAtomWithIdx(id).GetNeighbors()
                                if x.GetAtomicNum() == 8
                            ]
                        )
                    return oxygen_atom_id
            print(" no substructure found")
            return []

        bb = list(stk_mol.get_building_blocks())[0]
        hydrogen_atom_id = get_catalytic_bb_hydrogen_atoms(
            bb, patterns_catalytic
        )
        bb = list(stk_mol.get_building_blocks())[2]
        oxygen_atom_id = get_inhibitor_bb_oxygen_atoms(bb, patterns_inhibitor)
        oxygen_id = []
        hydrogen_id = []
        for b in stk_mol.get_atom_infos():
            if b.get_building_block_id() == 2:
                if b.get_building_block_atom().get_id() in oxygen_atom_id:
                    oxygen_id.append(b.get_atom().get_id())
            if b.get_building_block_id() == 0:
                if b.get_building_block_atom().get_id() in hydrogen_atom_id:
                    hydrogen_id.append(b.get_atom().get_id())
        return oxygen_id, hydrogen_id

    def get_atom_position_to_induce_constain_for_4bb(self, stk_mol):
        patterns_catalytic = [
            "NC(N)=O",
            "NC(N)=S",
            "O=C1C(N)=C(N)C1=O",
            "S=C1C(N)=C(N)C1=S",
            "O=C1C(N)=C1N",
            "O=C(C(N)=C(N)C1=O)C1=O",
            "O=S(N)(N)=O",
        ]
        patterns_inhibitor = [
            "*[N+]([O-])=O",
            "*[N+]=O([O-])",
        ]
        pmatchers = initialise_matchers(patterns_catalytic)

        def get_catalytic_bb_hydrogen_atoms(bb, patterns):
            pmatchers = initialise_matchers(patterns)
            mol = bb.to_rdkit_mol()
            for pm in pmatchers:
                if pm.HasMatch(mol):
                    match_id = pm.GetMatch(mol)
                    atom_type_id = [
                        mol.GetAtomWithIdx(x).GetAtomicNum() for x in match_id
                    ]
                    hydrogen_atom_id = []
                    for id in match_id:
                        hydrogen_atom_id.extend(
                            [
                                x.GetIdx()
                                for x in mol.GetAtomWithIdx(id).GetNeighbors()
                                if x.GetAtomicNum() == 1
                            ]
                        )
                    return hydrogen_atom_id
            print(" no substructure found")
            return []

        def get_inhibitor_bb_oxygen_atoms(bb, patterns):
            pmatchers = initialise_matchers(patterns)
            mol = bb.to_rdkit_mol()
            for pm in pmatchers:
                if pm.HasMatch(mol):
                    match_id = pm.GetMatch(mol)
                    atom_type_id = [
                        mol.GetAtomWithIdx(x).GetAtomicNum() for x in match_id
                    ]
                    oxygen_atom_id = []
                    for id in match_id:
                        oxygen_atom_id.extend(
                            [
                                x.GetIdx()
                                for x in mol.GetAtomWithIdx(id).GetNeighbors()
                                if x.GetAtomicNum() == 8
                            ]
                        )
                    return oxygen_atom_id
            print(" no substructure found")
            return []

        bb = list(stk_mol.get_building_blocks())[2]
        hydrogen_atom_id = get_catalytic_bb_hydrogen_atoms(
            bb, patterns_catalytic
        )
        bb = list(stk_mol.get_building_blocks())[0]
        oxygen_atom_id = get_inhibitor_bb_oxygen_atoms(bb, patterns_inhibitor)
        oxygen_id = []
        hydrogen_id = []
        for b in stk_mol.get_atom_infos():
            if b.get_building_block_id() == 0:
                if b.get_building_block_atom().get_id() in oxygen_atom_id:
                    oxygen_id.append(b.get_atom().get_id())
            if b.get_building_block_id() == 2:
                if b.get_building_block_atom().get_id() in hydrogen_atom_id:
                    hydrogen_id.append(b.get_atom().get_id())
        return oxygen_id, hydrogen_id


    def generate_constained_initial_geom(self, stk_mol):
        if len(list(stk_mol.get_building_blocks())) == 3:
            oxygen_id, hydrogen_id = (
                self.get_atom_position_to_induce_constain_for_3bb(stk_mol)
            )
        elif len(list(stk_mol.get_building_blocks()))  == 4:
            oxygen_id, hydrogen_id = (
                self.get_atom_position_to_induce_constain_for_4bb(stk_mol)
            )
    
        params = AllChem.ETKDGv3()
        mol = self.get_constructed_molecule_isomer(stk_mol)
        bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
        for o_id, h_id in zip(oxygen_id, hydrogen_id):
            bm[o_id, h_id] = 1.9
        bm[hydrogen_id[0], hydrogen_id[1]] = 2.5

        DG.DoTriangleSmoothing(bm)
        params.useRandomCoords = True

        params.SetBoundsMat(bm)
        cids = rdDistGeom.EmbedMultipleConfs(mol, 50, params)
        # order conformer by distance between oxygen and hydrogen
        distance_list = []
        conformers = mol.GetConformers()

        for o_id, h_id in zip(oxygen_id, hydrogen_id):
            distance_list.append(
                np.array(
                    [
                        rdMolTransforms.GetBondLength(conf, o_id, h_id)
                        for conf in conformers
                    ]
                )
            )
        conf = mol.GetConformer()
        distance_list[0].argsort()
        for conf_num in distance_list[0].argsort():
            mol_opt = stk_mol.with_position_matrix(
                position_matrix=conformers[int(conf_num)].GetPositions()
            )
            if stk.Smiles().get_key(mol_opt) == stk.Smiles().get_key(stk_mol):
                return mol_opt
        print("no conformer found")
        return stk_mol


class Optimise_ON(Optimise):
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
        self.isomers = isomers

    def load_database(self):
        client = pymongo.MongoClient(self.client)

        self.stk_mol_database = stk.ConstructedMoleculeMongoDb(
            client,
            database=self.db_mol,
        )

    def optimisation_sequence(self, xtb_path, output_dir):
        xtb = stko.OptimizerSequence(
            ETKDG_constaint(self.isomers),
            stko.XTB(
                xtb_path=xtb_path,
                output_dir=output_dir,
                unlimited_memory=True,
                num_cores=24,
            ),
        )
        return xtb 
    
