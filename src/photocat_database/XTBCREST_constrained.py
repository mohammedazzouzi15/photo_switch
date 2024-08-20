"""script to run crest for the optimized molecules with the constrained dihedral angle
starting from 4 different structures
to run this script you need to provide the Inchikey of the photocatalyst
use the following command to run the script
python XTBCREST_constrained.py -IE "Inchikey_photo_cat_E" -IZ "Inchikey_photo_cat_Z.
"""


from stko import XTBCREST
import subprocess as sp
import shutil
import uuid
import logging
import os
import stk
import pymongo


class XTBCREST_constrained(XTBCREST):
    def __init__(self, isomer, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.isomer = isomer

    def _run_crest(self, xyz: str, out_file: str) -> None:
        """Run CREST along side GFN-xTB.

        Parameters:

            xyz:
                The name of the input structure ``.xyz`` file.

            out_file:
                The name of output file with xTB results.

        """
        # Modify the memory limit.

        if self._unlimited_memory:
            memory = "ulimit -s unlimited ;"
        else:
            memory = ""

        if self._solvent is not None:
            solvent = f"--{self._solvent_model} {self._solvent}"
        else:
            solvent = ""

        # Set optimization level and type.
        optimization = f"-opt {self._opt_level}"
        mdlen = "" if self._mdlen is None else f"-mdlen {self._mdlen} "
        keepdirs = "-keepdir" if self._keepdir is True else ""
        if self._speed_setting is not None:
            speed_settings = f"-{self._speed_setting}"
            ewin = ""
        else:
            speed_settings = ""
            ewin = f"-ewin {self._ewin} "
        cross = "-nocross" if self._cross is False else ""

        cmd = (
            f"{memory} {self._crest_path} {xyz} "
            f"--gfn2 -T {self._num_cores} "
            f"-I det_control.in"
        )

        with open(out_file, "w") as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(
                cmd,
                # stdin=sp.PIPE,
                stdout=f,
                # stderr=sp.PIPE,
                # Shell is required to run complex arguments.
                shell=True,
            )
        print(cmd)

    def _write_detailed_control(self, atoms_list_dihedral: list) -> None:
        atoms_list_dihedral = " ,".join(str(e) for e in atoms_list_dihedral)
        if self.isomer == "E":
            dihedral_angle = 180
        elif self.isomer == "Z":
            dihedral_angle = 0
        string = (
            f"$constrain \n"
            f"dihedral: {atoms_list_dihedral}, {dihedral_angle} \n"
            "$end"
        )

        with open("det_control.in", "w") as f:
            f.write(string)

    def optimize(self, mol):
        """Optimize `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns:
        -------
        mol : :class:`.Molecule`
            The optimized molecule.

        """
        from photocat_database.calculators.Optimise_ON import ETKDG_constaint

        if self._output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self._output_dir
        output_dir = os.path.abspath(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        init_dir = os.getcwd()
        calculator = ETKDG_constaint()
        atoms_list_dihedral, rdkit_mol = calculator.get_NNdoublebond_atoms(mol)
        os.chdir(output_dir)
        self._write_detailed_control(atoms_list_dihedral)
        try:
            mol, complete = self._run_optimization(mol)
        finally:
            os.chdir(init_dir)

        if not complete:
            logging.warning(f"CREST run is incomplete for {mol}.")

        return mol


def run_crest_from_isomer_opt(Inchikey_photo_cat, isomer="E", state="_ON"):
    client = "mongodb://localhost:27017/"
    client = pymongo.MongoClient(client)
    database = f"constructed_photocat_{isomer}{state}"
    db_polymer = stk.ConstructedMoleculeMongoDb(
        client,
        database=database,
    )
    collection = "crest"
    collection = client[database][collection]
    if (
        collection.find_one({"InChIKey": Inchikey_photo_cat})
        is not None
    ):
        print("Already run")
        crest_results = collection.find_one({"InChIKey": Inchikey_photo_cat})
        return crest_results["crest_conformers"]
    constructed_molecule = db_polymer.get({"InChIKey": Inchikey_photo_cat})
    db_folder = f"data/crest/{database}/"
    xtb_opt_output_dir = os.path.join(db_folder, "Crest")
    output_dir = os.path.join(xtb_opt_output_dir, Inchikey_photo_cat)
    client = "mongodb://localhost:27017/"
    XTBCREST_constraint = XTBCREST_constrained(
        isomer=isomer,
        crest_path="/media/mohammed/Work/bin/crest",
        xtb_path="xtb",
        num_cores=24,
        gfn_version=2,
        charge=0,
        unlimited_memory=True,
        output_dir=output_dir,
    )
    XTBCREST_constraint.optimize(constructed_molecule)
    conformer_file = os.path.join(output_dir, "crest_conformers.xyz")
    conformer_file = os.path.abspath(conformer_file)
    polymer_xtb_opt_calc = {
        "crest_conformers": conformer_file,
    }
    collection.update_many(
                filter={"InChIKey": Inchikey_photo_cat},
                update={"$set": polymer_xtb_opt_calc},
                upsert=True,
            )

    return conformer_file
    # add the optimized molecule to the database
   


def main(Inchikey_photo_cat_E, Inchikey_photo_cat_Z):
    conformer_file = run_crest_from_isomer_opt(Inchikey_photo_cat_E, isomer="E", state="_ON")
    conformer_file = run_crest_from_isomer_opt(Inchikey_photo_cat_Z, isomer="Z", state="_ON")
    conformer_file = run_crest_from_isomer_opt(Inchikey_photo_cat_E, isomer="E", state="")
    conformer_file = run_crest_from_isomer_opt(Inchikey_photo_cat_Z, isomer="Z", state="")



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-IE","--Inchikey_photo_cat_E", type=str)
    parser.add_argument("-IZ","--Inchikey_photo_cat_Z", type=str)
    args = parser.parse_args()
    main(args.Inchikey_photo_cat_E, args.Inchikey_photo_cat_Z)

