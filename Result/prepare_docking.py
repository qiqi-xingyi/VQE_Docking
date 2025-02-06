# --*-- conding:utf-8 --*--
# @Time : 11/14/24 9:26 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : prepare_docking.py

from files_tool import QuantumResult
import os
import shutil
from files_tool import DockingFilePreparer
from files_tool import Mol2Translator

def create_full_protein(path):
    quantum_result = QuantumResult(path)

    # Read xyz process_data
    quantum_result.read_xyz()

    # Adjust the coordinate scale to match the average Cα-Cα distance
    quantum_result.adjust_scale()

    # Write the Cα PDB file using the scaled coordinates
    quantum_result.write_ca_pdb()

    # Prepare the alignment file required by Modeller
    quantum_result.prepare_alignment()

    # Generate the full atomic model using Modeller
    quantum_result.generate_full_model()


def rename_and_move_file(current_path, new_name, target_directory):

    new_path = os.path.join(os.path.dirname(current_path), new_name)
    os.rename(current_path, new_path)
    print(f"File renamed to {new_name}")
    os.makedirs(target_directory, exist_ok=True)
    target_path = os.path.join(target_directory, new_name)
    shutil.move(new_path, target_path)
    print(f"File moved to {target_path}")

def create_docking_file(Q_res_path, AF3_res_path, id):

    preparer = DockingFilePreparer(Q_res_path)
    q_translated_pdb = f"./process_data/best_group/{id}/full_model_trans.pdb"
    preparer.prepare_pdbqt(translate=True, output_translated_file=q_translated_pdb)

    preparer = DockingFilePreparer(AF3_res_path)
    af_translated_pdb = f"process_data/best_group/{id}/fold_{id}/fold_model_trans.pdb"
    preparer.prepare_pdbqt(translate=True, output_translated_file=af_translated_pdb)
    preparer.prepare_pdbqt()


def process_single_protein(protein_id, base_dir="process_data/best_group"):

    xyz_file_path = os.path.join(base_dir, protein_id, f"{protein_id}.xyz")
    Q_res_path = os.path.join(base_dir, protein_id, "full_model.pdb")

    AF3_res_path = os.path.join(base_dir, protein_id, f"fold_{protein_id}", f"fold_{protein_id}_model_4.cif")

    ligand_dir = os.path.join(base_dir, protein_id, "PDBbind_data", protein_id)
    mol2_input = os.path.join(ligand_dir, f"{protein_id}_ligand.mol2")
    mol2_output = os.path.join(ligand_dir, f"{protein_id}_ligand_trans.mol2")

    create_full_protein(xyz_file_path)

    modeller_output = "protein_full.B99990001.pdb"
    rename_and_move_file(
        modeller_output,
        new_name="full_model.pdb",
        target_directory=os.path.join(base_dir, protein_id)
    )

    create_docking_file(Q_res_path, AF3_res_path, protein_id)

    translator = Mol2Translator(mol2_input, mol2_output)
    translator.prepare_translated_mol2()


def process_all_proteins(root_dir="process_data/best_group"):

    for subfolder in os.listdir(root_dir):
        folder_path = os.path.join(root_dir, subfolder)
        if os.path.isdir(folder_path):
            protein_id = subfolder
            print(f"\n[INFO] Now processing protein_id: {protein_id}")
            try:
                process_single_protein(protein_id, base_dir=root_dir)
            except Exception as e:
                print(f"[ERROR] Failed to process {protein_id}: {e}")

if __name__ == '__main__':

    process_all_proteins("process_data/best_group")