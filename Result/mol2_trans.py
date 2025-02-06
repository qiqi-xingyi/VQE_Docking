# --*-- conding:utf-8 --*--
# @Time : 11/15/24 4:39 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : mol2_trans.py

from files_tool import Mol2Translator

if __name__ == '__main__':

    protein_id = '1a9m'

    # Initialize the translator with input and output file paths
    translator = Mol2Translator(f"./process_data/best_group/{protein_id}/PDBbind_data/{protein_id}/{protein_id}_ligand.mol2",
                                f"./process_data/best_group/{protein_id}/PDBbind_data/{protein_id}/{protein_id}_ligand_trans.mol2")

    # Perform the translation
    translator.prepare_translated_mol2()
