# --*-- conding:utf-8 --*--
# @Time : 11/10/24 2:07 PM
# @Author : Yuqi Zhang=
# @Email : yzhan135@kent.edu
# @File : create_full_protein.py

from files_tool import QuantumResult

if __name__ == '__main__':

    xyz_file_path = 'process_data/6mu3/6mu3_L.xyz'

    quantum_result = QuantumResult(xyz_file_path)

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






