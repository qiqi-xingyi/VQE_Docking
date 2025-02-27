# --*-- conding:utf-8 --*--
# @Time : 1/14/25 8:36 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : docking_others.py

import random
from Autodock_tool import AutoDockDocking
import os
import glob

if __name__ == "__main__":
    docking_output_path = 'docking_output_low_energy_6mu3'
    seed_log_file = f"./process_data/{docking_output_path}/seed_log.txt"
    pdbqt_folder = "./process_data/create_structure/6mu3_L"  # 设置存放 .pdbqt 文件的文件夹路径

    protein_id = '6mu3'

    # 确保输出目录存在
    os.makedirs(os.path.dirname(seed_log_file), exist_ok=True)

    # 记录种子的日志文件
    with open(seed_log_file, "w") as seed_file:
        seed_file.write("PDBQT_File\tSeed\n")

    # 遍历 .pdbqt 文件夹中的所有文件
    pdbqt_files = glob.glob(os.path.join(pdbqt_folder, "*.pdbqt"))

    for pdbqt_file in pdbqt_files:
        file_name = os.path.basename(pdbqt_file)
        chain_id = os.path.splitext(file_name)[0]  # 获取文件名（去掉后缀）
        seed = 69652  # 设置随机种子

        print(f"Processing file: {file_name} - Using seed: {seed}")

        # 记录当前文件和种子到日志文件
        with open(seed_log_file, "a") as seed_file:
            seed_file.write(f"{file_name}\t{seed}\n")

        quantum_output_dir = f"./process_data/{docking_output_path}/{chain_id}"
        os.makedirs(quantum_output_dir, exist_ok=True)

        # Quantum docking
        docking = AutoDockDocking(
            pdbqt_file,
            f"./process_data/6mu3/PDBbind_data/{protein_id}/{protein_id}_ligand_trans.mol2",
            quantum_output_dir,
            f"docking_log_{chain_id}.txt",
            seed
        )
        docking.run_docking()
        print(f"Docking for {file_name} completed.")
        print("\n")
