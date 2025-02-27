# --*-- conding:utf-8 --*--
# @Time : 11/10/24 4:03 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : dock_scorer.py

import random
from Autodock_tool import AutoDockDocking
import os

if __name__ == "__main__":

    num_trials = 20
    docking_output_path = 'docking_output_4zb8'
    seed_log_file = f"./process_data/{docking_output_path}/seed_log.txt"

    chain_id = '4zb8'
    protein_id = '4zb8'

    os.makedirs(os.path.dirname(seed_log_file), exist_ok=True)

    with open(seed_log_file, "w") as seed_file:
        seed_file.write("Trial\tSeed\n")

    for trial in range(1, num_trials + 1):
        # allowed_seeds = [10, 20, 30, 40]
        # seed = random.choice(allowed_seeds)

        seed = random.randint(1, 100000)
        # seed = 2
        print(f"Trial {trial} - Using seed: {seed}")

        with open(seed_log_file, "a") as seed_file:
            seed_file.write(f"{trial}\t{seed}\n")

        print(f'Quantum Trial {trial}:')
        quantum_output_dir = f"./process_data/{docking_output_path}/quantum_trial_{trial}"
        os.makedirs(quantum_output_dir, exist_ok=True)

        # Quantum docking
        docking = AutoDockDocking(
            f"./process_data/{chain_id}/full_model.pdbqt",
            f"./process_data/{chain_id}/PDBbind_data/{protein_id}/{protein_id}_ligand_trans.mol2",
            quantum_output_dir,
            f"docking_log_trial_{trial}.txt",
            seed
        )
        docking.run_docking()
        print(f"Quantum Trial {trial} docking completed.")
        print("\n")


        print(f'AF2 Trial {trial}:')
        af2_output_dir = f"./process_data/{docking_output_path}/af2_trial_{trial}"
        os.makedirs(af2_output_dir, exist_ok=True)

        # AF2 docking
        docking = AutoDockDocking(
            f"./process_data/{chain_id}/alphafold_predicted/fold_model_3.pdbqt",
            f"./process_data/{chain_id}/PDBbind_data/{protein_id}/{protein_id}_ligand_trans.mol2",
            af2_output_dir,
            f"af2_docking_log_trial_{trial}.txt",
            seed
        )
        docking.run_docking()
        print(f"AF2 Trial {trial} docking completed.")
        print("\n")


