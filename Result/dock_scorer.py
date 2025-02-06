# --*-- conding:utf-8 --*--
# @Time : 11/10/24 4:03 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : dock_scorer.py

import random
from Autodock_tool import AutoDockDocking
import os

def docking_test(protein_id):

    print(f'Protein ID:{protein_id}')

    num_trials = 20
    docking_output_path = f'docking_output/docking_output_{protein_id}'
    seed_log_file = f"./process_data/{docking_output_path}/seed_log.txt"

    os.makedirs(os.path.dirname(seed_log_file), exist_ok=True)

    with open(seed_log_file, "w") as seed_file:
        seed_file.write("Trial\tSeed\n")

    for trial in range(1, num_trials + 1):
        # allowed_seeds = [10, 20, 30, 40]
        # seed = random.choice(allowed_seeds)

        seed = random.randint(1, 100000)
        print(f"Trial {trial} - Using seed: {seed}")

        with open(seed_log_file, "a") as seed_file:
            seed_file.write(f"{trial}\t{seed}\n")

        print(f'Quantum Trial {trial}:')
        quantum_output_dir = f"./process_data/{docking_output_path}/quantum_trial_{trial}"
        os.makedirs(quantum_output_dir, exist_ok=True)

        # Quantum docking
        docking = AutoDockDocking(
            f"./process_data/new_group/{protein_id}/full_model.pdbqt",
            f"./process_data/new_group/{protein_id}/PDBbind_data/{protein_id}/{protein_id}_ligand_trans.mol2",
            quantum_output_dir,
            f"docking_log_trial_{trial}.txt",
            seed
        )
        docking.run_docking()
        print(f"Quantum Trial {trial} docking completed.")
        print("\n")

        print(f'AF3 Trial {trial}:')
        af3_output_dir = f"./process_data/{docking_output_path}/af3_trial_{trial}"
        os.makedirs(af3_output_dir, exist_ok=True)

        # AF3 docking
        docking = AutoDockDocking(
            f"./process_data/new_group/{protein_id}/fold_{protein_id}/fold_{protein_id}_model_4.pdbqt",
            f"./process_data/new_group/{protein_id}/PDBbind_data/{protein_id}/{protein_id}_ligand_trans.mol2",
            af3_output_dir,
            f"af3_docking_log_trial_{trial}.txt",
            seed
        )
        docking.run_docking()
        print(f"AF3 Trial {trial} docking completed.")
        print("\n")

if __name__ == "__main__":

    # protein_ids = ['1a9m','1qin','2xxx','3b26','6ugp']
    protein_ids = ['3ans']

    for obj in protein_ids:
        docking_test(obj)
