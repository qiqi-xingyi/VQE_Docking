# --*-- conding:utf-8 --*--
# @Time : 1/17/25 11:14 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Main.py

import os
import time
from typing import List, Tuple

from Protein_Folding import Peptide
from Protein_Folding.interactions.miyazawa_jernigan_interaction import MiyazawaJerniganInteraction
from Protein_Folding.penalty_parameters import PenaltyParameters
from Protein_Folding.protein_folding_problem import ProteinFoldingProblem
from qiskit_ibm_runtime import QiskitRuntimeService
from Qiskit_VQE import VQE5, StateCalculator

def predict_protein_structure(
    main_chain_sequence: str,
    protein_id: str,
    service: QiskitRuntimeService,
    max_iter: int = 150
):
    """
    Predicts a protein structure using a quantum VQE workflow and saves the results.

    :param main_chain_sequence: The main chain amino acid sequence (single-letter representation).
    :param protein_id: The identifier/name for the protein (used for output directories and files).
    :param service: An instance of QiskitRuntimeService for submitting quantum jobs.
    :param max_iter: Maximum iteration count for VQE optimization, default is 150.
    """
    print(f"Starting prediction for protein: {protein_id}, sequence: {main_chain_sequence}")

    # Initialize peptide with an empty side-chain sequence
    side_chain_sequences = ['' for _ in main_chain_sequence]
    peptide = Peptide(main_chain_sequence, side_chain_sequences)

    # Define interactions and penalty parameters
    mj_interaction = MiyazawaJerniganInteraction()
    penalty_terms = PenaltyParameters(10, 10, 10)

    # Create the protein folding problem
    protein_folding_problem = ProteinFoldingProblem(peptide, mj_interaction, penalty_terms)
    hamiltonian = protein_folding_problem.qubit_op()

    qubit_count = hamiltonian.num_qubits + 5
    print(f"Number of qubits required: {qubit_count}")

    # Initialize the VQE solver
    vqe_instance = VQE5(
        service=service,
        hamiltonian=hamiltonian,
        min_qubit_num=qubit_count,
        maxiter=max_iter
    )

    # Run the VQE and obtain results
    energy_list, final_result, ansatz, top_results = vqe_instance.run_vqe()

    # Save the energy list to a file
    output_energy_path = os.path.join("Result", "process_data", "best_group", protein_id, "System_Energy")
    os.makedirs(output_energy_path, exist_ok=True)
    with open(os.path.join(output_energy_path, f"energy_list_{protein_id}.txt"), 'w') as file:
        file.writelines(f"{energy}\n" for energy in energy_list)

    # Calculate and save the probability distribution
    state_calculator = StateCalculator(service, qubit_count, ansatz)
    probability_distribution = state_calculator.get_probability_distribution(final_result)

    output_prob_path = os.path.join("Result", "process_data", "best_group", protein_id, "Prob_distribution")
    os.makedirs(output_prob_path, exist_ok=True)
    with open(os.path.join(output_prob_path, "prob_distribution.txt"), 'w') as file:
        for state, probability in probability_distribution.items():
            file.write(f"{state}: {probability}\n")

    # Interpret the probability distribution and save the predicted structure
    protein_result = protein_folding_problem.interpret(probability_distribution)
    output_dir = os.path.join("Result", "process_data", "best_group", protein_id)
    protein_result.save_xyz_file(name=protein_id, path=output_dir)
    print("Protein structure saved as .xyz file")

    # Save the top-ranked results
    for rank, (energy_val, best_params) in enumerate(top_results, start=1):
        print(f"Top {rank} best energy = {energy_val}")

        best_prob_dist = state_calculator.get_probability_distribution(best_params)
        best_protein_result = protein_folding_problem.interpret(best_prob_dist)
        best_protein_result.save_xyz_file(
            name=f"{protein_id}_top_{rank}",
            path=output_dir
        )
        print(f"Protein structure for top {rank} result has been saved.")

    print(f"Finished processing: {protein_id}\n")

def main():
    """Main function to predict protein structures for a list of proteins."""
    # Initialize the Qiskit Runtime Service
    service = QiskitRuntimeService(
        channel='ibm_quantum',
        instance='YOUR_INSTANCE_NAME',  # Replace with your real instance
        token='YOUR_API_TOKEN'          # Replace with your real token
    )

    # List of proteins to process
    protein_list: List[Tuple[str, str]] = [
        ("DGKMKGLAF", "1qin"),
        ("IHGIGGFI", "1a9m"),
        ("KSIVDSGTTNLR", "1fkn"),
        ("NNLGTIAKSGT", "3b26"),
        ("GAVEDGATMTFF", "2xxx"),
        ("DWGGM", "3ans"),
        ("YAGYS", "6mu3")
    ]

    # Log execution times for each protein
    log_file_path = "execution_time_log.txt"
    with open(log_file_path, 'w') as log_file:
        log_file.write("Protein_ID\tExecution_Time(s)\n")

        for sequence, protein_name in protein_list:
            start_time = time.time()

            predict_protein_structure(
                main_chain_sequence=sequence,
                protein_id=protein_name,
                service=service,
                max_iter=150
            )

            execution_time = time.time() - start_time
            log_file.write(f"{protein_name}\t{execution_time:.2f}\n")

if __name__ == '__main__':
    main()