# --*-- conding:utf-8 --*--
# @Time : 11/15/24 1:06 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : get_ave_res.py

import os

class VinaDockingResultParser:
    def __init__(self, log_file):
        self.log_file = log_file

    def parse_results(self, start_line=36, end_line=47):
        """
        Parse docking results from the log file, focusing on a specified line range.

        Parameters:
            start_line (int): The starting line of the docking results in the log file.
            end_line (int): The ending line of the docking results in the log file.

        Returns:
            List of tuples with (affinity, rmsd_lower_bound, rmsd_upper_bound) for each mode.
        """
        results = []
        with open(self.log_file, 'r') as file:
            lines = file.readlines()[start_line - 1:end_line]
            for line in lines:
                parts = line.split()

                # Skip header or non-numeric lines
                if len(parts) < 4 or not parts[0].isdigit():
                    continue

                affinity = float(parts[1])
                rmsd_lb = float(parts[2])
                rmsd_ub = float(parts[3])
                results.append((affinity, rmsd_lb, rmsd_ub))
        return results

    def calculate_averages(self, results):
        """
        Calculate average affinity and RMSD values.

        Parameters:
            results (list of tuples): Each tuple contains (affinity, rmsd.l.b., rmsd.u.b.)

        Returns:
            avg_affinity (float): Average affinity (kcal/mol)
            avg_rmsd_lb (float): Average RMSD lower bound (Å)
            avg_rmsd_ub (float): Average RMSD upper bound (Å)
        """
        total_affinity = sum(r[0] for r in results)
        total_rmsd_lb = sum(r[1] for r in results)
        total_rmsd_ub = sum(r[2] for r in results)
        num_modes = len(results)

        avg_affinity = total_affinity / num_modes
        avg_rmsd_lb = total_rmsd_lb / num_modes
        avg_rmsd_ub = total_rmsd_ub / num_modes

        return avg_affinity, avg_rmsd_lb, avg_rmsd_ub


def process_all_results(base_dir, output_summary_file):
    quantum_results = []
    af2_results = []

    # Ensure the output directory exists
    output_dir = os.path.dirname(output_summary_file)
    os.makedirs(output_dir, exist_ok=True)

    for subdir, _, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".txt"):
                log_file = os.path.join(subdir, file)
                parser = VinaDockingResultParser(log_file)
                try:
                    results = parser.parse_results()
                    avg_affinity, avg_rmsd_lb, avg_rmsd_ub = parser.calculate_averages(results)

                    # Save results based on type
                    if "quantum" in log_file:
                        quantum_results.append((avg_affinity, avg_rmsd_lb, avg_rmsd_ub))
                    elif "af2" in log_file:
                        af2_results.append((avg_affinity, avg_rmsd_lb, avg_rmsd_ub))
                except Exception as e:
                    print(f"Error processing file {log_file}: {e}")

    # Calculate overall averages
    def calculate_overall_average(results):
        if not results:
            return 0, 0, 0
        avg_affinity = sum(r[0] for r in results) / len(results)
        avg_rmsd_lb = sum(r[1] for r in results) / len(results)
        avg_rmsd_ub = sum(r[2] for r in results) / len(results)
        return avg_affinity, avg_rmsd_lb, avg_rmsd_ub

    quantum_overall_avg = calculate_overall_average(quantum_results)
    af2_overall_avg = calculate_overall_average(af2_results)

    # Save results to summary file
    with open(output_summary_file, "w") as summary_file:
        summary_file.write("Quantum Results:\n")
        for i, (affinity, rmsd_lb, rmsd_ub) in enumerate(quantum_results, 1):
            summary_file.write(f"Trial {i}: Affinity = {affinity:.4f}, RMSD Lower Bound = {rmsd_lb:.4f}, RMSD Upper Bound = {rmsd_ub:.4f}\n")
        summary_file.write(f"Overall Average: Affinity = {quantum_overall_avg[0]:.4f}, RMSD Lower Bound = {quantum_overall_avg[1]:.4f}, RMSD Upper Bound = {quantum_overall_avg[2]:.4f}\n\n")

        summary_file.write("AF2 Results:\n")
        for i, (affinity, rmsd_lb, rmsd_ub) in enumerate(af2_results, 1):
            summary_file.write(f"Trial {i}: Affinity = {affinity:.4f}, RMSD Lower Bound = {rmsd_lb:.4f}, RMSD Upper Bound = {rmsd_ub:.4f}\n")
        summary_file.write(f"Overall Average: Affinity = {af2_overall_avg[0]:.4f}, RMSD Lower Bound = {af2_overall_avg[1]:.4f}, RMSD Upper Bound = {af2_overall_avg[2]:.4f}\n\n")

    print(f"Summary saved to {output_summary_file}")


if __name__ == "__main__":
    docking_output_path = "docking_output_low_energy_4zb8"

    base_dir = f"./process_data/{docking_output_path}"
    output_summary_file = f"./process_data/{docking_output_path}/summary_results.txt"
    process_all_results(base_dir, output_summary_file)