# --*-- conding:utf-8 --*--
# @Time : 11/10/24 3:40 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : utils.py

import os
import subprocess


class AutoDockDocking:
    def __init__(self, receptor_pdbqt, ligand_mol2, output_dir="./process_data/docking_output_6mu3", log_file_name='docking_log.txt', seed=2):
        self.receptor_pdbqt = receptor_pdbqt
        self.ligand_mol2 = ligand_mol2
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        self.log_file_name = log_file_name
        self.seed = seed


    def run_docking(self):
        """Runs docking using AutoDock Vina and outputs results."""
        # Step 1: Convert ligand from MOL2 to PDBQT format
        ligand_pdbqt = os.path.join(self.output_dir, "ligand.pdbqt")
        self._convert_mol2_to_pdbqt(self.ligand_mol2, ligand_pdbqt)

        # Step 2: Calculate center of mass for the docking box
        center_x, center_y, center_z = self.calculate_center_of_mass(self.ligand_mol2)

        # Step 3: Run AutoDock Vina docking with calculated box center
        output_file = os.path.join(self.output_dir, "docking_output.pdbqt")
        log_file = os.path.join(self.output_dir, self.log_file_name)

        print("Running AutoDock Vina docking...")
        self._run_vina(self.receptor_pdbqt, ligand_pdbqt, output_file, log_file, center_x, center_y, center_z)

        print(f"Docking complete. Results saved in {output_file}")
        return output_file, log_file

    def _convert_mol2_to_pdbqt(self, mol2_file, pdbqt_file):
        """Convert MOL2 file to PDBQT format using Open Babel."""
        if not self._is_tool_available("obabel"):
            raise EnvironmentError("Open Babel is not installed or not in PATH.")

        command = ["obabel", mol2_file, "-O", pdbqt_file, "--partialcharge", "gasteiger", "-h"]
        subprocess.run(command, check=True)
        print(f"Converted {mol2_file} to PDBQT format as {pdbqt_file}.")


    def calculate_center_of_mass(self, mol2_file):
        """Calculate the center of mass of the ligand from a MOL2 file."""
        atom_coordinates = []
        with open(mol2_file, 'r') as file:
            atom_section = False
            for line in file:
                # Check for the start of the ATOM section
                if line.startswith('@<TRIPOS>ATOM'):
                    atom_section = True
                    continue
                elif line.startswith('@<TRIPOS>') and atom_section:
                    # End of ATOM section
                    break

                # Read atom coordinates in the ATOM section
                if atom_section:
                    parts = line.split()
                    if len(parts) >= 5:
                        x, y, z = map(float, parts[2:5])
                        atom_coordinates.append((x, y, z))

        # Calculate the center of mass
        if atom_coordinates:
            x_coords, y_coords, z_coords = zip(*atom_coordinates)
            center_x = sum(x_coords) / len(x_coords)
            center_y = sum(y_coords) / len(y_coords)
            center_z = sum(z_coords) / len(z_coords)
            print(f"Calculated center of mass: X={center_x:.3f}, Y={center_y:.3f}, Z={center_z:.3f}")
            return center_x, center_y, center_z
        else:
            raise ValueError("No atom coordinates found in the MOL2 file.")

    def _run_vina(self, receptor_pdbqt, ligand_pdbqt, output_file, log_file, center_x, center_y, center_z):
        """Run AutoDock Vina docking with calculated box center."""
        size_x, size_y, size_z = 18, 18, 18  # Adjust box size as needed
        command = [
            "vina", "--receptor", receptor_pdbqt,
            "--ligand", ligand_pdbqt,
            "--out", output_file,
            "--center_x", str(center_x), "--center_y", str(center_y), "--center_z", str(center_z),
            "--size_x", str(size_x), "--size_y", str(size_y), "--size_z", str(size_z),
            "--exhaustiveness", "16",
            "--seed",f"{self.seed}"
        ]

        # Run the command and capture output
        try:
            with open(log_file, "w") as log:
                subprocess.run(command, stdout=log, stderr=log, check=True)
        except subprocess.CalledProcessError as e:
            print("Error: AutoDock Vina failed with the following command:")
            print(" ".join(command))
            print(f"Exit status: {e.returncode}")

    def _is_tool_available(self, tool_name):
        """Check if a tool is available on the system."""
        return subprocess.call(f"type {tool_name}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

    def parse_docking_results(self, log_file):
        """Parse AutoDock Vina log file to extract docking scores."""
        scores = []
        with open(log_file, "r") as file:
            for line in file:
                if "REMARK VINA RESULT:" in line:
                    parts = line.strip().split()
                    score = float(parts[3])  # Extract score (binding affinity in kcal/mol)
                    scores.append(score)
        print(f"Extracted docking scores: {scores}")
        return scores
