# --*-- conding:utf-8 --*--
# @Time : 11/15/24 4:33 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Mol_trans.py

import os


class Mol2Translator:
    """
    A class to translate the geometric center of a molecule in a MOL2 file to the origin.
    """

    def __init__(self, input_mol2, output_mol2):
        """
        Initialize the Mol2Translator class.

        Parameters:
        - input_mol2 (str): Path to the input MOL2 file.
        - output_mol2 (str): Path to save the translated MOL2 file.
        """
        self.input_mol2 = input_mol2
        self.output_mol2 = output_mol2
        self.atoms = []  # List to store atom information

    def parse_mol2(self):
        """
        Parse the MOL2 file and extract atom information.
        """
        if not os.path.exists(self.input_mol2):
            raise FileNotFoundError(f"The file {self.input_mol2} does not exist.")

        with open(self.input_mol2, 'r') as file:
            lines = file.readlines()

        atom_section = False
        for line in lines:
            line = line.strip()
            if line.startswith("@<TRIPOS>ATOM"):
                atom_section = True
                continue
            elif line.startswith("@<TRIPOS>") and atom_section:
                # End of ATOM section
                break

            if atom_section:
                parts = line.split()
                if len(parts) < 6:
                    continue  # Skip malformed lines
                atom_id = parts[0]
                atom_name = parts[1]
                try:
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])
                except ValueError:
                    raise ValueError(f"Invalid coordinates for atom {atom_id} in file {self.input_mol2}.")
                atom_type = parts[5]
                # Capture any additional columns beyond the first 6
                additional = parts[6:] if len(parts) > 6 else []
                self.atoms.append({
                    'atom_id': atom_id,
                    'atom_name': atom_name,
                    'x': x,
                    'y': y,
                    'z': z,
                    'atom_type': atom_type,
                    'additional': additional
                })

        if not self.atoms:
            raise ValueError(f"No atoms found in the MOL2 file {self.input_mol2}.")

    def calculate_geometric_center(self):
        """
        Calculate the geometric center of the molecule.

        Returns:
        - tuple: (center_x, center_y, center_z)
        """
        sum_x = sum(atom['x'] for atom in self.atoms)
        sum_y = sum(atom['y'] for atom in self.atoms)
        sum_z = sum(atom['z'] for atom in self.atoms)
        num_atoms = len(self.atoms)

        center_x = sum_x / num_atoms
        center_y = sum_y / num_atoms
        center_z = sum_z / num_atoms

        print(f"Geometric center: X={center_x:.3f}, Y={center_y:.3f}, Z={center_z:.3f}")
        return (center_x, center_y, center_z)

    def translate_atoms(self, center):
        """
        Translate all atoms so that the geometric center is at the origin.

        Parameters:
        - center (tuple): The geometric center coordinates (center_x, center_y, center_z).
        """
        center_x, center_y, center_z = center
        for atom in self.atoms:
            atom['x'] -= center_x
            atom['y'] -= center_y
            atom['z'] -= center_z

    def write_translated_mol2(self):
        """
        Write the translated atoms to a new MOL2 file, preserving the original structure.
        """
        with open(self.input_mol2, 'r') as file:
            lines = file.readlines()

        with open(self.output_mol2, 'w') as file:
            atom_section = False
            for line in lines:
                stripped_line = line.strip()
                if stripped_line.startswith("@<TRIPOS>ATOM"):
                    atom_section = True
                    file.write(line)
                    continue
                elif stripped_line.startswith("@<TRIPOS>") and atom_section:
                    atom_section = False
                    file.write(line)
                    continue

                if atom_section:
                    parts = line.split()
                    if len(parts) < 6:
                        file.write(line)
                        continue  # Write malformed lines as is
                    atom_id = parts[0]
                    # Find the corresponding atom in self.atoms
                    atom = next((a for a in self.atoms if a['atom_id'] == atom_id), None)
                    if atom is None:
                        file.write(line)
                        continue  # Write lines without matching atoms as is
                    # Reconstruct the atom line with translated coordinates
                    # Preserve formatting: align columns as in original
                    # Assuming original formatting uses spaces, not tabs
                    new_line = f"{atom['atom_id']:>6} {atom['atom_name']:<10} {atom['x']:>8.3f} {atom['y']:>8.3f} {atom['z']:>8.3f} {atom['atom_type']}"
                    # Append any additional columns beyond the atom type
                    if atom['additional']:
                        additional = ' ' + ' '.join(atom['additional'])
                        new_line += additional
                    new_line += '\n'
                    file.write(new_line)
                else:
                    file.write(line)

        print(f"Translated MOL2 file saved as {self.output_mol2}.")

    def prepare_translated_mol2(self):
        """
        Perform the full preparation: parse, calculate geometric center, translate, and write new MOL2.
        """
        self.parse_mol2()
        center = self.calculate_geometric_center()
        self.translate_atoms(center)
        self.write_translated_mol2()

