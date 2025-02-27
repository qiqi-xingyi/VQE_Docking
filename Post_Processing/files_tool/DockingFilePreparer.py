# --*-- conding:utf-8 --*--
# @Time : 11/14/24 3:08 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : DockingFilePreparer.py

import os
import subprocess
import numpy as np
from Bio.PDB import PDBParser, PDBIO, MMCIFParser

class DockingFilePreparer:
    def __init__(self, input_file):
        self.input_file = input_file
        self.file_extension = os.path.splitext(input_file)[-1].lower()
        self.output_pdbqt = input_file.replace(self.file_extension, ".pdbqt")

    # def prepare_pdbqt(self):
    #     # Step 1: Check if Open Babel is installed
    #     if not self._is_tool_available("obabel"):
    #         raise EnvironmentError("Open Babel is not installed or not in PATH.")
    #
    #     # Step 2: Convert CIF or PDB to PDBQT with Open Babel
    #     if self.file_extension in [".cif", ".pdb"]:
    #         print(f"Converting {self.file_extension.upper()} to PDBQT with Open Babel...")
    #         self._convert_to_pdbqt_with_obabel()
    #         print(f"Conversion complete. Output file: {self.output_pdbqt}")
    #     else:
    #         raise ValueError("Unsupported file format. Please provide a .cif or .pdb file.")
    def prepare_pdbqt(self, translate=False, output_translated_file=None):
        """
        Prepare a PDBQT file.

        Parameters:
        - translate (bool): Whether to translate the molecule to the origin.
        - output_translated_file (str): Path to save the translated PDB file.
        """
        # Step 1: Check if Open Babel is available
        if not self._is_tool_available("obabel"):
            raise EnvironmentError("Open Babel is not installed or not in PATH.")

        # Step 2: If translation is requested, translate the molecule to the origin
        if translate:
            if not output_translated_file:
                # Generate a translated file path by appending '_translated.pdb'
                output_translated_file = self.input_file.replace(self.file_extension, "_translated.pdb")
            print(f"Translating molecule to the origin...")
            self.translate_to_origin(self.input_file, output_translated_file)
            input_for_conversion = output_translated_file
        else:
            input_for_conversion = self.input_file

        # Step 3: Convert CIF or PDB to PDBQT using Open Babel
        if self.file_extension in [".cif", ".pdb"]:
            print(f"Converting {self.file_extension.upper()} to PDBQT with Open Babel...")
            self._convert_to_pdbqt_with_obabel(input_for_conversion, self.output_pdbqt)
            print(f"Conversion complete. Output file: {self.output_pdbqt}")
        else:
            raise ValueError("Unsupported file format. Please provide a .cif or .pdb file.")

    def _is_tool_available(self, tool_name):
        """Check if a tool is available on the system."""
        return subprocess.call(f"type {tool_name}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

    def _convert_to_pdbqt_with_obabel(self, input_file, output_pdbqt):
        """Use Open Babel to add hydrogens, assign Gasteiger charges, and output as PDBQT."""
        command = [
            "obabel", input_file, "-O", output_pdbqt,
            "--partialcharge", "gasteiger", "-h", "-xr"
        ]
        subprocess.run(command, check=True)

    def calculate_center_of_mass(self, structure_file):
        """
        Calculate the center of mass of all atoms in a PDB or CIF file.

        Parameters:
        - structure_file (str): Path to the PDB or CIF file.

        Returns:
        - numpy.ndarray: Coordinates of the center of mass.
        """
        if self.file_extension == ".pdb":
            parser = PDBParser(QUIET=True)
        elif self.file_extension == ".cif":
            parser = MMCIFParser(QUIET=True)
        else:
            raise ValueError("Only .pdb and .cif file formats are supported.")

        try:
            structure = parser.get_structure('structure', structure_file)
        except Exception as e:
            raise ValueError(f"Failed to parse file {structure_file}: {e}")

        atom_coords = []
        for atom in structure.get_atoms():
            coord = atom.get_coord()
            atom_coords.append(coord)

        if not atom_coords:
            raise ValueError(f"No atom coordinates found in file {structure_file}.")

        atom_coords = np.array(atom_coords)
        center_of_mass = atom_coords.mean(axis=0)
        print(f"Center of mass: X={center_of_mass[0]:.3f}, Y={center_of_mass[1]:.3f}, Z={center_of_mass[2]:.3f}")
        return center_of_mass

    def translate_to_origin(self, input_file, output_file):
        """
        Translate the molecule so that its center of mass is at the origin.

        Parameters:
        - input_file (str): Path to the original PDB or CIF file.
        - output_file (str): Path to save the translated PDB file.
        """
        if self.file_extension == ".pdb":
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('structure', input_file)
        elif self.file_extension == ".cif":
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure('structure', input_file)
        else:
            raise ValueError("Only .pdb and .cif file formats are supported.")

        center_of_mass = self.calculate_center_of_mass(input_file)
        translation = -center_of_mass  # Calculate translation vector

        # Translate all atoms
        for atom in structure.get_atoms():
            atom.coord += translation

        # Save the translated structure
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_file)
        print(f"Molecule has been translated to the origin and saved as {output_file}")