# --*-- conding:utf-8 --*--
# @Time : 11/12/24 1:15 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Quantum_res.py

import os
import numpy as np
from modeller import environ, log
from modeller.automodel import automodel, assess
from modeller import alignment, model

class QuantumResult:

    def __init__(self, xyz_file):
        self.xyz_file = xyz_file
        self.sequence = []
        self.coordinates = []
        self.scaled_coordinates = []

        # Determine the directory of the XYZ file
        self.base_dir = os.path.dirname(os.path.abspath(self.xyz_file))

        print(self.base_dir)

        # Define filenames based on the XYZ file name and base directory
        self.ca_pdb_filename = os.path.join(self.base_dir, 'ca_model.pdb')
        self.alignment_filename = os.path.join(self.base_dir, 'protein.ali')
        self.full_model_filename = os.path.join(self.base_dir, 'full_model.pdb')

        self.chain_id = "A"


    def read_xyz(self):
        with open(self.xyz_file, 'r') as f:
            lines = f.readlines()[2:]  # pass the first and second line (without any info)
            for line in lines:
                parts = line.strip().split()
                if len(parts) == 4:
                    res_name = parts[0]
                    x, y, z = map(float, parts[1:])
                    self.sequence.append(res_name)
                    self.coordinates.append((x, y, z))
                else:
                    print(f"Warning: Line '{line.strip()}' is malformed and will be skipped.")

    def adjust_scale(self, target_distance=3.8):

        def calculate_average_distance(coords):
            distances = []
            for i in range(len(coords) - 1):
                x1, y1, z1 = coords[i]
                x2, y2, z2 = coords[i + 1]
                distance = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
                distances.append(distance)
            return np.mean(distances)

        current_avg_distance = calculate_average_distance(self.coordinates)
        scale_factor = target_distance / current_avg_distance
        self.scaled_coordinates = [(x * scale_factor, y * scale_factor, z * scale_factor)
                                   for x, y, z in self.coordinates]
        print(f"Scale factor applied: {scale_factor:.4f}")

    def write_ca_pdb(self):
        res_name_map = {
            'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP',
            'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY',
            'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
            'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
            'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
        }

        with open(self.ca_pdb_filename, 'w') as f:
            for i, (res_name, (x, y, z)) in enumerate(zip(self.sequence, self.scaled_coordinates), start=1):
                res_name_3 = res_name_map.get(res_name.upper(), 'UNK')
                f.write(f"ATOM  {i:5d}  CA  {res_name_3} {self.chain_id}{i:4d}    "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
            f.write("END\n")
        print(f"PDB file '{self.ca_pdb_filename}' has been written.")

    def prepare_alignment(self, sequence_name='protein_full'):
        sequence_str = ''.join(self.sequence) + '*'
        start_residue = 1  # Starting residue number
        end_residue = len(self.sequence)  # Ending residue number

        with open(self.alignment_filename, 'w') as f:
            # Template entry (from Cα model)
            f.write(f">P1;ca_model\n")
            # Header line with 9 colons (10 fields), correct field order
            f.write(f"structureX:ca_model:{start_residue}:{self.chain_id}:{end_residue}:{self.chain_id}::::\n")
            f.write(sequence_str + "\n\n")
            # Target entry (full model)
            f.write(f">P1;{sequence_name}\n")
            # Header line with 9 colons (10 fields), correct field order
            f.write(f"sequence:{sequence_name}:{start_residue}:{self.chain_id}:{end_residue}:{self.chain_id}::::\n")
            f.write(sequence_str + "\n")
        print(f"Alignment file '{self.alignment_filename}' has been prepared.")

    def generate_full_model(self):

        log.none()  # Set to log.verbose() for detailed output
        env = environ()

        # Set the directory where Modeller will look for atom process_data
        env.io.atom_files_directory = [self.base_dir]

        # Read the alignment file
        aln = alignment(env)
        mdl = model(env, file=self.ca_pdb_filename, model_segment=(f'FIRST:{self.chain_id}', f'LAST:{self.chain_id}'))
        aln.append_model(mdl, align_codes='ca_model', atom_files=self.ca_pdb_filename)
        aln.append(file=self.alignment_filename, align_codes='protein_full')
        aln.align2d()

        # Define the automodel class
        class MyModel(automodel):
            def special_patches(self, aln):
                self.rename_segments(segment_ids=['A'])

        # Create and build the model
        a = MyModel(
            env,
            alnfile=self.alignment_filename,
            knowns='ca_model',
            sequence='protein_full',
            assess_methods=(assess.DOPE, assess.GA341),
        )
        a.run_code = '9999'
        a.starting_model = 1
        a.ending_model = 1
        a.make()

        # generated_model = os.path.join(self.base_dir, 'A_Beta_42_final.pdb')
        # if os.path.exists(generated_model):
        #     os.rename(generated_model, self.full_model_filename)
        #     print(f"Full atomic model '{self.full_model_filename}' has been generated.")
        # else:
        #     print("Error: The expected output model file was not found.")



