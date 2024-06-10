#  gaussian_info_extractor.py
#  Created by Ayesh Madushanka
#  Created on: May 16, 2024
import os
import pandas as pd


class GaussianInfoExtractor:
    
    """
    Class to Extract Qm parameters from Gaussian16 input file

    Args:
    current_dir_path (str): path for optimized compounds(each are inside a folder named Zinc_id)
    current_dir (str): Specific molecule's Zinc_id (folder names)
    
    """
    
    def __init__(self, current_dir_path, current_dir):
        self.current_dir_path = current_dir_path
        self.current_dir = current_dir
        self.output_data = None
        

    # read output files
    def output_file_reader(self) -> str:
        os.chdir(self.current_dir_path)
        gout_name = f"{self.current_dir}_run.out"
        gout_path = os.path.join(self.current_dir_path, gout_name)
        try:
            with open(gout_path, 'r') as f:
                lines = f.readlines()
            self.output_data = (lines, gout_name)
        except FileNotFoundError:
            print(f"Gaussian out file {gout_name} not found in {self.current_dir_path}, skipping g_out data.")
            
        return self.output_data
            
    
    # find frequencies and report imaginary frequencies
    def imaginary_freq_checker(self) -> list:
        lines, out_file = self.output_file_reader()
        target_word = "Frequencies"
        for line in lines:
            row = line.strip()
            
            if target_word in row:
                freq_line = row
                values_str = freq_line.split()
                values_str = values_str[2:]
                values = [float(value.strip()) for value in values_str]
                for value in values:
                    
                    if value < 0:
                        return self.current_dir
                  
    
    def intial_coord_extract(self) -> list:
        ini_xyz = []
        xyz_file = f"{self.current_dir}.xyz"
        with open(xyz_file, "r") as file:
            first_line = file.readline().rstrip()
        length = int(first_line)
        ini_coord = []
        ini_coordinates = "Symbolic Z-matrix:"
        lines, gauss_out = self.output_file_reader()
        
        for i, line in enumerate(lines):
            row = line.strip()
            
            if ini_coordinates in row:
                
                for index in range(length):
                    
                    initial_xyz = lines[index+2+i].strip().split()
                    ini_xyz.append(initial_xyz)
        return ini_xyz, length


# Extract QM parameters ffrom gaussian output file using keywords
    def qm_info_collector(self) -> list:
        opt = []
        charges = []
        lines, gauss_out = self.output_file_reader()
        print(self.current_dir)
        start_reading = False
        start_line = '#P Geom=AllCheck Guess=TCheck SCRF=Check GenChk RB3LYP/6-31G(2df,p) Fr'
        opt_coordinates = 'Standard orientation:'
        elec_E = 'SCF Done'
        isotropic = 'Isotropic polarizability'
        HOMO = 'Alpha  occ. eigenvalues'
        Muliliken_charges = 'Mulliken charges:'
        electro = 'Electronic spatial extent (au):'
        vibrational = 'Zero-point vibrational energy'
        dipol = 'Dipole moment (field-independent basis, Debye):'
        Internal_energy = 'Sum of electronic and thermal Energies='
        start_reading = False
        initial_coordinates, length = self.intial_coord_extract()  
        for i, line in enumerate(lines):
            row = line.strip()
            
            if start_line in row:
                start_reading = True
                
            if not start_reading:
                continue
            
            if opt_coordinates in row: 
                for num1 in range(length): 
                    num1 = num1 + 5
                    num1 = i + num1 
                    xyz = lines[num1].strip()
                    xyz_coord = xyz.split()[3:6]
                    opt.append(xyz_coord)
                    
            elif elec_E in row:
                tot_ele_E_0K = row.split()[4]
                
            elif HOMO in row:
                last_HOMO_line = row.split()[-1]
                first_LUMO_line = lines[i+1].strip().split()[4]
                
            elif isotropic in row:
                polarizability = row.split()[5]
                
            elif Muliliken_charges in row:
                for num in range(length): 
                    num = num + 2
                    num = i + num 
                    Mili = lines[num].strip().split()[2]
                    charges.append(Mili)
                    
            elif dipol in row:
                electronic = lines[i-2].strip().split()[5]
                dipol_moment = lines[i + 1].strip().split()[7]

            elif vibrational in row:
                zero_point = lines[i+1].strip().split()[0]
                rotational = lines[i-2].strip().split()
                rotational1 = rotational[3]
                rotational2 = rotational[4]
                rotational3= rotational[5]
                
            elif Internal_energy in row:
                Internal_energy_295K = row.split()[6]
                Enthalpy = lines[i + 1].strip().split()[6]
                Free_E = lines[i + 2].strip().split()[7]
                cv_entropy = lines[i + 6].strip()
                cv = cv_entropy.split()[2]
                entropy = cv_entropy.split()[3]
                
        return [initial_coordinates,
                opt, 
                charges, 
                tot_ele_E_0K, 
                last_HOMO_line, 
                first_LUMO_line,
                polarizability,
                electronic,
                dipol_moment,
                zero_point,
                rotational1,
                rotational2,
                rotational3,
                Internal_energy_295K,
                Enthalpy,
                Free_E,
                cv,
                entropy
                ]
    
      
    # making CSVs 
    def making_csvs(self) -> list:
        n = self.qm_info_collector()
        qm_parameters = (n[3:18])
        init_coordinates = n[0]
        init_name = ['atomic num', 'ini_x', 'init_y', 'init_z']
        df_init = pd.DataFrame(init_coordinates, columns=init_name)
        final_name = ['final_x', 'final_y', 'final_z']
        df_final = pd.DataFrame(n[1], columns=final_name)
        charge_name = ['Charge']
        df_charges = pd.DataFrame(n[2], columns = charge_name)
        df_coordinates = pd.concat([df_init, df_final, df_charges], axis=1)
        return qm_parameters, df_coordinates