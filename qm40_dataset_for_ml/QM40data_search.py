# QM40data_search.py
#  Created by Ayesh Madushanka
#  Created on: June 4, 2024

import os
import pandas as pd
from rdkit import Chem


# help function generation
def help_func():
    print("Welcome to the QM40 dataset data seraching modular!!!!")
    print("QM40 dataset includes a main dataset which has 150K molecules") 
    print("with calculated QM parameters described in below")
    print("--------------------------------------------------------------")
    QM_parameters1 = """
    Property              Description
    1 Zinc_id             Specific Id for the molecules
    2 smile               SMILE representation of the molecule
    3 Internal_E(0K)      Internal energy at 0 K
    4 HOMO                Energy of HOMO
    5 LUMO                Energy of LUMO
    6 HL_gap              Energy difference of (HOMO - LUMO)
    7 Polarizability      Isotropic polarizability
    8 spatial_extent      Electronic spatial extent
    9 dipol_mom           Dipole moment
    10 ZPE                Zero point energy
    11 rot1               Rotational constant1
    12 rot2               Rotational constant2
    13 rot3               Rotational constant3
    14 Inter_E(298)       Internal energy at 298.15 K
    15 Enthalpy           Enthalpy at 298.15 K
    16 Free_E             Free energy at 298.15 K
    17 CV                 Heat capacity at 298.15 K
    18 Entropy            Entropy at 298.15 K
    """ 
    print(QM_parameters1)
    print("Any of the property name can be used to find the QM parameters")
    print("of a perticular molecule")
    print("\n")
    print("If you need to access LmodA or optimized geometry of a molecule")
    print("Use get_geo_Lmod_data function with additional key word 'lmod'")
    print("lmod = 'lmod': LmodA dataset")
    print("lmod = 'geo': Optimized Geometry dataset with Mulliken charges")


class QM40DataExtractor:
    def __init__(self, smile_or_ZN_id, QM_parameter, lmod = None):
        self.smile_or_ZN_id = smile_or_ZN_id
        self.QM_parameter = QM_parameter
        self.lmod = lmod
        
    
    # Search QM parameter in the main dataset which include 150K molecules only providing smile or ID 
    def search_QM_parameters(self) -> None:
        hint_column = 0
        path = '/scratch/users/amahamadakalapuwage/testing_QM30/generate_dataset/target_folders/CSV_files/25_30_final.csv'
        df_main = pd.read_csv(path)
        
        if self.smile_or_ZN_id.startswith("ZINC"):
            hint_column = "Zinc_id"
        else:
            hint_column = "smile"     
        try:
            info = df_main.loc[df_main[hint_column] == self.smile_or_ZN_id]
            if info.empty:
                print(f"provided '{self.smile_or_ZN_id}' is not correct Zinc_id or smile")
            else:
                print(float(info[self.QM_parameter].values))   
        except KeyError:
            print(f"Error: There is no QM parameter called '{self.QM_parameter}' in QM40") 


    # To access geometry and LmodA datasets only providing smile or ID
    def get_geo_Lmod_data(self) -> None: 
        ZINC_ID = 0
        smile = 0
        atom_range = 0
        final_df = 0
        main_dataset_path = '/scratch/users/amahamadakalapuwage/testing_QM30/generate_dataset/target_folders/CSV_files/25_30_final.csv'
        folder_path = '/scratch/users/amahamadakalapuwage/testing_QM30/generate_dataset/target_folders'
        df_main = pd.read_csv(main_dataset_path)
        
        if self.smile_or_ZN_id .startswith("ZINC"):
            ZINC_ID = self.smile_or_ZN_id 
            get_smile = df_main[df_main['Zinc_id'] == self.smile_or_ZN_id]
            smile_string = get_smile['smile'].values
            smile = smile_string[0]
        else:
            hint_column = "smile"
            hint_row = df_main[df_main[hint_column] == self.smile_or_ZN_id]
            ZINC_ID = hint_row['Zinc_id'].values[0]
            smile = self.smile_or_ZN_id 
        try:
            if ZINC_ID is None:
                print(f"provided '{self.smile_or_ZN_id }' is not correct Zinc_id or smile")
            else:
                print(smile)
                mol = Chem.MolFromSmiles(smile)
                number_of_atoms = mol.GetNumAtoms()
                if 10 <= number_of_atoms < 15:
                    atom_range = '10_15'
                elif 15 <= number_of_atoms < 20:
                    atom_range = '15_20'        
                elif 20 <= number_of_atoms < 25:
                    atom_range = '20_25'
                elif 25 <= number_of_atoms < 30:
                    atom_range = '25_30'
                elif 30 <= number_of_atoms < 35:
                    atom_range = '30_35'
                elif 35 <= number_of_atoms < 40:
                    atom_range = '35_40'
                else:
                    print("provided Zinc_id or smile is incorrect")
        except KeyError:
            print(f"provided Zinc_id or smile is not available in QM40 dataset")   
        
        list_dir = [name for name in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, name))]
        
        for current_dir in list_dir:
            if current_dir.startswith(atom_range):
                current_mol_path = os.path.join(folder_path, current_dir)
                sub_dir = [name1 for name1 in os.listdir(current_mol_path) if os.path.isdir(os.path.join(current_mol_path, name1))]
                
                for current_sub_dir in sub_dir:
                    current_sub_dir_path = os.path.join(current_mol_path, current_sub_dir)
                    ZINC_dir = [name2 for name2 in os.listdir(current_sub_dir_path) if os.path.isdir(os.path.join(current_sub_dir_path, name2))]
                    
                    for current_ZN_dir in ZINC_dir:
                        if ZINC_ID in current_ZN_dir:
                            target_dir_path = os.path.join(current_sub_dir_path, current_ZN_dir)
                            if self.lmod == 'lmod':
                                target_file_name = f"{current_ZN_dir}_lmod.csv"
                                file_name = target_file_name
                            
                            else:
                                target_file_name = f"{current_ZN_dir}_xyz.csv"
                                file_name = target_file_name
                                
                            target_file_path = os.path.join(target_dir_path, file_name)
                            target_df = pd.read_csv(target_file_path)
                            if "Unnamed: 0" in target_df.columns:
                                target_df = target_df.drop("Unnamed: 0", axis=1)
                                final_df = target_df
                            
                            final_df = target_df
                            print(target_df)
