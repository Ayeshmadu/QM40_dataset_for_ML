#  Lmod_run.py
#  Created by Ayesh Madushanka
#  Created on: May 16 2024
import os
import sys
import shutil
import pandas as pd
from qm40_dataset_for_ml import utils as ut
from qm40_dataset_for_ml import gaussian_info_extractor
#%%

class RunLmodCalc:
    
    """
    1. Class to Run Local vibrational mode calculations.
    2. Generate a csv file from extracted QM information. (any name can be given)
    3. Generate the geometry and LModA csv files for each molecules

    Args:
    current_dir_path (str): path to the folder which has optimized compounds as separate folders
    dataset (str): csv file with the information of smile and Zinc_id (optional)
    
    """


    def __init__(self, current_dir_path, dataset=None):
            self.current_dir_path = current_dir_path
            current_path = os.getcwd()
            target_folder_path = os.path.join(current_path, "Lmod_geom_results")
            if not os.path.exists(target_folder_path):
                os.makedirs(target_folder_path)
            self.target_folder_path = target_folder_path
            self.dataset = dataset
            self.dataset_name = "QM_parameters.csv"
            
            
    # running g16 results and Lmod calculations to extract g16 info, coordinates and Lmod info
    def run_lmod_collect_g16_info(self) -> list:
        list_dir = [name for name in os.listdir(self.current_dir_path) if os.path.isdir(os.path.join(self.current_dir_path, name))]
        qm_parameter_collector = []
        smiles = []
        img_freq = []
        for current_folder in list_dir:
            current_folder_path = os.path.join(self.current_dir_path, current_folder)
            
            if ut.incomplete_running_checker(current_folder_path, current_folder):
                
                if self.dataset is not None:
                    smile = ut.get_smile(self.dataset, current_folder)
                    smiles.append(smile)    
                ut.lmod_clac(current_folder_path, current_folder, self.target_folder_path)
                g16 = gaussian_info_extractor.GaussianInfoExtractor(current_folder_path, current_folder)
                imaginary_freq = g16.imaginary_freq_checker()
                
                if imaginary_freq is not None:
                    img_freq.append(imaginary_freq)
                    
                parameters, coordinates = g16.making_csvs()
                qm_parameter_collector.append(parameters)
                coordinate_name = f"{current_folder_path}_xyz.csv"
                coordinates.to_csv(coordinate_name)
                os.chdir(self.target_folder_path)
                source_path = os.path.join(current_folder_path, coordinate_name)
                target_fol_path = os.path.join(self.target_folder_path, current_folder)
                shutil.copy2(source_path, target_fol_path)
        return qm_parameter_collector, smiles, img_freq


    # generate main dataset
    def final_dataset_generator(self) -> None:
        qm_parameters, smiles, img_freq = self.run_lmod_collect_g16_info()
        print(smiles)
        column_names = ['Internal_E(0K)', 'HOMO', 'LUMO', 'Polarizability', 'spatial extent', 'dipol_mom', 'ZPE', 'rot1', 'rot2', 'rot3', 'Inter_E(298)', 'Enthalpy', 'Free_E', 'CV', 'Entropy']
        df_init = pd.DataFrame(qm_parameters, columns=column_names)
  
        if self.dataset is not None:
            smile_df = pd.DataFrame(smiles)
            df_init = pd.concat([smile_df[0], df_init], axis=1)
            df_init = df_init.rename(columns={0: 'smile'})
            df_init.to_csv(self.dataset_name, index=False)
            print(f"Final dataset saved as: {self.dataset_name}")
        
        else:
            df_init.to_csv('QM_parameters.csv', index=False)
            print(f"Final dataset saved as: {self.dataset_name}")
            
        if len(img_freq) != 0:
            df_img_freq = pd.DataFrame(img_freq, columns="imaginary_freq")
            df_img_freq.to_csv('imaginary_frequency_detect.csv') 

# %%
