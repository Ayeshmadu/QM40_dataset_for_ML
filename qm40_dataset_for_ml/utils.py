# utils.py
#  Created by Ayesh Madushanka
#  Created on: May 16, 2024
import numpy as np
import pandas as pd
import os
import sys
import subprocess
import shutil

from rdkit.Chem import AllChem
from rdkit import Chem


# save dataframes as csvs
def save_df(datafrm: str, new_name: str):
    datafrm.to_csv(new_name, index=False)


#Function to keep smiles with specific atoms
def is_valid_smile(smile: str) -> bool: 
    smiles_atoms = ["C", "F", "O", "N", "S", "Cl"]
    mol = Chem.MolFromSmiles(smile)
    if mol is None:  
        return False
    for atom in mol.GetAtoms():
        atom_symbol = atom.GetSymbol()
        if atom_symbol not in smiles_atoms:
            return False
    return True


# Function to keep specific amount of heavy atoms in a smile
def heavy_atom_count(smiles: str, num_atoms: int=30) -> pd.DataFrame:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    return mol.GetNumAtoms() < num_atoms


# call valid smile function
def call_is_valid_smile(dataset: str, smile_column: str = "Column1") -> pd.DataFrame:
    df = pd.read_csv(dataset)
    filtered_df = df[df[smile_column].apply(is_valid_smile)]
    return filtered_df


# call heavy atom count function
def call_heavy_atom_count(dataset: str, smile_column: str = "Column1", threshold: int = 30) -> pd.DataFrame:
    df = pd.read_csv(dataset)
    filtered_df = df[df[smile_column].apply(heavy_atom_count, num_atoms=threshold)]
    return filtered_df


# Separate all the pdbs into 5 sub-directories
def separate_5_folders(source_dir_path: str, dest_dir_path: str):
    source_dir = source_dir_path
    dest_dir = dest_dir_path
    num_new_folders = 5
    subfolders = os.listdir(source_dir)
    subfolders_per_folder = len(subfolders) // num_new_folders 
    for i in range(num_new_folders):
        new_folder_path = os.path.join(dest_dir, f'folder_{i}')
        os.makedirs(new_folder_path, exist_ok=True)
        start_idx = i * subfolders_per_folder
        end_idx = start_idx + subfolders_per_folder
        for subfolder in subfolders[start_idx:end_idx]:
            source_subfolder = os.path.join(source_dir, subfolder)
            shutil.copytree(source_subfolder, os.path.join(new_folder_path, subfolder))


def atom_composition(mol):
    smiles_atoms = ["C", "F", "O", "N", "S", "Cl"]
    for atom in mol.GetAtoms():
        atom_symbol = atom.GetSymbol()
        if atom_symbol not in smiles_atoms:
            return False
    return True


def smile_sceening(dataset: str) -> list:
    zetoto10 = []
    tento15 = []
    fifteento20 = []
    twentyto25 = []
    twenty5to30 = []
    thirtyto35 = [] 
    df = pd.read_csv(dataset)
    
    for index, row in df.iterrows():
        Zinc_id = row['Column2']
        smile = row['Column1']
        mol = Chem.MolFromSmiles(smile)
        if mol is None:  
            return False
        if Chem.GetFormalCharge(mol) == 0:
            number_of_atoms = mol.GetNumAtoms()
            if atom_composition(mol):
                if 0 <= number_of_atoms < 10:
                    zetoto10.append([Zinc_id, smile])
                elif 10 <= number_of_atoms < 15:
                    tento15.append([Zinc_id, smile])
                elif 15 <= number_of_atoms < 20:
                    fifteento20.append([Zinc_id, smile])
                elif 20 <= number_of_atoms < 25:
                    twentyto25.append([Zinc_id, smile])
                elif 25 <= number_of_atoms < 30:
                    twenty5to30.append([Zinc_id, smile])
                elif 30 <= number_of_atoms < 35:
                    thirtyto35.append([Zinc_id, smile])
    return zetoto10, tento15, fifteento20, twentyto25, twenty5to30, thirtyto35

        
# creating PDB from smiles -> .xyz
def PDBfromSmiles(data: str, current_path: str) -> None:
    data = pd.read_csv(data) 
    for index, row in data.iterrows():
        id = row['Index']
        smile = row['smile']
        mol = Chem.MolFromSmiles(smile)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        file_name = f"{id}.pdb"
        folder_name = f"{id}"
        folder_path = os.path.join(current_path, folder_name)
        os.makedirs(folder_path, exist_ok=True)
        pdb_path = os.path.join(folder_path, file_name)
        Chem.MolToPDBFile(mol, pdb_path)
        xyz_file_name = f"{id}.xyz"
        xyz_filepath = os.path.join(folder_path, xyz_file_name)
        subprocess.call(["obabel", "-ipdb", pdb_path, "-oxyz", "-O", xyz_filepath])
   
    
# run XTB (Semi empherical QM calculation to optimize pdbs)    
def run_XTB(current_path: str) -> None:
    list_dir = [name for name in os.listdir(current_path) if os.path.isdir(name)]
    for current_dir in list_dir:
        workpath = os.path.join(current_path, current_dir)
        files = os.listdir(workpath)
        for file in files:
            if file.lower().endswith('.xyz'):
                original_name, extension_ori = os.path.splitext(file)
                new_file_name = f"{original_name}_xtb{extension_ori}"
                file = os.path.join(workpath, file)
                xtb_opt_file_path = os.path.join(workpath, "xtbopt.xyz")
                new_file_name_path = os.path.join(workpath, new_file_name)
                os.chdir(workpath)
                xtb_command = f"xtb {file} --opt"
                subprocess.run(xtb_command,shell=True)
                removable_files = ['charges',  'wbo', 'xtbopt.log', 'xtbrestart', 'xtbtopo.mol']
                try:
                    os.rename(xtb_opt_file_path, new_file_name_path)
                    for r_file in removable_files:
                        os.remove(r_file)
                except FileNotFoundError:
                    continue
                
                os.chdir(current_path)


#creating gaussian input file from .xyz file
def input_file_creator(current_dir_path: str) -> None:
    list_dir = [name for name in os.listdir(current_dir_path) if os.path.isdir(os.path.join(current_dir_path, name))]
    for current_dir in list_dir:
        command = []
        workpath = os.path.join(current_dir_path, current_dir)
        print(workpath)
        os.chdir(workpath)

        # run the python commands
        command.append('%chk=' + current_dir +'.chk \n')
        command.append('%NProcShared=16 \n')
        command.append('%mem=48GB \n') 
        command.append('\n')
        command.append('#p B3LYP/6-31G(2df,p) Opt Freq \n')
        command.append('\n')
        command.append( current_dir +'\n')
        command.append('\n')
        command.append('0 1 \n')
        
        xyz_name = f"{current_dir}_xtb.xyz"
        xyz_path = os.path.join(workpath, xyz_name)
        try:
            with open(xyz_path, "r") as f:
                lines = f.readlines()[2:]
                for single_line in lines:
                    command.append(single_line)    
            command.append('\n')
        except FileNotFoundError:
            print(f"xyz file {xyz_name} not found in {workpath}, skipping xyz data.")
            
        # Creating the .inp file
        Input_file = current_dir + ".inp"
        sbfile = open(Input_file,"w")

        for i in command :
            sbfile.write(i)
        sbfile.close()
        os.chdir(current_dir_path)

 
# sbatch submission for Gaussian
def sbatch_create_submission(current_dir_path: str) -> None:
    list_dir = [name for name in os.listdir(current_dir_path) if os.path.isdir(os.path.join(current_dir_path, name))] 
    for current_dir in list_dir:
        command = []
        workpath = os.path.join(current_dir_path, current_dir)
        os.chdir(workpath)
        command.append('module load gaussian/g16c \n')
        command.append('\n')
        command.append('gaussian_input_file='+ current_dir +'.inp \n')
        command.append('\n')
        command.append('job_directory="${SLURM_JOB_NAME}_${SLURM_JOB_ID}" \n')
        command.append('mkdir ${job_directory} \n')
        command.append('cp ${gaussian_input_file} ${gaussian_checkpoint} ${job_directory}/ \n')
        command.append('cd ${job_directory} \n')
        command.append('g16 < ${gaussian_input_file}')
        command.append('\n')
        sbatch = current_dir + "_run.sbatch"
        sbfile = open(sbatch,"w")
        sbfile.write('#!/bin/bash \n')
        sbfile.write('#SBATCH -J ' + current_dir + ' \n') 
        sbfile.write('#SBATCH -o ' + current_dir + '_run.out \n')
        sbfile.write('#SBATCH -p  standard-s \n')
        sbfile.write('#SBATCH -c 16 \n')
        sbfile.write('#SBATCH -N 1 \n')
        sbfile.write('#SBATCH --mem=50GB \n')
        sbfile.write('\n')

        for i in command :
            sbfile.write(i)
        sbfile.close()
        cmd = "sbatch " + sbatch
        os.system(cmd)
        os.chdir(current_dir_path)


# Extract the smile     
def get_smile(dataset: str, target_id: str) -> str:
    df = pd.read_csv(dataset)
    smile = df[df['Index'] == target_id]['smile'].values
    return smile 


# local vibrational mode calculations
def lmod_clac(path: str, current_dir: str, target_folder: str) -> None:
        subfolders = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
        subfolder_path = os.path.join(path, subfolders[0])
        os.chdir(subfolder_path)
        gen_fchk = f"formchk {current_dir}.chk {current_dir}.fchk"
        os.system("module load gaussian/g16c && " + gen_fchk)
        pdb_input_name = os.path.join(path, current_dir)
        Obabel_command = f"obabel -i pdb {pdb_input_name}.pdb -o mol -O ddd.mol"
        subprocess.run(Obabel_command, shell=True)
        command = creating_lmod_input(current_dir)
        Input_file = f"job.inp"
        sbfile = open(Input_file,"w")
        
        for i in command :
            sbfile.write(i)
        sbfile.close()
        lmod_command = f"python3.10 /work/group/catco/LModeA-github/LmodeA/lmodea-301/src/driver/LModeA.py job.inp"
        subprocess.run(lmod_command, shell=True)
        lmod_out_name = f"job.out"
        gout_path = os.path.join(subfolder_path, lmod_out_name)
        lmod_info = []
        target_word = "Local mode properties:" 
        try:
            with open(gout_path, 'r') as f:
                lines = f.readlines()
                output_data = (lines, lmod_out_name)
        except FileNotFoundError:
                print(f"Lmod out file {lmod_out_name} not found in {subfolder_path}, skipping lmod_out data.")
                
        for i, line in enumerate(lines):
            row = line.strip()
            if target_word in row:  
                for lmod_line in lines[i+4:i+100]: 
                    if not int(lmod_line.strip()[12:15]) == 0:
                        break
                    
                    lmod_split = lmod_line.split()
                    atom1 = lmod_split[1]
                    atom2 = lmod_split[2]
                    name = lmod_split[5]
                    lmodfc = lmod_split[7]
                    name_tag = "".join(char for char in name if char.isalpha())
                    lmod_info.append([atom1, atom2, name, name_tag, lmodfc])
                    
        column_names = ['Atom1', 'Atom2', 'Name', 'Name_tag', 'lmod']
        lmod_parametars = pd.DataFrame(lmod_info, columns = column_names)
        lmod_file_name = f"{current_dir}_lmod.csv"
        lmod_parametars.to_csv(lmod_file_name, index=False)
        tar_fol_path = os.path.join(target_folder, current_dir)
        if not os.path.exists(tar_fol_path):
            os.mkdir(tar_fol_path)
        lmod_file_path = os.path.join(subfolder_path, lmod_file_name)
        shutil.copy2(lmod_file_path, tar_fol_path)
        os.chdir(path)
   
   
# Local vibrational mode input file generator
def creating_lmod_input(current_dir: str) -> list:
    command = []
    command.append('$contrl \n')
    command.append('  qcprog = "gaussian" \n')
    command.append('  iprint = 1 \n') 
    command.append('  isymm = 1 \n')
    command.append('  iredun = 2 \n')
    command.append('  iacs = 2 \n')
    command.append('  nstep = 50 \n')
    command.append('$end \n')
    command.append('\n')
    command.append('\n')
    command.append('$qcdata \n')
    command.append('  fchk = "' + current_dir + '.fchk"\n')
    command.append('$end \n')
    command.append('\n')
    command.append('\n')
    command.append('$LocMod \n')
    command.append('  connect = "ddd.mol" \n')
    command.append('$End \n')
    command.append('\n')
    return command

# Incomplete running checker
def incomplete_running_checker(curr_path: str, current_file: str) -> bool:
    current_gout = f"{current_file}_run.out" 
    curr_file = os.path.join(curr_path, current_gout)
    
    with open(curr_file, 'r') as f:
        lines = f.readlines()
        
    for i, line in enumerate(lines):
        last_line = lines[-1].strip()
        target_phrase = "Normal termination of Gaussian"
        if target_phrase not in last_line:
            return False
    return True
