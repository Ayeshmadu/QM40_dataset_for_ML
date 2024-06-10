# QM40_dataset_for_ML


[![image](https://img.shields.io/pypi/v/QM40_dataset_for_ML.svg)](https://pypi.python.org/pypi/QM40_dataset_for_ML)
[![image](https://img.shields.io/conda/vn/conda-forge/QM40_dataset_for_ML.svg)](https://anaconda.org/conda-forge/QM40_dataset_for_ML)


**QM40 is a QMx type of dataset which includes 150K molecules optimized from B3LYP/6-31G(2df,p) level of theory in the Gaussian16 with QM parameters, optimized coordinates, Mulliken charges and Local vibrational mode parameters as a quantitative measurer of the bond strengths.**


-   Free software: MIT License
-   Documentation: https://Ayeshmadu.github.io/QM40_dataset_for_ML
    

## Features

-   Categorize smiles according to their heavy atom count.
-   Screen smiles with specific atoms.
-   Convert smiles to PDB and XYZ files.
-   Semi-empirical level of QM calculation (XTB).
-   Automated Gaussian16 input file generator.
-   Automated sbatch file generator for HPC.
-   Local vibrational mode (LmodA) calculations.
-   QM parameter, geometry, Mulliken charges, LmodA data extraction from Gaussian output files.
-   Extracted data converted into CSV files.

## Installation
```
pip install QM40-dataset-for-ML
```
## Dependencies

QM40_dataset_for_ML's Python dependencies are listed in its `requirements.txt` file. 

#### Linux
