# QM40data_search.py
#  Created by Ayesh Madushanka
#  Created on: June 4, 2024

import os
import pandas as pd
from rdkit import Chem


# help function generation
def help_func():
    print("Welcome to the QM40 dataset data seraching modular!!!!")
    print("QM40 dataset includes a main dataset which has 163K molecules")
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
    print(
        "If you need to access bond (with Lmod) or optimized geometry (with charges) of a molecule"
    )
    print("Use get_geo_Lmod_data function with additional key word 'lmod'")
    print("lmod_or_xyz = 'bond': bond dataset")
    print("lmod_or_xyz = 'xyz': Optimized Geometry dataset with Mulliken charges")


class QM40DataExtractor:
    def __init__(self, smile_or_ZN_id, QM_parameter=None, lmod_or_xyz=None):
        self.smile_or_ZN_id = smile_or_ZN_id
        self.QM_parameter = QM_parameter
        self.lmod_or_xyz = lmod_or_xyz
        self.current_path = os.getcwd()
        self.id_smile = self._search_identify()

    # Search QM parameter in the main dataset which include 150K molecules only providing smile or ID
    def _search_identify(self) -> None:
        hint_column = 0
        if self.smile_or_ZN_id.startswith("ZINC"):
            hint_column = "Zinc_id"
        else:
            hint_column = "smile"
        return hint_column

    def QM_parameter_serach(self) -> None:
        QM40_main_path = os.path.join(self.current_path, "QM40_main.csv")
        df_main = pd.read_csv(QM40_main_path)
        try:
            info = df_main.loc[df_main[self.id_smile] == self.smile_or_ZN_id]
            if info.empty:
                print(
                    f"provided '{self.smile_or_ZN_id}' is not correct Zinc_id or smile"
                )
            else:
                print(float(info[self.QM_parameter].values))
        except KeyError:
            print(
                f"Error: There is no QM parameter called '{self.QM_parameter}' in QM40"
            )

    # To access geometry and bond datasets only providing smile or ID
    def get_geo_Lmod_data(self) -> None:
        if self.lmod_or_xyz == "bond":
            bond_path = os.path.join(self.current_path, "QM40_bond.csv")
            df_bond = pd.read_csv(bond_path)
            molecule = df_bond[df_bond[self.id_smile] == self.smile_or_ZN_id]
            print(molecule)

        elif self.lmod_or_xyz == "xyz":
            xyz_path = os.path.join(self.current_path, "QM40_xyz.csv")
            df_xyz = pd.read_csv(xyz_path)
            molecule = df_xyz[df_xyz[self.id_smile] == self.smile_or_ZN_id]
            print(molecule)

        else:
            print(
                "Bond dataset and xyz dataset can only be accessed using 'bond' and 'xyz' keywords, respectively."
            )
