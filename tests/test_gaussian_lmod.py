# Since Gaussian16 and LmodA packages are commercial software and Closed Source Software respectively,
# this test model won't work on other systems.
# However we offer our Local vibrational mode calculation package upon request for academic porposes.

import os
import sys
sys.path.append("/scratch/users/amahamadakalapuwage/testing_QM30/my_template/QM40_dataset_for_ML/qm40_dataset_for_ml/")
import unittest
from Lmod_run import RunLmodCalc
import utils as ut


class TestQM40DatasetUtils(unittest.TestCase):
    """Tests for `qm40_utils` package."""
    
    
    def test_G16_LmodA(self):
        """Test smile to PDB conversion """
        current_dir = os.getcwd()
        current_dir1 = os.path.join(current_dir, 'g16_out')
        output_dir = os.path.join(current_dir, 'Lmod_geom_results')
        g16 = RunLmodCalc(current_dir1)
        actual_answer = g16.final_dataset_generator()

        # Assert 3 csv files called Qm_parameters.csv, opt_geometry and LmodA existence after G16 
        # and Local vibrational mode calculations
        self.assertTrue(os.path.isfile(os.path.join(output_dir, "QM_parameters.csv")))
        self.assertTrue(os.path.isfile(os.path.join(output_dir, "ZINC001729833512", "ZINC001729833512_lmod.csv")))
        self.assertTrue(os.path.isfile(os.path.join(output_dir, "ZINC001729833512", "ZINC001729833512_xyz.csv")))
    

# %%
