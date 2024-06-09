#%%

"""Tests for `qm40_dataset_for_ml` package."""

import os
import sys
sys.path.append("/scratch/users/amahamadakalapuwage/testing_QM30/my_template/QM40_dataset_for_ML/qm40_dataset_for_ml/")
import unittest
import utils as ut


class TestQM40DatasetUtils(unittest.TestCase):
    """Tests for `qm40_utils` package."""

    def test_is_valid_smile(self):
        """Test module takes only "C", "F", "O", "N", "S", "Cl" """
        smile = "CCO"
        expected_answer = True
        actual_answer = ut.is_valid_smile(smile)
        self.assertEqual(expected_answer, actual_answer)


    def test_heavy_atom_count(self):
        """Test number of atoms module """
        smile = "CCO"
        expected_answer = True
        actual_answer = ut.heavy_atom_count(smile, 15)
        self.assertEqual(expected_answer, actual_answer)
        
        
    def test_PDBfromSmiles(self):
        """Test smile to PDB conversion """
        current_dir = os.getcwd()
        actual_answer = ut.PDBfromSmiles('smile_pdb_check.csv', current_dir)
        
        # Assert folder existence (check for both ZINC001299846328 and ZINC001299987058)
        self.assertTrue(os.path.isdir(os.path.join(current_dir, "ZINC001299846328")))
        self.assertTrue(os.path.isdir(os.path.join(current_dir, "ZINC001299987058")))

        # Assert pdb and xyz existence within folder ZINC001299846328
        self.assertTrue(os.path.isfile(os.path.join(current_dir, "ZINC001299846328", "ZINC001299846328.pdb")))
        self.assertTrue(os.path.isfile(os.path.join(current_dir, "ZINC001299846328", "ZINC001299846328.xyz")))
        
        # Assert pdb and xyz existence within folder ZINC001299846328
        self.assertTrue(os.path.isfile(os.path.join(current_dir, "ZINC001299987058", "ZINC001299987058.pdb")))
        self.assertTrue(os.path.isfile(os.path.join(current_dir, "ZINC001299987058", "ZINC001299987058.xyz")))
        
        
    def test_run_XTB(self):
        """Test smile to PDB conversion """
        current_dir = os.getcwd()
        actual_answer = ut.run_XTB(current_dir)

        # Assert file existence called optimized xyz file in both ZINC001299846328 and ZINC001299987058
        self.assertTrue(os.path.isfile(os.path.join(current_dir, "ZINC001299846328", "ZINC001299846328_xtb.xyz")))
        self.assertTrue(os.path.isfile(os.path.join(current_dir, "ZINC001299987058", "ZINC001299987058_xtb.xyz")))
    


# %%
