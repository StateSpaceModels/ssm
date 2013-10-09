from Cmodel import Cmodel
import unittest
import copy
import json
import os

class TestCmodel(unittest.TestCase):

    def setUp(self):
        model = json.load(open(os.path.join('..' ,'example', 'foo', 'datapackages', 'model-seb-sir', 'datapackage.json')))
        self.m = Cmodel(model)

    def test_change_user_input(self):
        x = self.m.change_user_input('r0*2*correct_rate(v)')
        self.assertEqual(x, ['r0', '*', '2', '*', 'correct_rate', '(', 'v', ')'])

    def test_par_sv(self):
        self.assertEqual(self.m.par_sv, ['I_nyc', 'I_paris', 'S_nyc', 'S_paris'])

    def test_remainder(self):
        self.assertEqual(self.m.remainder, ['R_nyc', 'R_paris'])

    def test_par_proc(self):
        self.assertEqual(self.m.par_proc, ['r0_nyc', 'r0_paris', 'v', 'vol'])

    def test_par_obs(self):
        self.assertEqual(self.m.par_obs, ['phi', 'rep_all_CDC_inc', 'rep_all_google_inc', 'rep_nyc_CDC_inc', 'rep_paris_CDC_prev'])

    def test_par_forced(self):
        self.assertEqual(self.m.par_forced, ['N_nyc', 'N_paris', 'mu_b_nyc', 'mu_b_paris', 'mu_d_nyc', 'mu_d_paris', 'prop_all_CDC_inc', 'prop_all_google_inc', 'prop_nyc_CDC_inc', 'prop_paris_CDC_prev'])

    def test_par_diff(self):
        self.assertEqual(self.m.par_diff, ['diff__r0_nyc', 'diff__r0_paris'])

    def test_par_inc(self):
        self.assertEqual(self.m.par_inc, ['Inc_in_nyc', 'Inc_out'])

    def test_par_noise(self):
        self.assertEqual(self.m.par_noise, ['sto'])

    def test_par_other(self):
        self.assertEqual(self.m.par_other, [])


if __name__ == '__main__':
    unittest.main()
