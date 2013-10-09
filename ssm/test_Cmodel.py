from Cmodel import Cmodel
import unittest
import copy
import json
import os

class TestCmodel(unittest.TestCase):

    def setUp(self):
        model = json.load(open(os.path.join('..' ,'example', 'model', 'datapackage.json')))
        self.m = Cmodel(model) 
        
    def test_change_user_input(self):
        x = self.m.change_user_input('r0*2*correct_rate(v)')
        self.assertEqual(x, ['r0', '*', '2', '*', 'correct_rate', '(', 'v', ')'])

    def test_par_sv(self):
        self.assertEqual(set(self.m.par_sv), set(['S','R']))

    def test_remaider(self):
        self.assertEqual(self.m.remainder, 'I')

    def test_par_proc(self):
        self.assertEqual(set(self.m.par_proc), set(['r0', 'v', 'l', 'e', 'd', 'sto', 'alpha', 'vol', 'g']))

    def test_par_obs(self):
        self.assertEqual(set(self.m.par_obs), set(['rep','phi']))

    def test_par_fixed(self):
        self.assertEqual(set(self.m.par_fixed), set(['N', 'prop', 'mu_b', 'mu_d']))

    def test_par_fixed_obs(self):
        self.assertEqual(self.m.par_fixed_obs, set(['prop']))

    def test_drift_par_proc(self):
        self.assertEqual(set(self.m.drift_par_proc), set(['r0']))

    def test_vol_par_proc(self):
        self.assertEqual(set(self.m.vol_par_proc), set(['vol']))

    def test_drift_par_obs(self):
        self.assertEqual(self.m.drift_par_obs, [])

    def test_vol_par_obs(self):
        self.assertEqual(self.m.vol_par_obs, [])

    def test_drift_var(self):
        self.assertEqual(set(self.m.drift_var), set(['drift__par_proc__r0']))

    def test_proc_model(self):

        expected = [
            {'from': 'U', 'to': 'S',  'rate': 'mu_b*N'},
            {'from': 'S', 'to': 'E',  'rate': 'r0/N*v*(1.0+e*sin_t(d))*(N-S-R)', "tag": 'transmission', 'white_noise': {'name': 'white_noise__0', 'sd': 'sto'}}, #change  
            {'from': 'E', 'to': 'I', 'rate': '(1-alpha)*correct_rate(l)'},
            {'from': 'E', 'to': 'U',  'rate': 'alpha*correct_rate(l)'},
            {'from': 'E', 'to': 'U',  'rate': 'mu_d'},            
            {'from': 'S', 'to': 'U',  'rate': 'mu_d'},
            {'from': 'R', 'to': 'S',  'rate': 'g'},
            {'from': 'I', 'to': 'R', 'rate': '((1-alpha)*correct_rate(v))*(N-S-R)'}, #change
            {'from': 'I', 'to': 'U',  'rate': '(alpha*correct_rate(v)+mu_d)*(N-S-R)'} #change
        ]

        self.assertEqual(self.m.proc_model, expected)

    def test_obs_var(self):
        self.assertEqual(self.m.obs_var, ['Inc_out', 'Inc_in', 'Inc_weird', 'Prev', 'SI'])

    def test_obs_var_def(self):
        expected = [
            [{'from': 'I', 'to': 'R', 'rate': '((1-alpha)*correct_rate(v))*(N-S-R)'}, {'from': 'E', 'to': 'U', 'rate': 'mu_d'}],
            [{'from': 'S', 'to': 'E', 'rate': 'r0/N*v*(1.0+e*sin_t(d))*(N-S-R)', 'white_noise': {'name': 'white_noise__0', 'sd': 'sto'}}],
            [{'from': 'R', 'to': 'S', 'rate': 'g'}],
            ['I'],
            ['S', 'I']
        ]

        self.assertEqual(self.m.obs_var_def, expected)

if __name__ == '__main__':
    unittest.main()
