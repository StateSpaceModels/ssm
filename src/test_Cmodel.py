from Cmodel import Cmodel
import unittest
import copy
import os

class TestCmodel(unittest.TestCase):

    def setUp(self):
        self.m = Cmodel(os.path.join('..' ,'examples', 'foo', 'package.json'), "sir")

    def test_change_user_input(self):
        x = self.m.change_user_input('r0*2*correct_rate(v)')
        self.assertEqual(x, ['r0', '*', '2', '*', 'correct_rate', '(', 'v', ')'])

    def test_par_sv(self):
        self.assertEqual(self.m.par_sv, ['I_nyc', 'I_paris', 'S_nyc', 'S_paris'])

    def test_remainder(self):
        self.assertEqual(self.m.remainder, ['R_nyc', 'R_paris'])

    def test_par_proc(self):
        self.assertEqual(self.m.par_proc, ['r0_nyc', 'r0_paris', 'v'])

    def test_par_disp(self):
        self.assertEqual(self.m.par_disp, ['vol'])

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

    def test_make_C_term(self):
        terms = [
            {'x': 'mu_b_paris*(1.0+v*sin((v/N_paris+(mu_b_paris)))) + r0_paris', #input
             'h': 'mu_b_paris*(v*sin(mu_b_paris+v/N_paris)+1.0)+r0_paris', #expected human output
             'c': 'gsl_spline_eval(calc->spline[ORDER_mu_b_paris],t,calc->acc[ORDER_mu_b_paris])*(gsl_vector_get(par,ORDER_v)*sin(gsl_spline_eval(calc->spline[ORDER_mu_b_paris],t,calc->acc[ORDER_mu_b_paris])+gsl_vector_get(par,ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_paris],t,calc->acc[ORDER_N_paris]))+1.0)+diffed[ORDER_diff__r0_paris]'}, #expected C output

            {'x': 'N_paris-S_paris-I_paris+S_paris+I_paris',
             'h': 'N_paris',
             'c': 'gsl_spline_eval(calc->spline[ORDER_N_paris],t,calc->acc[ORDER_N_paris])'},

            {'x': 'rep_all_CDC_inc*(1.0-rep_all_CDC_inc)*prop_all_CDC_inc*x + (rep_all_CDC_inc*phi*prop_all_CDC_inc*x)**2',
             'h': 'pow(phi,2)*pow(prop_all_CDC_inc,2)*pow(rep_all_CDC_inc,2)*pow(x,2)+prop_all_CDC_inc*rep_all_CDC_inc*x*(-rep_all_CDC_inc+1.0)',
             'c': 'pow(gsl_vector_get(par,ORDER_phi),2)*pow(gsl_spline_eval(calc->spline[ORDER_prop_all_CDC_inc],t,calc->acc[ORDER_prop_all_CDC_inc]),2)*pow(gsl_vector_get(par,ORDER_rep_all_CDC_inc),2)*pow(x,2)+gsl_spline_eval(calc->spline[ORDER_prop_all_CDC_inc],t,calc->acc[ORDER_prop_all_CDC_inc])*gsl_vector_get(par,ORDER_rep_all_CDC_inc)*x*(-gsl_vector_get(par,ORDER_rep_all_CDC_inc)+1.0)'},
        ]

        for t in terms:
            self.assertEqual(self.m.make_C_term(t['x'], False, human=True), t['h'])
            self.assertEqual(self.m.make_C_term(t['x'], False, human=False), t['c'])

    def test_make_C_term_derivate(self):
        x = 'sin(2*M_PI*(t/ONE_YEAR +r0_paris))'
        h = '2*M_PI*cos(2*M_PI*(r0_paris+t/ONE_YEAR))'
        c = '2*M_PI*cos(2*M_PI*(diffed[ORDER_diff__r0_paris]+t/ONE_YEAR))'

        self.assertEqual(self.m.make_C_term(x, False, human=True, derivate='r0_paris'), h)
        self.assertEqual(self.m.make_C_term(x, False, human=False, derivate='r0_paris'), c)

    def test_make_C_term_skip_correct_rate(self):
        #correct_rate is only skipped for C code

        x = 'mu_b_paris*(1.0+correct_rate(v)*sin((correct_rate(v)/N_paris+(mu_b_paris)))) + r0_paris'
        c = 'gsl_spline_eval(calc->spline[ORDER_mu_b_paris],t,calc->acc[ORDER_mu_b_paris])*((gsl_vector_get(par,ORDER_v))*sin(gsl_spline_eval(calc->spline[ORDER_mu_b_paris],t,calc->acc[ORDER_mu_b_paris])+(gsl_vector_get(par,ORDER_v))/gsl_spline_eval(calc->spline[ORDER_N_paris],t,calc->acc[ORDER_N_paris]))+1.0)+diffed[ORDER_diff__r0_paris]'

        self.assertEqual(self.m.make_C_term(x, True, human=False), c)

        #    def test_make_C_term_extra_terms(self):
        #terms = [
        #   {'x': 'terms_forcing(v)', 
        #   'c': 'terms_forcing(gsl_vector_get(par,ORDER_v),t,p_data,cac)'}, 
        #
        #  {'x': 'heaviside(t-v)', 
        #   'c': 'heaviside(t-gsl_vector_get(par,ORDER_v))'}, 
        #
        #            {'x': 'heaviside((t-v))', 
        #   'c': 'heaviside(t-gsl_vector_get(par,ORDER_v))'}, 
        #
        #  {'x': 'ramp(t-v)', 
        #   'c': 'ramp(t-gsl_vector_get(par,ORDER_v))'}, 

        #            {'x': 'correct_rate(v)', 
        #   'c': 'ssm_correct_rate(gsl_vector_get(par,ORDER_v),dt)'}, 
        #]
            
        #for t in terms:
        #   self.assertEqual(self.m_diff.make_C_term(t['x'], False, human=False), t['c'])

        #    def test_cache_special_function_C(self):
        #
        #        caches = map(lambda x: self.m_diff.make_C_term(x, True), ['sin(2*M_PI*(t/ONE_YEAR +r0))', 'sin(2*M_PI*(t/ONE_YEAR +r0))', 'terms_forcing(v)*sin(2*M_PI*(t/ONE_YEAR +r0))*terms_forcing(v)'])
        #sf = self.m_diff.cache_special_function_C(caches, prefix='_sf[cac]')

        #        self.assertEqual(sf, ['sin(2*M_PI*(drifted[ORDER_drift__par_proc__r0][cac]+t/ONE_YEAR))', 'terms_forcing(par[ORDER_v][routers[ORDER_v]->map[cac]],t,p_data,cac)'])
        #self.assertEqual(caches, ['_sf[cac][0]', '_sf[cac][0]', 'pow(_sf[cac][1],2)*_sf[cac][0]'])



if __name__ == '__main__':
    unittest.main()
