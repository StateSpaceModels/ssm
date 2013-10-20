from Ccoder import Ccoder
import unittest
import copy
import json
import os

class TestCcoder(unittest.TestCase):

    def setUp(self):
        model = json.load(open(os.path.join('..' ,'examples', 'noise', 'datapackages', 'model-jdureau-noise', 'datapackage.json')))
        self.m_noise = Ccoder(model)

        model = json.load(open(os.path.join('..' ,'examples', 'noise', 'datapackages', 'model-jdureau-noise', 'datapackage.json')))
        del model["resources"][1]["data"][2]["white_noise"] 
        del model["resources"][1]["data"][3]["white_noise"]
        model['resources'].append({})
        model['resources'][4]['name']='sde'
        model['resources'][4]['data']={}
        model['resources'][4]['data']['drift']=[]
        model['resources'][4]['data']['drift'].append({})
        model['resources'][4]['data']['drift'][0]['name']='r0_paris'
        model['resources'][4]['data']['drift'][0]['f']=0.0
        model['resources'][4]['data']['drift'][0]['transformation']="log(r0_paris)"
        model['resources'][4]['data']['drift'].append({})
        model['resources'][4]['data']['drift'][1]['name']='r0_nyc'
        model['resources'][4]['data']['drift'][1]['f']=0.0
        model['resources'][4]['data']['drift'][1]['transformation']="log(r0_nyc)"
        model['resources'][4]['data']['dispersion']=[['vol',0],[0,'vol']]
        model['resources'][3]['data'][7]['name']= 'sto'
        diff_model = model
        self.m_diff = Ccoder(model)

        model = json.load(open(os.path.join('..' ,'examples', 'noise', 'datapackages', 'model-jdureau-noise', 'datapackage.json')))
        del model["resources"][1]["data"][2]["white_noise"] 
        del model["resources"][1]["data"][3]["white_noise"]
        model["resources"][1]["data"][0]["white_noise"] = {"name":"noise_SI", "sd": "sto"}
        model["resources"][1]["data"][1]["white_noise"] = {"name":"noise_SI2", "sd": "sto"}
        self.m_noise2 = Ccoder(model)

        model = json.load(open(os.path.join('..' ,'examples', 'noise', 'datapackages', 'model-jdureau-noise', 'datapackage.json')))
        del model["resources"][1]["data"][2]["white_noise"] 
        del model["resources"][1]["data"][3]["white_noise"]
        model["resources"][1]["data"][4]["white_noise"] = {"name":"noise_SI", "sd": "sto"}
        model["resources"][1]["data"][5]["white_noise"] = {"name":"noise_SI2", "sd": "sto"}
        self.m_noise3 = Ccoder(model)

        model = json.load(open(os.path.join('..' ,'examples', 'noise', 'datapackages', 'model-jdureau-noise', 'datapackage.json')))
        del model["resources"][1]["data"][2]["white_noise"] 
        del model["resources"][1]["data"][3]["white_noise"]
        model["resources"][1]["data"][8]["white_noise"] = {"name":"noise_SI", "sd": "sto"}
        model["resources"][1]["data"][9]["white_noise"] = {"name":"noise_SI2", "sd": "sto"}
        self.m_noise4 = Ccoder(model)

        model = json.load(open(os.path.join('..' ,'examples', 'noise', 'datapackages', 'model-jdureau-noise', 'datapackage.json')))
        del model["resources"][1]["data"][2]["white_noise"] 
        del model["resources"][1]["data"][3]["white_noise"]
        model["resources"][1]["data"][10]["white_noise"] = {"name":"noise_SI", "sd": "sto"}
        model["resources"][1]["data"][11]["white_noise"] = {"name":"noise_SI2", "sd": "sto"}
        self.m_noise5 = Ccoder(model)

        model = json.load(open(os.path.join('..' ,'examples', 'noise', 'datapackages', 'model-jdureau-noise', 'datapackage.json')))
        model["resources"][1]["data"][4]["white_noise"] = {"name":"noise_SI", "sd": "sto"}
        model["resources"][1]["data"][5]["white_noise"] = {"name":"noise_SI2", "sd": "sto"}
        self.m_noise6 = Ccoder(model)

        model = json.load(open(os.path.join('..' ,'examples', 'noise', 'datapackages', 'model-jdureau-noise', 'datapackage.json')))
        model["resources"][1]["data"][4]["white_noise"] = {"name":"noise_SI23", "sd": "sto"}
        model["resources"][1]["data"][5]["white_noise"] = {"name":"noise_SI24", "sd": "sto"}
        self.m_noise7 = Ccoder(model)

        model = diff_model
        model["resources"][1]["data"].append({"from": "R_paris",   "to": "I_paris",   "rate": "correct_rate(v)",            "description":"testing"})
        model["resources"][1]["data"].append({"from": "R_nyc",   "to": "I_nyc",   "rate": "correct_rate(v)",                "description":"testing"})
        self.m_diff2 = Ccoder(model)

        
    def test_make_C_term(self):
        terms = [
            {'x': 'mu_b_paris*(1.0+v*sin((v/N_paris+(mu_b_paris)))) + r0_paris', #input
             'h': 'mu_b_paris*(v*sin(mu_b_paris+v/N_paris)+1.0)+r0_paris', #expected human output
             'c': 'gsl_spline_eval(calc->spline[ORDER_mu_b_paris],t,calc->acc[ORDER_mu_b_paris])*(gsl_vector_get(par, ORDER_v)*sin(gsl_spline_eval(calc->spline[ORDER_mu_b_paris],t,calc->acc[ORDER_mu_b_paris])+gsl_vector_get(par, ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_paris],t,calc->acc[ORDER_N_paris]))+1.0)+diffed[ORDER_diff__r0_paris]'}, #expected C output

            {'x': 'N_paris-S_paris-I_paris+S_paris+I_paris',
             'h': 'N_paris',
             'c': 'gsl_spline_eval(calc->spline[ORDER_N_paris],t,calc->acc[ORDER_N_paris])'},

            {'x': 'rep_all_CDC_inc*(1.0-rep_all_CDC_inc)*prop_all_CDC_inc*x + (rep_all_CDC_inc*phi*prop_all_CDC_inc*x)**2',
             'h': 'pow(phi,2)*pow(prop_all_CDC_inc,2)*pow(rep_all_CDC_inc,2)*pow(x,2)+prop_all_CDC_inc*rep_all_CDC_inc*x*(-rep_all_CDC_inc+1.0)',
             'c': 'pow(gsl_vector_get(par, ORDER_phi),2)*pow(gsl_spline_eval(calc->spline[ORDER_prop_all_CDC_inc],t,calc->acc[ORDER_prop_all_CDC_inc]),2)*pow(gsl_vector_get(par, ORDER_rep_all_CDC_inc),2)*pow(x,2)+gsl_spline_eval(calc->spline[ORDER_prop_all_CDC_inc],t,calc->acc[ORDER_prop_all_CDC_inc])*gsl_vector_get(par, ORDER_rep_all_CDC_inc)*x*(-gsl_vector_get(par, ORDER_rep_all_CDC_inc)+1.0)'},
        ]

        for t in terms:
            self.assertEqual(self.m_diff.make_C_term(t['x'], False, human=True), t['h'])
            self.assertEqual(self.m_diff.make_C_term(t['x'], False, human=False), t['c'])

    def test_make_C_term_derivate(self):
        x = 'sin(2*M_PI*(t/ONE_YEAR +r0_paris))'
        h = '2*M_PI*cos(2*M_PI*(r0_paris+t/ONE_YEAR))'
        c = '2*M_PI*cos(2*M_PI*(diffed[ORDER_diff__r0_paris]+t/ONE_YEAR))'

        self.assertEqual(self.m_diff.make_C_term(x, False, human=True, derivate='r0_paris'), h)
        self.assertEqual(self.m_diff.make_C_term(x, False, human=False, derivate='r0_paris'), c)

    def test_make_C_term_skip_correct_rate(self):
        #correct_rate is only skipped for C code

        x = 'mu_b_paris*(1.0+correct_rate(v)*sin((correct_rate(v)/N_paris+(mu_b_paris)))) + r0_paris'
        c = 'gsl_spline_eval(calc->spline[ORDER_mu_b_paris],t,calc->acc[ORDER_mu_b_paris])*((gsl_vector_get(par, ORDER_v))*sin(gsl_spline_eval(calc->spline[ORDER_mu_b_paris],t,calc->acc[ORDER_mu_b_paris])+(gsl_vector_get(par, ORDER_v))/gsl_spline_eval(calc->spline[ORDER_N_paris],t,calc->acc[ORDER_N_paris]))+1.0)+diffed[ORDER_diff__r0_paris]'

        self.assertEqual(self.m_diff.make_C_term(x, True, human=False), c)

        #    def test_make_C_term_extra_terms(self):
        #terms = [
        #   {'x': 'terms_forcing(v)', 
        #   'c': 'terms_forcing(gsl_vector_get(par, ORDER_v),t,p_data,cac)'}, 
        #
        #  {'x': 'heaviside(t-v)', 
        #   'c': 'heaviside(t-gsl_vector_get(par, ORDER_v))'}, 
        #
        #            {'x': 'heaviside((t-v))', 
        #   'c': 'heaviside(t-gsl_vector_get(par, ORDER_v))'}, 
        #
        #  {'x': 'ramp(t-v)', 
        #   'c': 'ramp(t-gsl_vector_get(par, ORDER_v))'}, 

        #            {'x': 'correct_rate(v)', 
        #   'c': 'ssm_correct_rate(gsl_vector_get(par, ORDER_v),dt)'}, 
        #]
            
        #for t in terms:
        #   self.assertEqual(self.m_diff.make_C_term(t['x'], False, human=False), t['c'])

        #    def test_cache_special_function_C(self):
        #
        #        caches = map(lambda x: self.m_diff.make_C_term(x, True), ['sin(2*M_PI*(t/ONE_YEAR +r0))', 'sin(2*M_PI*(t/ONE_YEAR +r0))', 'terms_forcing(v)*sin(2*M_PI*(t/ONE_YEAR +r0))*terms_forcing(v)'])
        #sf = self.m_diff.cache_special_function_C(caches, prefix='_sf[cac]')

        #        self.assertEqual(sf, ['sin(2*M_PI*(drifted[ORDER_drift__par_proc__r0][cac]+t/ONE_YEAR))', 'terms_forcing(par[ORDER_v][routers[ORDER_v]->map[cac]],t,p_data,cac)'])
        #self.assertEqual(caches, ['_sf[cac][0]', '_sf[cac][0]', 'pow(_sf[cac][1],2)*_sf[cac][0]'])

    def test_eval_Q(self):
        calc_Q = self.m_noise.eval_Q()

        # testing env sto only
        # order: I_nyc, I_paris, S_nyc, S_paris, inc_all , inc_nyc
        term_paris = '((((r0_paris/N_paris*v*I_paris)*S_paris)*((sto)**2))*((r0_paris/N_paris*v*I_paris)*S_paris)))'
        term_nyc = '((((r0_nyc/N_nyc*v*I_nyc)*S_nyc)*((sto)**2))*((r0_nyc/N_nyc*v*I_nyc)*S_nyc)))'
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][0][0],'((1)*'+term_nyc+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][0][1],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][0][2],'((1)*'+term_nyc+'*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][0][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][0][4],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][0][5],'((1)*'+term_nyc+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][1][0],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][1][1],'((1)*'+term_paris+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][1][2],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][1][3],'((1)*'+term_paris+'*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][1][4],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][1][5],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][2][0],'((-1)*'+term_nyc+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][2][1],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][2][2],'((-1)*'+term_nyc+'*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][2][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][2][4],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][2][5],'((-1)*'+term_nyc+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][3][0],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][3][1],'((-1)*'+term_paris+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][3][2],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][3][3],'((-1)*'+term_paris+'*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][3][4],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][3][5],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][4][0],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][4][1],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][4][2],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][4][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][4][4],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][4][5],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][5][0],'((1)*'+term_nyc+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][5][1],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][5][2],'((1)*'+term_nyc+'*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][5][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][5][4],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][5][5],'((1)*'+term_nyc+'*(1)')

        # testing dem sto only
        term1_p = '(mu_b_paris*N_paris))'
        term2_p = '((r0_paris/N_paris*v*I_paris)*S_paris))'
        term3_p = '((correct_rate(v))*I_paris))'
        term4_p = '((mu_d_paris)*S_paris))'
        term5_p = '((mu_d_paris)*I_paris))'
        term1_n = '(mu_b_nyc*N_nyc))'
        term2_n = '((r0_nyc/N_nyc*v*I_nyc)*S_nyc))'
        term3_n = '((correct_rate(v))*I_nyc))'
        term4_n = '((mu_d_nyc)*S_nyc))'
        term5_n = '((mu_d_nyc)*I_nyc))'
        # order: I_nyc, I_paris, S_nyc, S_paris, inc_all , inc_nyc
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][0][0],'((1)*'+term2_n+'*(1) + ((-1)*'+term3_n+'*(-1) + ((-1)*'+term5_n+'*(-1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][0][1],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][0][2],'((1)*'+term2_n+'*(-1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][0][3],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][0][4],'((-1)*'+term3_n+'*(1) + ((-1)*'+term5_n+'*(1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][0][5],'((1)*'+term2_n+'*(1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][1][0],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][1][1],'((1)*'+term2_p+'*(1) + ((-1)*'+term3_p+'*(-1) + ((-1)*'+term5_p+'*(-1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][1][2],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][1][3],'((1)*'+term2_p+'*(-1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][1][4],'((-1)*'+term3_p+'*(1) + ((-1)*'+term5_p+'*(1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][1][5],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][2][0],'((-1)*'+term2_n+'*(1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][2][1],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][2][2],'((1)*'+term1_n+'*(1) + ((-1)*'+term2_n+'*(-1) + ((-1)*'+term4_n+'*(-1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][2][3],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][2][4],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][2][5],'((-1)*'+term2_n+'*(1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][3][0],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][3][1],'((-1)*'+term2_p+'*(1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][3][2],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][3][3],'((1)*'+term1_p+'*(1) + ((-1)*'+term2_p+'*(-1) + ((-1)*'+term4_p+'*(-1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][3][4],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][3][5],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][4][0],'((1)*'+term3_n+'*(-1) + ((1)*'+term5_n+'*(-1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][4][1],'((1)*'+term3_p+'*(-1) + ((1)*'+term5_p+'*(-1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][4][2],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][4][3],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][4][4],'((1)*'+term3_p+'*(1) + ((1)*'+term3_n+'*(1) + ((1)*'+term5_p+'*(1) + ((1)*'+term5_n+'*(1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][4][5],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][5][0],'((1)*'+term2_n+'*(1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][5][1],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][5][2],'((1)*'+term2_n+'*(-1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][5][3],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][5][4],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q_cm"][5][5],'((1)*'+term2_n+'*(1)')
        

    def test_eval_Q_tricky_cases(self):

        calc_Q = self.m_noise2.eval_Q()
        # testing env sto only for m_noise2 : WN on U->S
        term_p = '(((mu_b_paris*N_paris)*((sto)**2))*(mu_b_paris*N_paris)))'
        term_n = '(((mu_b_nyc*N_nyc)*((sto)**2))*(mu_b_nyc*N_nyc)))'
        for i in range(5):
            for j in range(5):
                if i==2 and j == 2:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],'((1)*'+term_n+'*(1)')
                elif i==3 and j == 3:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],'((1)*'+term_p+'*(1)')
                else:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],0)


        calc_Q = self.m_noise3.eval_Q()
        # testing env sto only for m_noise3 : WN on I->R
        term_p = '((((correct_rate(v))*I_paris)*((sto)**2))*((correct_rate(v))*I_paris)))'
        term_n = '((((correct_rate(v))*I_nyc)*((sto)**2))*((correct_rate(v))*I_nyc)))'
        for i in range(5):
            for j in range(5):
                if i==0 and j == 0:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],'((-1)*'+term_n+'*(-1)')
                elif i==1 and j == 1:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],'((-1)*'+term_p+'*(-1)')
                elif i==0  and j == 4:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],'((-1)*'+term_n+'*(1)')
                elif i==4  and j == 0:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],'((1)*'+term_n+'*(-1)')
                elif i==1  and j == 4:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],'((-1)*'+term_p+'*(1)')
                elif i==4  and j == 1:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],'((1)*'+term_p+'*(-1)')
                elif i==4  and j == 4:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],'((1)*'+term_p+'*(1) + ((1)*'+term_n+'*(1)')
                else:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],0)
        

        calc_Q = self.m_noise4.eval_Q()
        # testing env sto only for m_noise4 : WN on I->U
        term_p = '((((mu_d_paris)*I_paris)*((sto)**2))*((mu_d_paris)*I_paris)))'
        term_n = '((((mu_d_nyc)*I_nyc)*((sto)**2))*((mu_d_nyc)*I_nyc)))'
        for i in range(5):
            for j in range(5):
                if i==0 and j == 0:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],'((-1)*'+term_n+'*(-1)')
                elif i==1 and j == 1:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],'((-1)*'+term_p+'*(-1)')
                elif i==1  and j == 4:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],'((-1)*'+term_p+'*(1)')
                elif i==4  and j == 1:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],'((1)*'+term_p+'*(-1)')
                elif i==0  and j == 4:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],'((-1)*'+term_n+'*(1)')
                elif i==4  and j == 0:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],'((1)*'+term_n+'*(-1)')
                elif i==4  and j == 4:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],'((1)*'+term_p+'*(1) + ((1)*'+term_n+'*(1)')
                else:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],0)

                    

        calc_Q = self.m_noise5.eval_Q()
        # testing env sto only for m_noise5 : WN on R->U
        for i in range(5):
            for j in range(5):
                self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][i][j],0)

        calc_Q = self.m_noise6.eval_Q()
        # testing env sto only for m_noise6 : correlated WN on I->R and S->I
        term1_p = '((r0_paris/N_paris*v*I_paris)*S_paris)'
        term2_p = '((correct_rate(v))*I_paris)'
        term1_n = '((r0_nyc/N_nyc*v*I_nyc)*S_nyc)'
        term2_n = '((correct_rate(v))*I_nyc)'
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][0][0],'((1)*(('+term1_n+'*((sto)**2))*'+term1_n+') + (-1)*(('+term2_n+'*((sto)**2))*'+term1_n+'))*(1) + ((1)*(('+term1_n+'*((sto)**2))*'+term2_n+') + (-1)*(('+term2_n+'*((sto)**2))*'+term2_n+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][0][1],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][0][2],'((1)*(('+term1_n+'*((sto)**2))*'+term1_n+') + (-1)*(('+term2_n+'*((sto)**2))*'+term1_n+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][0][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][0][4],'((1)*(('+term1_n+'*((sto)**2))*'+term2_n+') + (-1)*(('+term2_n+'*((sto)**2))*'+term2_n+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][0][5],'((1)*(('+term1_n+'*((sto)**2))*'+term1_n+') + (-1)*(('+term2_n+'*((sto)**2))*'+term1_n+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][1][0],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][1][1],'((1)*(('+term1_p+'*((sto)**2))*'+term1_p+') + (-1)*(('+term2_p+'*((sto)**2))*'+term1_p+'))*(1) + ((1)*(('+term1_p+'*((sto)**2))*'+term2_p+') + (-1)*(('+term2_p+'*((sto)**2))*'+term2_p+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][1][2],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][1][3],'((1)*(('+term1_p+'*((sto)**2))*'+term1_p+') + (-1)*(('+term2_p+'*((sto)**2))*'+term1_p+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][1][4],'((1)*(('+term1_p+'*((sto)**2))*'+term2_p+') + (-1)*(('+term2_p+'*((sto)**2))*'+term2_p+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][1][5],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][2][0],'((-1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1) + ((-1)*(('+term1_n+'*((sto)**2))*'+term2_n+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][2][1],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][2][2],'((-1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][2][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][2][4],'((-1)*(('+term1_n+'*((sto)**2))*'+term2_n+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][2][5],'((-1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][3][0],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][3][1],'((-1)*(('+term1_p+'*((sto)**2))*'+term1_p+'))*(1) + ((-1)*(('+term1_p+'*((sto)**2))*'+term2_p+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][3][2],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][3][3],'((-1)*(('+term1_p+'*((sto)**2))*'+term1_p+'))*(-1)')
        
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][3][4],'((-1)*(('+term1_p+'*((sto)**2))*'+term2_p+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][3][5],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][4][0],'((1)*(('+term2_n+'*((sto)**2))*'+term1_n+'))*(1) + ((1)*(('+term2_n+'*((sto)**2))*'+term2_n+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][4][1],'((1)*(('+term2_p+'*((sto)**2))*'+term1_p+'))*(1) + ((1)*(('+term2_p+'*((sto)**2))*'+term2_p+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][4][2],'((1)*(('+term2_n+'*((sto)**2))*'+term1_n+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][4][3],'((1)*(('+term2_p+'*((sto)**2))*'+term1_p+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][4][4],'((1)*(('+term2_p+'*((sto)**2))*'+term2_p+'))*(1) + ((1)*(('+term2_n+'*((sto)**2))*'+term2_n+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][4][5],'((1)*(('+term2_n+'*((sto)**2))*'+term1_n+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][5][0],'((1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1) + ((1)*(('+term1_n+'*((sto)**2))*'+term2_n+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][5][1],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][5][2],'((1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][5][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][5][4],'((1)*(('+term1_n+'*((sto)**2))*'+term2_n+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][5][5],'((1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1)')
        
        calc_Q = self.m_noise7.eval_Q()
        # testing env sto only for m_noise7 : uncorrelated WN on I->R and S->I
        term1_p = '((r0_paris/N_paris*v*I_paris)*S_paris)'
        term2_p = '((correct_rate(v))*I_paris)'
        term1_n = '((r0_nyc/N_nyc*v*I_nyc)*S_nyc)'
        term2_n = '((correct_rate(v))*I_nyc)'
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][0][0],'((1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1) + ((-1)*(('+term2_n+'*((sto)**2))*'+term2_n+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][0][1],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][0][2],'((1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][0][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][0][4],'((-1)*(('+term2_n+'*((sto)**2))*'+term2_n+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][0][5],'((1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][1][0],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][1][1],'((1)*(('+term1_p+'*((sto)**2))*'+term1_p+'))*(1) + ((-1)*(('+term2_p+'*((sto)**2))*'+term2_p+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][1][2],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][1][3],'((1)*(('+term1_p+'*((sto)**2))*'+term1_p+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][1][4],'((-1)*(('+term2_p+'*((sto)**2))*'+term2_p+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][1][5],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][2][0],'((-1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][2][1],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][2][2],'((-1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][2][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][2][4],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][2][5],'((-1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][3][0],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][3][1],'((-1)*(('+term1_p+'*((sto)**2))*'+term1_p+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][3][2],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][3][3],'((-1)*(('+term1_p+'*((sto)**2))*'+term1_p+'))*(-1)')
        
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][3][4],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][3][5],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][4][0],'((1)*(('+term2_n+'*((sto)**2))*'+term2_n+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][4][1],'((1)*(('+term2_p+'*((sto)**2))*'+term2_p+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][4][2],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][4][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][4][4],'((1)*(('+term2_p+'*((sto)**2))*'+term2_p+'))*(1) + ((1)*(('+term2_n+'*((sto)**2))*'+term2_n+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][4][5],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][5][0],'((1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][5][1],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][5][2],'((1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][5][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][5][4],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q_cm"][5][5],'((1)*(('+term1_n+'*((sto)**2))*'+term1_n+'))*(1)')

    def test_jac(self):
        step_ode_sde = self.m_noise.step_ode_sde()
        jac = self.m_diff.jac(step_ode_sde['sf'])

        # testing jac
        # I ode - ((v)*I) - ((mu_d)*I) + ((r0/N*v*I)*S)
        self.assertEqual(jac["caches"][jac["jac"][0][0]],"-gsl_spline_eval(calc->spline[ORDER_mu_d_nyc],t,calc->acc[ORDER_mu_d_nyc])-(gsl_vector_get(par, ORDER_v))+X[ORDER_S_nyc]*diffed[ORDER_diff__r0_nyc]*gsl_vector_get(par, ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_nyc],t,calc->acc[ORDER_N_nyc])")
        self.assertEqual(jac["caches"][jac["jac"][1][1]],"-gsl_spline_eval(calc->spline[ORDER_mu_d_paris],t,calc->acc[ORDER_mu_d_paris])-(gsl_vector_get(par, ORDER_v))+X[ORDER_S_paris]*diffed[ORDER_diff__r0_paris]*gsl_vector_get(par, ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_paris],t,calc->acc[ORDER_N_paris])")
        self.assertEqual(jac["caches"][jac["jac_diff"][1][0]["value"]],"0")
        self.assertEqual(jac["caches"][jac["jac_diff"][0][1]["value"]],"0")
        
        # S ode - ((r0/N*v*I)*S) - ((mu_d)*S) + (mu_b*N)
        self.assertEqual(jac["caches"][jac["jac"][2][2]],"-X[ORDER_I_nyc]*diffed[ORDER_diff__r0_nyc]*gsl_vector_get(par, ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_nyc],t,calc->acc[ORDER_N_nyc])-gsl_spline_eval(calc->spline[ORDER_mu_d_nyc],t,calc->acc[ORDER_mu_d_nyc])")
        self.assertEqual(jac["caches"][jac["jac"][3][3]],"-X[ORDER_I_paris]*diffed[ORDER_diff__r0_paris]*gsl_vector_get(par, ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_paris],t,calc->acc[ORDER_N_paris])-gsl_spline_eval(calc->spline[ORDER_mu_d_paris],t,calc->acc[ORDER_mu_d_paris])")
        self.assertEqual(jac["caches"][jac["jac"][2][0]],"-X[ORDER_S_nyc]*diffed[ORDER_diff__r0_nyc]*gsl_vector_get(par, ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_nyc],t,calc->acc[ORDER_N_nyc])")
        self.assertEqual(jac["caches"][jac["jac"][3][1]],"-X[ORDER_S_paris]*diffed[ORDER_diff__r0_paris]*gsl_vector_get(par, ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_paris],t,calc->acc[ORDER_N_paris])")
        
        
        # testing jac_obs
        # all_inc
        self.assertEqual(jac["caches"][jac["jac_obs"][0][0]],"gsl_spline_eval(calc->spline[ORDER_mu_d_nyc],t,calc->acc[ORDER_mu_d_nyc])+(gsl_vector_get(par, ORDER_v))")
        self.assertEqual(jac["caches"][jac["jac_obs"][0][1]],"gsl_spline_eval(calc->spline[ORDER_mu_d_paris],t,calc->acc[ORDER_mu_d_paris])+(gsl_vector_get(par, ORDER_v))")
        self.assertEqual(jac["caches"][jac["jac_obs_diff"][0][0]["value"]],"0")
        self.assertEqual(jac["caches"][jac["jac_obs_diff"][0][1]["value"]],"0")
        # nyc_inc
        self.assertEqual(jac["caches"][jac["jac_obs"][1][0]],"X[ORDER_S_nyc]*diffed[ORDER_diff__r0_nyc]*gsl_vector_get(par, ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_nyc],t,calc->acc[ORDER_N_nyc])")
        self.assertEqual(jac["caches"][jac["jac_obs"][1][1]],"0")
        self.assertEqual(jac["caches"][jac["jac_obs_diff"][1][0]["value"]],"X[ORDER_I_nyc]*X[ORDER_S_nyc]*gsl_vector_get(par, ORDER_v)/gsl_spline_eval(calc->spline[ORDER_N_nyc],t,calc->acc[ORDER_N_nyc])")
        self.assertEqual(jac["caches"][jac["jac_obs_diff"][1][1]["value"]],"0")
        

if __name__ == '__main__':
    unittest.main()

