import unittest
import os
import subprocess
import shutil
import json
import numpy
from numpy import genfromtxt
from scipy import stats
from Builder import Builder
import math
import csv

class TestNoiseResults(unittest.TestCase):
      @classmethod
      def setUpClass(cls):
            print("")
            print("Testing classic numerical results on noise example")
            # copy noise from the examples and build it
            b = Builder(os.path.join(os.getenv("HOME"), 'ssm_test_model'), os.path.join('..' ,'example', 'noise', 'datapackages', 'model-jdureau-noise', 'datapackage.json'))
            b.prepare()
            b.code()
            b.write_data()

            os.chdir('/Users/dureaujoseph/ssm_test_model/C/templates')
            os.system('make clean')
            os.system('make')
            os.system('make install')
            os.chdir('/Users/dureaujoseph/ssm_test_model/')

      
      def setUp(self):
      # Things that need to be done before tests
            os.chdir('/Users/dureaujoseph/ssm_test_model/')

      @classmethod
      def tearDownClass(cls):
            os.chdir('/Users/dureaujoseph/ssm_test_model/')
            #            shutil.rmtree('/Users/dureaujoseph/ssm_test_model/')

      def test_kalman_map(self):
            os.system('./ksimplex --prior -M 1000 -c < /Users/dureaujoseph/ssm/example/noise/datapackage.json')
            tab = genfromtxt('trace_0.csv',delimiter=',',names=True)
            self.assertAlmostEqual(tab[374][5],-327.031)

class TestTransfsAndPMCMC(unittest.TestCase):
      @classmethod
      def setUpClass(cls):
            print("")
            print("Testing transformations and pMCMC by sampling from priors")
            # copy noise from the examples and build it
            #b = Builder(os.path.join(os.getenv("HOME"), 'ssm_test_model'), os.path.join('..' ,'example', 'noise', 'datapackages', 'model-jdureau-noise', 'datapackage.json'))
            #b.prepare()
            #b.code()
            #b.write_data()

            #os.chdir('/Users/dureaujoseph/ssm_test_model/C/templates')
            #os.system('make clean')
            #os.system('make')
            #os.system('make install')
            #os.chdir('/Users/dureaujoseph/ssm_test_model/')
            
      def setUp(self):
      # Things that need to be done before tests
            os.chdir('/Users/dureaujoseph/ssm_test_model/')

      @classmethod
      def tearDownClass(cls):
            os.chdir('/Users/dureaujoseph/ssm_test_model/')
            #shutil.rmtree('/Users/dureaujoseph/ssm_test_model/')


      def test_prior_unif_transf_log(self):
            os.chdir('/Users/dureaujoseph/ssm/example')
            shutil.copytree('noise','noise_test')
            f = open('/Users/dureaujoseph/ssm/example/noise_test/datapackage.json')
            j = json.load(f)
            j['resources'][0]['data']['r0_paris'] = 10
            with open('/Users/dureaujoseph/ssm/example/noise_test/datapackage.json','w') as outfile:
                  json.dump(j,outfile)
            f = open('/Users/dureaujoseph/ssm/example/noise_test/datapackages/model-jdureau-noise/datapackage.json')
            j = json.load(f)
            j['resources'][3]['data'][4]['prior']['upper']=10.5
            j['resources'][3]['data'][4]['prior']['lower']=9.5
            with open('/Users/dureaujoseph/ssm/example/noise_test/datapackages/model-jdureau-noise/datapackage.json','w') as outfile:
                  json.dump(j,outfile)


            os.chdir('/Users/dureaujoseph/ssm/ssm/')
            b = Builder(os.path.join(os.getenv("HOME"), 'ssm_test_model'), os.path.join('..' ,'example', 'noise_test', 'datapackages', 'model-jdureau-noise', 'datapackage.json'))
            b.prepare()
            b.code()
            b.write_data()

            os.chdir('/Users/dureaujoseph/ssm_test_model/C/templates')
            os.system('make clean')
            os.system('make')
            os.system('make install')
            
            os.chdir('/Users/dureaujoseph/ssm_test_model/')
            
            os.system('./pmcmc ode  -C 5000 -W 5000 -O 0 -M 20000 -c -a < /Users/dureaujoseph/ssm/example/noise_test/datapackage.json')
            
            shutil.rmtree('/Users/dureaujoseph/ssm/example/noise_test')
            test = self.call_test_unif('r0_paris')
            self.assertEqual(test,'1')



      def test_prior_normal_transf_log(self):
            os.chdir('/Users/dureaujoseph/ssm/example')
            shutil.copytree('noise','noise_test')
            f = open('/Users/dureaujoseph/ssm/example/noise_test/datapackage.json')
            j = json.load(f)
            j['resources'][0]['data']['r0_paris'] = 10
            with open('/Users/dureaujoseph/ssm/example/noise_test/datapackage.json','w') as outfile:
                  json.dump(j,outfile)
            f = open('/Users/dureaujoseph/ssm/example/noise_test/datapackages/model-jdureau-noise/datapackage.json')
            j = json.load(f)
            j['resources'][3]['data'][4]['prior']['distribution']='normal'
            j['resources'][3]['data'][4]['prior']['lower']=0
            j['resources'][3]['data'][4]['prior']['mean']=10
            j['resources'][3]['data'][4]['prior']['sd']=1
            del j['resources'][3]['data'][4]['prior']['upper']
            del j['resources'][3]['data'][4]['prior']['lower']
            with open('/Users/dureaujoseph/ssm/example/noise_test/datapackages/model-jdureau-noise/datapackage.json','w') as outfile:
                  json.dump(j,outfile)

            os.chdir('/Users/dureaujoseph/ssm/ssm/')
            b = Builder(os.path.join(os.getenv("HOME"), 'ssm_test_model'), os.path.join('..' ,'example', 'noise_test', 'datapackages', 'model-jdureau-noise', 'datapackage.json'))
            b.prepare()
            b.code()
            b.write_data()

            os.chdir('/Users/dureaujoseph/ssm_test_model/C/templates')
            os.system('make clean')
            os.system('make')
            os.system('make install')

            os.chdir('/Users/dureaujoseph/ssm_test_model/')
            
            os.system('./pmcmc ode  -C 5000 -W 5000 -O 0 -M 20000 -c -a < /Users/dureaujoseph/ssm/example/noise_test/datapackage.json')

            shutil.rmtree('/Users/dureaujoseph/ssm/example/noise_test')
            test = self.call_test_normal('r0_paris')
            self.assertEqual(test,'1')





      def test_prior_normal_and_unif_transf_log(self):
            os.chdir('/Users/dureaujoseph/ssm/example')
            shutil.copytree('noise','noise_test')
            f = open('/Users/dureaujoseph/ssm/example/noise_test/datapackage.json')
            j = json.load(f)
            j['resources'][0]['data']['r0_paris'] = 10
            j['resources'][0]['data']['r0_nyc'] = 10
            with open('/Users/dureaujoseph/ssm/example/noise_test/datapackage.json','w') as outfile:
                  json.dump(j,outfile)
            f = open('/Users/dureaujoseph/ssm/example/noise_test/datapackages/model-jdureau-noise/datapackage.json')
            j = json.load(f)
            j['resources'][3]['data'][4]['prior']['distribution']='normal'
            j['resources'][3]['data'][4]['prior']['lower']=0
            j['resources'][3]['data'][4]['prior']['mean']=10
            j['resources'][3]['data'][4]['prior']['sd']=1
            del j['resources'][3]['data'][4]['prior']['upper']
            del j['resources'][3]['data'][4]['prior']['lower']
            j['resources'][3]['data'][5]['prior']['lower']=9.5
            j['resources'][3]['data'][5]['prior']['upper']=10.5
            with open('/Users/dureaujoseph/ssm/example/noise_test/datapackages/model-jdureau-noise/datapackage.json','w') as outfile:
                  json.dump(j,outfile)

            os.chdir('/Users/dureaujoseph/ssm/ssm/src/C')
            os.system('make clean')
            os.system('make')
            os.system('make install')

            os.chdir('/Users/dureaujoseph/ssm/ssm/')
            b = Builder(os.path.join(os.getenv("HOME"), 'ssm_test_model'), os.path.join('..' ,'example', 'noise_test', 'datapackages', 'model-jdureau-noise', 'datapackage.json'))
            b.prepare()
            b.code()
            b.write_data()

            os.chdir('/Users/dureaujoseph/ssm_test_model/C/templates')
            os.system('make clean')
            os.system('make')
            os.system('make install')

            os.chdir('/Users/dureaujoseph/ssm_test_model/')
            
            os.system('./pmcmc ode  -C 1000 -W 100000 -O 0 -M 100000 -c -a < /Users/dureaujoseph/ssm/example/noise_test/datapackage.json')

            shutil.rmtree('/Users/dureaujoseph/ssm/example/noise_test')
            test1 = self.call_test_normal('r0_paris')
            test2 = self.call_test_unif('r0_nyc')
            print test1
            print test2
            self.assertEqual(int(test1)*int(test2),1)



      def test_prior_normal_transf_logit_ab(self):
            os.chdir('/Users/dureaujoseph/ssm/example')
            shutil.copytree('noise','noise_test')
            f = open('/Users/dureaujoseph/ssm/example/noise_test/datapackage.json')
            j = json.load(f)
            j['resources'][0]['data']['r0_paris'] = 10
            with open('/Users/dureaujoseph/ssm/example/noise_test/datapackage.json','w') as outfile:
                  json.dump(j,outfile)
            f = open('/Users/dureaujoseph/ssm/example/noise_test/datapackages/model-jdureau-noise/datapackage.json')
            j = json.load(f)
            j['resources'][3]['data'][4]['prior']['distribution']='normal'
            j['resources'][3]['data'][4]['prior']['lower']=0
            j['resources'][3]['data'][4]['prior']['upper']=20
            j['resources'][3]['data'][4]['prior']['mean']=10
            j['resources'][3]['data'][4]['prior']['sd']=1
            with open('/Users/dureaujoseph/ssm/example/noise_test/datapackages/model-jdureau-noise/datapackage.json','w') as outfile:
                  json.dump(j,outfile)

            os.chdir('/Users/dureaujoseph/ssm/ssm/')
            b = Builder(os.path.join(os.getenv("HOME"), 'ssm_test_model'), os.path.join('..' ,'example', 'noise_test', 'datapackages', 'model-jdureau-noise', 'datapackage.json'))
            b.prepare()
            b.code()
            b.write_data()

            os.chdir('/Users/dureaujoseph/ssm_test_model/C/templates')
            os.system('make clean')
            os.system('make')
            os.system('make install')

            os.chdir('/Users/dureaujoseph/ssm_test_model/')
            
            os.system('./pmcmc ode  -C 5000 -W 5000 -O 0 -M 20000 -c -a < /Users/dureaujoseph/ssm/example/noise_test/datapackage.json')

            shutil.rmtree('/Users/dureaujoseph/ssm/example/noise_test')
            test = self.call_test_normal('r0_paris')
            self.assertEqual(test,'1') 


      def test_prior_normal_transf_identity(self):
            os.chdir('/Users/dureaujoseph/ssm/example')
            shutil.copytree('noise','noise_test')
            f = open('/Users/dureaujoseph/ssm/example/noise_test/datapackage.json')
            j = json.load(f)
            j['resources'][0]['data']['r0_paris'] = 10
            with open('/Users/dureaujoseph/ssm/example/noise_test/datapackage.json','w') as outfile:
                  json.dump(j,outfile)
            f = open('/Users/dureaujoseph/ssm/example/noise_test/datapackages/model-jdureau-noise/datapackage.json')
            j = json.load(f)
            j['resources'][3]['data'][4]['prior']['distribution']='normal'
            j['resources'][3]['data'][4]['prior']['mean']=10
            j['resources'][3]['data'][4]['prior']['sd']=1
            del j['resources'][3]['data'][4]['prior']['upper']
            del j['resources'][3]['data'][4]['prior']['lower']
            with open('/Users/dureaujoseph/ssm/example/noise_test/datapackages/model-jdureau-noise/datapackage.json','w') as outfile:
                  json.dump(j,outfile)


            os.chdir('/Users/dureaujoseph/ssm/ssm/')
            b = Builder(os.path.join(os.getenv("HOME"), 'ssm_test_model'), os.path.join('..' ,'example', 'noise_test', 'datapackages', 'model-jdureau-noise', 'datapackage.json'))
            b.prepare()
            b.code()
            b.write_data()

            os.chdir('/Users/dureaujoseph/ssm_test_model/C/templates')
            os.system('make clean')
            os.system('make')
            os.system('make install')

            os.chdir('/Users/dureaujoseph/ssm_test_model/')
            
            os.system('./pmcmc ode  -C 5000 -W 5000 -O 0 -M 20000 -c -a < /Users/dureaujoseph/ssm/example/noise_test/datapackage.json')

            shutil.rmtree('/Users/dureaujoseph/ssm/example/noise_test')
            test = self.call_test_normal('r0_paris')
            self.assertEqual(test,'1') 

      def call_test_unif(self, varname):
            shutil.copyfile(Root+'/TestsR/test_unif.r','/Users/dureaujoseph/ssm_test_model/test_unif.r')
            os.chdir('/Users/dureaujoseph/ssm_test_model/')
            os.system('R --vanilla < test_unif.r ' + varname  + '> /dev/null 2>&1')
            f = open("outfile.txt","r")
            x = f.readlines()
            return x[0]

      def call_test_normal(self,varname):
            shutil.copyfile(Root+'/TestsR/test_normal.r','/Users/dureaujoseph/ssm_test_model/test_normal.r')
            os.chdir('/Users/dureaujoseph/ssm_test_model/')
            os.system('R --vanilla < test_normal.r ' + varname  + '> /dev/null 2>&1')
            f = open("outfile.txt","r")
            x = f.readlines()
            return x[0]

class TestKalmanOnDiffusions(unittest.TestCase):
      @classmethod
      def setUpClass(cls):
            print("")
            print("Testing EKF on linear examples")

            # copy noise from the examples, add a diffusing parameter for tests, and build it
            
            os.chdir('/Users/dureaujoseph/ssm/example')
            shutil.copytree('noise','noise_test')
            f = open('/Users/dureaujoseph/ssm/example/noise_test/datapackages/model-jdureau-noise/datapackage.json')
            j = json.load(f)
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][23]['name'] = 'test_par'
            j['resources'][3]['data'][23]['prior'] = {}
            j['resources'][3]['data'][23]['prior']['distribution']='fixed'
            j['resources'][3]['data'][23]['prior']['value']=0
            j['resources'][3]['data'][24]['name'] = 'test_vol'
            j['resources'][3]['data'][24]['prior'] = {}
            j['resources'][3]['data'][24]['prior']['distribution']='fixed'
            j['resources'][3]['data'][24]['prior']['value']= 0.1428571 #0.3779645 # 1/sqrt(7)
            j['resources'].append({})
            j['resources'][4]['name']='sde'
            j['resources'][4]['data']={}
            j['resources'][4]['data']['drift']=[]
            j['resources'][4]['data']['drift'].append({})
            j['resources'][4]['data']['drift'][0]['name']='test_par'
            j['resources'][4]['data']['drift'][0]['f']=0.0
            j['resources'][4]['data']['dispersion']=[['test_vol']]
            with open('/Users/dureaujoseph/ssm/example/noise_test/datapackages/model-jdureau-noise/datapackage.json','w') as outfile:
                  json.dump(j,outfile)

            os.chdir('/Users/dureaujoseph/ssm/ssm/src/C/')
            os.system('make clean')
            os.system('make')
            os.system('make install')
            
            os.chdir('/Users/dureaujoseph/ssm/ssm/')
            b = Builder(os.path.join(os.getenv("HOME"), 'ssm_test_model'), os.path.join('..' ,'example', 'noise_test', 'datapackages', 'model-jdureau-noise', 'datapackage.json'))
            b.prepare()
            b.code()
            b.write_data()

            os.chdir('/Users/dureaujoseph/ssm_test_model/C/templates')
            os.system('make clean')
            os.system('make')
            os.system('make install')


            os.chdir('/Users/dureaujoseph/ssm_test_model/')

            
      
      def setUp(self):
      # Things that need to be done before tests
            os.chdir('/Users/dureaujoseph/ssm_test_model/')

      @classmethod
      def tearDownClass(cls):
            os.chdir('/Users/dureaujoseph/ssm_test_model/')
            shutil.rmtree('/Users/dureaujoseph/ssm/example/noise_test')
            
      def test_1step(self):

            os.system('./kalman -O 2  -c -x < /Users/dureaujoseph/ssm/example/noise_test/datapackage.json')
            tab = genfromtxt('hat_0.csv',delimiter=',',names=True)
            self.assertAlmostEqual(tab['lower_test_par'][0],-1.96/math.sqrt(7),5)

      def test_10step(self):
            os.system('./kalman -O 10  -c -x < /Users/dureaujoseph/ssm/example/noise_test/datapackage.json')
            tab = genfromtxt('hat_0.csv',delimiter=',',names=True)
            self.assertAlmostEqual(tab['lower_test_par'][9]/math.sqrt(10),-1.96/math.sqrt(7),5)

class TestSMCSDEagainstKalman(unittest.TestCase):
      @classmethod
      def setUpClass(cls):
            print("")
            print("Testing SMC against EKF on linear examples")

            os.chdir('/Users/dureaujoseph/ssm/example')
            shutil.copytree('noise','noise_test')
            f = open('/Users/dureaujoseph/ssm/example/noise_test/datapackages/model-jdureau-noise/datapackage.json')
            j = json.load(f)
            j['resources'][0]["data"] = []
            j['resources'][0]["data"].append({"name":"paris","composition":["S","I","I2","E1","E2"]})
            j['resources'][1]["data"] = []
            j['resources'][1]["data"].append({"from":"S", "to":"I", "rate":"1000/S", "white_noise": {"name":"noise_SI", "sd": "vol"}})
            j['resources'][1]["data"].append({"from":"I", "to":"I2", "rate":"1000/I", "white_noise": {"name":"noise_SI", "sd": "vol"}})
            j['resources'][1]["data"].append({"from":"I2", "to":"S", "rate":"1000/I2", "white_noise": {"name":"noise_SI2", "sd": "vol2"}})
            j['resources'][1]["data"].append({"from":"E1", "to":"E2", "rate":"3/N", "tracked": ["all_inc_out"]})

            j['resources'][2]["data"] = [j['resources'][2]["data"][0]]

            j['resources'][3]["data"] = []
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][0]['name'] = 'r0'
            j['resources'][3]['data'][0]['prior'] = {}
            j['resources'][3]['data'][0]['prior']['distribution']='uniform'
            j['resources'][3]['data'][0]['prior']['lower']=50
            j['resources'][3]['data'][0]['prior']['upper']=150
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][1]['name'] = 'S'
            j['resources'][3]['data'][1]['prior'] = {}
            j['resources'][3]['data'][1]['prior']['distribution']='fixed'
            j['resources'][3]['data'][1]['prior']['value']=200000
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][2]['name'] = 'I'
            j['resources'][3]['data'][2]['prior'] = {}
            j['resources'][3]['data'][2]['prior']['distribution']='fixed'
            j['resources'][3]['data'][2]['prior']['value']=200000
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][3]['name'] = 'I2'
            j['resources'][3]['data'][3]['prior'] = {}
            j['resources'][3]['data'][3]['prior']['distribution']='fixed'
            j['resources'][3]['data'][3]['prior']['value']=200000
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][4]['name'] = 'E1'
            j['resources'][3]['data'][4]['prior'] = {}
            j['resources'][3]['data'][4]['prior']['distribution']='fixed'
            j['resources'][3]['data'][4]['prior']['value']=200000
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][5]['name'] = 'vol'
            j['resources'][3]['data'][5]['prior'] = {}
            j['resources'][3]['data'][5]['prior']['distribution']='fixed'
            j['resources'][3]['data'][5]['prior']['value']=0.01428571
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][6]['name'] = 'vol2'
            j['resources'][3]['data'][6]['prior'] = {}
            j['resources'][3]['data'][6]['prior']['distribution']='fixed'
            j['resources'][3]['data'][6]['prior']['value']=0.02857143
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][7]['name'] = 'rep_all_CDC_inc'
            j['resources'][3]['data'][7]['prior'] = {}
            j['resources'][3]['data'][7]['prior']['distribution']='fixed'
            j['resources'][3]['data'][7]['prior']['value']=0.2
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][8]['name'] = 'prop_all_CDC_inc'
            j['resources'][3]['data'][8]['prior'] = {}
            j['resources'][3]['data'][8]['prior']['distribution']='fixed'
            j['resources'][3]['data'][8]['prior']['value']=0.2
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][9]['name'] = 'phi'
            j['resources'][3]['data'][9]['prior'] = {}
            j['resources'][3]['data'][9]['prior']['distribution']='fixed'
            j['resources'][3]['data'][9]['prior']['value']=0.2
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][10]['name'] = 'E2'
            j['resources'][3]['data'][10]['prior'] = {}
            j['resources'][3]['data'][10]['prior']['distribution']='fixed'
            j['resources'][3]['data'][10]['prior']['value']=200000
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][11]['name'] = 'N'
            j['resources'][3]['data'][11]['prior'] = {}
            j['resources'][3]['data'][11]['prior']['distribution']='fixed'
            j['resources'][3]['data'][11]['prior']['value']=0.2


            j['resources'].append({})
            j['resources'][4]['name']='sde'
            j['resources'][4]['data']={}
            j['resources'][4]['data']['drift']=[]
            j['resources'][4]['data']['drift'].append({})
            j['resources'][4]['data']['drift'][0]['name']='r0'
            j['resources'][4]['data']['drift'][0]['f']=0.0
            j['resources'][4]['data']['dispersion']=[['vol']]
            with open('/Users/dureaujoseph/ssm/example/noise_test/datapackages/model-jdureau-noise/datapackage.json','w') as outfile:
                  json.dump(j,outfile)

            f = open('/Users/dureaujoseph/ssm/example/noise_test/datapackage.json')
            j = json.load(f)
            j['resources'][0]["data"]={}
            j['resources'][0]["data"]["r0"] = 100
            j['resources'][1]["data"]={}
            j['resources'][1]["data"]["r0"] = {"r0":0.02}
            with open('/Users/dureaujoseph/ssm/example/noise_test/datapackage.json','w') as outfile:
                  json.dump(j,outfile)

            os.chdir('/Users/dureaujoseph/ssm/ssm/src/C/')
            os.system('make clean')
            os.system('make')
            os.system('make install')
            
            os.chdir('/Users/dureaujoseph/ssm/ssm/')
            b = Builder(os.path.join(os.getenv("HOME"), 'ssm_test_model'), os.path.join('..' ,'example', 'noise_test', 'datapackages', 'model-jdureau-noise', 'datapackage.json'))
            b.prepare()
            b.code()
            b.write_data()

            os.chdir('/Users/dureaujoseph/ssm_test_model/C/templates')
            os.system('make clean')
            os.system('make')
            os.system('make install')


            os.chdir('/Users/dureaujoseph/ssm_test_model/')

      def setUp(self):
      # Things that need to be done before tests
            os.chdir('/Users/dureaujoseph/ssm_test_model/')

      @classmethod
      def tearDownClass(cls):
            os.chdir('/Users/dureaujoseph/ssm_test_model/')
            shutil.rmtree('/Users/dureaujoseph/ssm/example/noise_test')
                        
      def test_only_env_sto(self):
            os.system('./kalman -O 2  -c -x --no_dem_sto < /Users/dureaujoseph/ssm/example/noise_test/datapackage.json')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 500
            os.system('./smc sde --no_dem_sto -t -O 2 -x -J ' + str(nparts) + ' -I 1 -N 4 --dt 0.0001  < /Users/dureaujoseph/ssm/example/noise_test/datapackage.json')
            tab1 = genfromtxt('hat_1.csv',delimiter=',',names=True)
            
            meanSMC0a = tab1[1][13]
            meanEKF0a = tab0[1][13]
            q975SMC0a = tab1[1][15]
            q975EKF0a = tab0[1][15]
            meanSMC0b = tab1[1][10]
            meanEKF0b = tab0[1][10]
            q975SMC0b = tab1[1][12]
            q975EKF0b = tab0[1][12]
            meanSMC0r0 = tab1[1][19]
            meanEKF0r0 = tab0[1][19]
            q975SMC0r0 = tab1[1][21]
            q975EKF0r0 = tab0[1][21]


            # Tests on 97.5% quantiles
            # based on on CLT for empirical quantile given in "Statistics and Data Analysis for Financial Engineering, Rupert 2011"
            self.assertTrue(abs((q975SMC0r0-q975EKF0r0)/(math.sqrt(0.975*0.025)/(stats.norm.pdf(q975EKF0r0-meanEKF0r0,loc=0,scale=(q975EKF0r0-meanEKF0r0)/1.96)*math.sqrt(nparts))))<1.96)
            self.assertTrue(abs((q975SMC0a-q975EKF0a)/(math.sqrt(0.975*0.025)/(stats.norm.pdf(q975EKF0a-meanEKF0a,loc=0,scale=(q975EKF0a-meanEKF0a)/1.96)*math.sqrt(nparts))))<1.96)
            self.assertTrue(abs((q975SMC0b-q975EKF0b)/(math.sqrt(0.975*0.025)/(stats.norm.pdf(q975EKF0b-meanEKF0b,loc=0,scale=(q975EKF0b-meanEKF0b)/1.96)*math.sqrt(nparts))))<1.96)


class TestpMCMCsmoothing(unittest.TestCase):
      @classmethod
      def setUpClass(cls):
            print("")
            print("Testing smoothing density provided by pMCMC")
            os.chdir('/Users/dureaujoseph/ssm/example')
            shutil.copytree('noise','noise_test')
            f = open('/Users/dureaujoseph/ssm/example/noise_test/datapackages/model-jdureau-noise/datapackage.json')
            j = json.load(f)
            j['resources'][0]["data"] = []
            j['resources'][0]["data"].append({"name":"paris","composition":["S","I","I2","E1","E2"]})
            j['resources'][1]["data"] = []
            j['resources'][1]["data"].append({"from":"S", "to":"I", "rate":"1000/S", "white_noise": {"name":"noise_SI", "sd": "vol"}})
            j['resources'][1]["data"].append({"from":"I", "to":"I2", "rate":"1000/I", "white_noise": {"name":"noise_SI", "sd": "vol"}})
            j['resources'][1]["data"].append({"from":"I2", "to":"S", "rate":"1000/I2", "white_noise": {"name":"noise_SI2", "sd": "vol2"}})
            j['resources'][1]["data"].append({"from":"E1", "to":"E2", "rate":"3/N", "tracked": ["all_inc_out"]})

            j['resources'][2]["data"] = [j['resources'][2]["data"][0]]

            j['resources'][3]["data"] = []
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][0]['name'] = 'r0'
            j['resources'][3]['data'][0]['prior'] = {}
            j['resources'][3]['data'][0]['prior']['distribution']='fixed'
            j['resources'][3]['data'][0]['prior']['value']=100
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][1]['name'] = 'S'
            j['resources'][3]['data'][1]['prior'] = {}
            j['resources'][3]['data'][1]['prior']['distribution']='fixed'
            j['resources'][3]['data'][1]['prior']['value']=200000
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][2]['name'] = 'I'
            j['resources'][3]['data'][2]['prior'] = {}
            j['resources'][3]['data'][2]['prior']['distribution']='fixed'
            j['resources'][3]['data'][2]['prior']['value']=200000
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][3]['name'] = 'I2'
            j['resources'][3]['data'][3]['prior'] = {}
            j['resources'][3]['data'][3]['prior']['distribution']='fixed'
            j['resources'][3]['data'][3]['prior']['value']=200000
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][4]['name'] = 'E1'
            j['resources'][3]['data'][4]['prior'] = {}
            j['resources'][3]['data'][4]['prior']['distribution']='fixed'
            j['resources'][3]['data'][4]['prior']['value']=200000
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][5]['name'] = 'vol'
            j['resources'][3]['data'][5]['prior'] = {}
            j['resources'][3]['data'][5]['prior']['distribution']='fixed'
            j['resources'][3]['data'][5]['prior']['value']=0.01428571
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][6]['name'] = 'vol2'
            j['resources'][3]['data'][6]['prior'] = {}
            j['resources'][3]['data'][6]['prior']['distribution']='fixed'
            j['resources'][3]['data'][6]['prior']['value']=0.02857143
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][7]['name'] = 'rep_all_CDC_inc'
            j['resources'][3]['data'][7]['prior'] = {}
            j['resources'][3]['data'][7]['prior']['distribution']='fixed'
            j['resources'][3]['data'][7]['prior']['value']=0.2
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][8]['name'] = 'prop_all_CDC_inc'
            j['resources'][3]['data'][8]['prior'] = {}
            j['resources'][3]['data'][8]['prior']['distribution']='fixed'
            j['resources'][3]['data'][8]['prior']['value']=0.2
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][9]['name'] = 'phi'
            j['resources'][3]['data'][9]['prior'] = {}
            j['resources'][3]['data'][9]['prior']['distribution']='fixed'
            j['resources'][3]['data'][9]['prior']['value']=0.2
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][10]['name'] = 'E2'
            j['resources'][3]['data'][10]['prior'] = {}
            j['resources'][3]['data'][10]['prior']['distribution']='uniform'
            j['resources'][3]['data'][10]['prior']['lower']=199995
            j['resources'][3]['data'][10]['prior']['upper']=200010
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][11]['name'] = 'N'
            j['resources'][3]['data'][11]['prior'] = {}
            j['resources'][3]['data'][11]['prior']['distribution']='fixed'
            j['resources'][3]['data'][11]['prior']['value']=0.2


            j['resources'].append({})
            j['resources'][4]['name']='sde'
            j['resources'][4]['data']={}
            j['resources'][4]['data']['drift']=[]
            j['resources'][4]['data']['drift'].append({})
            j['resources'][4]['data']['drift'][0]['name']='r0'
            j['resources'][4]['data']['drift'][0]['f']=0.0
            j['resources'][4]['data']['dispersion']=[['vol']]
            with open('/Users/dureaujoseph/ssm/example/noise_test/datapackages/model-jdureau-noise/datapackage.json','w') as outfile:
                  json.dump(j,outfile)

            f = open('/Users/dureaujoseph/ssm/example/noise_test/datapackage.json')
            j = json.load(f)
            j['resources'][0]["data"]={}
            j['resources'][0]["data"]["E2"] = 200000
            j['resources'][1]["data"]={}
            j['resources'][1]["data"]["E2"] = {"E2":0.02}
            with open('/Users/dureaujoseph/ssm/example/noise_test/datapackage.json','w') as outfile:
                  json.dump(j,outfile)

            os.chdir('/Users/dureaujoseph/ssm/ssm/src/C/')
            os.system('make clean')
            os.system('make')
            os.system('make install')
            
            os.chdir('/Users/dureaujoseph/ssm/ssm/')
            b = Builder(os.path.join(os.getenv("HOME"), 'ssm_test_model'), os.path.join('..' ,'example', 'noise_test', 'datapackages', 'model-jdureau-noise', 'datapackage.json'))
            b.prepare()
            b.code()
            b.write_data()

            os.chdir('/Users/dureaujoseph/ssm_test_model/C/templates')
            os.system('make clean')
            os.system('make')
            os.system('make install')


            os.chdir('/Users/dureaujoseph/ssm_test_model/')
            
      def setUp(self):
      # Things that need to be done before tests            
            os.chdir('/Users/dureaujoseph/ssm_test_model/')

      @classmethod
      def tearDownClass(cls):
            shutil.rmtree('/Users/dureaujoseph/ssm/example/noise_test')

      def test_particle_genealogy(self):
            # We assess that the exploration of the genealogy to reconstruct sampled traj is correct by checking the trajectory
            # of the random walk through the distribution of its increments
            
            os.system('./kalman -O 2  -c -x --no_dem_sto < /Users/dureaujoseph/ssm/example/noise_test/datapackage.json')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 100
            nbiters = 2
            
            os.system('./pmcmc sde -t -x --no_dem_sto -J ' + str(nparts) + ' -M ' + str(nbiters) + ' -I 1 -N 4 -T '  + str(nbiters) + '  < /Users/dureaujoseph/ssm/example/noise_test/datapackage.json')          
            tab1 = genfromtxt('X_1.csv',delimiter=',',names=True)

            x = []
            for i in range(1,53):
                  # Compute and normalise increments of sampled path
                  x.append((float(tab1[52 - i][7])-float(tab1[52 - i + 1][7]))/(math.sqrt(7)*0.1/7))          

            self.assertTrue(abs((stats.scoreatpercentile(x,50)-0)/(math.sqrt(0.5*0.5)/(stats.norm.pdf(0,loc=0,scale=1)*math.sqrt(len(x)))))<1.96)

      def test_correct_timing(self):
            # We check that timing has not been broken when playing with the indexes reconstructing the sampled path in pmcmc

            os.system('./kalman -O 3  -c -x --no_dem_sto < /Users/dureaujoseph/ssm/example/noise_test/datapackage.json')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 2
            nbiters = 500
            os.system('./pmcmc sde -t -x --no_dem_sto -O 4 -J ' + str(nparts) + ' -M ' + str(nbiters) + ' -I 2 -N 4 -T '  + str(nbiters) + '  < /Users/dureaujoseph/ssm/example/noise_test/datapackage.json')          
            tab1 = genfromtxt('X_2.csv',delimiter=',',names=True)

            ### t = 1 (first obs)
            x = []
            for i in range(1,nbiters):
                  x.append(float(tab1[3*i - 2][7])) 
                  
            meanSMC0r0 = numpy.mean(x)
            meanEKF0r0 = tab0[1][19]
            q975SMC0r0 = stats.scoreatpercentile(x,97.5)
            q975EKF0r0 = tab0[1][21]

            # Tests on 97.5% quantiles
            # based on on CLT for empirical quantile given in "Statistics and Data Analysis for Financial Engineering, Rupert 2011"
            self.assertTrue(abs((q975SMC0r0-q975EKF0r0)/(math.sqrt(0.975*0.025)/(stats.norm.pdf(q975EKF0r0-meanEKF0r0,loc=0,scale=(q975EKF0r0-meanEKF0r0)/1.96)*math.sqrt(nparts))))<1.96)

            ### t = 2 (second obs)
            x = []
            for i in range(1,nbiters):
                  x.append(float(tab1[3*i - 3][7]))

            meanSMC0r0 = numpy.mean(x)
            meanEKF0r0 = tab0[2][19]
            q975SMC0r0 = stats.scoreatpercentile(x,97.5)
            q975EKF0r0 = tab0[2][21]

            # Tests on 97.5% quantiles
            # based on on CLT for empirical quantile given in "Statistics and Data Analysis for Financial Engineering, Rupert 2011"
            self.assertTrue(abs((q975SMC0r0-q975EKF0r0)/(math.sqrt(0.975*0.025)/(stats.norm.pdf(q975EKF0r0-meanEKF0r0,loc=0,scale=(q975EKF0r0-meanEKF0r0)/1.96)*math.sqrt(nparts))))<1.96)

      
            

class TestpMCMCsmoothingWithNaNs(unittest.TestCase):
      @classmethod
      def setUpClass(cls):
            print("")
            print("Testing smoothing density provided by pMCMC, in presence of NaNs")
            os.chdir('/Users/dureaujoseph/ssm/example')
            shutil.copytree('noise','noise_test')
            f = open('/Users/dureaujoseph/ssm/example/noise_test/datapackages/model-jdureau-noise/datapackage.json')
            j = json.load(f)
            j['resources'][0]["data"] = []
            j['resources'][0]["data"].append({"name":"paris","composition":["S","I","I2","E1","E2"]})
            j['resources'][1]["data"] = []
            j['resources'][1]["data"].append({"from":"S", "to":"I", "rate":"1000/S", "white_noise": {"name":"noise_SI", "sd": "vol"}})
            j['resources'][1]["data"].append({"from":"I", "to":"I2", "rate":"1000/I", "white_noise": {"name":"noise_SI", "sd": "vol"}})
            j['resources'][1]["data"].append({"from":"I2", "to":"S", "rate":"1000/I2", "white_noise": {"name":"noise_SI2", "sd": "vol2"}})
            j['resources'][1]["data"].append({"from":"E1", "to":"E2", "rate":"3/N", "tracked": ["all_inc_out"]})

            j['resources'][2]["data"] = [j['resources'][2]["data"][0]]

            j['resources'][3]["data"] = []
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][0]['name'] = 'r0'
            j['resources'][3]['data'][0]['prior'] = {}
            j['resources'][3]['data'][0]['prior']['distribution']='fixed'
            j['resources'][3]['data'][0]['prior']['value']=100
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][1]['name'] = 'S'
            j['resources'][3]['data'][1]['prior'] = {}
            j['resources'][3]['data'][1]['prior']['distribution']='fixed'
            j['resources'][3]['data'][1]['prior']['value']=200000
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][2]['name'] = 'I'
            j['resources'][3]['data'][2]['prior'] = {}
            j['resources'][3]['data'][2]['prior']['distribution']='fixed'
            j['resources'][3]['data'][2]['prior']['value']=200000
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][3]['name'] = 'I2'
            j['resources'][3]['data'][3]['prior'] = {}
            j['resources'][3]['data'][3]['prior']['distribution']='fixed'
            j['resources'][3]['data'][3]['prior']['value']=200000
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][4]['name'] = 'E1'
            j['resources'][3]['data'][4]['prior'] = {}
            j['resources'][3]['data'][4]['prior']['distribution']='fixed'
            j['resources'][3]['data'][4]['prior']['value']=200000
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][5]['name'] = 'vol'
            j['resources'][3]['data'][5]['prior'] = {}
            j['resources'][3]['data'][5]['prior']['distribution']='fixed'
            j['resources'][3]['data'][5]['prior']['value']=0.01428571
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][6]['name'] = 'vol2'
            j['resources'][3]['data'][6]['prior'] = {}
            j['resources'][3]['data'][6]['prior']['distribution']='fixed'
            j['resources'][3]['data'][6]['prior']['value']=0.02857143
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][7]['name'] = 'rep_all_CDC_inc'
            j['resources'][3]['data'][7]['prior'] = {}
            j['resources'][3]['data'][7]['prior']['distribution']='fixed'
            j['resources'][3]['data'][7]['prior']['value']=0.2
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][8]['name'] = 'prop_all_CDC_inc'
            j['resources'][3]['data'][8]['prior'] = {}
            j['resources'][3]['data'][8]['prior']['distribution']='fixed'
            j['resources'][3]['data'][8]['prior']['value']=0.2
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][9]['name'] = 'phi'
            j['resources'][3]['data'][9]['prior'] = {}
            j['resources'][3]['data'][9]['prior']['distribution']='fixed'
            j['resources'][3]['data'][9]['prior']['value']=0.2
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][10]['name'] = 'E2'
            j['resources'][3]['data'][10]['prior'] = {}
            j['resources'][3]['data'][10]['prior']['distribution']='uniform'
            j['resources'][3]['data'][10]['prior']['lower']=199995
            j['resources'][3]['data'][10]['prior']['upper']=200010
            j['resources'][3]['data'].append({})
            j['resources'][3]['data'][11]['name'] = 'N'
            j['resources'][3]['data'][11]['prior'] = {}
            j['resources'][3]['data'][11]['prior']['distribution']='fixed'
            j['resources'][3]['data'][11]['prior']['value']=0.2


            j['resources'].append({})
            j['resources'][4]['name']='sde'
            j['resources'][4]['data']={}
            j['resources'][4]['data']['drift']=[]
            j['resources'][4]['data']['drift'].append({})
            j['resources'][4]['data']['drift'][0]['name']='r0'
            j['resources'][4]['data']['drift'][0]['f']=0.0
            j['resources'][4]['data']['dispersion']=[['vol']]
            with open('/Users/dureaujoseph/ssm/example/noise_test/datapackages/model-jdureau-noise/datapackage.json','w') as outfile:
                  json.dump(j,outfile)

            f = open('/Users/dureaujoseph/ssm/example/noise_test/datapackage.json')
            j = json.load(f)
            j['resources'][0]["data"]={}
            j['resources'][0]["data"]["E2"] = 200000
            j['resources'][1]["data"]={}
            j['resources'][1]["data"]["E2"] = {"E2":0.02}
            with open('/Users/dureaujoseph/ssm/example/noise_test/datapackage.json','w') as outfile:
                  json.dump(j,outfile)

            os.chdir('/Users/dureaujoseph/ssm/example/noise_test/datapackages/model-jdureau-noise/datapackages/data-jdureau-test/data')
            shutil.copy('data.csv','data_temp.csv')

            with open('data_temp.csv','rU') as csvfile:
                  with open('data.csv','wb') as csvfile2:
                        file = csv.reader(csvfile)
                        file2 = csv.writer(csvfile2)

                        ind = 0
                        for row in file:
                              if ind>0:
                                    i = 0
                                    for x in row:
                                          if i > 0:
                                                csvfile2.write(',')
                                          if ind > 2 & i > 0:
                                                csvfile2.write('null')
                                          else :
                                                csvfile2.write(str(x))
                                          i = i+1
                                    csvfile2.write('\n')
                              else:
                                    i = 0
                                    for x in row:
                                          if i > 0:
                                                csvfile2.write(',')
                                          csvfile2.write(x)
                                          i = i+1
                                    csvfile2.write('\n')
                              ind = ind + 1

                  

            os.chdir('/Users/dureaujoseph/ssm/ssm/src/C/')
            os.system('make clean')
            os.system('make')
            os.system('make install')
            
            os.chdir('/Users/dureaujoseph/ssm/ssm/')
            b = Builder(os.path.join(os.getenv("HOME"), 'ssm_test_model'), os.path.join('..' ,'example', 'noise_test', 'datapackages', 'model-jdureau-noise', 'datapackage.json'))
            b.prepare()
            b.code()
            b.write_data()

            os.chdir('/Users/dureaujoseph/ssm_test_model/C/templates')
            os.system('make clean')
            os.system('make')
            os.system('make install')


            os.chdir('/Users/dureaujoseph/ssm_test_model/')


      def setUp(self):
            # Things that need to be done before tests            
            os.chdir('/Users/dureaujoseph/ssm_test_model/')

      @classmethod
      def tearDownClass(cls):
            shutil.rmtree('/Users/dureaujoseph/ssm/example/noise_test')

      def test_particle_genealogy(self):
            # We assess that the exploration of the genealogy to reconstruct sampled traj is correct by checking the trajectory
            # of the random walk through the distribution of its increments
            
            os.system('./kalman -O 2  -c -x --no_dem_sto < /Users/dureaujoseph/ssm/example/noise_test/datapackage.json')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 100
            nbiters = 2
            
            os.system('./pmcmc sde -t -x --no_dem_sto -J ' + str(nparts) + ' -M ' + str(nbiters) + ' -I 1 -N 4 -T '  + str(nbiters) + '  < /Users/dureaujoseph/ssm/example/noise_test/datapackage.json')          
            tab1 = genfromtxt('X_1.csv',delimiter=',',names=True)

            x = []
            for i in range(1,53):
                  # Compute and normalise increments of sampled path
                  x.append((float(tab1[52 - i][7])-float(tab1[52 - i + 1][7]))/(math.sqrt(7)*0.1/7))          

            self.assertTrue(abs((stats.scoreatpercentile(x,50)-0)/(math.sqrt(0.5*0.5)/(stats.norm.pdf(0,loc=0,scale=1)*math.sqrt(len(x)))))<1.96)

      def test_correct_timing(self):
            # We check that timing has not been broken when playing with the indexes reconstructing the sampled path in pmcmc

            os.system('./kalman -O 3  -c -x --no_dem_sto < /Users/dureaujoseph/ssm/example/noise_test/datapackage.json')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 2
            nbiters = 500
            os.system('./pmcmc sde -t -x --no_dem_sto -O 4 -J ' + str(nparts) + ' -M ' + str(nbiters) + ' -I 2 -N 4 -T '  + str(nbiters) + '  < /Users/dureaujoseph/ssm/example/noise_test/datapackage.json')          
            tab1 = genfromtxt('X_2.csv',delimiter=',',names=True)

            ### t = 1 (first obs)
            x = []
            for i in range(1,nbiters):
                  x.append(float(tab1[3*i - 2][7])) 
                  
            meanSMC0r0 = numpy.mean(x)
            meanEKF0r0 = tab0[1][19]
            q975SMC0r0 = stats.scoreatpercentile(x,97.5)
            q975EKF0r0 = tab0[1][21]

            # Tests on 97.5% quantiles
            # based on on CLT for empirical quantile given in "Statistics and Data Analysis for Financial Engineering, Rupert 2011"
            self.assertTrue(abs((q975SMC0r0-q975EKF0r0)/(math.sqrt(0.975*0.025)/(stats.norm.pdf(q975EKF0r0-meanEKF0r0,loc=0,scale=(q975EKF0r0-meanEKF0r0)/1.96)*math.sqrt(nparts))))<1.96)

            ### t = 2 (second obs)
            x = []
            for i in range(1,nbiters):
                  x.append(float(tab1[3*i - 3][7]))

            meanSMC0r0 = numpy.mean(x)
            meanEKF0r0 = tab0[2][19]
            q975SMC0r0 = stats.scoreatpercentile(x,97.5)
            q975EKF0r0 = tab0[2][21]

            # Tests on 97.5% quantiles
            # based on on CLT for empirical quantile given in "Statistics and Data Analysis for Financial Engineering, Rupert 2011"
            self.assertTrue(abs((q975SMC0r0-q975EKF0r0)/(math.sqrt(0.975*0.025)/(stats.norm.pdf(q975EKF0r0-meanEKF0r0,loc=0,scale=(q975EKF0r0-meanEKF0r0)/1.96)*math.sqrt(nparts))))<1.96)

      


def suite_TestNoiseResults():
      suite = unittest.TestSuite()
      suite.addTest(unittest.makeSuite(TestNoiseResults))
      return suite
            
def suite_TestTransfsAndPMCMC():
      suite = unittest.TestSuite()
      suite.addTest(unittest.makeSuite(TestTransfsAndPMCMC))
      return suite

def suite_TestKalmanOnDiffusions():
      suite = unittest.TestSuite()
      suite.addTest(unittest.makeSuite(TestKalmanOnDiffusions))
      return suite

def suite_SMCSDEagainstKalman():
      suite = unittest.TestSuite()
      suite.addTest(unittest.makeSuite(TestSMCSDEagainstKalman))
      return suite

def suite_pMCMCsmoothing():
      suite = unittest.TestSuite()
      suite.addTest(unittest.makeSuite(TestpMCMCsmoothing))
      return suite

if __name__ == '__main__' :

      run_NoiseResults = 1
      run_TransfsAndPMCMC = 0
      run_KalmanOnDiffusions = 0
      run_SMCSDEagainstKalman = 0
      run_pMCMCsmoothing = 0
      run_pMCMCsmoothingWithNaNs = 0

      Root = os.getcwd()

      suite = unittest.TestSuite()
      if run_NoiseResults:
            suite.addTest(unittest.makeSuite(TestNoiseResults))

      if run_TransfsAndPMCMC:
            suite.addTest(unittest.makeSuite(TestTransfsAndPMCMC))

      if run_KalmanOnDiffusions:
            suite.addTest(unittest.makeSuite(TestKalmanOnDiffusions))

      if run_SMCSDEagainstKalman:
            suite.addTest(unittest.makeSuite(TestSMCSDEagainstKalman))

      if run_pMCMCsmoothing:
            suite.addTest(unittest.makeSuite(TestpMCMCsmoothing))
            
      if run_pMCMCsmoothingWithNaNs:
            suite.addTest(unittest.makeSuite(TestpMCMCsmoothingWithNaNs))

      unittest.TextTestRunner().run(suite)
