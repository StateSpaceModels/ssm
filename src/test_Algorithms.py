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

            os.chdir(Root)
            os.system(Root + '/../bin/ssm install  ' + Root + '/../examples/noise/package.json')
            
      def setUp(self):
      # Things that need to be done before tests
            os.chdir(Root + '/ssm_model')

      @classmethod
      def tearDownClass(cls):
            #os.chdir(Root + '/ssm_model')
            shutil.rmtree(Root + '/ssm_model')

      def test_ode(self):
            os.system('./simplex -M 10000 < ../../examples/noise/package.json | ./simplex -M 10000 | ./simplex -M 100000 > theta.json')
            os.system('./smc ode -J 1  --trace < theta.json')
            tab = genfromtxt('trace_0.csv',delimiter=',',names=True).tolist()
            self.assertAlmostEqual(tab[5],-824.653)
            os.system('./smc ode -J 10  --trace < theta.json')
            tab = genfromtxt('trace_0.csv',delimiter=',',names=True).tolist()
            self.assertAlmostEqual(tab[5],-824.653)

      def test_kalman_map(self):
            os.system('./ksimplex -M 1000 --trace < ' + Root + '/../examples/noise/package.json')
            tab = genfromtxt('trace_0.csv',delimiter=',',names=True)
            self.assertAlmostEqual(tab[214][5],-464.806)

class TestTransfsAndPMCMC(unittest.TestCase):
      @classmethod
      def setUpClass(cls):
            print("")
            print("Testing transformations and pMCMC by sampling from priors")
            
      def setUp(self):
      # Things that need to be done before tests
            os.chdir(Root)

      @classmethod
      def tearDownClass(cls):
            os.chdir(Root)


      def test_prior_unif_transf_log(self):
            os.chdir(Root + '/../examples')
            shutil.copytree('noise','noise_test')
            f = open(Root + '/../examples/noise_test/package.json')
            j = json.load(f)
            j['resources'][8]['data']['r0_paris'] = 10
            with open(Root + '/../examples/noise_test/package.json','w') as outfile:
                  json.dump(j,outfile)
            f = open(Root + '/../examples/noise_test/package.json')
            j = json.load(f)
            j['resources'][2]['data']['upper']=10.5
            j['resources'][2]['data']['lower']=9.5
            with open(Root + '/../examples/noise_test/package.json','w') as outfile:
                  json.dump(j,outfile)

            os.chdir(Root)
            os.system(Root + '/../bin/ssm install  ' + Root + '/../examples/noise_test/package.json')

            os.chdir(Root + '/ssm_model')            
            os.system('./pmcmc ode  -C 5000 -W 5000 -O 0 -M 20000 -c -a < ' + Root + '/../examples/noise_test/package.json')
            
            test = self.call_test_unif('r0_paris')
            shutil.rmtree(Root + '/ssm_model')
            shutil.rmtree(Root + '/../examples/noise_test')

            self.assertEqual(test,'1')



      def test_prior_normal_transf_log(self):
            os.chdir(Root + '/../examples')
            shutil.copytree('noise','noise_test')
            f = open(Root + '/../examples/noise_test/package.json')
            j = json.load(f)
            j['resources'][8]['data']['r0_paris'] = 10
            with open(Root + '/../examples/noise_test/package.json','w') as outfile:
                  json.dump(j,outfile)
            f = open(Root + '/../examples/noise_test/package.json')
            j = json.load(f)
            j['resources'][2]['data']['distribution']='normal'
            j['resources'][2]['data']['lower']=0
            j['resources'][2]['data']['mean']=10
            j['resources'][2]['data']['sd']=1
            del j['resources'][2]['data']['upper']
            del j['resources'][2]['data']['lower']
            with open(Root + '/../examples/noise_test/package.json','w') as outfile:
                  json.dump(j,outfile)

            os.chdir(Root)
            os.system(Root + '/../bin/ssm install  ' + Root + '/../examples/noise_test/package.json')

            os.chdir(Root + '/ssm_model')            
            os.system('./pmcmc ode  -C 5000 -W 5000 -O 0 -M 20000 -c -a < '+ Root + '/../examples/noise_test/package.json')

            test = self.call_test_normal('r0_paris')
            shutil.rmtree(Root + '/ssm_model')
            shutil.rmtree(Root + '/../examples/noise_test')

            self.assertEqual(test,'1')





      def test_prior_normal_and_unif_transf_log(self):
            os.chdir(Root + '/../examples')
            shutil.copytree('noise','noise_test')
            f = open(Root + '/../examples/noise_test/package.json')
            j = json.load(f)
            j['resources'][8]['data']['r0_paris'] = 10
            j['resources'][8]['data']['r0_nyc'] = 10
            with open(Root + '/../examples/noise_test/package.json','w') as outfile:
                  json.dump(j,outfile)
            f = open(Root + '/../examples/noise_test/package.json')
            j = json.load(f)
            j['resources'][2]['data']['distribution']='normal'
            j['resources'][2]['data']['lower']=0
            j['resources'][2]['data']['mean']=10
            j['resources'][2]['data']['sd']=1
            del j['resources'][2]['data']['upper']
            del j['resources'][2]['data']['lower']
            j['resources'][3]['data']['lower']=9.5
            j['resources'][3]['data']['lower']=10.5
            with open(Root + '/../examples/noise_test/package.json','w') as outfile:
                  json.dump(j,outfile)            

            os.chdir(Root)
            os.system(Root + '/../bin/ssm install  ' + Root + '/../examples/noise_test/package.json')

            os.chdir(Root + '/ssm_model')            
            
            os.system('./pmcmc ode  -C 1000 -W 100000 -O 0 -M 100000 -c -a < ' + Root + '/../examples/noise_test/package.json')

            test1 = self.call_test_normal('r0_paris')
            test2 = self.call_test_unif('r0_nyc')
            shutil.rmtree(Root + '/ssm_model')
            shutil.rmtree(Root + '/../examples/noise_test')

            self.assertEqual(int(test1)*int(test2),1)



      def test_prior_normal_transf_logit_ab(self):
            os.chdir(Root + '/../examples')
            shutil.copytree('noise','noise_test')
            f = open(Root + '/../examples/noise_test/package.json')
            j = json.load(f)
            j['resources'][8]['data']['r0_paris'] = 10
            with open(Root + '/../examples/noise_test/package.json','w') as outfile:
                  json.dump(j,outfile)
            f = open(Root + '/../examples/noise_test/package.json')
            j = json.load(f)
            j['resources'][2]['data']['distribution']='normal'
            j['resources'][2]['data']['lower']=0
            j['resources'][2]['data']['upper']=20
            j['resources'][2]['data']['mean']=10
            j['resources'][2]['data']['sd']=1
            with open(Root + '/../examples/noise_test/package.json','w') as outfile:
                  json.dump(j,outfile)   

            os.chdir(Root)
            os.system(Root + '/../bin/ssm install  ' + Root + '/../examples/noise_test/package.json')

            os.chdir(Root + '/ssm_model')            
            
            os.system('./pmcmc ode  -C 5000 -W 5000 -O 0 -M 20000 -c -a < ' + Root + '/../examples/noise_test/package.json')

            test = self.call_test_normal('r0_paris')
            shutil.rmtree(Root + '/ssm_model')
            shutil.rmtree(Root + '/../examples/noise_test')

            self.assertEqual(test,'1') 


      def test_prior_normal_transf_identity(self):
            os.chdir(Root + '/../examples')
            shutil.copytree('noise','noise_test')
            f = open(Root + '/../examples/noise_test/package.json')
            j = json.load(f)
            j['resources'][8]['data']['r0_paris'] = 10
            with open(Root + '/../examples/noise_test/package.json','w') as outfile:
                  json.dump(j,outfile)
            f = open(Root + '/../examples/noise_test/package.json')
            j = json.load(f)
            j['resources'][2]['data']['distribution']='normal'
            j['resources'][2]['data']['mean']=10
            j['resources'][2]['data']['sd']=1
            del j['resources'][2]['data']['upper']
            del j['resources'][2]['data']['lower']

            with open(Root + '/../examples/noise_test/package.json','w') as outfile:
                  json.dump(j,outfile)

            os.chdir(Root)
            os.system(Root + '/../bin/ssm install  ' + Root + '/../examples/noise_test/package.json')
            os.chdir(Root + '/ssm_model')
            
            os.system('./pmcmc ode  -C 5000 -W 5000 -O 0 -M 20000 -c -a < ' + Root + '/../examples/noise_test/package.json')

            test = self.call_test_normal('r0_paris')
            shutil.rmtree(Root + '/ssm_model')
            shutil.rmtree(Root + '/../examples/noise_test')

            self.assertEqual(test,'1') 

      def call_test_unif(self, varname):
            shutil.copyfile(Root+'/TestsR/test_unif.r',Root + '/ssm_model/test_unif.r')
            os.chdir(Root + '/ssm_model')
            os.system('R --vanilla < test_unif.r ' + varname  + '> /dev/null 2>&1')
            f = open("outfile.txt","r")
            x = f.readlines()
            return x[0]

      def call_test_normal(self,varname):
            shutil.copyfile(Root+'/TestsR/test_normal.r',Root + '/ssm_model/test_normal.r')
            os.chdir(Root + '/ssm_model')
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
            
            os.chdir(Root + '/../examples')
            shutil.copytree('noise','noise_test')
            f = open(Root + '/../examples/noise_test/package.json')
            j = json.load(f)
            j['model']['inputs'].insert(0,{})
            j['model']['inputs'].insert(0,{})
            j['model']['inputs'][0]['name'] = 'test_par'
            j['model']['inputs'][0]['description'] = ''
            j['model']['inputs'][0]['data'] = { "resource": "test_par"}
            j['model']['inputs'][1]['name'] = 'test_vol'
            j['model']['inputs'][1]['description'] = ''
            j['model']['inputs'][1]['data'] = { "resource": "test_vol"}
            j['resources'].insert(0,{}) 
            j['resources'].insert(0,{})
            j['resources'][0]['name'] = 'test_par'
            j['resources'][0]['description'] = ''
            j['resources'][0]['data'] = { "distribution":"fixed", "value" : 0}
            j['resources'][1]['name'] = 'test_vol'
            j['resources'][1]['description'] = ''
            j['resources'][1]['data'] = { "distribution":"fixed", "value" : 0.1428571}            
            j['model']['sde'] = {}
            j['model']['sde']['drift']=[]
            j['model']['sde']['drift'].append({})
            j['model']['sde']['drift'][0]['name']='test_par'
            j['model']['sde']['drift'][0]['f']=0.0
            j['model']['sde']['dispersion']=[['test_vol']]

            with open(Root + '/../examples/noise_test/package.json','w') as outfile:
                  json.dump(j,outfile)

            os.chdir(Root)
            os.system(Root + '/../bin/ssm install  ' + Root + '/../examples/noise_test/package.json')
            
            os.chdir(Root + '/ssm_model')
            
      
      def setUp(self):
      # Things that need to be done before tests
            os.chdir(Root + '/ssm_model')

      @classmethod
      def tearDownClass(cls):
            os.chdir(Root + '/ssm_model/')
            shutil.rmtree(Root + '/../examples/noise_test')
            
      def test_1step(self):

            os.system('./kalman -O 2  -c -x < ' + Root + '/../examples/noise_test/package.json')
            tab = genfromtxt('hat_0.csv',delimiter=',',names=True)
            self.assertAlmostEqual(tab['lower_test_par'][0],-1.96/math.sqrt(7),5)

      def test_10step(self):
            os.system('./kalman -O 10  -c -x < ' + Root + '/../examples/noise_test/package.json')
            tab = genfromtxt('hat_0.csv',delimiter=',',names=True)
            self.assertAlmostEqual(tab['lower_test_par'][9]/math.sqrt(10),-1.96/math.sqrt(7),5)

class TestSMCSDEagainstKalman(unittest.TestCase):
      @classmethod
      def setUpClass(cls):
            print("")
            print("Testing SMC against EKF on linear examples")

            os.chdir(Root + '/../examples')
            shutil.copytree('noise','noise_test')
            f = open(Root + '/../examples/noise_test/package.json')
            j = json.load(f)
            j["model"]["populations"] = []
            j["model"]["populations"].append({"name":"paris","composition":["S","I","I2","E1","E2"]})
            j["model"]["reactions"] = []
            j["model"]["reactions"].append({"from":"S", "to":"I", "rate":"1000/S", "white_noise": {"name":"noise_SI", "sd": "vol"}})
            j["model"]["reactions"].append({"from":"I", "to":"I2", "rate":"1000/I", "white_noise": {"name":"noise_SI", "sd": "vol"}})
            j["model"]["reactions"].append({"from":"I2", "to":"S", "rate":"1000/I2", "white_noise": {"name":"noise_SI2", "sd": "vol2"}})
            j["model"]["reactions"].append({"from":"E1", "to":"E2", "rate":"3/N", "tracked": ["all_inc_out"]})

            j["model"]["observations"] = [j["model"]["observations"][0]]

            j["model"]["inputs"] = []
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'r0'
            j["model"]["inputs"][0]['data'] = {'resource':'r0'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'S'
            j["model"]["inputs"][0]['data'] = {'resource':'S'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'I'
            j["model"]["inputs"][0]['data'] = {'resource':'I'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'I2'
            j["model"]["inputs"][0]['data'] = {'resource':'I2'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'E1'
            j["model"]["inputs"][0]['data'] = {'resource':'E1'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'vol'
            j["model"]["inputs"][0]['data'] = {'resource':'vol'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'vol2'
            j["model"]["inputs"][0]['data'] = {'resource':'vol2'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'rep_all_CDC_inc'
            j["model"]["inputs"][0]['data'] = {'resource':'rep_all_CDC_inc'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'prop_all_CDC_inc'
            j["model"]["inputs"][0]['data'] = {'resource':'prop_all_CDC_inc'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'phi'
            j["model"]["inputs"][0]['data'] = {'resource':'phi'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'E2'
            j["model"]["inputs"][0]['data'] = {'resource':'E2'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'N'
            j["model"]["inputs"][0]['data'] = {'resource':'N'}

            j["resources"]=[]
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'r0'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'uniform'
            j["resources"][0]['data']['lower'] = 50
            j["resources"][0]['data']['upper'] = 150
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'S'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 200000
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'I'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 200000
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'I2'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 200000
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'E1'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 200000
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'vol'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 0.01428571
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'vol2'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 0.02857143
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'rep_all_CDC_inc'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 0.2
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'prop_all_CDC_inc'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 0.2
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'phi'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 0.2
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'E2'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 200000
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'N'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 0.2

            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'values'
            j["resources"][0]['data'] = { "r0": 100 }

            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'covariance'
            j["resources"][0]['data'] = { "r0": { "r0": 0.02 }}

            j["model"]["sde"] = {}
            j["model"]["sde"]["drift"] = []
            j["model"]["sde"]["drift"].append({ "name": "r0", "f": 0.0})
            j["model"]["sde"]["dispersion"] = [['vol']]
            
            with open(Root + '/../examples/noise_test/package.json','w') as outfile:
                  json.dump(j,outfile)

            os.chdir(Root)
            os.system(Root + '/../bin/ssm install  ' + Root + '/../examples/noise_test/package.json')

            os.chdir(Root + '/ssm_model')

      def setUp(self):
      # Things that need to be done before tests
            os.chdir(Root + '/ssm_model')

      @classmethod
      def tearDownClass(cls):
            os.chdir(Root + '/ssm_model/')
            shutil.rmtree(Root + '/../examples/noise_test')
                        
      def test_only_env_sto(self):
            os.system('./kalman -O 2  -c -x --no_dem_sto < ' + Root + '/../examples/noise_test/package.json')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 500
            os.system('./smc sde --no_dem_sto -t -O 2 -x -J ' + str(nparts) + ' -I 1 -N 4 --dt 0.0001  < ' + Root + '/../examples/noise_test/package.json')
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
            os.chdir(Root + '/../examples')
            shutil.copytree('noise','noise_test')
            f = open(Root + '/../examples/noise_test/package.json')
            j = json.load(f)
            j["model"]["populations"] = []
            j["model"]["populations"].append({"name":"paris","composition":["S","I","I2","E1","E2"]})
            j["model"]["reactions"] = []
            j["model"]["reactions"].append({"from":"S", "to":"I", "rate":"1000/S", "white_noise": {"name":"noise_SI", "sd": "vol"}})
            j["model"]["reactions"].append({"from":"I", "to":"I2", "rate":"1000/I", "white_noise": {"name":"noise_SI", "sd": "vol"}})
            j["model"]["reactions"].append({"from":"I2", "to":"S", "rate":"1000/I2", "white_noise": {"name":"noise_SI2", "sd": "vol2"}})
            j["model"]["reactions"].append({"from":"E1", "to":"E2", "rate":"3/N", "tracked": ["all_inc_out"]})

            j["model"]["observations"] = [j["model"]["observations"][0]]

            j["model"]["inputs"] = []
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'r0'
            j["model"]["inputs"][0]['data'] = {'resource':'r0'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'S'
            j["model"]["inputs"][0]['data'] = {'resource':'S'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'I'
            j["model"]["inputs"][0]['data'] = {'resource':'I'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'I2'
            j["model"]["inputs"][0]['data'] = {'resource':'I2'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'E1'
            j["model"]["inputs"][0]['data'] = {'resource':'E1'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'vol'
            j["model"]["inputs"][0]['data'] = {'resource':'vol'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'vol2'
            j["model"]["inputs"][0]['data'] = {'resource':'vol2'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'rep_all_CDC_inc'
            j["model"]["inputs"][0]['data'] = {'resource':'rep_all_CDC_inc'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'prop_all_CDC_inc'
            j["model"]["inputs"][0]['data'] = {'resource':'prop_all_CDC_inc'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'phi'
            j["model"]["inputs"][0]['data'] = {'resource':'phi'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'E2'
            j["model"]["inputs"][0]['data'] = {'resource':'E2'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'N'
            j["model"]["inputs"][0]['data'] = {'resource':'N'}

            j["resources"]=[]
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'r0'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 100
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'S'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 200000
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'I'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 200000
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'I2'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 200000
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'E1'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 200000
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'vol'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 0.01428571
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'vol2'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 0.02857143
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'rep_all_CDC_inc'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 0.2
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'prop_all_CDC_inc'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 0.2
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'phi'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 0.2
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'E2'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'uniform'
            j["resources"][0]['data']['lower'] = 200010
            j["resources"][0]['data']['upper'] = 199995
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'N'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 0.2

            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'values'
            j["resources"][0]['data'] = { "E2": 200000 }

            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'covariance'
            j["resources"][0]['data'] = { "E2": { "E2": 0.02 }}

            j["model"]["sde"] = {}
            j["model"]["sde"]["drift"] = []
            j["model"]["sde"]["drift"].append({ "name": "r0", "f": 0.0})
            j["model"]["sde"]["dispersion"] = [['vol']]

            with open(Root + '/../examples/noise_test/package.json','w') as outfile:
                  json.dump(j,outfile)


            os.chdir(Root)
            os.system(Root + '/../bin/ssm install  ' + Root + '/../examples/noise_test/package.json')

            os.chdir(Root + '/ssm_model')
            
      def setUp(self):
      # Things that need to be done before tests            
            os.chdir(Root + '/ssm_model')

      @classmethod
      def tearDownClass(cls):
            shutil.rmtree(Root + '/../examples/noise_test')

      def test_particle_genealogy(self):
            # We assess that the exploration of the genealogy to reconstruct sampled traj is correct by checking the trajectory
            # of the random walk through the distribution of its increments
            
            os.system('./kalman -O 2  -c -x --no_dem_sto < ' + Root + '/../examples/noise_test/package.json')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 100
            nbiters = 2
            
            os.system('./pmcmc sde -t -x --no_dem_sto -J ' + str(nparts) + ' -M ' + str(nbiters) + ' -I 1 -N 4 -T '  + str(nbiters) + '  < ' + Root + '/../examples/noise_test/package.json')          
            tab1 = genfromtxt('X_1.csv',delimiter=',',names=True)

            x = []
            for i in range(1,53):
                  # Compute and normalise increments of sampled path
                  x.append((float(tab1[52 - i][7])-float(tab1[52 - i + 1][7]))/(math.sqrt(7)*0.1/7))          

            self.assertTrue(abs((stats.scoreatpercentile(x,50)-0)/(math.sqrt(0.5*0.5)/(stats.norm.pdf(0,loc=0,scale=1)*math.sqrt(len(x)))))<1.96)

      def test_correct_timing(self):
            # We check that timing has not been broken when playing with the indexes reconstructing the sampled path in pmcmc

            os.system('./kalman -O 3  -c -x --no_dem_sto < ' + Root + '/../examples/noise_test/package.json')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 2
            nbiters = 500
            os.system('./pmcmc sde -t -x --no_dem_sto -O 4 -J ' + str(nparts) + ' -M ' + str(nbiters) + ' -I 2 -N 4 -T '  + str(nbiters) + '  < ' + Root + '/../examples/noise_test/package.json')          
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
            os.chdir(Root + '/../examples')
            shutil.copytree('noise','noise_test')
            f = open(Root + '/../examples/noise_test/package.json')
            j = json.load(f)
            j["model"]["populations"] = []
            j["model"]["populations"].append({"name":"paris","composition":["S","I","I2","E1","E2"]})
            j["model"]["reactions"] = []
            j["model"]["reactions"].append({"from":"S", "to":"I", "rate":"1000/S", "white_noise": {"name":"noise_SI", "sd": "vol"}})
            j["model"]["reactions"].append({"from":"I", "to":"I2", "rate":"1000/I", "white_noise": {"name":"noise_SI", "sd": "vol"}})
            j["model"]["reactions"].append({"from":"I2", "to":"S", "rate":"1000/I2", "white_noise": {"name":"noise_SI2", "sd": "vol2"}})
            j["model"]["reactions"].append({"from":"E1", "to":"E2", "rate":"3/N", "tracked": ["all_inc_out"]})

            j["model"]["observations"] = [j["model"]["observations"][0]]

            j["model"]["inputs"] = []
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'r0'
            j["model"]["inputs"][0]['data'] = {'resource':'r0'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'S'
            j["model"]["inputs"][0]['data'] = {'resource':'S'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'I'
            j["model"]["inputs"][0]['data'] = {'resource':'I'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'I2'
            j["model"]["inputs"][0]['data'] = {'resource':'I2'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'E1'
            j["model"]["inputs"][0]['data'] = {'resource':'E1'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'vol'
            j["model"]["inputs"][0]['data'] = {'resource':'vol'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'vol2'
            j["model"]["inputs"][0]['data'] = {'resource':'vol2'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'rep_all_CDC_inc'
            j["model"]["inputs"][0]['data'] = {'resource':'rep_all_CDC_inc'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'prop_all_CDC_inc'
            j["model"]["inputs"][0]['data'] = {'resource':'prop_all_CDC_inc'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'phi'
            j["model"]["inputs"][0]['data'] = {'resource':'phi'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'E2'
            j["model"]["inputs"][0]['data'] = {'resource':'E2'}
            j["model"]["inputs"].insert(0,{})
            j["model"]["inputs"][0]['name'] = 'N'
            j["model"]["inputs"][0]['data'] = {'resource':'N'}

            j["resources"]=[]
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'r0'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 100
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'S'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 200000
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'I'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 200000
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'I2'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 200000
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'E1'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 200000
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'vol'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 0.01428571
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'vol2'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 0.02857143
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'rep_all_CDC_inc'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 0.2
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'prop_all_CDC_inc'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 0.2
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'phi'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 0.2
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'E2'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'uniform'
            j["resources"][0]['data']['lower'] = 200010
            j["resources"][0]['data']['upper'] = 199995
            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'N'
            j["resources"][0]['data'] = {}
            j["resources"][0]['data']['distribution'] = 'fixed'
            j["resources"][0]['data']['value'] = 0.2

            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'values'
            j["resources"][0]['data'] = { "E2": 200000 }

            j["resources"].insert(0,{})
            j["resources"][0]['name'] = 'covariance'
            j["resources"][0]['data'] = { "E2": { "E2": 0.02 }}

            j["model"]["sde"] = {}
            j["model"]["sde"]["drift"] = []
            j["model"]["sde"]["drift"].append({ "name": "r0", "f": 0.0})
            j["model"]["sde"]["dispersion"] = [['vol']]
            
            with open(Root + '/../examples/noise_test/package.json','w') as outfile:
                  json.dump(j,outfile)


            os.chdir(Root + '/../examples/noise_test/node_modules/noise-data/data')
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

            os.chdir(Root)
            os.system(Root + '/../bin/ssm install  ' + Root + '/../examples/noise_test/package.json')

            os.chdir(Root + '/ssm_model')


      def setUp(self):
            # Things that need to be done before tests            
            os.chdir(Root + '/ssm_model')

      @classmethod
      def tearDownClass(cls):
            os.chdir(Root + '/ssm_model')
            #            shutil.rmtree(Root + '/../examples/noise_test')

      def test_particle_genealogy(self):
            # We assess that the exploration of the genealogy to reconstruct sampled traj is correct by checking the trajectory
            # of the random walk through the distribution of its increments
            
            os.system('./kalman -O 2  -c -x --no_dem_sto < ' + Root + '/../examples/noise_test/package.json')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 100
            nbiters = 2
            
            os.system('./pmcmc sde -t -x --no_dem_sto -J ' + str(nparts) + ' -M ' + str(nbiters) + ' -I 2 -N 4 -T '  + str(nbiters) + '  < ' + Root + '/../examples/noise_test/package.json')          
            tab1 = genfromtxt('X_2.csv',delimiter=',',names=True)

            x = []
            for i in range(1,53):
                  # Compute and normalise increments of sampled path
                  x.append((float(tab1[52 - i][7])-float(tab1[52 - i + 1][7]))/(math.sqrt(7)*0.1/7))          

            self.assertTrue(abs((stats.scoreatpercentile(x,50)-0)/(math.sqrt(0.5*0.5)/(stats.norm.pdf(0,loc=0,scale=1)*math.sqrt(len(x)))))<1.96)

      def test_correct_timing(self):
            # We check that timing has not been broken when playing with the indexes reconstructing the sampled path in pmcmc

            os.system('./kalman -O 3  -c -x --no_dem_sto < ' + Root + '/../examples/noise_test/package.json')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 2
            nbiters = 500
            os.system('./pmcmc sde -t -x --no_dem_sto -O 4 -J ' + str(nparts) + ' -M ' + str(nbiters) + ' -I 3 -N 4 -T '  + str(nbiters) + '  < ' + Root + '/../examples/noise_test/package.json')          
            tab1 = genfromtxt('X_3.csv',delimiter=',',names=True)

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
      run_TransfsAndPMCMC = 1
      run_KalmanOnDiffusions = 1
      run_SMCSDEagainstKalman = 1
      run_pMCMCsmoothing = 1
      run_pMCMCsmoothingWithNaNs = 1

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
