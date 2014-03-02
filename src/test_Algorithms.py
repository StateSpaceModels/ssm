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

            os.system('cp ' + Root + '/../examples/noise/ssm.json .')
            os.system('cp ' + Root + '/../examples/noise/theta.json .')
            os.system('cp -r ' + Root + '/../examples/noise/data .')
            os.system('ssm')

      def setUp(self):
      # Things that need to be done before tests
            os.chdir(Root + '/bin')

      @classmethod
      def tearDownClass(cls):
            #os.chdir(Root + '/bin')
            shutil.rmtree(Root + '/bin')

      def test_ode(self):
            os.system('./simplex -M 10000 < ../../examples/noise/theta.json | ./simplex -M 10000 | ./simplex -M 100000 > theta.json')
            os.system('./smc ode -J 1  --trace < theta.json')
            tab = genfromtxt('trace_0.csv',delimiter=',',names=True).tolist()
            self.assertAlmostEqual(tab[5],-824.598)
            os.system('./smc ode -J 10  --trace < theta.json')
            tab = genfromtxt('trace_0.csv',delimiter=',',names=True).tolist()
            self.assertAlmostEqual(tab[5],-824.598)

      def test_kalman_map(self):
            os.system('./ksimplex -M 1000 --trace < ' + Root + '/../examples/noise/theta.json')
            tab = genfromtxt('trace_0.csv',delimiter=',',names=True)
            self.assertAlmostEqual(tab[450][5],-508.607)

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
            f = open(Root + '/../examples/noise_test/theta.json')
            j = json.load(f)
            j['resources'][0]['data']['r0_paris'] = 10
            with open(Root + '/../examples/noise_test/theta.json','w') as outfile:
                  json.dump(j,outfile)
            f = open(Root + '/../examples/noise_test/data/r0_paris.json')
            j = json.load(f)
            j['distributionParameter'][0]['value']=9.5
            j['distributionParameter'][1]['value']=10.5
            with open(Root + '/../examples/noise_test/data/r0_paris.json','w') as outfile:
                  json.dump(j,outfile)

            os.chdir(Root)
            os.system('cp ' + Root + '/../examples/noise_test/ssm.json .')
            os.system('cp ' + Root + '/../examples/noise_test/theta.json .')
            os.system('cp -r ' + Root + '/../examples/noise_test/data .')
            os.system('ssm')


            os.chdir(Root + '/bin')
            os.system('./pmcmc ode  -C 5000 -W 5000 -O 0 -M 20000 -c -a < ' + Root + '/../examples/noise_test/theta.json')

            test = self.call_test_unif('r0_paris')
            shutil.rmtree(Root + '/bin')
            shutil.rmtree(Root + '/../examples/noise_test')

            self.assertEqual(test,'1')



      def test_prior_normal_transf_log(self):
            os.chdir(Root + '/../examples')
            shutil.copytree('noise','noise_test')
            f = open(Root + '/../examples/noise_test/theta.json')
            j = json.load(f)
            j['resources'][0]['data']['r0_paris'] = 10
            with open(Root + '/../examples/noise_test/theta.json','w') as outfile:
                  json.dump(j,outfile)
            f = open(Root + '/../examples/noise_test/data/r0_paris.json')
            j = json.load(f)
            j['distributionParameter'][0]['value']=0
            j['name']='normal'
            j['distributionParameter'][1] = { "name": "mean", "value": 10 }
            j['distributionParameter'].append({ "name": "sd", "value": 1 })
            with open(Root + '/../examples/noise_test/data/r0_paris.json','w') as outfile:
                  json.dump(j,outfile)

            os.chdir(Root)
            os.system('cp ' + Root + '/../examples/noise_test/ssm.json .')
            os.system('cp ' + Root + '/../examples/noise_test/theta.json .')
            os.system('cp -r ' + Root + '/../examples/noise_test/data .')
            os.system('ssm')

            os.chdir(Root + '/bin')
            os.system('./pmcmc ode  -C 5000 -W 5000 -O 0 -M 20000 -c -a < '+ Root + '/../examples/noise_test/theta.json')

            test = self.call_test_normal('r0_paris')
            shutil.rmtree(Root + '/bin')
            shutil.rmtree(Root + '/../examples/noise_test')

            self.assertEqual(test,'1')





      def test_prior_normal_and_unif_transf_log(self):
            os.chdir(Root + '/../examples')
            shutil.copytree('noise','noise_test')
            f = open(Root + '/../examples/noise_test/theta.json')
            j = json.load(f)
            j['resources'][0]['data']['r0_paris'] = 10
            j['resources'][0]['data']['r0_nyc'] = 10
            with open(Root + '/../examples/noise_test/theta.json','w') as outfile:
                  json.dump(j,outfile)

            f = open(Root + '/../examples/noise_test/data/r0_paris.json')
            j = json.load(f)
            j['distributionParameter'][0]['value']=0
            j['name']='normal'
            j['distributionParameter'][1] = { "name": "mean", "value": 10 }
            j['distributionParameter'].append({ "name": "sd", "value": 1 })
            with open(Root + '/../examples/noise_test/data/r0_paris.json','w') as outfile:
                  json.dump(j,outfile)

            f = open(Root + '/../examples/noise_test/data/r0_nyc.json')
            j = json.load(f)
            j['distributionParameter'][0]['value']=9.5
            j['distributionParameter'][1]['value']=10.5
            with open(Root + '/../examples/noise_test/data/r0_nyc.json','w') as outfile:
                  json.dump(j,outfile)

            os.chdir(Root)
            os.system('cp ' + Root + '/../examples/noise_test/ssm.json .')
            os.system('cp ' + Root + '/../examples/noise_test/theta.json .')
            os.system('cp -r ' + Root + '/../examples/noise_test/data .')
            os.system('ssm')

            os.chdir(Root + '/bin')

            os.system('./pmcmc ode  -C 1000 -W 100000 -O 0 -M 100000 -c -a < ' + Root + '/../examples/noise_test/theta.json')

            test1 = self.call_test_normal('r0_paris')
            test2 = self.call_test_unif('r0_nyc')
            shutil.rmtree(Root + '/bin')
            shutil.rmtree(Root + '/../examples/noise_test')

            self.assertEqual(int(test1)*int(test2),1)



      def test_prior_normal_transf_logit_ab(self):
            os.chdir(Root + '/../examples')
            shutil.copytree('noise','noise_test')
            f = open(Root + '/../examples/noise_test/theta.json')
            j = json.load(f)
            j['resources'][0]['data']['r0_paris'] = 10
            with open(Root + '/../examples/noise_test/theta.json','w') as outfile:
                  json.dump(j,outfile)

            f = open(Root + '/../examples/noise_test/data/r0_paris.json')
            j = json.load(f)
            j['distributionParameter'][0]['value']=0
            j['distributionParameter'][1]['value']=20
            j['name']='normal'
            j['distributionParameter'].append({ "name": "mean", "value": 10 })
            j['distributionParameter'].append({ "name": "sd", "value": 1 })
            with open(Root + '/../examples/noise_test/data/r0_paris.json','w') as outfile:
                  json.dump(j,outfile)


            os.chdir(Root)
            os.system('cp ' + Root + '/../examples/noise_test/ssm.json .')
            os.system('cp ' + Root + '/../examples/noise_test/theta.json .')
            os.system('cp -r ' + Root + '/../examples/noise_test/data .')
            os.system('ssm')


            os.chdir(Root + '/bin')

            os.system('./pmcmc ode  -C 5000 -W 5000 -O 0 -M 20000 -c -a < ' + Root + '/../examples/noise_test/theta.json')

            test = self.call_test_normal('r0_paris')
            shutil.rmtree(Root + '/bin')
            shutil.rmtree(Root + '/../examples/noise_test')

            self.assertEqual(test,'1')


      def test_prior_normal_transf_identity(self):
            os.chdir(Root + '/../examples')
            shutil.copytree('noise','noise_test')
            f = open(Root + '/../examples/noise_test/theta.json')
            j = json.load(f)
            j['resources'][0]['data']['r0_paris'] = 10
            with open(Root + '/../examples/noise_test/theta.json','w') as outfile:
                  json.dump(j,outfile)
            f = open(Root + '/../examples/noise_test/data/r0_paris.json')
            j = json.load(f)
            j['name']='normal'
            j['distributionParameter'][0] = { "name": "mean", "value": 10 }
            j['distributionParameter'][1] = { "name": "sd", "value": 1 }
            with open(Root + '/../examples/noise_test/data/r0_paris.json','w') as outfile:
                  json.dump(j,outfile)

            os.chdir(Root)
            os.system('cp ' + Root + '/../examples/noise_test/ssm.json .')
            os.system('cp ' + Root + '/../examples/noise_test/theta.json .')
            os.system('cp -r ' + Root + '/../examples/noise_test/data .')
            os.system('ssm')


            os.chdir(Root + '/bin')

            os.system('./pmcmc ode  -C 5000 -W 5000 -O 0 -M 20000 -c -a < ' + Root + '/../examples/noise_test/theta.json')

            test = self.call_test_normal('r0_paris')
            shutil.rmtree(Root + '/bin')
            shutil.rmtree(Root + '/../examples/noise_test')

            self.assertEqual(test,'1')

      def call_test_unif(self, varname):
            shutil.copyfile(Root+'/TestsR/test_unif.r',Root + '/bin/test_unif.r')
            os.chdir(Root + '/bin')
            os.system('R --vanilla < test_unif.r ' + varname  + '> /dev/null 2>&1')
            f = open("outfile.txt","r")
            x = f.readlines()
            return x[0]

      def call_test_normal(self,varname):
            shutil.copyfile(Root+'/TestsR/test_normal.r',Root + '/bin/test_normal.r')
            os.chdir(Root + '/bin')
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
            f = open(Root + '/../examples/noise_test/ssm.json')
            j = json.load(f)
            j['inputs'].insert(0,{})
            j['inputs'].insert(0,{})
            j['inputs'][0]['name'] = 'test_par'
            j['inputs'][0]['description'] = ''
            j['inputs'][0]['require'] = { "name": "test_par", "path": "data/test_par.json"}
            j['inputs'][1]['name'] = 'test_vol'
            j['inputs'][1]['description'] = ''
            j['inputs'][1]['require'] = { "name": "test_vol", "path": "data/test_vol.json"}
            j['sde'] = {}
            j['sde']['drift']=[]
            j['sde']['drift'].append({})
            j['sde']['drift'][0]['name']='test_par'
            j['sde']['drift'][0]['f']=0.0
            j['sde']['dispersion']=[['test_vol']]
            with open(Root + '/../examples/noise_test/ssm.json','w') as outfile:
                  json.dump(j,outfile)

            os.system('cp ' + Root + '/../examples/noise_test/data/sto.json ' + Root + '/../examples/noise_test/data/test_par.json')
            f = open(Root + '/../examples/noise_test/data/test_par.json')
            j = json.load(f)
            j['name'] = "dirac"
            j['distributionParameter'] = [{ "value" : 0 }]
            with open(Root + '/../examples/noise_test/data/test_par.json','w') as outfile:
              json.dump(j,outfile)
            os.system('cp ' + Root + '/../examples/noise_test/data/r0_paris.json ' + Root + '/../examples/noise_test/data/test_vol.json')
            f = open(Root + '/../examples/noise_test/data/test_vol.json')
            j = json.load(f)
            j['name'] = "dirac"
            j['distributionParameter'] = [{ "value" : 0.1428571 }]
            with open(Root + '/../examples/noise_test/data/test_vol.json','w') as outfile:
              json.dump(j,outfile)

            os.chdir(Root)
            os.system('cp ' + Root + '/../examples/noise_test/ssm.json .')
            os.system('cp ' + Root + '/../examples/noise_test/theta.json .')
            os.system('cp -r ' + Root + '/../examples/noise_test/data .')
            os.system('ssm')

            os.chdir(Root + '/bin')


      def setUp(self):
      # Things that need to be done before tests
            os.chdir(Root + '/bin')

      @classmethod
      def tearDownClass(cls):
            os.chdir(Root + '/bin/')
            shutil.rmtree(Root + '/../examples/noise_test')

      def test_1step(self):

            os.system('./kalman -O 2  -c -x < ' + Root + '/../examples/noise_test/theta.json')
            tab = genfromtxt('hat_0.csv',delimiter=',',names=True)
            self.assertAlmostEqual(tab['lower_test_par'][0],-1.96/math.sqrt(7),5)

      def test_10step(self):
            os.system('./kalman -O 10  -c -x < ' + Root + '/../examples/noise_test/theta.json')
            tab = genfromtxt('hat_0.csv',delimiter=',',names=True)
            self.assertAlmostEqual(tab['lower_test_par'][9]/math.sqrt(10),-1.96/math.sqrt(7),5)

class TestSMCSDEagainstKalman(unittest.TestCase):
      @classmethod
      def setUpClass(cls):
            print("")
            print("Testing SMC against EKF on linear examples")

            os.chdir(Root + '/../examples')
#            shutil.copytree('noise','noise_test')
            f = open(Root + '/../examples/noise_test/ssm.json')
            j = json.load(f)
            j["populations"] = []
            j["populations"].append({"name":"paris","composition":["S","I","I2","E1","E2"]})
            j["reactions"] = []
            j["reactions"].append({"from":"S", "to":"I", "rate":"1000/S", "white_noise": {"name":"noise_SI", "sd": "vol"}})
            j["reactions"].append({"from":"I", "to":"I2", "rate":"1000/I", "white_noise": {"name":"noise_SI", "sd": "vol"}})
            j["reactions"].append({"from":"I2", "to":"S", "rate":"1000/I2", "white_noise": {"name":"noise_SI2", "sd": "vol2"}})
            j["reactions"].append({"from":"E1", "to":"E2", "rate":"3/N", "accumulators": ["all_inc_out"]})

            j["observations"] = [j["observations"][0]]

            j["data"] = [j["data"][0]]

            j["inputs"] = []
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'r0'
            j["inputs"][0]['require'] = {'name':'r0', 'path': 'data/r0.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'S'
            j["inputs"][0]['require'] = {'name':'S', 'path': 'data/S.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'I'
            j["inputs"][0]['require'] = {'name':'I', 'path': 'data/I.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'I2'
            j["inputs"][0]['require'] = {'name':'I2', 'path': 'data/I2.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'E1'
            j["inputs"][0]['require'] = {'name':'E1', 'path': 'data/E1.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'vol'
            j["inputs"][0]['require'] = {'name':'vol', 'path': 'data/vol.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'vol2'
            j["inputs"][0]['require'] = {'name':'vol2', 'path': 'data/vol2.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'rep_all_CDC_inc'
            j["inputs"][0]['require'] = {'name':'rep_all_CDC_inc', 'path': 'data/rep_all_CDC_inc.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'prop_all_CDC_inc'
            j["inputs"][0]['require'] = {'name':'prop_all_CDC_inc', 'path': 'data/prop_all_CDC_inc.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'phi'
            j["inputs"][0]['require'] = {'name':'phi', 'path': 'data/phi.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'E2'
            j["inputs"][0]['require'] = {'name':'E2', 'path': 'data/E2.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'N'
            j["inputs"][0]['require'] = {'name':'N', 'path': 'data/N.json'}

            j["sde"] = {}
            j["sde"]["drift"] = []
            j["sde"]["drift"].append({ "name": "r0", "f": 0.0})
            j["sde"]["dispersion"] = [['vol']]

            with open(Root + '/../examples/noise_test/ssm.json','w') as outfile:
                  json.dump(j,outfile)


            j={}
            j['name'] = 'uniform'
            j['distributionParameter'] = []
            j['distributionParameter'].append({ 'name': 'lower', 'value': 50 })
            j['distributionParameter'].append({ 'name': 'upper', 'value': 150 })
            with open(Root + '/../examples/noise_test/data/r0.json','w') as outfile:
                  json.dump(j,outfile)
            j={}
            j['name'] = 'dirac'
            j['distributionParameter'] = []
            j['distributionParameter'].append({ 'value': 200000 })
            with open(Root + '/../examples/noise_test/data/S.json','w') as outfile:
                  json.dump(j,outfile)
            with open(Root + '/../examples/noise_test/data/I.json','w') as outfile:
                  json.dump(j,outfile)
            with open(Root + '/../examples/noise_test/data/I2.json','w') as outfile:
                  json.dump(j,outfile)
            with open(Root + '/../examples/noise_test/data/E1.json','w') as outfile:
                  json.dump(j,outfile)
            with open(Root + '/../examples/noise_test/data/E2.json','w') as outfile:
                  json.dump(j,outfile)
            j={}
            j['name'] = 'dirac'
            j['distributionParameter'] = []
            j['distributionParameter'].append({ 'value': 0.01428571 })
            with open(Root + '/../examples/noise_test/data/vol.json','w') as outfile:
                  json.dump(j,outfile)
            j={}
            j['name'] = 'dirac'
            j['distributionParameter'] = []
            j['distributionParameter'].append({ 'value': 0.02857143 })
            with open(Root + '/../examples/noise_test/data/vol2.json','w') as outfile:
                  json.dump(j,outfile)
            j={}
            j['name'] = 'dirac'
            j['distributionParameter'] = []
            j['distributionParameter'].append({ 'value': 0.2 })
            with open(Root + '/../examples/noise_test/data/rep_all_CDC_inv.json','w') as outfile:
                  json.dump(j,outfile)
            with open(Root + '/../examples/noise_test/data/prop_all_CDC_inc.json','w') as outfile:
                  json.dump(j,outfile)
            with open(Root + '/../examples/noise_test/data/phi.json','w') as outfile:
                  json.dump(j,outfile)
            with open(Root + '/../examples/noise_test/data/N.json','w') as outfile:
                  json.dump(j,outfile)

            j = { 'resources': [] }
            j['resources'].append({ 'name': 'values', 'data': {'r0': 100 }})
            j['resources'].append({ 'name': 'covariance', 'data': {'r0': { 'r0': 0.02 }}})


            with open(Root + '/../examples/noise_test/theta.json','w') as outfile:
                  json.dump(j,outfile)

            os.chdir(Root)
            os.system('cp ' + Root + '/../examples/noise_test/ssm.json .')
            os.system('cp ' + Root + '/../examples/noise_test/theta.json .')
            os.system('cp -r ' + Root + '/../examples/noise_test/data .')
            os.system('ssm')

            os.chdir(Root + '/bin')

      def setUp(self):
      # Things that need to be done before tests
            os.chdir(Root + '/bin')

      @classmethod
      def tearDownClass(cls):
            os.chdir(Root + '/bin/')
            shutil.rmtree(Root + '/../examples/noise_test')

      def test_only_env_sto(self):
            os.system('./kalman -O 2  -c -x --no_dem_sto < ' + Root + '/../examples/noise_test/theta.json')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 500
            os.system('./smc sde --no_dem_sto -t -O 2 -x -J ' + str(nparts) + ' -I 1 -N 4 --dt 0.0001  < ' + Root + '/../examples/noise_test/theta.json')
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

            f = open(Root + '/../examples/noise_test/ssm.json')
            j = json.load(f)
            j["populations"] = []
            j["populations"].append({"name":"paris","composition":["S","I","I2","E1","E2"]})
            j["reactions"] = []
            j["reactions"].append({"from":"S", "to":"I", "rate":"1000/S", "white_noise": {"name":"noise_SI", "sd": "vol"}})
            j["reactions"].append({"from":"I", "to":"I2", "rate":"1000/I", "white_noise": {"name":"noise_SI", "sd": "vol"}})
            j["reactions"].append({"from":"I2", "to":"S", "rate":"1000/I2", "white_noise": {"name":"noise_SI2", "sd": "vol2"}})
            j["reactions"].append({"from":"E1", "to":"E2", "rate":"3/N", "accumulators": ["all_inc_out"]})

            j["observations"] = [j["observations"][0]]

            j["data"] = [j["data"][0]]

            j["inputs"] = []
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'r0'
            j["inputs"][0]['require'] = {'name':'r0', 'path': 'data/r0.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'S'
            j["inputs"][0]['require'] = {'name':'S', 'path': 'data/S.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'I'
            j["inputs"][0]['require'] = {'name':'I', 'path': 'data/I.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'I2'
            j["inputs"][0]['require'] = {'name':'I2', 'path': 'data/I2.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'E1'
            j["inputs"][0]['require'] = {'name':'E1', 'path': 'data/E1.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'vol'
            j["inputs"][0]['require'] = {'name':'vol', 'path': 'data/vol.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'vol2'
            j["inputs"][0]['require'] = {'name':'vol2', 'path': 'data/vol2.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'rep_all_CDC_inc'
            j["inputs"][0]['require'] = {'name':'rep_all_CDC_inc', 'path': 'data/rep_all_CDC_inc.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'prop_all_CDC_inc'
            j["inputs"][0]['require'] = {'name':'prop_all_CDC_inc', 'path': 'data/prop_all_CDC_inc.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'phi'
            j["inputs"][0]['require'] = {'name':'phi', 'path': 'data/phi.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'E2'
            j["inputs"][0]['require'] = {'name':'E2', 'path': 'data/E2.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'N'
            j["inputs"][0]['require'] = {'name':'N', 'path': 'data/N.json'}

            j["sde"] = {}
            j["sde"]["drift"] = []
            j["sde"]["drift"].append({ "name": "r0", "f": 0.0})
            j["sde"]["dispersion"] = [['vol']]

            with open(Root + '/../examples/noise_test/ssm.json','w') as outfile:
                  json.dump(j,outfile)


            j={}
            j['name'] = 'dirac'
            j['distributionParameter'] = []
            j['distributionParameter'].append({ 'value': 100 })
            with open(Root + '/../examples/noise_test/data/r0.json','w') as outfile:
                  json.dump(j,outfile)
            j={}
            j['name'] = 'dirac'
            j['distributionParameter'] = []
            j['distributionParameter'].append({ 'value': 200000 })
            with open(Root + '/../examples/noise_test/data/S.json','w') as outfile:
                  json.dump(j,outfile)
            with open(Root + '/../examples/noise_test/data/I.json','w') as outfile:
                  json.dump(j,outfile)
            with open(Root + '/../examples/noise_test/data/I2.json','w') as outfile:
                  json.dump(j,outfile)
            with open(Root + '/../examples/noise_test/data/E1.json','w') as outfile:
                  json.dump(j,outfile)

            j={}
            j['name'] = 'uniform'
            j['distributionParameter'] = []
            j['distributionParameter'].append({ 'name': 'lower', 'value': 199995 })
            j['distributionParameter'].append({ 'name': 'upper', 'value': 200010 })
            with open(Root + '/../examples/noise_test/data/E2.json','w') as outfile:
                  json.dump(j,outfile)
            j={}
            j['name'] = 'dirac'
            j['distributionParameter'] = []
            j['distributionParameter'].append({ 'value': 0.01428571 })
            with open(Root + '/../examples/noise_test/data/vol.json','w') as outfile:
                  json.dump(j,outfile)
            j={}
            j['name'] = 'dirac'
            j['distributionParameter'] = []
            j['distributionParameter'].append({ 'value': 0.02857143 })
            with open(Root + '/../examples/noise_test/data/vol2.json','w') as outfile:
                  json.dump(j,outfile)
            j={}
            j['name'] = 'dirac'
            j['distributionParameter'] = []
            j['distributionParameter'].append({ 'value': 0.2 })
            with open(Root + '/../examples/noise_test/data/rep_all_CDC_inv.json','w') as outfile:
                  json.dump(j,outfile)
            with open(Root + '/../examples/noise_test/data/prop_all_CDC_inc.json','w') as outfile:
                  json.dump(j,outfile)
            with open(Root + '/../examples/noise_test/data/phi.json','w') as outfile:
                  json.dump(j,outfile)
            with open(Root + '/../examples/noise_test/data/N.json','w') as outfile:
                  json.dump(j,outfile)

            j = { 'resources': [] }
            j['resources'].append({ 'name': 'values', 'data': {'E2': 200000 }})
            j['resources'].append({ 'name': 'covariance', 'data': {'E2': { 'E2': 0.02 }}})

            with open(Root + '/../examples/noise_test/theta.json','w') as outfile:
                  json.dump(j,outfile)


            os.chdir(Root)
            os.system('cp ' + Root + '/../examples/noise_test/ssm.json .')
            os.system('cp ' + Root + '/../examples/noise_test/theta.json .')
            os.system('cp -r ' + Root + '/../examples/noise_test/data .')
            os.system('ssm')

            os.chdir(Root + '/bin')

      def setUp(self):
      # Things that need to be done before tests
            os.chdir(Root + '/bin')

      @classmethod
      def tearDownClass(cls):
            shutil.rmtree(Root + '/../examples/noise_test')

      def test_particle_genealogy(self):
            # We assess that the exploration of the genealogy to reconstruct sampled traj is correct by checking the trajectory
            # of the random walk through the distribution of its increments

            os.system('./kalman -O 2  -c -x --no_dem_sto < ' + Root + '/../examples/noise_test/theta.json')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 100
            nbiters = 2

            os.system('./pmcmc sde -t -x --no_dem_sto -J ' + str(nparts) + ' -M ' + str(nbiters) + ' -I 1 -N 4 -T '  + str(nbiters) + '  < ' + Root + '/../examples/noise_test/theta.json')
            tab1 = genfromtxt('X_1.csv',delimiter=',',names=True)

            x = []
            for i in range(1,53):
                  # Compute and normalise increments of sampled path
                  x.append((float(tab1[52 - i][7])-float(tab1[52 - i + 1][7]))/(math.sqrt(7)*0.1/7))

            self.assertTrue(abs((stats.scoreatpercentile(x,50)-0)/(math.sqrt(0.5*0.5)/(stats.norm.pdf(0,loc=0,scale=1)*math.sqrt(len(x)))))<1.96)

      def test_correct_timing(self):
            # We check that timing has not been broken when playing with the indexes reconstructing the sampled path in pmcmc

            os.system('./kalman -O 3  -c -x --no_dem_sto < ' + Root + '/../examples/noise_test/theta.json')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 2
            nbiters = 500
            os.system('./pmcmc sde -t -x --no_dem_sto -O 4 -J ' + str(nparts) + ' -M ' + str(nbiters) + ' -I 2 -N 4 -T '  + str(nbiters) + '  < ' + Root + '/../examples/noise_test/theta.json')
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

            f = open(Root + '/../examples/noise_test/ssm.json')
            j = json.load(f)
            j["populations"] = []
            j["populations"].append({"name":"paris","composition":["S","I","I2","E1","E2"]})
            j["reactions"] = []
            j["reactions"].append({"from":"S", "to":"I", "rate":"1000/S", "white_noise": {"name":"noise_SI", "sd": "vol"}})
            j["reactions"].append({"from":"I", "to":"I2", "rate":"1000/I", "white_noise": {"name":"noise_SI", "sd": "vol"}})
            j["reactions"].append({"from":"I2", "to":"S", "rate":"1000/I2", "white_noise": {"name":"noise_SI2", "sd": "vol2"}})
            j["reactions"].append({"from":"E1", "to":"E2", "rate":"3/N", "accumulators": ["all_inc_out"]})

            j["observations"] = [j["observations"][0]]

            j["data"] = [j["data"][0]]

            j["inputs"] = []
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'r0'
            j["inputs"][0]['require'] = {'name':'r0', 'path': 'data/r0.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'S'
            j["inputs"][0]['require'] = {'name':'S', 'path': 'data/S.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'I'
            j["inputs"][0]['require'] = {'name':'I', 'path': 'data/I.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'I2'
            j["inputs"][0]['require'] = {'name':'I2', 'path': 'data/I2.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'E1'
            j["inputs"][0]['require'] = {'name':'E1', 'path': 'data/E1.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'vol'
            j["inputs"][0]['require'] = {'name':'vol', 'path': 'data/vol.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'vol2'
            j["inputs"][0]['require'] = {'name':'vol2', 'path': 'data/vol2.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'rep_all_CDC_inc'
            j["inputs"][0]['require'] = {'name':'rep_all_CDC_inc', 'path': 'data/rep_all_CDC_inc.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'prop_all_CDC_inc'
            j["inputs"][0]['require'] = {'name':'prop_all_CDC_inc', 'path': 'data/prop_all_CDC_inc.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'phi'
            j["inputs"][0]['require'] = {'name':'phi', 'path': 'data/phi.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'E2'
            j["inputs"][0]['require'] = {'name':'E2', 'path': 'data/E2.json'}
            j["inputs"].insert(0,{})
            j["inputs"][0]['name'] = 'N'
            j["inputs"][0]['require'] = {'name':'N', 'path': 'data/N.json'}

            j["sde"] = {}
            j["sde"]["drift"] = []
            j["sde"]["drift"].append({ "name": "r0", "f": 0.0})
            j["sde"]["dispersion"] = [['vol']]

            with open(Root + '/../examples/noise_test/ssm.json','w') as outfile:
                  json.dump(j,outfile)


            j={}
            j['name'] = 'dirac'
            j['distributionParameter'] = []
            j['distributionParameter'].append({ 'value': 100 })
            with open(Root + '/../examples/noise_test/data/r0.json','w') as outfile:
                  json.dump(j,outfile)
            j={}
            j['name'] = 'dirac'
            j['distributionParameter'] = []
            j['distributionParameter'].append({ 'value': 200000 })
            with open(Root + '/../examples/noise_test/data/S.json','w') as outfile:
                  json.dump(j,outfile)
            with open(Root + '/../examples/noise_test/data/I.json','w') as outfile:
                  json.dump(j,outfile)
            with open(Root + '/../examples/noise_test/data/I2.json','w') as outfile:
                  json.dump(j,outfile)
            with open(Root + '/../examples/noise_test/data/E1.json','w') as outfile:
                  json.dump(j,outfile)

            j={}
            j['name'] = 'uniform'
            j['distributionParameter'] = []
            j['distributionParameter'].append({ 'name': 'lower', 'value': 199995 })
            j['distributionParameter'].append({ 'name': 'upper', 'value': 200010 })
            with open(Root + '/../examples/noise_test/data/E2.json','w') as outfile:
                  json.dump(j,outfile)
            j={}
            j['name'] = 'dirac'
            j['distributionParameter'] = []
            j['distributionParameter'].append({ 'value': 0.01428571 })
            with open(Root + '/../examples/noise_test/data/vol.json','w') as outfile:
                  json.dump(j,outfile)
            j={}
            j['name'] = 'dirac'
            j['distributionParameter'] = []
            j['distributionParameter'].append({ 'value': 0.02857143 })
            with open(Root + '/../examples/noise_test/data/vol2.json','w') as outfile:
                  json.dump(j,outfile)
            j={}
            j['name'] = 'dirac'
            j['distributionParameter'] = []
            j['distributionParameter'].append({ 'value': 0.2 })
            with open(Root + '/../examples/noise_test/data/rep_all_CDC_inv.json','w') as outfile:
                  json.dump(j,outfile)
            with open(Root + '/../examples/noise_test/data/prop_all_CDC_inc.json','w') as outfile:
                  json.dump(j,outfile)
            with open(Root + '/../examples/noise_test/data/phi.json','w') as outfile:
                  json.dump(j,outfile)
            with open(Root + '/../examples/noise_test/data/N.json','w') as outfile:
                  json.dump(j,outfile)

            j = { 'resources': [] }
            j['resources'].append({ 'name': 'values', 'data': {'E2': 200000 }})
            j['resources'].append({ 'name': 'covariance', 'data': {'E2': { 'E2': 0.02 }}})

            with open(Root + '/../examples/noise_test/theta.json','w') as outfile:
                  json.dump(j,outfile)


            os.chdir(Root + '/../examples/noise_test/data')
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
            os.system('cp ' + Root + '/../examples/noise_test/ssm.json .')
            os.system('cp ' + Root + '/../examples/noise_test/theta.json .')
            os.system('cp -r ' + Root + '/../examples/noise_test/data .')
            os.system('ssm')


            os.chdir(Root + '/bin')


      def setUp(self):
            # Things that need to be done before tests
            os.chdir(Root + '/bin')

      @classmethod
      def tearDownClass(cls):
            os.chdir(Root + '/bin')
            shutil.rmtree(Root + '/../examples/noise_test')

      def test_particle_genealogy(self):
            # We assess that the exploration of the genealogy to reconstruct sampled traj is correct by checking the trajectory
            # of the random walk through the distribution of its increments

            os.system('./kalman -O 2  -c -x --no_dem_sto < ' + Root + '/../examples/noise_test/theta.json')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 1000
            nbiters = 2

            os.system('./pmcmc sde -t -x --no_dem_sto -J ' + str(nparts) + ' -M ' + str(nbiters) + ' -I 2 -N 4 -T '  + str(nbiters) + '  < ' + Root + '/../examples/noise_test/theta.json')
            tab1 = genfromtxt('X_2.csv',delimiter=',',names=True)

            x = []
            for i in range(1,53):
                  # Compute and normalise increments of sampled path
                  x.append((float(tab1[52 - i][7])-float(tab1[52 - i + 1][7]))/(math.sqrt(7)*0.1/7))

            print(abs((stats.scoreatpercentile(x,50)-0)/(math.sqrt(0.5*0.5)/(stats.norm.pdf(0,loc=0,scale=1)*math.sqrt(len(x))))))
            self.assertTrue(abs((stats.scoreatpercentile(x,50)-0)/(math.sqrt(0.5*0.5)/(stats.norm.pdf(0,loc=0,scale=1)*math.sqrt(len(x)))))<1.96)

      def test_correct_timing(self):
            # We check that timing has not been broken when playing with the indexes reconstructing the sampled path in pmcmc

            os.system('./kalman -O 3  -c -x --no_dem_sto < ' + Root + '/../examples/noise_test/theta.json')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 2
            nbiters = 500
            os.system('./pmcmc sde -t -x --no_dem_sto -O 4 -J ' + str(nparts) + ' -M ' + str(nbiters) + ' -I 3 -N 4 -T '  + str(nbiters) + '  < ' + Root + '/../examples/noise_test/theta.json')
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

      run_NoiseResults = 0
      run_TransfsAndPMCMC = 0
      run_KalmanOnDiffusions = 0
      run_SMCSDEagainstKalman = 0
      run_pMCMCsmoothing = 0
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
