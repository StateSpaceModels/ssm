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
            os.system('./smc < /Users/dureaujoseph/ssm/example/foo/datapackage.json')

      
      def setUp(self):
      # Things that need to be done before tests
            os.system('./smc < /Users/dureaujoseph/ssm/example/foo/datapackage.json')

      @classmethod
      def tearDownClass(cls):
            os.system('./smc < /Users/dureaujoseph/ssm/example/foo/datapackage.json')
            #shutil.rmtree('/Users/dureaujoseph/ssm_test_model/')

      def test_kalman_map(self):
            os.system('./ksimplex --prior -M 100 -c < /Users/dureaujoseph/ssm/example/noise/datapackage.json')
            tab = genfromtxt('trace_0.csv',delimiter=',',names=True)
            self.assertAlmostEqual(tab[35][5],-729.491)

class TestTransfsAndPMCMC(unittest.TestCase):
      @classmethod
      def setUpClass(cls):
            print("")
            print("Testing transformations and pMCMC by sampling from priors")
            # copy noise from the examples and build it
            shutil.copytree(Root + '/../example/noise',Root + '/noise')     
            os.chdir(Root+'/noise')
            os.system('plom build -t theta.json --local'  + "> /dev/null 2>&1")
            
      def setUp(self):
      # Things that need to be done before tests
            os.chdir(Root+'/noise/model')

      @classmethod
      def tearDownClass(cls):
            shutil.rmtree(Root + '/noise')

            
      def test_prior_unif_transf_log(self):
            os.system('plom pipe theta.json -S r0:city1__all:distribution:uniform,r0:city1__all:guess:10,r0:city1__all:min:9.5,r0:city1__all:max:10.5 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full --quiet')
            os.system('plom pipe theta.json -T -C -S r0:city1__all:distribution:uniform,r0:city1__all:guess:10,r0:city1__all:min:9.5,r0:city1__all:max:10.5 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full --quiet')
            test = self.call_test_unif('r0.city1__all')
            self.assertEqual(test,'1')

      def test_prior_normal_transf_log(self):
            os.system('plom pipe theta.json -S r0:city1__all:distribution:normal,r0:city1__all:guess:10,r0:city1__all:min:8.04,r0:city1__all:max:11.96 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full --quiet')
            os.system('plom pipe theta.json -T -C -S r0:city1__all:distribution:normal,r0:city1__all:guess:10,r0:city1__all:min:8.04,r0:city1__all:max:11.96 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full --quiet')
            test = self.call_test_normal('r0.city1__all')
            self.assertEqual(test,'1')

      def test_prior_normal_and_unif_transf_log(self):
            os.system('plom pipe theta.json -S r0:city1__all:distribution:normal,r0:city1__all:guess:10,r0:city1__all:min:8.04,r0:city1__all:max:11.96,r0:city2__all:distribution:uniform,r0:city2__all:guess:10,r0:city2__all:min:9.5,r0:city2__all:max:10.5 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full --quiet')
            os.system('plom pipe theta.json -T -C -S r0:city1__all:distribution:normal,r0:city1__all:guess:10,r0:city1__all:min:8.04,r0:city1__all:max:11.96,r0:city2__all:distribution:uniform,r0:city2__all:guess:10,r0:city2__all:min:9.5,r0:city2__all:max:10.5 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full --quiet')
            test1 = self.call_test_normal('r0.city1__all')
            test2 = self.call_test_unif('r0.city2__all')
            self.assertEqual(int(test1)*int(test2),1)

      def test_prior_unif_transf_logit_ab(self):
            os.system('plom pipe theta.json -S r0:all:transformation:logit_ab,r0:city1__all:distribution:uniform,r0:city1__all:guess:10,r0:city1__all:min:9.5,r0:city1__all:max:10.5 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full --quiet')
            os.system('plom pipe theta.json -T -C -S r0:all:transformation:logit_ab,r0:city1__all:distribution:uniform,r0:city1__all:guess:10,r0:city1__all:min:9.5,r0:city1__all:max:10.5 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full --quiet')
            test = self.call_test_unif('r0.city1__all')
            self.assertEqual(test,'1')

      def test_prior_unif_transf_identity(self):
            os.system('plom pipe theta.json -S r0:all:transformation:identity,r0:city1__all:distribution:uniform,r0:city1__all:guess:10,r0:city1__all:min:9.5,r0:city1__all:max:10.5 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full --quiet')
            os.system('plom pipe theta.json -T -C -S r0:all:transformation:identity,r0:city1__all:distribution:uniform,r0:city1__all:guess:10,r0:city1__all:min:9.5,r0:city1__all:max:10.5 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full --quiet')
            test = self.call_test_unif('r0.city1__all')
            self.assertEqual(test,'1')

      def test_prior_normal_transf_identity(self):
            os.system('plom pipe theta.json -S r0:all:transformation:identity,r0:city1__all:distribution:normal,r0:city1__all:guess:10,r0:city1__all:min:8.04,r0:city1__all:max:11.96 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full --quiet')
            os.system('plom pipe theta.json -T -C -S r0:all:transformation:identity,r0:city1__all:distribution:normal,r0:city1__all:guess:10,r0:city1__all:min:8.04,r0:city1__all:max:11.96 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full --quiet')
            test = self.call_test_normal('r0.city1__all')
            self.assertEqual(test,'1')

      def call_test_unif(self,varname):
            shutil.copyfile(Root+'/TestsR/test_unif.r',Root+'/noise/model/test_R0_unif.r')
            os.system('R --vanilla < test_R0_unif.r ' + varname  + "> /dev/null 2>&1")
            f = open(Root+"/noise/model/outfile.txt","r")
            x = f.readlines()
            return x[0]

      def call_test_normal(self,varname):
            shutil.copyfile(Root+'/TestsR/test_normal.r',Root+'/noise/model/test_R0_normal.r')
            os.system('R --vanilla < test_R0_normal.r ' + varname  + "> /dev/null 2>&1")
            f = open(Root+"/noise/model/outfile.txt","r")
            x = f.readlines()
            return x[0]


class TestKalmanOnDiffusions(unittest.TestCase):
      @classmethod
      def setUpClass(cls):
            print("")
            print("Testing EKF on linear examples")
            # copy noise from the examples, add a diffusing parameter for tests, and build it
            shutil.copytree(Root + '/../example/noise',Root + '/noise_test_diff')     
            os.chdir(Root+'/noise_test_diff')
            p = json.load(open('process.json'))
            t = json.load(open('theta.json'))
            p["parameter"].append({'id':'test_par','comment':'will be freely diffusing'})
            p["parameter"].append({'id':'test_vol','comment':'will be equal to 1'})
            p["diffusion"]=[]
            p["diffusion"].append({'parameter':'test_par','volatility':'test_vol','drift':0.0})
            t["parameter"]['test_par']={'min':1,'max':1,'guess':1,'sd_transf':0.0}
            t["parameter"]['test_vol']={'min':1,'max':1,'guess':1,'sd_transf':0.0,'unit':'W'}
            json.dump(p,open('process.json','w'))
            json.dump(t,open('theta.json','w'))
            os.system('plom build -t theta.json --local'  + "> /dev/null 2>&1")
      
      def setUp(self):
      # Things that need to be done before tests
            os.chdir(Root+'/noise_test_diff/model')

      @classmethod
      def tearDownClass(cls):
            shutil.rmtree(Root + '/noise_test_diff')
            
      def test_1step(self):
            os.system('plom pipe theta.json | ./kalman -o 2 --quiet')
            tab = genfromtxt('hat_0.csv',delimiter=',',names=True)
            self.assertAlmostEqual(tab['low95difftest_parall'][1],-1.96/math.sqrt(7),5)

      def test_10step(self):
            os.system('plom pipe theta.json | ./kalman -o 10 --quiet')
            tab = genfromtxt('hat_0.csv',delimiter=',',names=True)
            self.assertAlmostEqual(tab['low95difftest_parall'][10]/math.sqrt(10),-1.96/math.sqrt(7),5)

class TestSMCSDEagainstKalman(unittest.TestCase):
      @classmethod
      def setUpClass(cls):
            print("")
            print("Testing SMC against EKF on linear examples")
            shutil.copytree(Root + '/../example/linear', Root + '/linear')     
            os.chdir(Root+'/linear')
            os.system('plom build -t theta.json --local'  + "> /dev/null 2>&1")
            
      def setUp(self):
      # Things that need to be done before tests
            os.chdir(Root+'/linear/model')

      @classmethod
      def tearDownClass(cls):
            shutil.rmtree(Root + '/linear')

      def test_only_env_sto(self):
            
            os.system('plom pipe theta.json |  ./kalman --no_dem_sto --traj -o 2 --quiet')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 500
            os.system('plom pipe theta.json |  ./smc sde --no_dem_sto --traj -o 2 -J ' + str(nparts) + ' -i 1 -N 4 --DT 0.0001 --quiet'  + "> /dev/null 2>&1")
            tab1 = genfromtxt('hat_1.csv',delimiter=',',names=True)

            meanSMC0a = tab1[1][2]
            meanEKF0a = tab0[1][2]
            q975SMC0a = tab1[1][3]
            q975EKF0a = tab0[1][3]
            meanSMC0b = tab1[1][8]
            meanEKF0b = tab0[1][8]
            q975SMC0b = tab1[1][9]
            q975EKF0b = tab0[1][9]
            meanSMC0r0 = math.log(tab1[1][20])
            meanEKF0r0 = tab0[1][20]
            q975SMC0r0 = math.log(tab1[1][21])
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
            shutil.copytree(Root + '/../example/linear',Root + '/linear')     
            os.chdir(Root+'/linear')
            os.system('plom build -t theta.json --local'  + "> /dev/null 2>&1")
            
      def setUp(self):
      # Things that need to be done before tests
            os.chdir(Root+'/linear/model')

      @classmethod
      def tearDownClass(cls):
            shutil.rmtree(Root + '/linear')

      def test_particle_genealogy(self):
            # We assess that the exploration of the genealogy to reconstruct sampled traj is correct by checking the trajectory
            # of the random walk through the distribution of its increments
            
            os.system('plom pipe theta.json |  ./kalman --no_dem_sto --traj -o 2 --quiet')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 100
            nbiters = 1
            os.system('plom pipe theta.json |  ./pmcmc sde --full --no_dem_sto -J ' + str(nparts) + ' -M ' + str(nbiters) + ' -i 1 -N 4 -n '  + str(nbiters) + ' --quiet'  + "> /dev/null 2>&1")          
            tab1 = genfromtxt('X_1.csv',delimiter=',',names=True)

            x = []
            for i in range(1,53):
                  # Compute and normalise increments of sampled path
                  x.append((math.log(float(tab1[52 - i][8]))-math.log(float(tab1[52 - i + 1][8])))/(math.sqrt(7)*0.1/7))          

            self.assertTrue(abs((stats.scoreatpercentile(x,50)-0)/(math.sqrt(0.5*0.5)/(stats.norm.pdf(0,loc=0,scale=1)*math.sqrt(len(x)))))<1.96)

            
      def test_correct_timing(self):
            # We check that timing has not been broken when playing with the indexes reconstructing the sampled path in pmcmc
            
            os.system('plom pipe theta.json |  ./kalman --no_dem_sto --traj -o 2 --quiet')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 2
            nbiters = 500
            os.system('plom pipe theta.json |  ./pmcmc sde --full --no_dem_sto -o 2 -J ' + str(nparts) + ' -M ' + str(nbiters) + ' -i 2 -N 4 --DT 0.0001 -n '  + str(nbiters) + ' --quiet'  + "> /dev/null 2>&1")          
            tab1 = genfromtxt('X_2.csv',delimiter=',',names=True)

            ### t = 1 (first obs)
            x = []
            for i in range(1,nbiters):
                  x.append(math.log(float(tab1[3*i - 2][8])))
                  
            meanSMC0r0 = numpy.mean(x)
            meanEKF0r0 = tab0[1][20]
            q975SMC0r0 = stats.scoreatpercentile(x,97.5)
            q975EKF0r0 = tab0[1][21]

            # Tests on 97.5% quantiles
            # based on on CLT for empirical quantile given in "Statistics and Data Analysis for Financial Engineering, Rupert 2011"
            self.assertTrue(abs((q975SMC0r0-q975EKF0r0)/(math.sqrt(0.975*0.025)/(stats.norm.pdf(q975EKF0r0-meanEKF0r0,loc=0,scale=(q975EKF0r0-meanEKF0r0)/1.96)*math.sqrt(nparts))))<1.96)

            ### t = 2 (second obs)
            x = []
            for i in range(1,nbiters):
                  x.append(math.log(float(tab1[3*i - 3][8])))

            meanSMC0r0 = numpy.mean(x)
            meanEKF0r0 = tab0[2][20]
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
            shutil.copytree(Root + '/../example/linear',Root + '/linear')     
            os.chdir(Root+'/linear')
            os.system('mv data/dataNaNs.csv data/data.csv')
            os.system('plom build -t theta.json --local'  + "> /dev/null 2>&1")
            
      def setUp(self):
      # Things that need to be done before tests
            os.chdir(Root+'/linear/model')

      @classmethod
      def tearDownClass(cls):
            shutil.rmtree(Root + '/linear')

      def test_particle_genealogy(self):
            # We assess that the exploration of the genealogy to reconstruct sampled traj is correct by checking the trajectory
            # of the random walk through the distribution of its increments
            
            os.system('plom pipe theta.json |  ./kalman --no_dem_sto --traj -o 2 --quiet')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 100
            nbiters = 1
            os.system('plom pipe theta.json |  ./pmcmc sde --full --no_dem_sto -J ' + str(nparts) + ' -M ' + str(nbiters) + ' -i 1 -N 4 -n '  + str(nbiters) + ' --quiet'  + "> /dev/null 2>&1")          
            tab1 = genfromtxt('X_1.csv',delimiter=',',names=True)

            x = []
            for i in range(1,53):
                  # Compute and normalise increments of sampled path
                  x.append((math.log(float(tab1[52 - i][8]))-math.log(float(tab1[52 - i + 1][8])))/(math.sqrt(7)*0.1/7))          

            self.assertTrue(abs((stats.scoreatpercentile(x,50)-0)/(math.sqrt(0.5*0.5)/(stats.norm.pdf(0,loc=0,scale=1)*math.sqrt(len(x)))))<1.96)

            
      def test_correct_timing(self):
            # We check that timing has not been broken when playing with the indexes reconstructing the sampled path in pmcmc
            
            os.system('plom pipe theta.json |  ./kalman --no_dem_sto --traj -o 2 --quiet')
            tab0 = genfromtxt('hat_0.csv',delimiter=',',names=True)
            nparts = 2
            nbiters = 500
            os.system('plom pipe theta.json |  ./pmcmc sde --full --no_dem_sto -o 2 -J ' + str(nparts) + ' -M ' + str(nbiters) + ' -i 2 -N 4 --DT 0.0001 -n '  + str(nbiters) + ' --quiet'  + "> /dev/null 2>&1")          
            tab1 = genfromtxt('X_2.csv',delimiter=',',names=True)

            ### t = 1 (first obs)
            x = []
            for i in range(1,nbiters):
                  x.append(math.log(float(tab1[3*i - 2][8])))
                  
            meanSMC0r0 = numpy.mean(x)
            meanEKF0r0 = tab0[1][20]
            q975SMC0r0 = stats.scoreatpercentile(x,97.5)
            q975EKF0r0 = tab0[1][21]

            # Tests on 97.5% quantiles
            # based on on CLT for empirical quantile given in "Statistics and Data Analysis for Financial Engineering, Rupert 2011"
            self.assertTrue(abs((q975SMC0r0-q975EKF0r0)/(math.sqrt(0.975*0.025)/(stats.norm.pdf(q975EKF0r0-meanEKF0r0,loc=0,scale=(q975EKF0r0-meanEKF0r0)/1.96)*math.sqrt(nparts))))<1.96)

            ### t = 2 (second obs)
            x = []
            for i in range(1,nbiters):
                  x.append(math.log(float(tab1[3*i - 3][8])))

            meanSMC0r0 = numpy.mean(x)
            meanEKF0r0 = tab0[2][20]
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
