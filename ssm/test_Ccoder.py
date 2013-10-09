from Cmodel import Cmodel
import unittest
import copy
import json
import os

class TestCcoder(unittest.TestCase):

    def setUp(self):

        model = json.load(open(os.path.join('..' ,'noise', 'model', 'datapackage.json')))
        self.m_noise = Cmodel(model)

        model = json.load(open(os.path.join('..' ,'diff', 'model', 'datapackage.json')))
        self.m_diff = Cmodel(model)

        model = json.load(open(os.path.join('..' ,'noise', 'model', 'datapackage.json')))
        del model["resources"][1]["data"][2]["white_noise"] 
        del model["resources"][1]["data"][3]["white_noise"]
        model["resources"][1]["data"][0]["white_noise"] = {"id":"noise_SI", "sd": "sto"}
        model["resources"][1]["data"][1]["white_noise"] = {"id":"noise_SI", "sd": "sto"}
        self.m_noise2 = Cmodel(model)

        model = json.load(open(os.path.join('..' ,'noise', 'model', 'datapackage.json')))
        del model["resources"][1]["data"][2]["white_noise"] 
        del model["resources"][1]["data"][3]["white_noise"]
        model["resources"][1]["data"][4]["white_noise"] = {"id":"noise_SI", "sd": "sto"}
        model["resources"][1]["data"][5]["white_noise"] = {"id":"noise_SI", "sd": "sto"}
        self.m_noise3 = Cmodel(model)

        model = json.load(open(os.path.join('..' ,'noise', 'model', 'datapackage.json')))
        del model["resources"][1]["data"][2]["white_noise"] 
        del model["resources"][1]["data"][3]["white_noise"]
        model["resources"][1]["data"][8]["white_noise"] = {"id":"noise_SI", "sd": "sto"}
        model["resources"][1]["data"][9]["white_noise"] = {"id":"noise_SI", "sd": "sto"}
        self.m_noise4 = Cmodel(model)

        model = json.load(open(os.path.join('..' ,'noise', 'model', 'datapackage.json')))
        model["resources"][1]["data"][4]["white_noise"] = {"id":"noise_SI", "sd": "sto"}
        model["resources"][1]["data"][5]["white_noise"] = {"id":"noise_SI", "sd": "sto"}
        self.m_noise6 = Cmodel(model)

        model = json.load(open(os.path.join('..' ,'noise', 'model', 'datapackage.json')))
        model["resources"][1]["data"][4]["white_noise"] = {"id":"noise_SI2", "sd": "sto"}
        model["resources"][1]["data"][5]["white_noise"] = {"id":"noise_SI2", "sd": "sto"}
        self.m_noise7 = Cmodel(model)

        model = json.load(open(os.path.join('..' ,'diff', 'model', 'datapackage.json')))
        model["resources"][1]["data"].append({"from": "R_paris",   "to": "I_paris",   "rate": "correct_rate(v)",            "description":"testing"})
        model["resources"][1]["data"].append({"from": "R_nyc",   "to": "I_nyc",   "rate": "correct_rate(v)",                "description":"testing"})
        self.m_diff2 = Cmodel(model)

        
    def test_make_C_term(self):
        print 'test init'





if __name__ == '__main__':
    unittest.main()

