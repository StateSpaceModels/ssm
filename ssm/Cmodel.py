##########################################################################
#    This file is part of ssm.
#
#    ssm is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ssm is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with ssm.  If not, see
#    <http://www.gnu.org/licenses/>.
#########################################################################

import copy
import sys


class Cmodel:

    """
    parse a JSON model description
    """


    def __init__(self, model,  **kwargs):
        self.model = model

        self.op = set(['+', '-', '*', '/', '^', ',', '(', ')']) ##!!!CAN'T contain square bracket '[' ']'
        self.reserved = set(['U', 'x'])
        self.special_functions = set(['terms_forcing', 'heaviside', 'ramp', 'slowstep', 'sin', 'cos', 'correct_rate'])

        ###########################################################################

        self.remainder = [x['remainder']['id'] for x in self.get_resource('populations') if 'remainder' in x]
        self.ur = ['U'] + self.remainder

        #par_fixed
        par_fixed = [x['id'] for x in self.get_resource('parameters') if 'prior' in x and 'path' in x['prior']]
        self.par_fixed = sorted(par_fixed)


        #par_sv and par_inc (incidence)
        par_sv = set()
        par_inc = set()


        for r in self.get_resource('reactions'):
            if r['from'] not in self.ur:
                par_sv.add(r['from'])
            if r['to'] not in self.ur:
                par_sv.add(r['to'])

            if "tracked" in r:
                for inc in r['tracked']:
                    par_inc.add(inc)

        self.par_sv = sorted(list(par_sv))
        self.par_inc = sorted(list(par_inc))

        #par proc
        par_proc = set()
        for r in self.get_resource('reactions'):
            el =  self.change_user_input(r['rate'])
            for e in el:
                if e not in self.op and e not in self.special_functions and e not in self.par_sv and e not in self.par_fixed:
                    try:
                        float(e)
                    except ValueError:
                        par_proc.add(e)

        self.par_proc = sorted(list(par_proc))


        #par_obs
        par_obs = set();
        priors = [x['id'] for x in self.get_resource('parameters')]
        for o in self.get_resource('observations'):
            for p in [o['pdf']['mean'], o['pdf']['var']]:
                el =  self.change_user_input(p)
                for e in el:
                    if e not in self.op and e not in self.special_functions and e not in self.par_sv and e not in self.par_proc and e not in self.par_fixed and e not in self.par_inc:
                        try:
                            float(e)
                        except ValueError:
                                par_obs.add(e)

        self.par_obs = sorted(list(par_obs))

        #par_diff (state variable for diffusions)
        par_diff = []

        for x in self.get_resource('sde')['skeleton']:
            par_diff.append(x['id'])

        self.par_diff = ['diff__' + x for x in sorted(par_diff)]

        #par_vol (volatilites)
        self.par_vol = sorted(set([x for subl in self.get_resource('sde')['sigma'] for x in subl if x != 0]))


        #map prior id to id
        self.map_prior_id2id = {}
        for p in self.get_resource('parameters'):
            if 'prior' in p and 'id' in p['prior']:
                self.map_prior_id2id[p['prior']['id']] = p['id']        


        # proc_model
        self.proc_model = copy.deepcopy(self.get_resource('reactions'))


        # obs_model
        self.obs_model = copy.deepcopy(self.get_resource('observations'))

        #fix rates:
        #replace ^ by ** for sympy
        # We treat reaction starting from remainder as reaction
        # starting from U that is rate -> rate * from size. It results
        # in simpler code in Ccoder.py. We also replace remainder by
        # N- sum(par_sv) in the rates (and in the rates ONLY)

        remainder_def = {}
        for x in self.get_resource('populations'):
            if 'remainder' in x:
                remainder_def[x['remainder']['id']] = '({0}-{1})'.format(x['remainder']['pop_size'], '-'.join([s for s in x['composition'] if s != x['remainder']['id']]))

        resolve_remainder = lambda x: remainder_def[x] if x in self.remainder else x

        for i, m in enumerate(self.proc_model):
            self.proc_model[i]['rate'] = self.proc_model[i]['rate'].replace('^', '**')

            if self.proc_model[i]['from'] in self.remainder:
                self.proc_model[i]['rate'] = '({0})*{1}'.format(self.proc_model[i]['rate'], self.proc_model[i]['from'])

            self.proc_model[i]['rate'] = ''.join(map(resolve_remainder, self.change_user_input(m['rate'])))


        for i, m in enumerate(self.obs_model):
            for x in m['pdf']:
                if x != "distribution":
                    self.obs_model[i]['pdf'][x] = self.obs_model[i]['pdf'][x].replace('^', '**')
                    self.obs_model[i]['pdf'][x] = ''.join(map(resolve_remainder, self.change_user_input(self.obs_model[i]['pdf'][x])))

        ## incidence def
        self.par_inc_def = []
        for inc in self.par_inc:
            self.par_inc_def.append([x for x in self.proc_model if "tracked" in x and inc in x['tracked'] ])



    def change_user_input(self, reaction):
        """transform the reaction in smtg that we can parse in a programming language:
        example: change_user_input('r0*2*correct_rate(v)') -> ['r0', '*', '2', '*', 'correct_rate', '(', 'v', ')']"""

        myreaction=reaction.replace(' ','') ##get rid of whitespaces
        mylist=[]
        mystring=''

        for i in range(len(myreaction)):

            if myreaction[i] in self.op :
                if len(mystring)>0:
                    mylist.append(mystring)
                    mystring=''
                mylist.append(myreaction[i])
            else:
                mystring += myreaction[i]

        if len(mystring)>0: ##the string doesn't end with an operator
            mylist.append(mystring)

        return mylist


    def get_resource(self, key):
        """get resource key from a datapackage.json"""

        r = [x for x in self.model['resources'] if x['name']==key]

        return r and r[0]['data']




if __name__=="__main__":

    import json
    import os

    model = json.load(open(os.path.join('..' ,'example', 'model', 'datapackage.json')))
    m = Cmodel(model)

