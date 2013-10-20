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
        self.special_functions = set(['terms_forcing', 'heaviside', 'ramp', 'slowstep', 'sin', 'cos', 'correct_rate', 'sqrt'])

        ###########################################################################

        self.remainder = sorted([x['remainder']['name'] for x in self.get_resource('populations') if 'remainder' in x])
        self.ur = ['U'] + self.remainder

        parameters = self.get_resource('parameters')
        sde = self.get_resource('sde')
        reactions = self.get_resource('reactions')
        observations = self.get_resource('observations')

        #par_forced (covariates)
        par_forced = [x['name'] for x in parameters if 'schema' in x]
        self.par_forced = sorted(par_forced)

        #par_sv and par_inc (incidence)
        par_sv = set()
        par_inc = set()

        for r in reactions:
            if r['from'] not in self.ur:
                par_sv.add(r['from'])
            if r['to'] not in self.ur:
                par_sv.add(r['to'])

            if "tracked" in r:
                for inc in r['tracked']:
                    par_inc.add(inc)

        self.par_sv = sorted(list(par_sv))
        self.par_inc = sorted(list(par_inc))

        #par proc and par_noise
        par_proc = set()
        par_noise = set()
        self.white_noise = []
        for r in reactions:
            el =  self.change_user_input(r['rate'])
            for e in el:
                if e not in self.op and e not in self.special_functions and e not in self.par_sv and e not in self.par_forced:
                    try:
                        float(e)
                    except ValueError:
                        par_proc.add(e)

            if 'white_noise' in r:
                par_noise.add(r['white_noise']['sd'])
                if r['white_noise']['name'] not in [y['name'] for y in self.white_noise]:
                    self.white_noise.append(r['white_noise'])

        self.par_noise = sorted(list(par_noise))
        self.par_proc = sorted(list(par_proc))


        #par_obs
        par_obs = set();
        priors = [x['name'] for x in parameters]
        for o in observations:
            for p in [o['pdf']['mean'], o['pdf']['sd']]:
                el =  self.change_user_input(p)
                for e in el:
                    if e not in self.op and e not in self.special_functions and e not in self.par_sv and e not in self.par_noise and e not in self.par_proc and e not in self.par_forced and e not in self.par_inc:
                        try:
                            float(e)
                        except ValueError:
                            par_obs.add(e)

        self.par_obs = sorted(list(par_obs))

        ##par_disp (parameter involve **only** in dispertion (nowhere else)
        disp = [x for subl in sde['dispersion'] for x in subl if x != 0] if 'dispersion' in sde else []
        par_disp = set()
        for x in disp:
            el =  self.change_user_input(x)
            for e in el:
                if e not in self.op and e not in self.special_functions and e not in self.par_sv and e not in self.par_proc and e not in self.par_obs and e not in self.par_noise and e not in self.par_forced:
                    try:
                        float(e)
                    except ValueError:
                        par_disp.add(e)

        self.par_disp = sorted(list(par_disp))

        #par_diff (state variable for diffusions)
        par_diff = []
        if sde:
            for x in sde.get('drift', []):
                par_diff.append(x['name'])

        self.par_diff = ['diff__' + x for x in sorted(par_diff)]

        ##par_other
        par_ssm = self.par_sv + self.par_inc + self.remainder + self.par_diff + self.par_noise + self.par_proc +  self.par_obs + self.par_forced + self.par_disp
        self.par_other = sorted([x['name'] for x in parameters if x['name'] not in par_ssm])

        ##all parameters
        self.all_par = par_ssm + self.par_other + ['t']

        ##orders in nav->states and nav->parameters
        ## !!par_sv must be first in both order_states and order_parameters, remainder must be last in order_states
        self.order_states = {x:i for i,x in enumerate(self.par_sv + self.par_inc + self.par_diff + self.remainder)}
        self.order_parameters = {x:i for i,x in enumerate(self.par_sv + self.par_noise + self.par_proc + self.par_disp + self.par_obs + self.par_other)}

        #map prior name to name
        self.map_prior_name2name = {}
        self.map_name2prior_name = {}
        for p in parameters:
            if 'prior' in p and 'name' in p['prior']:
                self.map_prior_name2name[p['prior']['name']] = p['name']
                self.map_name2prior_name[p['name']] = p['prior']['name']
            else:
                self.map_name2prior_name[p['name']] = p['name']


        # proc_model
        self.proc_model = copy.deepcopy(reactions)

        # obs_model
        self.obs_model = copy.deepcopy(observations)

        #fix rates:
        #replace ^ by ** for sympy
        # We treat reaction starting from remainder as reaction
        # starting from U that is rate -> rate * from size. It results
        # in simpler code in Ccoder.py. We also replace remainder by
        # N- sum(par_sv) in the rates (and in the rates ONLY)

        remainder_def = {}
        for x in self.get_resource('populations'):
            if 'remainder' in x:
                remainder_def[x['remainder']['name']] = '({0}-{1})'.format(x['remainder']['pop_size'], '-'.join([s for s in x['composition'] if s != x['remainder']['name']]))

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

    model = json.load(open(os.path.join('..' ,'examples', 'foo', 'datapackages', 'model-seb-sir', 'datapackage.json')))
    m = Cmodel(model)
