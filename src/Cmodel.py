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
import os
import os.path
import json
from sympy import diff, Symbol, sympify, simplify
from sympy.solvers import solve
from sympy.printing import ccode

class ModelError(Exception):
    def __init__(self, value):
        self.value = value
        def __str__(self):
            return repr(self.value)

class Cmodel:

    """
    parse a JSON model description
    """

    def __init__(self, dpkgRoot, dpkg, **kwargs):
        self.dpkgRoot = os.path.abspath(unicode(dpkgRoot, 'utf8'))
        self.model = copy.deepcopy(dpkg)

        self.op = set(['+', '-', '*', '/', ',', '(', ')']) ##!!!CAN'T contain square bracket '[' ']'
        self.reserved = set(['U', 'x', 't', 'E', 'LN2', 'LN10','LOG2E', 'LOG10E', 'PI', 'SQRT1_2', 'SQRT2']) #JS Math Global Object
        self.special_functions = set(['terms_forcing', 'heaviside', 'ramp', 'slowstep', 'sigmoid', 'sin', 'cos', 'correct_rate', 'ssm_correct_rate', 'sqrt', 'pow', 'exp', 'log'])

        self.remainder = sorted([x['remainder']['name'] for x in self.model['populations'] if 'remainder' in x])
        self.ur = ['U'] + self.remainder

        #resolve links for priors (named parameters here) for every
        #inputs, replace the require hash ({path: name: }) by it's
        #corresponding resource. Note that if (and only if) the
        #require hash has a name property, transformations have to be
        #done in terms of this name.
        for i, p in enumerate(self.model['inputs']):
            if 'require' in p:
                if 'fields' not in p['require']: ##only for priors (in S|S|M priors HAVE TO BE in JSON/JSON-LD) and covariates in csv
                    try:
                        rpath = os.path.join(self.dpkgRoot, p['require']['path'])
                        resource = json.load(open(rpath))
                    except:
                        raise ModelError('invalid data for ' + p['name'] + ' could not load ' + p['require']['path'])

                    if 'name' not in p['require']:
                        p['require']['name'] = p['name']

                    #semantic to SSM transfo
                    prior = { 'distribution': resource['name'] if resource['name'] != 'dirac' else 'fixed' }
                    for x in resource['distributionParameter']:
                        prior[x['name'] if 'name' in x else 'value'] = x['value']

                    self.model['inputs'][i]['data'] = prior


        parameters = self.model['inputs']
        sde = self.model.get('sde', {})

        reactions = self.model['reactions']
        observations = self.model['observations']

        #par_forced (covariates)
        par_forced = [x['name'] for x in parameters if 'require' in x and 'fields' in x['require']]
        self.par_forced = sorted(par_forced)

        #par_sv and par_inc (incidence)
        par_sv = set()
        par_inc = set()

        for r in reactions:
            if r['from'] not in self.ur:
                par_sv.add(r['from'])
            if r['to'] not in self.ur:
                par_sv.add(r['to'])

            if "accumulators" in r:
                for inc in r['accumulators']:
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
                if e not in self.op and e not in self.reserved and e not in self.special_functions and e not in self.par_sv and e not in self.par_forced:
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
        for o in observations:
            if o['distribution'] == 'poisson':
                o['sd'] = 'sqrt('+o['mean']+')'

            for p in [o['mean'], o['sd']]:
                el =  self.change_user_input(p)
                for e in el:
                    if e not in self.op and e not in self.reserved and e not in self.special_functions and e not in self.par_sv and e not in self.par_noise and e not in self.par_proc and e not in self.par_forced and e not in self.par_inc:
                        try:
                            float(e)
                        except ValueError:
                            par_obs.add(e)

        self.par_obs = sorted(list(par_obs))

        ##par_disp (parameter involve **only** in dispersion (nowhere else))
        disp = [x for subl in sde['dispersion'] for x in subl if x != 0] if 'dispersion' in sde else []
        par_disp = set()
        for x in disp:
            el =  self.change_user_input(x)
            for e in el:
                if e not in self.op and e not in self.reserved and e not in self.special_functions and e not in self.par_sv and e not in self.par_proc and e not in self.par_obs and e not in self.par_noise and e not in self.par_forced:
                    try:
                        float(e)
                    except ValueError:
                        par_disp.add(e)

        self.par_disp = sorted(list(par_disp))

        ##par_diff (state variable for diffusions)
        par_diff = []
        if sde:
            for x in sde.get('drift', []):
                par_diff.append(x['name'])

        ##NOTE: do not sort par_diff as the dispersion matrix is read in the unsorted order
        self.par_diff = ['diff__' + x for x in par_diff]

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
            if 'require' in p and 'fields' not in p['require'] and 'name' in p['require']:
                self.map_prior_name2name[p['require']['name']] = p['name']
                self.map_name2prior_name[p['name']] = p['require']['name']
            else:
                self.map_name2prior_name[p['name']] = p['name']

        # proc_model
        self.proc_model = copy.deepcopy(reactions)

        # obs_model
        self.obs_model = copy.deepcopy(observations)

        #fix rates:
        #replace pow by ** for sympy
        # We treat reaction starting from remainder as reaction
        # starting from U that is rate -> rate * from size. It results
        # in simpler code in Ccoder.py. We also replace remainder by
        # N- sum(par_sv) in the rates (and in the rates ONLY)

        remainder_def = {}
        for x in self.model['populations']:
            if 'remainder' in x:
                remainder_def[x['remainder']['name']] = '({0}-{1})'.format(x['remainder']['pop_size'], '-'.join([s for s in x['composition'] if s != x['remainder']['name']]))

        resolve_remainder = lambda x: remainder_def[x] if x in self.remainder else x

        for i, m in enumerate(self.proc_model):
            self.proc_model[i]['rate'] = self.pow2star(self.proc_model[i]['rate'])

            if self.proc_model[i]['from'] in self.remainder:
                self.proc_model[i]['rate'] = '({0})*{1}'.format(self.proc_model[i]['rate'], self.proc_model[i]['from'])

            self.proc_model[i]['rate'] = ''.join(map(resolve_remainder, self.change_user_input(m['rate'])))


        for i, m in enumerate(self.obs_model):
            for x in m:
                if x != "distribution" and x!= 'name' and x !='start':
                    self.obs_model[i][x] = self.pow2star(self.obs_model[i][x])
                    self.obs_model[i][x] = ''.join(map(resolve_remainder, self.change_user_input(self.obs_model[i][x])))

        ## incidence def
        self.par_inc_def = []
        for inc in self.par_inc:
            self.par_inc_def.append([x for x in self.proc_model if "accumulators" in x and inc in x['accumulators'] ])


    def change_user_input(self, term):
        """transform the term in smtg that we can parse in a programming language:
        example: change_user_input('r0*2*correct_rate(v)') -> ['r0', '*', '2', '*', 'correct_rate', '(', 'v', ')']"""

        myterm=term.replace(' ','') ##get rid of whitespaces
        mylist=[]
        mystring=''

        for i in range(len(myterm)):

            if myterm[i] in self.op :
                if len(mystring)>0:
                    mylist.append(mystring)
                    mystring=''
                mylist.append(myterm[i])
            else:
                mystring += myterm[i]

        if len(mystring)>0: ##the string doesn't end with an operator
            mylist.append(mystring)

        return mylist


    def pow2star(self, term):
        """replace pow(a,b) by (a)**(b) so that Sympy works"""


        terms = self.change_user_input(term)

        start = 0
        for i, x in enumerate(terms):
            if x == 'pow':
                start = i
                break

        if not start and not terms[0] == 'pow':
            return term

        pos = 1 #counter for open parenthesis
        ind = start+2 #skip first parenthesis
        lhs = ''
        rhs = ''
        left = True
        while (ind < len(terms)):
            if terms[ind] == '(':
                pos += 1
            if terms[ind] == ')':
                pos -= 1

            if pos == 1 and terms[ind] == ',':
                left = False
            else:
                if left:
                    lhs += terms[ind]
                else:
                    if not (pos == 0 and terms[ind] == ')'):
                        rhs += terms[ind]
                    else:
                        break;
            ind += 1

        return self.pow2star(''.join(terms).replace('pow({0},{1})'.format(lhs, rhs), '({0})**({1})'.format(lhs, rhs)))




    def toC(self, term, no_correct_rate, force_par=False, xify=None, human=False, set_t0=False):

        if term == xify:
            term = 'x'

        if term in self.map_prior_name2name:
            term = self.map_prior_name2name[term]

        if term == 'correct_rate':
            return '' if no_correct_rate else 'ssm_correct_rate'

        if human:
            return term

        if term in self.par_sv or term in self.par_inc:
            if force_par:
                return 'gsl_vector_get(par,ORDER_{0})'.format(term)
            else:
                return 'X[ORDER_{0}]'.format(term)

        elif term in self.par_forced:
            return 'gsl_spline_eval(calc->spline[ORDER_{0}],{1},calc->acc[ORDER_{0}])'.format(term, '0.0' if set_t0 else 't')

        elif term in self.par_proc or term in self.par_noise or term in self.par_disp or term in self.par_obs or term in self.par_other:
            if ('diff__' + term) in self.par_diff:
                return 'diffed[ORDER_diff__{0}]'.format(term)
            else:
                return 'gsl_vector_get(par,ORDER_{0})'.format(term)

        else: ##r is an operator or x
            return term


    def generator_C(self, term, no_correct_rate, force_par=False, xify=None, human=False, set_t0=False):
        """add extra terms (for C code) at the end of special functions (support nested special functions)"""

        terms = self.change_user_input(term)

        ind = 0
        Cterm = ''
        stack = []
        while (ind < len(terms)):

            if terms[ind] in self.special_functions:
                myf = terms[ind]
                Cterm += self.toC(myf, no_correct_rate, force_par=force_par, xify=xify, human=human, set_t0=set_t0) + '('
                ind += 2 #skip first parenthesis
                stack.append({"f": myf, "pos": 1}) #pos: counter for open parenthesis
            else:
                if stack:
                    if terms[ind] == '(':
                        stack[-1]['pos'] += 1
                        Cterm += '('
                    elif terms[ind] == ')':
                        stack[-1]['pos'] -= 1
                        if stack[-1]['pos'] == 0:
                            ftoclose = stack.pop()
                            ##add extra terms (no whitespace)
                            if not human:
                                if ftoclose['f'] == 'correct_rate' and not no_correct_rate:
                                    Cterm += ',dt)'
                                else:
                                    Cterm += ')'
                            else:
                                Cterm += ')'

                        else:
                            Cterm += ')'

                    else:
                        Cterm += self.toC(terms[ind], no_correct_rate, force_par=force_par, xify=xify, human=human, set_t0=set_t0)

                    ind +=1
                else:
                    Cterm += self.toC(terms[ind], no_correct_rate, force_par=force_par, xify=xify, human=human, set_t0=set_t0)
                    ind += 1


        return Cterm


    def make_C_term(self, term, no_correct_rate, derivate=None, inverse=None, human=False, force_par=False, xify=None, set_t0=False):

        """transform a term into its ssm C expression OR the ssm C
        expression of its derivate, differentiating against the
        derivate (if derivate not None) OR compute inverse function
        """

        #prefix all the state variable and parameters by ssm___ to
        #avoid namespace collision with Sympy as QCOSINE letters are
        #used by SymPy

        myterm = self.change_user_input(term)
        safe = ''

        for r in myterm:
            if r in self.all_par:
                safe += 'ssm___' + r
            elif inverse and r == inverse:
                safe += 'ssm___' + r
            else:
                safe += r

        if derivate:
	    sy = Symbol(str('ssm___' + derivate)) if derivate != 'x' else Symbol(derivate)
	    pterm = diff(sympify(safe), sy)
        elif inverse:
            if inverse in myterm:
                sy = Symbol(str('ssm___' + inverse))
                pterm = solve(sympify(safe), sy)
                if not pterm:
                    raise SsmError("can't find a solution to " + term + "=0 solving for " + inverse)
                elif len(pterm)!=1:
                    raise SsmError("no unique solution for " + term + "=0 solving for " + inverse)
                else:
                    pterm = pterm[0]
            else:
                pterm = sympify(safe)

        else:
            pterm = sympify(safe)

        #remove the ssm___ prefix
        #term = ccode(simplify(pterm)).replace('ssm___', '') ##NOTE simplify is just too slow to be used...
        term = ccode(pterm).replace('ssm___', '')

        #make the ssm C expression
        return self.generator_C(term, no_correct_rate, force_par=force_par, xify=xify, human=human, set_t0=set_t0)


if __name__=="__main__":

    dpkgRoot = os.path.join('..' ,'examples', 'foo')
    dpkg = json.load(open(os.path.join(dpkgRoot, 'ssm.json')))
    m = Cmodel(dpkgRoot, dpkg)
