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
from Cmodel import Cmodel
from sympy import diff, Symbol, sympify, simplify
from sympy.solvers import solve
from sympy.printing import ccode
import copy

class SsmError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class Ccoder(Cmodel):
    """write the C code from the user input coming from the web interface..."""

    def __init__(self, model,  **kwargs):
        Cmodel.__init__(self, model,  **kwargs)


    def toC(self, term, no_correct_rate, force_par=False, xify=None, human=False, set_t0=False):

        if term == xify:
            term = 'x'

        if term in self.map_prior_id2id:
            term = self.map_prior_id2id[term]

        if term == 'correct_rate':
            return '' if no_correct_rate else 'ssm_correct_rate'

        if human:
            return term

        if term in self.par_sv or term in self.par_inc:
            if force_par:
                return 'gsl_vector_get(par, ORDER_{0})'.format(term)
            else:
                return 'X[ORDER_{0}]'.format(term)

        elif term in self.par_fixed:
            return 'gsl_spline_eval(calc->spline[ORDER_{0}],{1},calc->acc[ORDER_{0}])'.format(term, '0.0' if set_t0 else 't')

        elif term in self.par_proc or term in self.par_noise or term in self.par_obs or term in self.par_other:
            if ('diff__' + term) in self.par_diff:
                return 'diffed[ORDER_diff__{0}]'.format(term)
            else:
                return 'gsl_vector_get(par, ORDER_{0})'.format(term)

        else: ##r is an operator or x
            return term


    def generator_C(self, term, no_correct_rate, force_par=False, xify=None, human=False, set_t0=False):

        terms = self.change_user_input(term)

        ind = 0
        Cterm = ''
        while (ind < len(terms)):

            if terms[ind] in self.special_functions:
                myf = terms[ind]

                Cterm += self.toC(myf, no_correct_rate, force_par=force_par, xify=xify, human=human, set_t0=set_t0) + '('
                ind += 2 #skip first parenthesis

                pos = 1 #counter for open parenthesis
                while pos > 0:
                    if terms[ind] == '(':
                        pos += 1
                    if terms[ind] == ')':
                        pos -= 1

                    if pos >0:
                        Cterm += self.toC(terms[ind], no_correct_rate, force_par=force_par, xify=xify, human=human, set_t0=set_t0)
                        ind += 1

                ##add extra terms (no whitespace)
                if not human:
                    if myf == 'terms_forcing': ##TODO fix
                        Cterm += ',t,p_data,cac'
                    elif myf == 'correct_rate' and not no_correct_rate:
                        Cterm += ',dt'

                ##close bracket
                Cterm += terms[ind]
                ind += 1

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


    def parameters(self):
        """
        Everything needed to create ssm_parameter_t, ssm_state_t and load ssm_input_t
        """
        parameters = copy.deepcopy(self.get_resource('parameters'))

        for p in parameters:
            if 'transformation' in p:
                p['f_user2par'] = self.make_C_term(p['transformation'], True, force_par=True, xify=p['prior']['id'] if ('prior' in p and 'id' in p['prior']) else p['id'], set_t0=True)

                ## inverse of the transformation function
                ## if 'prior' in p and 'id' in p['prior']:
                ##     p['f_par2user'] = self.make_C_term(p['transformation']+ '-' + p['id'], True, inverse=p['prior']['id'], force_par=True, xify=p['id'], set_t0=True)

            if 'to_prior' in p:
                p['f_2prior'] = self.make_C_term(p['to_prior'], True)


        drifts = self.get_resource('sde')
        drifts = drifts and drifts['drift']
        #TODO support ode drifts += self.get_resource('ode')

        states = self.par_sv + self.par_inc
        pars = self.par_sv + self.par_noise + self.par_proc + self.par_obs + self.par_other

        #make C code for f_, f_inv f_der, f_der_inv
        for p in drifts:
            if 'transformation' in p:
                p['f'] = self.make_C_term(p['transformation'], True, human=True, xify=p['id'], set_t0=True)
                p['f_inv'] = self.make_C_term(p['transformation']+ '- x', True, inverse=p['id'], human=True, set_t0=True)
                p['f_der'] = self.make_C_term(p['transformation'], True, derivate=p['id'], human=True, xify=p['id'], set_t0=True)
                p['f_der_inv'] = self.make_C_term(p['f_inv'], True, derivate='x', human=True, set_t0=True)

            if p['id'] in self.order_parameters:
                p['offset_ic'] = self.order_parameters[p['id']]
                

        #sort parameters
        #start by making dict:
        pdict = {x['id']:x for x in parameters}
        sdict = {'diff__' + x['id']: x for x in drifts}

        f_remainders = {}
        f_remainders_par = {}
        f_remainders_var = {}
        for x in self.get_resource('populations'):
            if 'remainder' in x:
                rem = x['remainder']['id']
                eq = x['remainder']['pop_size'] + ' - ' + ' - '.join([r for r in x['composition'] if r != rem])
                f_remainders[rem] = self.make_C_term(eq, True)
                f_remainders_par[rem] = self.make_C_term(eq, True, force_par=True, set_t0=True)
                eq = ''
                for x_i in x['composition']:
                    for x_j in x['composition']:
                        if eq != '':
                            eq += ' + '
                        eq += 'gsl_matrix_get(&Ct.matrix,' + str(self.order_states[x_i]) +','  + str(self.order_states[x_j]) + ')';
                f_remainders_var[rem] = eq;

        # Initial compartment sizes in cases of no remainder
        ic = []
        for x in self.get_resource('populations'):
            if 'remainder' not in x:
                ic.append([self.make_C_term(eq, True, force_par=True, set_t0=True) for t in x['composition']])


        return {
            'parameters': parameters,
            'order_parameters': self.order_parameters,
            'order_states': self.order_states,
            'drifts': drifts,
            'par_sv': self.par_sv,
            'states': states,
            'remainders': self.remainder,
            'f_remainders': f_remainders,
            'f_remainders_par': f_remainders_par,
            'f_remainders_var': f_remainders_var,
            'ic': ic,
            'sde': [sdict[x] for x in self.par_diff],
            'pars': [pdict[x] for x in pars]
        }


    def observed(self):
        ##WARNING right now only the discretized normal is supported.
        ##TODO: generalization for different distribution

        obs = copy.deepcopy(self.obs_model)

        for x in obs:
            x['pdf']['mean'] = self.make_C_term(x['pdf']['sd'], True)
            x['pdf']['sd'] = self.make_C_term(x['pdf']['sd'], True)

        return {'observed': obs}


    def iterators(self):

        return {
            'state': {
                'sv': [self.order_states[x] for x in self.par_sv],
                'remainders': [self.order_states[x] for x in self.remainder],
                'inc': [self.order_states[x] for x in self.par_inc],
                'sv_inc': [self.order_states[x] for x in (self.par_sv + self.par_inc)],
                'diff': [self.order_states[x] for x in self.par_diff]                
            },
            'parameter': {
                'all': [self.order_parameters[x] for x in (self.par_sv + self.par_noise + self.par_proc + self.par_obs + self.par_other)],
                'noise': [self.order_parameters[x] for x in self.par_noise],
                'icsv': [self.order_parameters[x] for x in self.par_sv],
                'icdiff': [self.order_parameters[x.split('diff__')[1]] for x in self.par_diff]
            }
        }



    def orders(self):
        """
        #define and #undef
        """
        univ_rem = ['U']
        if self.remainder:
            univ_rem += self.remainder

        return {
            'var': [{'name': x, 'order': self.order_parameters[x]} for x in (self.par_sv + self.par_noise + self.par_proc + self.par_obs + self.par_other)],
            'diff': [{'name': x, 'order': o} for o,x in enumerate(self.par_diff) ],
            'inc': [{'name': x, 'order': self.order_states[x]} for x in self.par_inc],
            'covariates': [{'name': x, 'order': o} for o,x in enumerate(self.par_fixed)] ,
            'univ_rem': [{'name': x, 'order': len(self.par_sv)+o} for o,x in enumerate(univ_rem) ]
        }




    def cache_special_function_C(self, caches_C, sf=None, prefix='_sf'):
        """caches_C: List of cached expression in C
        caches_C is modified in place
        sf: an optional list of unique special function to be cached
        returns sf (created if sf input is None)
        """

        if not sf:
            sf = []
            for term in caches_C:
                if any([x in term for x in self.special_functions]):
                    terms = self.change_user_input(term)
                    ind = 0
                    while (ind < len(terms)):
                        if terms[ind] in self.special_functions:
                            f = terms[ind] + '('
                            ind += 2 #skip first parenthesis
                            pos = 1 #counter for open parenthesis
                            while pos > 0:
                                if terms[ind] == '(':
                                    pos += 1
                                if terms[ind] == ')':
                                    pos -= 1

                                f += terms[ind]
                                ind +=1

                            sf.append(f)
                        else:
                            ind += 1

            sf = list(set(sf))

        for i, term in enumerate(caches_C):
            if any([x in term for x in self.special_functions]):
                for s in sf:
                    caches_C[i] = caches_C[i].replace(s, prefix + '[{0}]'.format(sf.index(s)))

        return sf



    def alloc_psr(self):
        Clist = []
        univ = ['U']
        if self.remainder:
            univ += self.remainder

        for s in self.par_sv + univ:
            nbreac = len([r for r in self.proc_model if r['from']==s]) +1 ##+1 to stay in the same compartment or to declare smtg in case of no reaction (not super clean but makes C code easier...)
            Clist.append({'state':s, 'nb_reaction': nbreac})

        return Clist

    def step_psr(self):

        """
        prob and update for Poisson with stochastic rate step function

        prob general structure:

        sum=...;
        if(sum>0.0){
            prob[0]=(1-exp(-sum))*(rate/sum);
            ...
            prob[last]=1-sum(prob);
        }
        else{
            prob[0]=0;
            ...
            prob[last]=1;
        }
        we need the if statement to avoid division by 0
        """

        ###########
        ## prob  ##
        ###########
        proc_model = copy.deepcopy(self.proc_model) ##we are going to modify it...

        ##make the rates noisy (if needed e.g r
        rates = set()
        for r in proc_model:
            if r['from'] not in (['U'] + self.remainder):
                myrate = r['rate']
                if 'white_noise' in r:
                    myrate = '({0})*{1}'.format(myrate, r['white_noise']['id'])

                rates.add(myrate)

        rates = list(rates)
        caches = map(lambda x: self.make_C_term(x, False), rates)
        sf = self.cache_special_function_C(caches)

        for r in proc_model:
            if r['from'] not in (['U'] + self.remainder):
                myrate = r['rate']
                if 'white_noise' in r:
                    myrate = '({0})*{1}'.format(myrate, r['white_noise']['id'])

                r['ind_cache'] = rates.index(myrate)

        Ccode=''

        for s in self.par_sv:
            myexit = [r for r in proc_model if r['from'] == s]
            exitlist=[]

            if len(myexit)>0:

                for e in myexit:
                    exitlist.append('_r[{0}]*dt'.format(e['ind_cache']))

                Csum= 'sum = ' + '+'.join(exitlist) + ';\n'
                Ccode += Csum+ 'if(sum>0.0){\none_minus_exp_sum = (1.0-exp(-sum));\n'
                Cprob=''
                sumprob='1.0'
                for reacnb in range(len(exitlist)):
                    Cprob += 'calc->prob[ORDER_{0}][{1}] = one_minus_exp_sum*(({2})/sum);\n'.format(s, reacnb, exitlist[reacnb])
                    sumprob += ' - calc->prob[ORDER_{0}][{1}]'.format(s, reacnb)

                Cprob += 'calc->prob[ORDER_{0}][{1}] = '.format(s,len(exitlist)) + sumprob + ';\n'
                Ccode += Cprob+ '}\n'
                Ccode +='else{\n'

                Celse=''
                for reacnb in range(len(exitlist)):
                    Celse += 'calc->prob[ORDER_{0}][{1}] = 0.0;\n'.format(s, reacnb)

                Celse += 'calc->prob[ORDER_{0}][{1}] = 1.0;\n'.format(s,len(exitlist))+'}\n\n'

                Ccode += Celse

        ############
        ## update ##
        ############

        incDict = dict([(x,'') for x in self.par_sv])

        for s in self.par_sv: ##stay in the same compartment
            myexit = [r for r in self.proc_model if r['from'] == s]
            if len(myexit)>0: ##only if you can exit from this compartment in this case the remaining has a sense
                incDict[s] += 'calc->inc[ORDER_{0}][{1}]'.format(s, len(myexit))
            else:
                incDict[s] += 'X[ORDER_{0}]'.format(s)

        for s in self.par_sv: #come in from other compartments
            myinput = [r for r in self.proc_model if r['from'] == s]
            for nbreac in range(len(myinput)):
                if myinput[nbreac]['to'] not in (['U'] + self.remainder): ##we exclude deaths or transitions to remainder in the update
                    incDict[myinput[nbreac]['to']] += ' + calc->inc[ORDER_{0}][{1}]'.format(myinput[nbreac]['from'], nbreac)


        ##we add flow from (['U'] + self.remainder) (Poisson term). We want to cache those flow so that the incidences can be computed
        poisson = []
        for s in (['U'] + self.remainder):
            reac_from_univ = [r for r in self.proc_model if (r['from'] == s and (r['to'] not in (['U'] + self.remainder)) )]
            for nbreac in range(len(reac_from_univ)):
                myrate = self.make_C_term(reac_from_univ[nbreac]['rate'], False)
                if 'white_noise' in reac_from_univ[nbreac]:
                    myrate = '({0})*{1}'.format(myrate, reac_from_univ[nbreac]['white_noise']['id'])

                poisson.append('calc->inc[ORDER_{0}][{1}] = gsl_ran_poisson(calc->randgsl, ({2})*dt)'.format(s, nbreac, myrate))
                incDict[reac_from_univ[nbreac]['to']] += ' + calc->inc[ORDER_{0}][{1}]'.format(s, nbreac)

        Cstring=''
        for s in self.par_sv:
            Cstring += 'X[ORDER_{0}] = {1};\n'.format(s, incDict[s])


        return {'code': Ccode, 'caches': caches, 'sf': sf, 'poisson': poisson, 'update_code': Cstring}


    def step_psr_inc(self):
        """generate C code to compute the dynamic of the observed
        **incidence** in case of stochastic models (euler multinomial)
        and put in into

        Clist = [{'right_hand_side': }]

        """
        Clist = []

        for i in range(len(self.par_inc_def)):
            right_hand_side=''

            for j in range(len(self.par_inc_def[i])):
                id_out = [self.proc_model.index(r) for r in self.proc_model if ((r['from'] == self.par_inc_def[i][j]['from']) and (r['to'] == self.par_inc_def[i][j]['to']) and (r['rate'] == self.par_inc_def[i][j]['rate']))]
                for o in id_out:
                    myexit = [r for r in self.proc_model if r['from']==self.proc_model[o]['from']]
                    right_hand_side += ' + calc->inc[ORDER_{0}][{1}]'.format(self.par_inc_def[i][j]['from'], myexit.index(self.proc_model[o]))

            Clist.append({'index': i, 'right_hand_side':right_hand_side})

        return Clist


    def step_psr_multinomial(self):
        draw = []
        for s in self.par_sv:
            nbexit = len([r for r in self.proc_model if r['from']==s])
            if nbexit>0:
                draw.append({'state': s, 'nb_exit': nbexit+1}) ##+1 to stay in the compartment

        return draw



    def step_ode_sde(self):
        """
        Generates ODE and SDEs
        note: sf are used in self.jac() for Lyapunov exp computations
        """

        proc_model = copy.deepcopy(self.proc_model) ##we are going to modify it...

        odeDict = dict([(x, []) for x in self.par_sv])

        rates = list(set(r['rate'] for r in proc_model))

        caches = map(lambda x: self.make_C_term(x, True), rates)
        sf = self.cache_special_function_C(caches, prefix='_sf')

        for i, r in enumerate(proc_model):
            r['ind_cache'] = rates.index(r['rate'])
            r['ind_dem_sto'] = i


        def get_rhs_term(sign, cached, reaction):
            if 'white_noise' in reaction:
                noise_id = reaction['white_noise']['id']
                noise_sd = self.toC(reaction['white_noise']['sd'], False)
            else:
                noise_id = None
                noise_sd= None

            return {'sign': sign, 'term': cached, 'noise_id': noise_id, 'noise_sd': noise_sd, 'ind_dem_sto': reaction['ind_dem_sto']}


        ################################
        ##Dynamic of the state variables
        ################################

        ##outputs
        for r in proc_model:
            if r['from'] not in (['U'] + self.remainder):
                cached = '_r[{0}]*X[ORDER_{1}]'.format(r['ind_cache'], r['from'])
                odeDict[r['from']].append(get_rhs_term('-', cached, r))

        ##inputs
        for r in proc_model:
            if r['to'] not in (['U'] + self.remainder):
                if r['from'] not in (['U'] + self.remainder):
                    cached = '_r[{0}]*X[ORDER_{1}]'.format(r['ind_cache'], r['from'])
                else:
                    cached= '_r[{0}]'.format(r['ind_cache'])

                odeDict[r['to']].append(get_rhs_term('+', cached, r))


        #######################################
        ##Dynamic of the observed **incidence**
        #######################################

        obs_list = []

        for i in range(len(self.par_inc_def)):
            eq = []

            if isinstance(self.par_inc_def[i][0], dict): ##incidence
                for j in range(len(self.par_inc_def[i])):
                    id_out = [proc_model.index(r) for r in proc_model if ((r['from'] == self.par_inc_def[i][j]['from']) and (r['to'] == self.par_inc_def[i][j]['to']) and (r['rate'] == self.par_inc_def[i][j]['rate'])) ]
                    for o in id_out:
                        reaction = proc_model[o]
                        if self.par_inc_def[i][j]['from'] in (['U'] + self.remainder):
                            cached = '_r[{0}]'.format(reaction['ind_cache'])
                        else:
                            cached = '_r[{0}]*X[ORDER_{1}]'.format(reaction['ind_cache'], self.par_inc_def[i][j]['from'])

                        eq.append(get_rhs_term('+', cached, reaction))

                obs_list.append({'index':i, 'eq': eq})


        ##############################################################################################################
        ##we create the ODE and  4 versions of the SDE system (no_dem_sto, no_white_noise, no_dem_sto_no_white_noise and full)
        ##############################################################################################################
        unique_noises_id = [x['id'] for x in self.white_noise]
        dem_sto_id = ['dem_sto__' +str(i) for i, x in enumerate(self.proc_model)]

        def eq_dem_env(eq_list):
            eq = ''  #deter skeleton
            dem = '' #demographic stochasticity
            env = '' #env stochasticity

            for x in eq_list:
                eq += ' {0} ({1})'.format(x['sign'], x['term'])

                #dem sto
                dem += '{0} sqrt(({1}))*dem_sto__{2}'.format(x['sign'], x['term'], x['ind_dem_sto'])

                #env sto
                if x['noise_id']:
                    env += '{0} ({1})*{2}*{3}'.format(x['sign'], x['term'], x['noise_sd'], x['noise_id'])

            return (eq, dem, env)


        func = {'no_dem_sto': {'proc': {'system':[], 'noises': unique_noises_id},
                               'obs': []},
                'no_white_noise': {'proc': {'system':[], 'noises': dem_sto_id},
                               'obs': []},
                'full': {'proc': {'system':[], 'noises': dem_sto_id + unique_noises_id},
                         'obs': []},
                'no_dem_sto_no_white_noise': {'proc':{'system':[], 'noises':[]},
                                          'obs':[]},
                'ode': {'proc':{'system':[], 'noises':[]},
                        'obs':[]}}


        #state variables
        for i, s in enumerate(self.par_sv):

            eq, dem, env = eq_dem_env(odeDict[s])
            if env:
                env = '+ ' + env

            #TODO get rid of the 'dt' for Euler Maruyama (should be handled on the C side as it is the case for sqrt(dt) for the stochastic part)'
            func['ode']['proc']['system'].append({'index': i, 'eq': eq})
            func['no_dem_sto_no_white_noise']['proc']['system'].append({'index': i, 'eq': '({0})*dt'.format(eq)})
            func['no_dem_sto']['proc']['system'].append({'index': i, 'eq': '({0})*dt {1}'.format(eq, env)})
            func['no_white_noise']['proc']['system'].append({'index': i, 'eq': '({0})*dt + {1}'.format(eq, dem)})
            func['full']['proc']['system'].append({'index': i, 'eq': '({0})*dt + {1} {2}'.format(eq, dem, env)})

        #observed incidence
        for myobs in obs_list:

            eq, dem, env = eq_dem_env(myobs['eq'])
            if env:
                env = ' + ' + env

            #TODO get rid of the 'dt' for Euler Maruyama (should be handled on the C side as it is the case for sqrt(dt) for the stochastic part)'
            func['ode']['obs'].append({'index': myobs['index'], 'eq': eq})
            func['no_dem_sto_no_white_noise']['obs'].append({'index': myobs['index'], 'eq': '({0})*dt'.format(eq)})
            func['no_dem_sto']['obs'].append({'index': myobs['index'], 'eq': '({0})*dt {1}'.format(eq, env)})
            func['no_white_noise']['obs'].append({'index': myobs['index'], 'eq': '({0})*dt {1}'.format(eq, dem)})
            func['full']['obs'].append({'index': myobs['index'], 'eq': '({0})*dt + {1} {2}'.format(eq, dem, env)})


        return {'func': func, 'caches': caches, 'sf': sf}


    def compute_diff(self):

        sde = self.get_resource('sde')
        if sde and 'dispersion' in sde:
            dispersion = sde['dispersion']
            diff.terms = []
            diff.n_browns = len(dispersion[0])
            for x in dispersion:
                term = ''
                for i, y in enumerate(x):
                    if y:
                        term += (' + ' if term else '') + self.make_C_term(y, True) + ' * _w[{0}]'.format(i)
                diff.terms.append(term)
            return diff

        else:
            return []





    def jac(self, sf_jac_only):
        """compute jacobian matrix of the process model (including
        observed variable) using Sympy


        sf_jac_only: list of cached special function generated by
        self.print_ode() used to get the index of caches_C for the
        jacobian matrix of simulation methods

        """

        my_model = copy.deepcopy(self.proc_model)
        odeDict = dict([(x,'') for x in self.par_sv])

        ##############################
        ###   Build odeDict
        ##############################

        ##outputs
        for r in my_model:
            if r['from'] not in (['U'] + self.remainder):
                rate= ' - (({0})*{1})'.format(r['rate'], r['from'])
                odeDict[r['from']] += rate

        ##inputs
        for r in my_model:
            if r['to'] not in (['U'] + self.remainder):
                if r['from'] not in (['U'] + self.remainder):
                    rate= ' + (({0})*{1})'.format(r['rate'], r['from'])
                    odeDict[r['to']] += rate
                else:
                    rate= ' + ({0})'.format(r['rate'])
                    odeDict[r['to']] += rate

        ##observed equations
        obsList = []


        for i in range(len(self.par_inc_def)):
            eq = ''
            for j in range(len(self.par_inc_def[i])):
                reaction = self.par_inc_def[i][j]
                if reaction['from'] in (['U'] + self.remainder):
                    eq += ' + ({0})'.format(reaction['rate'])
                else:
                    eq += ' + (({0})*{1})'.format(reaction['rate'], reaction['from'])

            obsList.append(eq)

        ####################
        ### Jacobian
        ####################

        ##derive process model equations (odeDict) per par_sv
        caches = []
        caches_jac_only = []

        jac = []
        jac_only = []
        jac_diff = []


        for s in range(len(self.par_sv)):
            jac.append([])
            jac_only.append([])

            if self.par_diff:
                jac_diff.append([])

            for sy in self.par_sv:
                Cterm = self.make_C_term(odeDict[self.par_sv[s]], True, derivate=sy)
                jac[s].append(Cterm)
                jac_only[s].append(Cterm)
                caches.append(Cterm)
                caches_jac_only.append(Cterm)

            #see doc of kalman.c diff_derivative()
            for sy in self.par_diff:
                Cterm = self.make_C_term(odeDict[self.par_sv[s]], True, derivate=sy.split('diff__')[1])
                jac_diff[s].append({'value': Cterm,
                                    'der': self.make_C_term(sy, True),
                                    'name': sy,
                                    'order': self.order_states[sy]})
                caches.append(Cterm)

        ##derive observation equations (obsList) per par_sv
        jac_obs = []
        jac_obs_diff = []

        for o in range(len(obsList)):
            jac_obs.append([])
            if self.par_diff:
                jac_obs_diff.append([])

            for sy in self.par_sv:
                Cterm = self.make_C_term(obsList[o], True, derivate=sy)
                jac_obs[o].append(Cterm)
                caches.append(Cterm)

            #see doc of kalman.c diff_derivative()
            for sy in self.par_diff:
                Cterm = self.make_C_term(obsList[o], True, derivate=sy)
                jac_obs_diff[o].append({'value': Cterm,
                                        'der': self.make_C_term(sy, True),
                                        'name': sy,
                                        'order': self.order_states[sy]})
                caches.append(Cterm)


        ##cache rates and remove duplicates
        caches = list(set(caches))
        caches_jac_only = list(set(caches_jac_only))

        ##replace with index of caches (will be _r[index] in C)
        for s in range(len(self.par_sv)):
            for i in range(len(self.par_sv)):
                Cterm = jac[s][i]
                jac[s][i] = caches.index(Cterm)
                jac_only[s][i] = caches_jac_only.index(Cterm)

            for i in range(len(self.par_diff)):
                jac_diff[s][i]['value'] = caches.index(jac_diff[s][i]['value'])


        for o in range(len(obsList)):
            for i in range(len(self.par_sv)):
                jac_obs[o][i] = caches.index(jac_obs[o][i])

            for i in range(len(self.par_diff)):
                jac_obs_diff[o][i]['value'] = caches.index(jac_obs_diff[o][i]['value'])


        ##special function that have to be cached (caches is transformed by self.cache_special_function_)
        sf = self.cache_special_function_C(caches, prefix='_sf')
        ##for jac_only (used for Lyapunov exp computations only, sf is shared with the one of print_ode. We just update caches_jac_only)
        self.cache_special_function_C(caches_jac_only, sf=sf_jac_only, prefix='_sf')

        return {'jac_only': jac_only,
                'jac': jac,
                'jac_obs': jac_obs,
                'jac_diff': jac_diff,
                'jac_obs_diff': jac_obs_diff,
                'caches': caches,
                'sf': sf,
                'caches_jac_only': caches_jac_only}


    def Ht(self):
        """compute jacobian matrix of the mean of the obs process (assumed to be Gaussian) using Sympy"""

        proc_model = copy.deepcopy(self.proc_model) ##we are going to modify it...
        obs = copy.deepcopy(self.obs_model)
        N_REAC = len(proc_model)
        N_PAR_SV = len(self.par_sv)
        N_PAR_INC = len(self.par_inc)
        N_DIFF = len(self.par_diff)

        Ht_sv = []
        Ht_inc = []
        Ht_diff = []

        ## Derivatives of observed means against state variables
        for s in range(len(self.par_sv)):
            Ht_sv.append([])
            for x in obs:
                Cterm = self.make_C_term(x['pdf']['mean'], True, derivate=self.par_sv[s])
                Ht_sv[s].append(Cterm)

        ## Derivatives of observed means against incidence variables
        for s in range(len(self.par_inc)):
            Ht_inc.append([])
            for x in obs:
                Cterm = self.make_C_term(x['pdf']['mean'], True, derivate=self.par_inc[s])
                Ht_inc[s].append(Cterm)

        ## Derivatives of observed means against diffusing variables
        for s in range(len(self.par_diff)):
            Ht_diff.append([])
            for x in obs:
                Cterm = self.make_C_term(x['pdf']['mean'], True, derivate=self.par_diff[s])
                Ht_diff[s].append(Cterm)

        return {'Ht_sv': Ht_sv,
                'Ht_inc': Ht_inc,
                'Ht_diff': Ht_diff}


    def h_grads(self):
        """compute the gradients of the observation functions using Sympy in order to compute the prediction variance through first-order Taylor expansions"""
        obs = copy.deepcopy(self.obs_model)
        h_grads = {}

        for x in obs:
            term = {}
            term['id'] = x['id']
            term['grads'] = []
            for s in (self.par_sv + self.par_inc + self.par_diff):
                Cterm = self.make_C_term(x['pdf']['mean'], True, derivate=s if 'diff__' not in s else s.split('diff__')[1])
                if Cterm!='0':
                    grad = {}
                    grad['Cterm'] = Cterm
                    grad['ind'] = self.order_states[s]
                    term['grads'].append(grad)

            h_grads[x['id']]=term


        return {'h_grads': h_grads}



    def eval_Q(self, debug = False):
        """

        The construction of Qsv is based on three levels:
         - states: state variables and observations (s)
         - reactions (r)
         - noise terms (n)

        At the reaction level, Qr is a two-blocks diagonal matrix: Qr_dem and Qr_env.
        Qr_dem corresponds to demographic noise and has reaction rates on the diagonal.
        Qr_env corresponds to white noises. It is built from Qn through Lr.
        Qn is a diagonal matrix which diagonal terms correspond to squarred amplitude of white noises.
        The stoechiometric matrices L are used to switch from one level to another:
              Qr_env = Lr Qn Lr'  and Qs = Ls Qr Ls'

        In particular, Lr has reaction rates in term (i,j) if reaction i is concerned by white noise j.
        Ls has +1 or -1 in term (i,j) if reaction j goes to or leaves from state i, and O's everywhere else.

        Note: we assume only one environmental noise term per reaction


        """
        proc_model = copy.deepcopy(self.proc_model) ##we are going to modify it...
        N_REAC = len(proc_model)
        N_PAR_SV = len(self.par_sv)
        N_PAR_INC = len(self.par_inc)
        N_DIFF = len(self.par_diff)

        unique_noises_names = [x['id'] for x in self.white_noise]
        N_ENV_STO_UNIQUE = len(unique_noises_names)

        ##add sd and order properties to noisy reactions
        N_ENV_STO = 0
        for reaction in proc_model:
            if 'white_noise' in reaction:
                reaction['order_env_sto_unique'] = unique_noises_names.index(reaction['white_noise']['id'])
                reaction['order_env_sto'] = N_ENV_STO
                N_ENV_STO += 1


        s = N_REAC + N_ENV_STO ## number of noise terms (potentially non-independent)
        ##for demographic stochasticity, one independent noise term per reaction

        Ls = [[0]*s for x in range(N_PAR_SV + N_PAR_INC)]
        Qs = [[0]*(N_PAR_SV + N_PAR_INC) for x in range(N_PAR_SV + N_PAR_INC)]
        Qr = [[0]*s for x in range(s)]
        Qr_dem = [[0]*s for x in range(N_REAC)]
        Qr_sto = [[0]*s for x in range(N_ENV_STO)]
        Lr = [[0]*N_ENV_STO_UNIQUE for x in range(N_ENV_STO)]
        Qn = [[0]*N_ENV_STO_UNIQUE for x in range(N_ENV_STO_UNIQUE)]


        ###########################################
        #    Create Ls and Qr_dem                 #
        ###########################################

        #state variables
        for B_dem_ind, r in enumerate(proc_model):
            is_noise = 'white_noise' in r
            if is_noise:
                B_sto_ind = N_REAC + r['order_env_sto']

            if r['from'] not in (['U'] + self.remainder):
                i = self.par_sv.index(r['from'])
                Ls[i][B_dem_ind] -= 1 ##demographic stochasticity
                if is_noise:
                    Ls[i][B_sto_ind] -= 1 ##env stochasticity

                Qc_term = '({0})*{1}'.format(r['rate'], r['from'])
            else:
                Qc_term = r['rate']

            if r['to'] not in (['U'] + self.remainder):
                i = self.par_sv.index(r['to'])
                Ls[i][B_dem_ind] += 1
                if is_noise:
                    Ls[i][B_sto_ind] += 1

            Qr_dem[B_dem_ind][B_dem_ind] =  Qc_term

        # incidence variables
        for i in range(len(self.par_inc_def)): #(for every incidence variable)
            for B_dem_ind, r in enumerate(proc_model):
                # for every incidence
                for inc in self.par_inc_def[i]:
                    # if it involves incidence
                    if (r['from'] == inc['from']) and (r['to'] == inc['to']) and (r['rate'] == inc['rate']):
                        Ls[N_PAR_SV + i][B_dem_ind] += 1
                        if 'white_noise' in r:
                            B_sto_ind = N_REAC + r['order_env_sto']
                            Ls[N_PAR_SV + i][B_sto_ind] += 1


        ############################
        ## Create Qr_env = Lr Qn Lr'
        ############################
        for r in proc_model:
            if 'white_noise' in r:
                if r['from'] not in (['U'] + self.remainder):
                    Qn_term = '({0})*{1}'.format(r['rate'], r['from'])
                else:
                    Qn_term = r['rate']

                Lr[r['order_env_sto']][r['order_env_sto_unique']] = Qn_term
                Qn[r['order_env_sto_unique']][r['order_env_sto_unique']] = '({0})**2'.format(r['white_noise']['sd'])



        def matrix_product(A, B):
            if not A or not B:
                return []

            res = [[0]*len(B[0]) for x in range(len(A))]

            for i in range(len(A)):
                for j in range(len(B[0])):
                    for k in range(len(B)):
                        if (A[i][k] and B[k][j]):
                            term = ('({0})*({1})').format(A[i][k], B[k][j])

                            if res[i][j]:
                               res[i][j] = res[i][j] + ' + {0}'.format(term)
                            else:
                               res[i][j] = term


            return res


        Qr_env = matrix_product(Lr, Qn)
        Qr_env = matrix_product(Qr_env, zip(*Lr))

        for i in range(N_ENV_STO):
            for j in range(N_ENV_STO):
                Qr[N_REAC+i][N_REAC+j] = Qr_env[i][j]

        #we fill Qr with Qc_dem and Qc_env
        for i in range(N_REAC):
            for j in range(N_REAC):
                Qr[i][j] = Qr_dem[i][j]


        #we split Ls into Ls_dem and Ls_env
        Ls_dem = [[0]*N_REAC for x in range(N_PAR_SV + N_PAR_INC)]
        for i in range(N_PAR_SV + N_PAR_INC):
            for j in range(N_REAC):
                Ls_dem[i][j] = Ls[i][j]

        Ls_env = [[0]*N_ENV_STO for x in range(N_PAR_SV + N_PAR_INC)]
        for i in range(N_PAR_SV + N_PAR_INC):
            for j in range(N_ENV_STO):
                Ls_env[i][j] = Ls[i][N_REAC + j]

        ############################
        ## Create Q_sde
        ############################

        sde = self.get_resource('sde')
        if sde and 'dispersion' in sde:
            dispersion = sde['dispersion']
            # Q_sde = dispersion * dispersion'
            Q_sde = matrix_product(dispersion, zip(*dispersion))


        #####################################################################################
        ##we create 4 versions of Q (no_dem_sto, no_env_sto, no_dem_sto_no_env_sto and full)
        #####################################################################################

        Qs = matrix_product(Ls, Qr)
        Qs = matrix_product(Qs, zip(*Ls))

        Qs_dem = matrix_product(Ls_dem, Qr_dem)
        Qs_dem = matrix_product(Qs_dem, zip(*Ls_dem))

        Qs_env = matrix_product(Ls_env, Qr_env)
        Qs_env = matrix_product(Qs_env, zip(*Ls_env))

        calc_Q = {'no_dem_sto': {'Q_proc':[],
                                 'Q_inc':[],
                                 'Q_cm': Qs_env,
                                 'Q_sde': []},
                  'no_env_sto': {'Q_proc':[],
                                 'Q_inc':[],
                                 'Q_cm': Qs_dem,
                                 'Q_sde': []},
                  'full': {'Q_proc':[],
                           'Q_inc':[],
                           'Q_cm': Qs,
                           'Q_sde': []},
                  'no_dem_sto_no_env_sto':{'Q_proc':[],
                                           'Q_inc':[],
                                           'Q_cm': [],
                                           'Q_sde': []}}

        if debug:
            for k in calc_Q:

                print '\n\nNon null term of Q_'+ k
                print "sv:"
                for i, x in enumerate(self.par_sv):
                    print i, x

                print "obs:"
                for i, x in enumerate(self.par_inc_def):
                    print N_PAR_SV+ i, x

                for i in range(len(calc_Q[k]['Q_cm'])):
                    for j in range(i+1):
                        if calc_Q[k]['Q_cm'][i][j]:
                            print '----------'
                            #print Q[i][j]
                            print 'Q_cm[{0}][{1}]: '.format(i, j),  self.make_C_term(calc_Q[k]['Q_cm'][i][j], True, human=True)
                            if i != j:
                                print 'Q_cm[{0}][{1}] == Q_cm[{1}][{0}]: '.format(i, j), self.make_C_term(calc_Q[k]['Q_cm'][i][j], True, human=True) == self.make_C_term(calc_Q[k]['Q_cm'][j][i], True, human=True)


        #convert in a version easy to template in C
        #Note that we only template the lower triangle (Q is symmetrical)
        for k, tpl in calc_Q.iteritems():
            if tpl['Q_cm']:
                for i  in range(len(tpl['Q_cm'])):
                    for j in range(i+1):
                        if tpl['Q_cm'][i][j]:
                            if i< N_PAR_SV and j < N_PAR_SV:
                                tpl['Q_proc'].append({'i': i, 'j': j, 'term': self.make_C_term(tpl['Q_cm'][i][j], True)})
                            else:
                                tpl['Q_inc'].append({'i': {'is_inc': False, 'ind': i} if i < N_PAR_SV else {'is_inc': True, 'ind': i - N_PAR_SV},
                                                     'j': {'is_inc': False, 'ind': j} if j < N_PAR_SV else {'is_inc': True, 'ind': j - N_PAR_SV},
                                                     'term': self.make_C_term(tpl['Q_cm'][i][j], True)})
            if sde:
                for i in range(len(Q_sde)):
                    for j in range(i+1):
                        if Q_sde[i][j]:
                            tpl['Q_sde'].append({'i': i, 'j': j, 'term': self.make_C_term(Q_sde[i][j], True)})


        ##cache special functions
        for key in calc_Q:
            if calc_Q[key]['Q_cm']:

                optim_rates_proc = [x['term'] for x in calc_Q[key]['Q_proc']]
                optim_rates_inc = [x['term'] for x in calc_Q[key]['Q_inc']]
                optim_rates = optim_rates_proc + optim_rates_inc

                calc_Q[key]['sf'] = self.cache_special_function_C(optim_rates, prefix='_sf[cac]')

                for i in range(len(optim_rates_proc)):
                    calc_Q[key]['Q_proc'][i]['term'] = optim_rates[i]

                n_proc = len(optim_rates_proc)
                for i in range(len(optim_rates_inc)):
                    calc_Q[key]['Q_inc'][i]['term'] = optim_rates[n_proc + i]

            else:
                calc_Q[key]['sf'] = []




        return calc_Q


if __name__=="__main__":

    """test Ccoder"""

    import json
    import os

    model = json.load(open(os.path.join('..' ,'example', 'model', 'datapackage.json')))
    m = Ccoder(model)
