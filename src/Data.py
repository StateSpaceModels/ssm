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
import os
import os.path
import json
import dateutil.parser
import datetime
import sys
import csv
from Ccoder import Ccoder

class DataError(Exception):
    def __init__(self, value):
        self.value = value
        def __str__(self):
            return repr(self.value)

class Data(Ccoder):

    def __init__(self, path_rendered, path, model_name, **kwargs):
        Ccoder.__init__(self, path, model_name, **kwargs)

        self.starts = [dateutil.parser.parse(x['start']).replace(tzinfo=None) for x in self.obs_model]
        self.t0 =  min(self.starts)

        try:        
            #the .data.json contains all the data and metadata comming from the data package and it's dependencies'
            self._data = json.load(open(os.path.join(path_rendered, self.model['name'],'.data.json')))
        except ValueError, IOError:
            raise DataError('could not process .data.json')
        
    def get_data(self, data_inputs_name):

        try:
            data = [x for x in self._data if x['name'] == data_inputs_name][0]['data']
        except IndexError:
            raise DataError('invalid resource name')

        return data


    def prepare_data(self):
        ##TODO pad begining in case t0 are different (so that reset zero is respected!!)

        obs = copy.deepcopy(self.obs_model)

        obs_id = [x['name'] for x in obs]

        dateset = set()
        data = {}
        for i, x in enumerate(obs):
            try:                
                linked_data = [d for d in self.model['data'] if d['name'] == x['name']][0]                
            except:
                raise DataError('invalid data link')                          

            date_name = linked_data['data'][0]['field']
            x_name = linked_data['data'][1]['field']
            
            oobs = copy.deepcopy(x)
            data[x['name']] = x
            data[x['name']]['order'] = i
            data[x['name']]['ind_inc_reset'] = [self.order_states[s] for s in self.get_inc_reset(oobs)]
            data[x['name']]['dict'] = {d[date_name]:d[x_name] for d in self.get_data(x['name'])}

            if 'transformation' in x:                
                on = linked_data['data'][1]['name'] #!!!! the hash need a name property (in addition to datapackage resource and field)
                f = eval('lambda {0}: {1}'.format(on, x['transformation']))
            else:
                f = lambda v: v

            data[x['name']]['f'] = f

            dateset |= set(data[x['name']]['dict'].keys())

        dates = list(dateset)
        dates.sort()

        data_C = []
        for d in dates:
            pd = dateutil.parser.parse(d).replace(tzinfo=None)
            row = {
                'date': pd.isoformat(),
                'observed': [],
                'values': [],
                'reset': [],
                'time': (pd-self.t0).days
            }

            for x in obs_id:
                if d in data[x]['dict']:
                    row['reset'].extend(data[x]['ind_inc_reset'])
                    if data[x]['dict'][d] is not None:
                        row['observed'].append(data[x]['order'])
                        row['values'].append(data[x]['f'](data[x]['dict'][d]))

            row['reset'] = list(set(row['reset']))
            data_C.append(row)

        return data_C

    def prepare_covariates(self):

        parameters = {x['name']:x for x in self.model['inputs']}

        data_C = []

        for p in self.par_forced:
            if 'transformation' in parameters[p]:
                f = eval('lambda {0}: {1}'.format(parameters[p]['data'][1]['name'], parameters[p]['transformation']))
            else:
                f = lambda x: x

            data =  self.get_data(p)
            date_name = parameters[p]['data'][0]['field']
            name = parameters[p]['data'][1]['field']
            x = []; y = []
            for d in data:
                if d[name] != None:
                    x.append((dateutil.parser.parse(d[date_name]).replace(tzinfo=None)-self.t0).days)
                    y.append(f(d[name]))

            data_C.append({'name': p, 'x': x, 'y': y})

        return data_C


if __name__=="__main__":

    d = Data(os.path.join('..' ,'examples', 'foo', 'ssm_model'), os.path.join('..' ,'examples', 'foo', 'package.json'), "sir")
    print d.prepare_covariates()
    print d.prepare_data()[29]
