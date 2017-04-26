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
import dateutil.tz
import datetime
import sys
import csv
from Ccoder import Ccoder

def parse_date(x):
    """
    parse a ISO8601 string into a datetime and convert it into a date
    TODO support for datetime (necessitate work on C side)...

    """

    return dateutil.parser.parse(x).date()

class DataError(Exception):
    def __init__(self, value):
        self.value = value
        def __str__(self):
            return repr(self.value)

class Data(Ccoder):

    def __init__(self, path_rendered, dpkgRoot, dpkg, **kwargs):
        path_rendered = unicode(path_rendered, 'utf8')
        Ccoder.__init__(self, dpkgRoot, dpkg, **kwargs)

        self.starts = [parse_date(x['start']) for x in self.obs_model]
        self.t0 =  min(self.starts)

        try:        
            #the .data.json contains all the data and metadata comming from the data package and it's dependencies'
            self._data = json.load(open(os.path.join(path_rendered, '.data.json')))
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

            date_name = linked_data['require']['fields'][0]
            x_name = linked_data['require']['fields'][1]

            oobs = copy.deepcopy(x)
            data[x['name']] = x
            data[x['name']]['order'] = i
            data[x['name']]['ind_inc_reset'] = [self.order_states[s] for s in self.get_inc_reset(oobs)]
            data[x['name']]['dict'] = {d[date_name]:d[x_name] for d in self.get_data(x['name'])}

            if 'transformation' in x:                
                if 'name' not in linked_data['require']:
                    raise DataError('the require hash need a name property (the transformation has to be done in terms of this name)')
                on = linked_data['require']['name']
                f = eval('lambda {0}: {1}'.format(on, x['transformation']))
            else:
                f = lambda v: v

            data[x['name']]['f'] = f

            dateset |= set(data[x['name']]['dict'].keys())

        dates = list(dateset)
        dates.sort()

        data_C = []
    
        for i, d in enumerate(dates):
            pd = parse_date(d)
            row = {
                'date': pd.isoformat(),
                'observed': [],
                'values': [],
                'reset': [],
                'time': (pd-self.t0).days
            }

            if i > 0:
                for x in obs_id:
                    if dates[i-1] in data[x]['dict']:
                        row['reset'].extend(data[x]['ind_inc_reset'])
            
            for x in obs_id:
                if d in data[x]['dict']:        
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
                f = eval('lambda {0}: {1}'.format(parameters[p]['require']['fields'][1], parameters[p]['transformation']))
            else:
                f = lambda x: x

            data =  self.get_data(p)
            date_name = parameters[p]['require']['fields'][0]
            name = parameters[p]['require']['fields'][1]
            x = []; y = []
            for d in data:
                if d[name] != None:
                    x.append((parse_date(d[date_name])-self.t0).days)
                    y.append(f(d[name]))

            data_C.append({'name': p, 'x': x, 'y': y})

        return data_C


if __name__=="__main__":

    dpkgRoot = os.path.join('..' ,'examples', 'foo')
    dpkg = json.load(open(os.path.join(dpkgRoot, 'package.json')))
    d = Data(os.path.join('..' ,'examples', 'foo', 'bin'), dpkgRoot, dpkg)
    print d.prepare_covariates()
    print d.prepare_data()[29]
    
