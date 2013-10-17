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
import datetime
import sys
import csv
from Cmodel import Cmodel

class DataError(Exception):
    def __init__(self, value):
        self.value = value
        def __str__(self):
            return repr(self.value)


class Data(Cmodel):

    def __init__(self, path,  **kwargs):
        self.path = os.path.abspath(path)
        model = json.load(open(self.path))
        Cmodel.__init__(self, model,  **kwargs)
        self.root = os.path.dirname(self.path)

        self.starts = [datetime.datetime.strptime(x['start'], "%Y-%m-%d").date() for x in self.obs_model]
        self.t0 =  min(self.starts)


    def cast(self, row):
        for k, v in row.iteritems():
            if k == 'date':
                row[k] = datetime.datetime.strptime(v, "%Y-%m-%d").date()
            else :
                try:
                    row[k] = float(v)
                except ValueError:
                    row[k] = None

        return row


    def get_field(self, resource, field, root=None):
        root = root or self.root
            
        try:
            f = [x for x in resource['schema']['fields'] if x['name'] == field][0]
        except IndexError:
            raise DataError("invalid field: " + field)            

        if 'foreignkey' in f:
            root = os.path.join(root, 'datapackages', f['foreignkey']['datapackage'])

            try:
                datapackage = json.load(open(os.path.join(root, 'datapackage.json')))                
                resource = [x for x in datapackage['resources'] if x['name'] == f['foreignkey']['resource']][0]
            except:
                print "Unexpected error:", sys.exc_info()[0]
                raise DataError("invalid foreignkey for " + f['name'])
        
            return self.get_field(resource, f['foreignkey']['field'], root=root)

        elif 'path' in resource:
            with open(os.path.join(root, resource['path']), 'rb') as f:
                reader = csv.DictReader(f)
                res = [self.cast({field: x[field]}) for x in reader]
                if not res:
                    raise DataError("invalid field for " + f['foreignkey']['field'])

                return res
            
        elif 'data' in resource:
            try:
                res = [{f['foreignkey']['field']: x[f['foreignkey']['field']]} for x in resource['data']]
            except:
                raise DataError("invalid field for " + f['foreignkey']['field'])                
                
            return res


        else:
            raise DataError('could not get data ' + resource['name'] + ' ' + fied)
        

    def get_data(self, resource):

        data = []

        for f in resource['schema']['fields']:
            data.append(self.get_field(resource, f['name']))
    
        sizes = [len(x) for x in data]
        if(len(set(sizes))) != 1:
            raise DataError("invalid tabular data for " + resource['name'])

        res = []
        for i in range(len(data[0])):
            row = {}
            for j in range(len(data)):
                row.update(data[j][i])
            res.append(row)

        return res


    def get_inc_reset(self, pdf):
        inc = set()
        for x in pdf:
            if x != "distribution":
                for e in self.change_user_input(pdf[x]):
                    if e in self.par_inc:
                        inc.add(e)

        return inc


    def prepare_data(self):

        ##TODO pad begining in case t0 are different (so that reset zero is respected!!)

        obs = copy.deepcopy(self.obs_model)

        obs_id = [x['name'] for x in obs]

        dateset = set()
        data = {}
        for i, x in enumerate(obs):
            data[x['name']] = x
            data[x['name']]['order'] = i
            data[x['name']]['ind_inc_reset'] = [self.order_states[s] for s in self.get_inc_reset(x['pdf'])]
            data[x['name']]['dict'] = {d['date']:d[x['name']] for d in self.get_data(x)}

            if 'transformation' in x:
                if 'schema' in x:
                    on = [y['name'] for y in x['schema']['fields'] if y['name'] != 'date'][0]
                else:
                    on = x['name']

                f = eval('lambda {0}: {1}'.format(on, x['transformation']))
            else:
                f = lambda v: v

            data[x['name']]['f'] = f

            dateset |= set(data[x['name']]['dict'].keys())


        dates = list(dateset)
        dates.sort()

        data_C = []
        for d in dates:
            row = {
                'date': d.isoformat(),
                'observed': [],
                'values': [],
                'reset': [],
                'time': (d-self.t0).days
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

        parameters = {x['name']:x for x in self.get_resource('parameters')}

        data_C = []

        for p in self.par_forced:
            if 'transformation' in parameters[p]:
                x = [x['name'] for x in p['schema']['fields'] if x['name'] != 'date'][0]
                f = eval('lambda {0}: {1}'.format(x, parameters[p]['transformation']))
            else:
                f = lambda x: x

            data =  self.get_data(parameters[p])
            name = [x for x in data[0].keys() if x!= 'date'][0]
            x = []; y = []
            for d in data:
                if d[name] != None:
                    x.append((d['date']-self.t0).days)
                    y.append(f(d[name]))

            data_C.append({'name': p, 'x': x, 'y': y})

        return data_C


if __name__=="__main__":

    d = Data(os.path.join('..' ,'example', 'foo', 'datapackages', 'model-seb-sir', 'datapackage.json'))
    print d.prepare_covariates()
    print d.prepare_data()[29]
