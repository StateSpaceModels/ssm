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

        obs = self.get_resource('observations')
        self.starts = [datetime.datetime.strptime(x['start'], "%Y-%m-%d").date() for x in obs]
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



    def get_data(self, path, root=None, name=None):

        root = root or self.root

        plist = os.path.normpath(path).split(os.sep)
        if plist[0] == 'datapackages':
            try:
                datapackage =  json.load(open(os.path.join(*([root] + plist[0:2] + ['datapackage.json']))))
            except:
                print "Unexpected error:", sys.exc_info()[0]
                raise DataError("invalid link: " + path)



            r = [x for x in datapackage['resources'] if x['name']==plist[2]]
            if not r:
                raise DataError("invalid link: " + path)

            resource = r[0]
            if 'data' in resource:
                if len(plist) == 3:
                    return [self.cast(x) for x in resource['data']]
                else:
                    return [self.cast({'data': x['date'], plist[3]: x[plist[3]]}) for x in resource['data']]

            elif 'path' in resource:
                return self.get_data(resource['path'], root = os.path.join(*([root] + plist[0:2])), name = None if len(plist) < 4 else plist[3])

        else:
            with open(os.path.join(root, path), 'rb') as f:
                reader = csv.DictReader(f)

                if name:
                    return [self.cast({'date': x['date'], name: x[name]}) for x in reader]
                else:
                    return [self.cast(x) for x in reader]


    def prepare_data(self):

        ##TODO pad begining in case t0 are different (so that reset zero is respected!!)

        obs = copy.deepcopy(self.get_resource('observations'))

        obs_id = [x['id'] for x in obs]

        dateset = set()
        data = {}
        for i, x in enumerate(obs):
            data[x['id']] = x
            data[x['id']]['order'] = i
            data[x['id']]['data']['dict'] = {d['date']:d[x['id']] for d in self.get_data(x['data']['path'])}

            if 'transformation' in data[x['id']]:
                f = eval('lambda {0}: {1}'.format(data[x['id']]['data']['id'], data[x['id']]['transformation']))
            else:
                f = lambda x: x

            data[x['id']]['data']['f'] = f

            dateset |= set(data[x['id']]['data']['dict'].keys())


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
                if d in data[x]['data']['dict']:
                    row['reset'].append(data[x]['order'])
                    if data[x]['data']['dict'][d] is not None:
                        row['observed'].append(data[x]['order'])
                        row['values'].append(data[x]['data']['f'](data[x]['data']['dict'][d]))

            data_C.append(row)

        return data_C


    def prepare_covariates(self):

        parameters = {x['id']:x for x in self.get_resource('parameters')}

        data_C = []

        for p in self.par_forced:
            if 'transformation' in parameters[p]:
                f = eval('lambda {0}: {1}'.format(parameters[p]['prior']['id'], parameters[p]['transformation']))
            else:
                f = lambda x: x

            data =  self.get_data(parameters[p]['prior']['path'])
            name = [x for x in data[0].keys() if x!= 'date'][0]
            x = []; y = []
            for d in data:
                if d[name] != None:
                    x.append((d['date']-self.t0).days)
                    y.append(f(d[name]))

            data_C.append({'id': p, 'x': x, 'y': y})

        return data_C


if __name__=="__main__":

    d = Data(os.path.join('..' ,'example', 'foo', 'datapackages', 'model-seb-sir', 'datapackage.json'))
    print d.prepare_covariates()
