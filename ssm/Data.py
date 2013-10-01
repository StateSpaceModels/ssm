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

import os
import os.path
import json
import datetime
import sys
import csv

class DataError(Exception):
    def __init__(self, value):
        self.value = value
        def __str__(self):
            return repr(self.value)


class Data:

    def __init__(self, path,  **kwargs):
        self.path = os.path.abspath(path)
        self.model = json.load(open(self.path))
        self.root = os.path.dirname(self.path)


    def cast(self, row):
        for k, v in row.iteritems():
            if k == 'date':
                row[k] = datetime.datetime.strptime(v, "%Y-%m-%d").date()
            else :
                row[k] = float(v)

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
                    return resource['data']
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



if __name__=="__main__":

    d = Data(os.path.join('..' ,'example', 'model', 'datapackage.json'))
    print d.get_data("datapackages/data-sballest-noise/N/nyc")
