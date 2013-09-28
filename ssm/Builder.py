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
import tarfile
import shutil

from Ccoder import Ccoder

from jinja2 import Environment, FileSystemLoader

class Builder(Ccoder):
    """build a model"""

    def __init__(self, path_rendered, model,  **kwargs):
        Ccoder.__init__(self, model,  **kwargs)

        self.path_rendered = path_rendered
        self.env = Environment(loader=FileSystemLoader(os.path.join(self.path_rendered, 'C', 'templates')))

    def prepare(self, path_templates=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src', 'C'), replace=True):
        """
        copy templates to path_rendered
        """

        ##this function is called only when a new user has created or edited a model whose name is unique (primary key) so it is the only one able to recreate a model...
        if replace:
            if os.path.exists(self.path_rendered):
                shutil.rmtree(self.path_rendered)

        #copy templates to uploads/rendered/user_name/model_id
        if not os.path.exists(self.path_rendered):
            shutil.copytree(path_templates, os.path.join(self.path_rendered, 'C'))

    def archive(self, replace=True):
        """make a tarball"""

        tar = tarfile.open(os.path.join(os.path.dirname(self.path_rendered), os.path.basename(self.path_rendered)+'.tar.gz'), "w:gz")
        tar.add(self.path_rendered, arcname=os.path.basename(self.path_rendered))
        tar.close()

        if replace:
            if os.path.exists(self.path_rendered):
                shutil.rmtree(self.path_rendered)

    def render(self, prefix, data):

        template = self.env.get_template(prefix + '_template.c')
        with open(os.path.join(self.path_rendered, 'C', 'templates', prefix + ".c"), "w") as f:
            f.write(template.render(data))
            os.remove(os.path.join(self.path_rendered, 'C', 'templates', prefix + '_template.c'))

    def code(self):
        """generate C code for MIF, Simplex, pMCMC, Kalman, simulation, ..."""
                
        parameters = self.parameters()
        self.render('transform', parameters)
        self.render('input', parameters)

        observed = self.observed()
        self.render('observed', observed)





if __name__=="__main__":

    import json
    import os
    
    model = json.load(open(os.path.join('..' ,'example', 'model', 'datapackage.json')))        
    b = Builder(os.path.join(os.getenv("HOME"), 'ssm_test_model'), model)
 
    b.prepare()
    b.code()
