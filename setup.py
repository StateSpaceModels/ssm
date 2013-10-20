#!/usr/bin/env python

from distutils.core import setup

setup(name='ssm',
      version='0.0.0',
      description='Inference for State Space Models',
      author='Sebastien Ballesteros',
      author_email='sebastien@plom.io',
      url='http://www.plom.io',
      packages=['ssm'],
      package_dir={'ssm': 'src'},
      package_data={'ssm': ['C/templates/*']}
)
