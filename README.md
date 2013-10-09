S|S|M
=====

_Pipable_ plug-and-play inference methods for time series analysis with *S*tate *S*pace *M*odels.

    cat theta.json | ./simplex -M 10000 | ./ksimplex -M 10000 > mle.json
    cat mle.json | ./kmcmc -M 100000 | ./pmcmc -J 1000 -M 500000

Overview
========

All the methods are implemented in plain C.  The C code contain
generic part (working with any models) and model specific part.  The
specific parts are templated using Python.

Installation
============

##Dependencies

C:
- [gsl](http://www.gnu.org/software/gsl/)
- [zmq](http://www.zeromq.org/)
- [jansson](http://www.digip.org/jansson/)

Python:
- [Python 2.7.x](www.python.org/)
- [Jinja2](http://jinja.pocoo.org/docs/)
- [SymPy](http://sympy.org/)

On OSX with [homebrew](http://mxcl.github.io/homebrew/) and [pip](https://pypi.python.org/pypi/pip):

    brew install jansson zmq gsl
    sudo pip install jinja2 sympy

On Ubuntu:

    apt-get update
    apt-get install -y python-software-properties python g++ make build-essential
    add-apt-repository -y ppa:chris-lea/zeromq
    apt-get update
    apt-get install -y libzmq-dev libjansson-dev python-sympy python-jinja2 libgsl0-dev
 

##Building the C libraries
in ssm/src/C:

    make
    make install

Note: this will install the libs in ```/usr/local```. Edit the variable
```PREFIX``` of the```Makefile``` to change this destination.

##Creating and installing the python package

At the root of the repo run:

    python setup.py sdist
    
in the package directory (that you will find after unpacking the tarball in ```dist/```):

    python setup.py install


License
=======

GPL version 3 or any later version.
