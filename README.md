S|S|M
=====

_Pipable_ plug-and-play inference methods for time series analysis
with *S*tate *S*pace *M*odels.

    cat theta.json | ./simplex -M 10000 | ./ksimplex -M 10000 > mle.json
    cat mle.json | ./kmcmc -M 100000 | ./pmcmc -J 1000 -M 500000


[![NPM](https://nodei.co/npm/ssm.png)](https://nodei.co/npm/ssm/)

All the methods are implemented in C. The C code contain generic part
(working with any models) and model specific part.  The specific parts
are templated using Python and [SymPy](http://sympy.org/) for symbolic
mathematics.

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
 

##Building the standalone C libraries

in src/C:

    make
    make install

Note: this will install the libs in ```/usr/local```. Edit the variable
```PREFIX``` of the```Makefile``` to change this destination.


##Installing the npm package

    npm install -g ssm


# Building a model from a data pacakge (package.json file)

From the command line:

    ssm build <package.json> [options]

This will build the binaries (in ```path_model_coded_in_C/```)


Tests
=====

C code (whith [clar](https://github.com/vmg/clar)):
In ```tests/```:

    make test

(see Makefile)

Python code:
In ```src/```:

    python -m unittest discover



License
=======

GPL version 3 or any later version.
