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

##Installing the required dependencies

C:
- [gsl](http://www.gnu.org/software/gsl/) (>= 1.15)
- [zmq](http://www.zeromq.org/) (3.2 release)
- [jansson](http://www.digip.org/jansson/) (>= 2.4)

Python:
- [Python 2.7.x](www.python.org/)
- [Jinja2](http://jinja.pocoo.org/docs/)
- [SymPy](http://sympy.org/)
- [dateutil](http://labix.org/python-dateutil)

On OSX with [homebrew](http://mxcl.github.io/homebrew/) and [pip](https://pypi.python.org/pypi/pip):

    brew install jansson zmq gsl
    sudo pip install jinja2 sympy python-dateutil

On Ubuntu:

    apt-get update
    apt-get install -y python-software-properties python g++ make build-essential
    add-apt-repository -y ppa:chris-lea/zeromq
    apt-get update
    apt-get install -y libzmq-dev libjansson-dev python-sympy python-jinja2 python-dateutil libgsl0-dev
 
##Installing S|S|M itself

    npm install -g ssm

Note: requires that all the C and python dependencies have been
installed _before_. It will build the standalone C libraries _and_
install them in ```/usr/local```. Edit the variable ```PREFIX``` of
the ```Makefile``` (in ```src/C```) to change this destination.

Pull requests are welcome for a .gyp file and windows support!

# Building models from a data package

From the command line:

    ssm build <package.json> [options]


Tests
=====

    npm test

Notes:
The C code is tested with [clar](https://github.com/vmg/clar)


License
=======

GPL version 3 or any later version.
