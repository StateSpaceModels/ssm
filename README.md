S|S|M
=====

Inference for time series analysis with *S*tate *S*pace *M*odels, like
playing with duplo blocks.

    cat theta.json | ./simplex -M 10000 | ./ksimplex -M 10000 > mle.json
    cat mle.json | ./kmcmc -M 100000 | ./pmcmc -J 1000 -M 500000 --trace > yeaaah.json

Maths, methods and algorithms
=============================

For more details on the modeling framework and on the algorithms 
available in SSM, see the [documentation](https://github.com/JDureau/ssm/raw/master/doc/doc.pdf).


Installation
============

All the methods are implemented in C. The C code contains generic parts
(working with any models) and model specific parts.  The specific parts
are templated using Python and [SymPy](http://sympy.org/) for symbolic
calculations. [JavaScript](https://brendaneich.com/brendaneich_content/uploads/CapitolJS.021.png)
is used to glue things together and add features on top of the C core.

## Installing the required dependencies

<a href="http://en.wikipedia.org/wiki/C_(programming_language)">C</a>:
- [gsl](http://www.gnu.org/software/gsl/) (>= 1.15)
- [zmq](http://www.zeromq.org/) (3.2 release)
- [jansson](http://www.digip.org/jansson/) (>= 2.4)

[Python 2.7.x](www.python.org/)
- [Jinja2](http://jinja.pocoo.org/docs/)
- [SymPy](http://sympy.org/)
- [dateutil](http://labix.org/python-dateutil)

[Node.js](http://nodejs.org/)

On Ubuntu:

    apt-get update
    apt-get install -y python-software-properties python g++ make build-essential
    add-apt-repository -y ppa:chris-lea/node.js
    add-apt-repository -y ppa:chris-lea/zeromq
    apt-get update
    apt-get install -y nodejs libzmq3-dev libjansson-dev python-sympy python-jinja2 python-dateutil libgsl0-dev

On OSX with [homebrew](http://brew.sh/) and [pip](https://pypi.python.org/pypi/pip):

    brew install jansson zmq gsl node
    sudo pip install jinja2 sympy python-dateutil
    
## Installing S|S|M itself

Download the library from [this link](http://5.135.176.187:3000/). 

Unzip it, and run the following in the terminal:

    cd ssm
    npm install -g

Note: this step requires that all the C and python dependencies have been
installed _before_ as it also builds the standalone C libraries.
We recommend _not_ to use ```sudo``` for this command.

If (and only if) you _have to_ use ```sudo``` to install package
globaly (```-g```) then proceed differently:

    git clone https://github.com/JDureau/ssm.git
    cd ssm
    npm install -g

Pull requests are welcome for a .gyp file and windows support!


## Last check

To make sure that everything is properly installed and ssm runs smoothly,
try the following:

    cd ssm/examples/tutorial
    ssm
    
If the model compiles successfully, you're good to go!


## Any problem with the install?

Check out the _Issues_ section of this repository (top right of this webpage), and if no issue solves your problem, post a new one describing it.

Usage
=====

What follows use [this example](https://github.com/JDureau/ssm/tree/master/examples/tutorial).
All the paths will be relative to this directory.

## Data and parameters (priors)

Data have to be in [CSV](http://tools.ietf.org/html/rfc4180) format
(following [RFC 4180](http://tools.ietf.org/html/rfc4180) ). Data MUST
contain unique headers AND a minimum of 2 columns, with one being
dates following the [ISO 8601](http://en.wikipedia.org/wiki/ISO_8601)
date definition (YYYY-MM-DD).

For instance:

    $ head data/data.csv
    
    "date","cases"
    "2012-08-02",5
    "2012-08-09",5
    "2012-08-16",6
    "2012-08-23",12
    "2012-08-30",null

Parameters (priors) have to be specified in [JSON](http://json.org/) or [JSON-LD](http://json-ld.org/) following:

    $ cat data/pr_v.json

```json
{
  "name": "normal",
  "distributionParameter" : [
    { "name" : "mean",  "value" : 12.5,   "unitCode": "DAY" },
    { "name" : "sd",    "value" : 3.8265, "unitCode": "DAY" },
    { "name" : "lower", "value" : 0,      "unitCode": "DAY" }
  ]
}
```

## Model

A model is described in [JSON](http://www.json.org/), typically in a
```ssm.json``` file.

S|S|M support any State Space Model built as system of ordinary or 
stochastic differential equations, a compartmental model, or a 
combination thereof.

The syntax to define a model is fully described as JSON
[schema](http://json-schema.org/)
[here](https://raw.github.com/JDureau/ssm/master/json-schema/model-schema.json).

### Link to the data

The first thing to do when writting a model is to _link_ it to the
data it explains.

    $ cat ssm.json | json data
```json
"data": [
  { 
    "name": "cases", 
    "require": { "path": "data/data.csv", "fields": ["date", "cases"] },
  }
]
```

The ```data.require``` property is a link pointing to a time-series.
A link is an object with 3 properties:
- ```path``` (mandatory), the path to the linked resource (in CSV or JSON)
- ```fields``` necessary only in case of resources containing data in
  CSV. In this later case, the first field must be the name of the
  column containing the dates of the time series and the second one
  the name of the column containing the actual values.
- ```name``` the name under which the resource should be imported.

Note that ```data``` itself can be a list so that multiple time-series
can be handled.

### Link to the priors and covariates

The same link objects are used to point to the resources that will be
used as priors or covariate of the model.

    $ cat ssm.json | json inputs
```json
"inputs": [
  {
    "name": "r0", 
    "description": "Basic reproduction number", 
    "require": { "name": "r0", "path": "data/r0.json" } 
  },
  { 
    "name": "v",
    "description": "Recovery rate",
    "require": { "name":  "pr_v", "path": "data/pr_v.json" },
    "transformation": "1/pr_v",
    "to_resource": "1/v" 
  },
  {
    "name": "S", 
    "description": "Number of susceptible",
    "require": { "name": "S", "path": "data/S.json" } 
  },
  { 
    "name": "I",
    "description": "Number of infectious", 
    "require": { "name": "I", "path": "data/I.json" } 
  },
  { 
    "name": "R", 
    "description": "Number of recovered",
    "require": { "name": "R", "path": "data/R.json" } 
  },
  { 
    "name": "rep",
    "description": "Reporting rate",
    "require": { "name": "rep", "path": "data/rep.json" } 
  }
],
```

Note that this linking stage also allows to include some
_transformations_ so that a relation can be established between your
model requirement and existing priors or covariates living in other
datapackages. For example ```v``` (a rate) is linked to a prior
expressed in duration: ```pr_v``` through an inverse transformation. 


### Process Model

The process model can be expressed as an ODE, an SDE or a compartmental
model defining a Poisson process (potentialy with stochastic rates).  
Let's take the example of a simple Susceptible-Infected-Recovered
compartmental model for population dynamics. The process model
contains the following properties:

the populations 

    $ cat ssm.json | json populations
```json
"populations": [
  {"name": "NYC", "composition": ["S", "I", "R"]}
]
```
and the reactions, defining the process model

    $ cat ssm.json | json reactions
```json
"reactions": [
  {"from": "S", "to": "I", "rate": "r0/(S+I+R)*v*I", "description": "infection", "accumulators": ["Inc"]},
  {"from": "I", "to": "R", "rate": "v", "description":"recovery"}
]
```
Note that the populations object is a list. Structured populatiols can be
defined by appending terms to the list.


An ```sde``` property can be added in case you want that some
parameters follow diffusions (see
[here](https://github.com/JDureau/ssm/blob/master/examples/foo/package.json)
for an example, and [here](http://arxiv.org/abs/1203.5950) for 
references). White environmental noise can also be added to the reaction
as in this [example](https://raw.github.com/JDureau/ssm/master/examples/noise/package.json)
(references [here](http://arxiv.org/abs/0802.0021)).

The ```accumulators``` property allows to defined new state variable
(here ```Inc```) that will accumulate the flow of the reaction they
label.  Accumulators state variables are reset to 0 for each data
point related to the accumulator state.

### Observation model

One observation model has to be defined per observed time-series.

    $ cat ssm.json | json observations
```json
"observations": [
  {
    "name": "cases",
    "start": "2012-07-26",
    "distribution": "discretized_normal",
    "mean": "rep * Inc",
    "sd": "sqrt(rep * ( 1.0 - rep ) * Inc )"
  }
]
```
SSM currently implements the `discretized_normal`, `poisson` and `binomial`observation processes. See the developer [documentation](doc/developers/dev_help.md) to contribute and add your favourite distribution. 

### Initial conditions

Finally, values of the parameters and the covariance matrix between
them need need to be defined in a separate JSON file typicaly named
```theta.json```. ```theta.json``` will be used as initial values for
inference algorithms:

    $ cat theta.json
```json
"resources": [
  {
    "name": "values",
    "description": "initial values for the parameters",
    "data": {
      "r0": 25.0,
      "pr_v": 11.0
    }
  },
  {
    "name": "covariance",
    "description": "covariance matrix",
    "data": {
      "r0": {"r0": 0.04, "pr_v": 0.01},
      "pr_v": {"pr_v": 0.02, "r0": 0.01}
    }
  }
]
```  
Only the diagonal terms are mandatory for the covariance matrix.


## Installing a model from a configuration file

At the root of a directory with a configuration file (```ssm.json```), run

    $ ssm [options]

This will build executables (in
```bin/```) for several inference and simulation methods
([MIF](http://www.pnas.org/content/103/49/18438),
[pMCMC](http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2009.00736.x/abstract),
[simplex](http://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method),
[SMC](http://en.wikipedia.org/wiki/Particle_filter),
[Kalman filters](http://en.wikipedia.org/wiki/Kalman_filter), ...)
customized to different implementation of you model
([ode](http://en.wikipedia.org/wiki/Ordinary_differential_equation),
[sde](http://en.wikipedia.org/wiki/Stochastic_differential_equation),
[poisson process with stochastic rates](http://arxiv.org/pdf/0802.0021.pdf),
...).


All the methods are directly ready for *parallel computing* (using
multiple cores of a machine _and_ leveraging a cluster of machines).

Run ```./method --help``` in ```bin/``` to get help and see the different
implementations and options supported by the method.
In the same way, help for the ```ssm``` command can be obtained with
```ssm --help```

## Inference like playing with duplo blocks

Everything that follows supposes that we are in ```bin/``` and that ```theta.json``` has been moved into ```bin/```.

Let's start by plotting the data

with [R](http://www.r-project.org/):

```R
data <- read.csv('../data/data.csv', na.strings='null')
plot(as.Date(data$date), data$cases, type='s')
```

Let's run a first simulation:

     $ cat theta.json | ./simul --traj

And add the simulated trajectory to our first plot

```R
traj <- read.csv('X_0.csv')
lines(as.Date(traj$date), traj$cases, type='s', col='red')
```

Let's infer the parameters to get a better fit

     $ cat theta.json | ./simplex -M 10000 --trace > mle.json

let's read the values found:

     $ cat mle.json | json resources | json -c "this.name=='values'"
```json
[
 {
   "name": "values",
   "data": {
     "pr_v": 19.379285906561037,
     "r0": 29.528755614881494
   }
 }
]
```

Let's plot the evolution of the parameters:

```R
trace <- read.csv('trace_0.csv')
layout(matrix(1:3,1,3))
plot(trace$index, trace$r0, type='l')
plot(trace$index, trace$pr_v, type='l')
plot(trace$index, trace$fitness, type='l')
```

Now let's redo a simulation with these values (```mle.json```):

     $ cat mle.json | ./simul --traj -v

and replot the results:

```R
plot(as.Date(data$date), data$cases, type='s')
traj <- read.csv('X_0.csv')
lines(as.Date(traj$date), traj$cases, type='s', col='red')
```

to realize that the fit is now much better.

And now in one line:

    $ cat theta.json | ./simplex -M 10000 --trace | ./simul --traj | json resources | json -c "this.name=='values'"
```json 
[
  {
    "name": "values",
    "data": {
      "r0": 29.528755614881494,
      "pr_v": 19.379285906561037
    }
  }
]
```

Let's get some posteriors and sample some trajectories by adding a
pmcmc at the end of our pipeline (we actualy add 2 of them to skip
the convergence of the mcmc algorithm).

     $ cat theta.json | ./simplex -M 10000 | ./pmcmc -M 10000 | ./pmcmc -M 100000 --trace --traj  | json resources | json -c 'this.name=="summary"'
```json
[
 {
   "name": "summary",
   "data": {
     "id": 0,
     "log_ltp": -186.70579009197556,
     "AICc": 363.94320971360844,
     "n_parameters": 2,
     "AIC": 363.6765430469418,
     "DIC": 363.6802334782078,
     "log_likelihood": -179.8382715234709,
     "sum_squares": null,
     "n_data": 48
   }
 }
]
```

Some posteriors plots (still with R)

```R
 trace <- read.csv('trace_0.csv')
 layout(matrix(1:2,1,2))
 hist(trace$r0)
 hist(trace$pr_v)
```

The sampled trajectories

```R
 traj <- read.csv('X_0.csv')
 plot(as.Date(data$date), data$cases, type='s')
 samples <- unique(traj$index)
 for(i in samples){
   lines(as.Date(traj$date[traj$index == i]), traj$cases[traj$index == i], type='s', col='red')
 }
```

## Be cautious

Always validate your results... SSM outputs are fully compatible with
[CODA](http://cran.r-project.org/web/packages/coda/index.html).

In addition to the diagnostic provided by
[CODA](http://cran.r-project.org/web/packages/coda/index.html), you
can run S|S|M algorithn with the ```--diag``` option to add some
diagnostic outputs.  For instance let's run a particle filter with
1000 particles (```--J```) with a stochastic version of our model
(```psr```) after a simplex:

    $ cat theta.json | ./simplex -M 10000 | ./smc psr -J 1000 --diag  --verbose

the ```--diag``` option give us access to the prediction residuals and
the effective sample size. Let's plot these quantities

```R
  diag <- read.csv('diag_0.csv')
  layout(matrix(1:3,3,1))

  #data vs prediction
  plot(as.Date(data$date), data$cases, type='p')
  lines(as.Date(diag$date), diag$pred_cases, type='p', col='red')

  #prediction residuals
  plot(as.Date(diag$date), diag$res_cases, type='p')
  abline(h=0, lty=2)

  #effective sample size
  plot(as.Date(diag$date), diag$ess, type='s')
```

## Parallel computing

Let's say that you want to run a particle filter of a stochastic
version of our previous model with 1000 particles on your 4 cores
machines (```--n_thread```). Also instead of plotting 1000
trajectories you just want a summary of the empirical confindence
envelopes (```--hat```).

    $ cat theta.json | ./smc psr -J 1000 --n_thread 4 --hat

Let's plot the trajectories

```R
  hat <- read.csv('hat_0.csv')
  plot(as.Date(hat$date), hat$mean_cases, type='s')
  lines(as.Date(hat$date), hat$lower_cases, type='s', lty=2)
  lines(as.Date(hat$date), hat$upper_cases, type='s', lty=2)
```

Your machine is not enough ? You can use several.  First let's
transform our ```smc``` into a _server_ that will dispatch some work to
several _workers_ (living on different machines).

    $ cat theta.json | ./smc psr -J 1000 --tcp

All the algorithm shipped with S|S|M can be transformed into servers
with the ```--tcp``` option.

Now let's start some workers giving them the address of the server.

    $ cat theta.json | ./worker psr smc --server 127.0.0.1 &
    $ cat theta.json | ./worker psr smc --server 127.0.0.1 &

Note that you can add workers at any time during a run.

# Plugins


## Piping to the future

S|S|M can also be used to perform predictions.

```ssm-predict``` allows to re-create initial conditions adapted to
the ```simul``` program from the trace and trajectories sampled from
the posterior distributions obtained after Bayesian methods
(```pmcmc```, ```kmcmc```).

You can install this plugin with

    npm install -g ssm-predict

And use it with

    $ ssm-predict theta.json X_0.csv trace_0.csv 2012-11-22 | ./simul --start 2012-11-22 --end 2013-12-25 --verbose --hat

We can plot the results of this prediction taking care to extend the
xlim on our first plot. For the prediction we ran ```simul``` with the
```--hat``` option that will output empirical credible envelop
instead of all the projected trajectories (as does ```--traj```).

```R
  data <- read.csv('../data/data.csv', na.strings='null')
  plot(as.Date(data$date), data$cases, type='s', xlim=c(min(as.Date(data$date)), as.Date('2013-12-25')))
  
  traj <- read.csv('X_0.csv') #from the previous run
  samples <- unique(traj$index)
  for(i in samples){
      lines(as.Date(traj$date[traj$index == i]), traj$cases[traj$index == i], type='s', col='red')
  }
      
  hat <- read.csv('hat_0.csv') #from the current run
  lines(as.Date(hat$date), hat$mean_cases, type='s' , col='blue')
  lines(as.Date(hat$date), hat$lower_cases, type='s', lty=2, col='blue')
  lines(as.Date(hat$date), hat$upper_cases, type='s', lty=2, col='blue')
```


License
=======

GPL version 3 or any later version.
