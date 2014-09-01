# Help documentation for SSM developers

This document provides some help for developers who wants to contribute to the SSM library by adding new functionalities that require to modify the source code.

## How to add a new observation process?

SSM originally implements the `discretized_normal` observation process. Adding a new observation distribution require to modify four files: 

* `src/C/templates/observed_template.c`
* `lib/validate.js`
* `src/Ccoder.py`
* `src/Cmodel.py`

Let's illustrate this by adding the `poisson` distribution, which is fully parametrized by its `mean` and would be specified in `ssm.json` as follows:

```json
"observations" : [
        {
            "name" : "incidence_observed",
            "start" : "1789-07-14",
            "distribution" : "poisson",
            "mean" : "reporting_rate * incidence"
        }
    ]
```

where `reporting_rate` is the rate at which new cases are reported and would be defined in the `inputs`.

### Modify `src/C/templates/observed_template.c`

This file contains the template for the likelihood (`f_likelihood_tpl_`) and random generator (`f_obs_ran_tpl_`) functions of the observation process. Both need to be modified by adding `elif` commands to switch between distributions and make calls to suitable functions of the [GSL](https://www.gnu.org/software/gsl/manual/html_node/The-Poisson-Distribution.html) library to perform the core calculation.

```C
static double f_likelihood_tpl_{{ x.name }}(double y, ssm_X_t *p_X, ssm_par_t *par, ssm_calc_t *calc, double t)
{
    double like;
    double *X = p_X->proj;

    {% if x.distribution == 'discretized_normal' %}
    
    double gsl_mu = {{ x.mean }};
    double gsl_sd = {{ x.sd }};
    if (y > 0.0) {
        like = gsl_cdf_gaussian_P(y + 0.5 - gsl_mu, gsl_sd) - gsl_cdf_gaussian_P(y - 0.5 - gsl_mu, gsl_sd);
    } else {
        like = gsl_cdf_gaussian_P(y + 0.5 - gsl_mu, gsl_sd);
    }

    {% elif x.distribution == 'poisson' %}

    double gsl_mu = {{ x.mean }};
    like = gsl_ran_poisson_pdf(rint(y), gsl_mu);

    {% endif %}

    return like;
}
```
```C
static double f_obs_ran_tpl_{{ x.name }}(ssm_X_t *p_X, ssm_par_t *par, ssm_calc_t *calc, double t)
{
    double *X = p_X->proj;

    {% if x.distribution == 'discretized_normal' %}
    
    double gsl_mu = {{ x.mean }};
    double gsl_sd = {{ x.sd }};
    double yobs = gsl_mu + gsl_ran_gaussian(calc->randgsl, gsl_sd);
    yobs = (yobs >0) ? yobs : 0.0;

    {% elif x.distribution == 'poisson' %}

    double gsl_mu = {{ x.mean }};
    double yobs = gsl_ran_poisson(calc->randgsl, gsl_mu);

    {% endif %}

    return yobs;
}
```
### Modify `lib/validate.js`

This file contains the function `module.exports` that checks that the model is semantically complete and correct, and throw error messages otherwise.

```javascript
model.observations.forEach(function(obs){

  if(obs.distribution == 'discretized_normal'){
    
    parseRate(obs.mean).forEach(function(t){
     if( (allowedTerms.indexOf(t) == -1) && (isNaN(parseFloat(t))) ){
       throw new Error(util.format("In 'observations', 'mean' %s, the term %s cannot be understood by SSM. Please define it.",obs.mean,t));
     };
   });
    
    parseRate(obs.sd).forEach(function(t){
     if( (allowedTerms.indexOf(t) == -1) && (isNaN(parseFloat(t))) ){
       throw new Error(util.format("In 'observations', 'sd' %s, the term %s cannot be understood by SSM. Please define it.",obs.sd,t));
     };
   });
    
  } else if (obs.distribution == 'poisson'){
   
    parseRate(obs.mean).forEach(function(t){
     if( (allowedTerms.indexOf(t) == -1) && (isNaN(parseFloat(t))) ){
       throw new Error(util.format("In 'observations', 'mean' %s, the term %s cannot be understood by SSM. Please define it.",obs.mean,t));
     };
   });
  };
});

model.observations.forEach(function(obs){
  
  if( (obs.distribution != 'discretized_normal') && (obs.distribution != 'poisson')){
    throw new Error("For the moment, SSM only supports 'discretized_normal' and 'poisson' distributions for observations.");
  }
});
```

### Modify `src/Ccoder.py`

This file contains the method `observed` which defines all the parameters needed by SSM to handle the distribution. In particular, several algorithms in SSM perform Gaussian approximations, which require us to define the `mean` and `sd` (standard deviation) parameters.

```python
def observed(self):

    obs = copy.deepcopy(self.obs_model)
    
    for x in obs:
        if x['distribution'] == 'discretized_normal':
            x['mean'] = self.make_C_term(x['mean'], True)
            x['sd'] = self.make_C_term(x['sd'], True)
        elif x['distribution'] == 'poisson':
            x['mean'] = self.make_C_term(x['mean'], True)
            x['sd'] = self.make_C_term(x['sd'], True)

    return {'observed': obs}
```
Note that all the parameters of the distribution need to appear in addition to `mean` and `sd`. For instance, if one wants to add a binomial distribution parametrized by its probability `p` and sample size `n`, these two parameters must also be defined here.

### Modify `src/Cmodel.py`

Finally, the appropriate `mean` and `sd` of the Gaussian approximation must be defined in this file. We can use the fact that when `mean > 1000`, the Poisson distribution can be approximated by a normal distribution with the same `mean` and `sd = sqrt(mean)`. Thus, we only need to compute the standard deviation:

```python
for o in observations:
            if o['distribution'] == 'poisson':
                o['sd'] = 'sqrt('+o['mean']+')'

            for p in [o['mean'], o['sd']]:
                el =  self.change_user_input(p)
                for e in el:
                    if e not in self.op and e not in self.reserved and e not in self.special_functions and e not in self.par_sv and e not in self.par_noise and e not in self.par_proc and e not in self.par_forced and e not in self.par_inc:
                        try:
                            float(e)
                        except ValueError:
                            par_obs.add(e)

        self.par_obs = sorted(list(par_obs))
```

