var fs = require('fs')
, tv4 = require("tv4")
, _ = require('underscore')
, util = require('util');

var op = [ '+', '-', '*', '/', ',', '(', ')' ]
, modelSchema = require('../json-schema/model-schema.json')
, funcsSSM = ['heaviside', 'ramp', 'correct_rate', 'sigmoid']
, jsmaths =  Object.getOwnPropertyNames(Math);

function parseRate(rate){
  rate = rate.replace(/\s+/g, '');

  var s = ''
  , l = [];

  for (var i = 0; i< rate.length; i++){
    if (op.indexOf(rate[i]) !== -1){
      if(s.length){
        l.push(s);
        s = '';
      }
      l.push(rate[i]);
    } else {
      s += rate[i];
    }
  }

  if (s.length){
    l.push(s);
  }

  return l;
};

/**
 * check that the model is semantically complete and correct
 */
 module.exports = function(model){

  if( !tv4.validate(model, modelSchema) ){
    // check that model complies with corresponding schema-json
    throw new Error(tv4.error.message);
  } 
  
  //
  // Compartmental model state variables
  //
  
  var stateVariables = [];
  model.populations.forEach(function(pop){
    stateVariables = stateVariables.concat(pop.composition);
  });
  
  var remainders = [], popSizes = [];
  model.populations.forEach(function(pop){
    if ('remainder' in pop){
      remainders.push(pop.remainder.name);
      popSizes.push(pop.remainder.pop_size);
    } 
  });
  
  var stateVariablesNoRem = _.difference(stateVariables, remainders);
  
  var parameters = model.inputs.map(function(input){return input.name;});
  
  var observations = model.observations.map(function(x){return x.name;});
  
  var data = model.data.map(function(x){return x.name;});
  
  var incidences = []
  model.reactions.forEach(function(r){
    if('accumulators' in r){
      r.accumulators.forEach(function(x){incidences.push(x)});
    }
  });
  incidences = _.uniq(incidences);

  // no ode for the moment
  if('ode' in model){
    throw new Error("For the moment, SSM does not support ode's. This will come sortly.");
  }
  
  // No repetition in populations. 
  if ( _.uniq(stateVariables).length != stateVariables.length ){
    throw new Error("In 'populations', state variable cannot be repeated.");
  }
  
  // all 'to' and 'from' are state variables
  model.reactions.forEach(function(r){
    if( (r.to != 'U') && (stateVariables.indexOf(r.to) == -1) ){
      throw new Error("In 'reactions', 'to' fields must designate state variables defined in 'populations'.");
    }
    if( (r.from != 'U') && (stateVariables.indexOf(r.from) == -1) ){
      throw new Error("In 'reactions', 'from' fields must designate state variables defined in 'populations'.");
    }
  });
  
  // to and from must be different
  model.reactions.forEach(function(r){
    if(r.from == r.to){
      throw new Error("In each reaction, 'to' and 'from' fields must be different.");
    }
  });
  
  // all state variables but remainder in parameters
  stateVariablesNoRem.forEach(function(s){
    if( parameters.indexOf(s) == -1 ){
      throw new Error(util.format("State variable %s is not defined in 'inputs'.",s));
    }
  });
  
  // Population sizes in parameters
  popSizes.forEach(function(s){
    if( parameters.indexOf(s) == -1 ){
      throw new Error(util.format("Population size %s is not defined in 'inputs'.",s));
    }
  });
  
  // no accumulators incidence from remainder
  model.reactions.forEach(function(r){
    if( (remainders.indexOf(r.from) != -1) && ('accumulators' in r) ){
      if( r.accumulators.length > 0){
       throw new Error(util.format("The incidence variable %s is ill-defined. As %s has been defined as a remainder state variable, reactions leaving from it will be ignored.",_.first(r.accumulators),r.from));
     } 
   }
 });
  
  // all t0's are the same.
  var t0s = model.observations.map(function(o){return o.start;});

  if ( _.uniq(t0s).length != 1 ){
    throw new Error("For the moment, SSM does not support observations starting on different dates. However, your resources can contain data starting from anytime before 'start' (SSM will ignore them).");
  }
  
  // SSM understands reaction rates
  var allowedTerms = op.concat(parameters, funcsSSM, jsmaths,'t'); //terms allowed in rates of the process model
  model.reactions.forEach(function(r){
    parseRate(r.rate).forEach(function(t){
      if( (allowedTerms.indexOf(t) == -1) &&  (isNaN(parseFloat(t))) ){
        throw new Error(util.format("In 'rate' %s, the term %s cannot be understood by SSM. Please define it.",r.rate,t));
      }
    });
    if('white_noise' in r){
      if('sd' in r.white_noise){
        parseRate(r.white_noise.sd).forEach(function(t){
          if( (allowedTerms.indexOf(t) == -1) && (isNaN(parseFloat(t))) ){
            throw new Error(util.format("In 'rate' %s, the term %s cannot be understood by SSM. Please define it.",r.rate,t));
          }
        });
      } else {
        throw new Error(util.format("Please specify the amplitude of the white noise %s through with an 'sd' field.",r.white_noise.name));
      }
    }
  });

  
  // observations exactly map to data
  if (!_.isEqual(observations,data)){
    throw new Error(util.format("There is no one-to-one mapping between the elements of 'data' and the ones of 'observations' (see %s).",_.first(_.difference(observations,data))));
  }
  
  // Observed variables either correspond to state variables or accumulators flows.
  var allowedTerms = op.concat(parameters, funcsSSM, incidences, jsmaths, 't'); //terms allowed in rates of the process model

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

    } else if (obs.distribution == 'binomial'){


      parseRate(obs.p).forEach(function(t){
        if( (allowedTerms.indexOf(t) == -1) && (isNaN(parseFloat(t))) ){
          throw new Error(util.format("In 'observations', 'p' %s, the term %s cannot be understood by SSM. Please define it.",obs.p,t));
        };
      });

      parseRate(obs.n).forEach(function(t){
        if( (allowedTerms.indexOf(t) == -1) && (isNaN(parseFloat(t))) ){
          throw new Error(util.format("In 'observations', 'n' %s, the term %s cannot be understood by SSM. Please define it.",obs.n,t));
        };
      });

    };

  });

  // SSM only supports discretised_normal and poisson observation distributions so far
  model.observations.forEach(function(obs){
    if( (obs.distribution != 'discretized_normal') && (obs.distribution != 'poisson') && (obs.distribution != 'binomial')){
      throw new Error("For the moment, SSM only supports 'discretized_normal', 'poisson' and 'binomial' distributions for observations.");
    }
  });

  if('sde' in model){
    // Zero SDE drifts
    model.sde.drift.forEach(function(d){
      if( d.f != 0.0 ){
       throw new Error("For the moment, SSM does not support non-zero drifts in SDE's.");
     }
   });
    
    // SDE transformations are understood by SSM
    allowedTerms = op.concat(parameters, funcsSSM, jsmaths,'t'); //terms allowed in rates of the process model
    model.sde.drift.forEach(function(r){
      if ('transformation' in r){
       parseRate(r.transformation).forEach(function(t){
         if( (allowedTerms.indexOf(t) == -1) && (isNaN(parseFloat(t))) ){
           throw new Error(util.format("In SDE's, 'transformation' %s, the term %s cannot be understood by SSM. Please define it.",r.transformation, t));
         }
       });
     }
   });
    
    // SDE dispersions are understood by SSM
    model.sde.dispersion.forEach(function(tmp){
      tmp.forEach(function(r){
       if (isNaN(parseFloat(r))){
         parseRate(r).forEach(function(t){
           if( (allowedTerms.indexOf(t) == -1) && (isNaN(parseFloat(t))) ){
             throw new Error(util.format("In SDE's, the dispersion term %s cannot be understood by SSM. Please define it.",t));
           }
         });
       }
     });
    });
    
    // Variables on which SDE's are defined must belong to parameters, but cannot be state variables.
    model.sde.drift.forEach(function(d){
      if ( parameters.indexOf(d.name) == -1 ){
       throw new Error(util.format("The variable %s on which you define an SDE is not defined in 'inputs'. Please define it.",d.name));
     }
     if ( (stateVariables.indexOf(d.name) != -1) || (popSizes.indexOf(d.name) != -1) ){
       throw new Error(util.format("SDE's cannot be defined on a state variable or population size (%s).",d.name));
     }
   });
  }

};
