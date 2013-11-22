var fs = require('fs')
  , tv4 = require("tv4")
  , _ = require('underscore')
  , util = require('util'); 

var op = ['+', '-', '*', '/', ',', '(', ')']
  , pkgSchema = require('./../json-schema/package-schema.json')
  , modelSchema = require('./../json-schema/model-schema.json')
  , priorSchema = require('./../json-schema/prior-schema.json')
  , funcsSSM = ['heaviside','ramp','sin','cos','tan','correct_rate']
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

module.exports = function(dpkg, emitter){
  
  if( !tv4.validate(dpkg, pkgSchema) ){
    // check that package.json complies with standard package.json requirements
    throw new Error(tv4.error.message);
    return;
  }

  if (!(_.contains(_.keys(dpkg), 'model'))){
    throw new Error("A 'model' object is required in your package.json to be installed with SSM");
    return;
  }

  var model = dpkg.model;
  if( !tv4.validate(model, modelSchema) ){
    // check that model complies with corresponding schema-json
    throw new Error(tv4.error.message);
    return;
  } 

    
  /**
   * check that the model is semantically complete and correct
   */


  //
  // Compartmental model state variables
  //
  
  var stateVariables = []
  model.populations.forEach(function(pop){pop.composition.map(function(x) {stateVariables.push(x)})});
  
  var remainders = []
  model.populations.forEach(function(pop){
    if (_.contains(_.keys(pop),'remainder')){
      remainders.push(pop.remainder.name);
    }
  });
  
  var stateVariablesNoRem = _.difference(stateVariables, remainders)
  
  var popSizes = []
  model.populations.forEach(function(pop){
    if (_.contains(_.keys(pop),'remainder')){
      popSizes.push(pop.remainder.pop_size);
    }
  });
  
  var parameters = []
  model.inputs.forEach(function(input){parameters.push(input.name)});
  
  var estimatedParameters = []
  model.inputs.forEach(function(input){
    if (_.contains(_.keys(input),'data')){
      if( !Array.isArray(input.data)){
	estimatedParameters.push(input.name)
      }
    }
  });
  
  var observations = []
  model.observations.forEach(function(obs){observations.push(obs.name)});
  
  var data = []
  model.data.forEach(function(obs){data.push(obs.name)});
  
  var incidences = []
  model.reactions.forEach(function(r){
    if(_.contains(_.keys(r),'tracked')){
      r.tracked.forEach(function(x){incidences.push(x)});
    }
  });
  incidences = _.uniq(incidences);
  
  var priors = []
  dpkg.resources.forEach(function(r){
    if( _.contains(_.keys(r),'name')){
      if(estimatedParameters.indexOf(r.name) != -1){
	priors.push(r)
      }
    } else {
      throw new Error("Every resource object needs to have a 'name' property'.");
    }
  });
  
  priors.forEach(function(p){
    if( !tv4.validate(p, priorSchema) ){
      // check that priors complies with corresponding schema-json
      throw new Error(tv4.error.message);
    }
  });
  
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
  
  // no tracked incidence from remainder
  model.reactions.forEach(function(r){
    if( (remainders.indexOf(r.from) != -1) && _.contains(_.keys(r),'tracked') ){
      if( r.tracked.length > 0){
	throw new Error(util.format("The incidence variable %s is ill-defined. As %s has been defined as a remainder state variable, reactions leaving from it will be ignored.",_.first(r.tracked),r.from));
      } 
    }
  });

  // all t0's are the same.
  t0s = []
  model.observations.forEach(function(o){
    t0s.push(o.start);
  });
  if ( _.uniq(t0s).length != 1 ){
    throw new Error("For the moment, SSM does not support observations starting on different dates. However, your resources can contain data starting from anytime before 'start' (SSM will ignore them).");
  }
  
  // SSM understands reaction rates
  var AllowedTerms = op.concat(parameters, funcsSSM, jsmaths,'t'); //terms allowed in rates of the process model
  model.reactions.forEach(function(r){
    parseRate(r.rate).forEach(function(t){
      if( (AllowedTerms.indexOf(t) == -1) && (t.indexOf('Math.') != 0) && (isNaN(parseFloat(t))) ){
	throw new Error(util.format("In 'rate' %s, the term %s cannot be understood by SSM. Please define it.",r.rate,t));
      }
    });
  });
  
  // observations exactly map to data
  if (!_.isEqual(observations,data)){
    throw new Error(util.format("There is no one-to-one mapping between the elements of 'data' and the ones of 'observations' (see %s).",_.first(_.difference(observations,data))));
  }
  
  // Observed variables either correspond to state variables or tracked flows.
  var AllowedTerms = op.concat(parameters, funcsSSM, incidences, jsmaths, 't'); //terms allowed in rates of the process model
  model.observations.forEach(function(obs){
    parseRate(obs.mean).forEach(function(t){
      if( (AllowedTerms.indexOf(t) == -1) && (t.indexOf('Math.') != 0) && (isNaN(parseFloat(t))) ){
	throw new Error(util.format("In 'observations', 'mean' %s, the term %s cannot be understood by SSM. Please define it.",obs.mean,t));
      };
    });
    parseRate(obs.sd).forEach(function(t){
      if( (AllowedTerms.indexOf(t) == -1) && (t.indexOf('Math.') != 0) && (isNaN(parseFloat(t))) ){
	throw new Error(util.format("In 'observations', 'sd' %s, the term %s cannot be understood by SSM. Please define it.",obs.sd,t));
      };
    });
  });
  
  // SSM only supports discretised_normal observation distribution so far
  model.observations.forEach(function(obs){
    if(obs.distribution != 'discretized_normal'){
      throw new Error("For the moment, SSM only supports 'discretized_normal' distributions for observations.");
    }
  });
  
  
  if(_.contains(_.keys(model),'sde')){
    // Zero SDE drifts
    model.sde.drift.forEach(function(d){
      if( d.f != 0.0 ){
	throw new Error("For the moment, SSM does not support non-zero drifts in SDE's.");
      }
    });
    
    // SDE transformations are understood by SSM
    var AllowedTerms = op.concat(parameters, funcsSSM, jsmaths,'t'); //terms allowed in rates of the process model
    model.sde.drift.forEach(function(r){
      if (_.contains(_.keys(r),'transformation')){
	parseRate(r.transformation).forEach(function(t){
	  if( (AllowedTerms.indexOf(t) == -1) && (t.indexOf('Math.') != 0) && (isNaN(parseFloat(t))) ){
	    throw new Error(util.format("In SDE's, 'transformation' %s, the term %s cannot be understood by SSM. Please define it.",r.transformation, t));
	  }
	});
      }
    });
    
    // SDE dispersions are understood by SSM
    var AllowedTerms = op.concat(parameters, funcsSSM, jsmaths,'t'); //terms allowed in rates of the process model
    model.sde.dispersion.forEach(function(tmp){
      tmp.forEach(function(r){
	if (isNaN(parseFloat(r))){
	  parseRate(r).forEach(function(t){
	    if( (AllowedTerms.indexOf(t) == -1) && (t.indexOf('Math.') != 0) && (isNaN(parseFloat(t))) ){
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




