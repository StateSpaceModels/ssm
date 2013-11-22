var fs = require('fs')
  , tv4 = require("tv4")
  , _ = require('underscore')
  , util = require('util'); 

var op = ['+', '-', '*', '/', ',', '(', ')']
  , schema = require('./schema.json')
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

module.exports = function(pathDpkg, emitter, callback){

  function fail(err){
    if(err){
      emitter.emit('error', err);
      process.exit(1);
    }
  };
  
  //var valid = tv4.validate(noise, schema);
  //if (valid){
  //  console.log(valid);
  //} else {
  //  console.log(tv4.error);
  //}


  try {
    /**
     * check that package.json is a proper JSON file
     */ 
    dpkg = JSON.parse(fs.readFileSync(pathDpkg));

    /**
     * check that package.json is a proper SSM file
     */ 
    tv4check = tv4.validate(dpkg, schema);
    if (tv4check) {

    /**
     * check that the model is semantically complete and correct
     */
      model = dpkg.model
      //
      // Compartmental model state variables
      //

      var StateVariables = []
      model.populations.forEach(function(pop){pop.composition.map(function(x) {StateVariables.push(x)})});

      var Remainders = []
      model.populations.forEach(function(pop){
	if (_.contains(_.keys(pop),'remainder')){
	  Remainders.push(pop.remainder.name);
        }
      });

      var StateVariablesNoRem = _.difference(StateVariables, Remainders)
      
      var PopSizes = []
      model.populations.forEach(function(pop){
	if (_.contains(_.keys(pop),'remainder')){
	  PopSizes.push(pop.remainder.pop_size);
        }
      });

      var Parameters = []
      model.inputs.forEach(function(input){Parameters.push(input.name)});

      var EstimatedParameters = []
      model.inputs.forEach(function(input){
	if (_.contains(_.keys(input),'data')){
	  if( !Array.isArray(input.data)){
	    EstimatedParameters.push(input.name)
	  }
	}
      });

      var Observations = []
      model.observations.forEach(function(obs){Observations.push(obs.name)});

      var Data = []
      model.data.forEach(function(obs){Data.push(obs.name)});

      var Incidences = []
      model.reactions.forEach(function(r){
	if(_.contains(_.keys(r),'tracked')){
	  r.tracked.forEach(function(x){Incidences.push(x)});
	}
      });
      Incidences = _.uniq(Incidences);

      // No repetition in populations. 
      if ( _.uniq(StateVariables).length != StateVariables.length ){
	return fail("In 'populations', state variable cannot be repeated.");
      }

      // all 'to' and 'from' are state variables
      model.reactions.forEach(function(r){
	if( (r.to != 'U') && (StateVariables.indexOf(r.to) == -1) ){
	  return fail("In 'reactions', 'to' fields must designate state variables defined in 'populations'.");
	}
	if( (r.from != 'U') && (StateVariables.indexOf(r.from) == -1) ){
	  return fail("In 'reactions', 'from' fields must designate state variables defined in 'populations'.");
	}
      });

      // to and from must be different
      model.reactions.forEach(function(r){
	if(r.from == r.to){
	  return fail("In each reaction, 'to' and 'from' fields must be different.");
	}
      });

      // all state variables but remainder in parameters
      StateVariablesNoRem.forEach(function(s){
	if( Parameters.indexOf(s) == -1 ){
	  return fail(util.format("State variable %s is not defined in 'inputs'.",s));
	}
      })

      // Population sizes in parameters
      PopSizes.forEach(function(s){
	if( Parameters.indexOf(s) == -1 ){
	  return fail(util.format("Population size %s is not defined in 'inputs'.",s));
	}
      })

      // no tracked incidence from remainder
      model.reactions.forEach(function(r){
	if( (Remainders.indexOf(r.from) != -1) && _.contains(_.keys(r),'tracked') ){
	  if( r.tracked.length > 0){
	    return fail(util.format("The incidence variable %s is ill-defined. As %s has been defined as a remainder state variable, reactions leaving from it will be ignored.",_.first(r.tracked),r.from));
	  } 
	}
      });

      // all t0's are the same.
      t0s = []
      model.observations.forEach(function(o){
	t0s.push(o.start);
      });
      if ( _.uniq(t0s).length != 1 ){
	return fail("For the moment, SSM does not support observations starting on different dates. However, your resources can contain data starting from anytime before 'start' (SSM will ignore them).");
      }

      // SSM understands reaction rates
      var AllowedTerms = op.concat(Parameters, funcsSSM, jsmaths,'t'); //terms allowed in rates of the process model
      model.reactions.forEach(function(r){
	parseRate(r.rate).forEach(function(t){
	  if( (AllowedTerms.indexOf(t) == -1) && (t.indexOf('Math.') != 0) && (isNaN(parseFloat(t))) ){
	    return fail(util.format("In 'rate' %s, the term %s cannot be understood by SSM. Please define it.",r.rate,t));
	  }
	});
      });

      // Observations exactly map to data
      if (!_.isEqual(Observations,Data)){
	return fail(util.format("There is no one-to-one mapping between the elements of 'data' and the ones of 'observations' (see %s).",_.first(_.difference(Observations,Data))));
      }

      // Observed variables either correspond to state variables or tracked flows.
      var AllowedTerms = op.concat(Parameters, funcsSSM, Incidences, jsmaths, 't'); //terms allowed in rates of the process model
      model.observations.forEach(function(obs){
	parseRate(obs.mean).forEach(function(t){
	  if( (AllowedTerms.indexOf(t) == -1) && (t.indexOf('Math.') != 0) && (isNaN(parseFloat(t))) ){
	    return fail(util.format("In 'observations', 'mean' %s, the term %s cannot be understood by SSM. Please define it.",obs.mean,t));
	  };
	});
	parseRate(obs.sd).forEach(function(t){
	  if( (AllowedTerms.indexOf(t) == -1) && (t.indexOf('Math.') != 0) && (isNaN(parseFloat(t))) ){
	    console.log(t);
	    console.log(parseFloat('r'));
	    return fail(util.format("In 'observations', 'sd' %s, the term %s cannot be understood by SSM. Please define it.",obs.sd,t));
	  };
	});
      });

      // SSM only supports discretised_normal observation distribution so far
      model.observations.forEach(function(obs){
	if(obs.distribution != 'discretized_normal'){
	  return fail("For the moment, SSM only supports 'discretized_normal' distributions for observations.");
	}
      });


      if(_.contains(_.keys(model),'sde')){
	// Zero SDE drifts
	model.sde.drift.forEach(function(d){
	  if( d.f != 0.0 ){
	    return fail("For the moment, SSM does not support non-zero drifts in SDE's.");
	  }
	});
	
	// SDE transformations are understood by SSM
	var AllowedTerms = op.concat(Parameters, funcsSSM, jsmaths,'t'); //terms allowed in rates of the process model
	model.sde.drift.forEach(function(r){
	  parseRate(r.transformation).forEach(function(t){
	    if( (AllowedTerms.indexOf(t) == -1) && (t.indexOf('Math.') != 0) && (isNaN(parseFloat(t))) ){
	      return fail(util.format("In SDE's, 'transformation' %s, the term %s cannot be understood by SSM. Please define it.",r.transformation, t));
	    }
	  });
	});
	
	// SDE dispersions are understood by SSM
	var AllowedTerms = op.concat(Parameters, funcsSSM, jsmaths,'t'); //terms allowed in rates of the process model
	model.sde.dispersion.forEach(function(tmp){
	  tmp.forEach(function(r){
	    if (isNaN(parseFloat(r))){
	      parseRate(r).forEach(function(t){
		if( (AllowedTerms.indexOf(t) == -1) && (t.indexOf('Math.') != 0) && (isNaN(parseFloat(t))) ){
		  return fail(util.format("In SDE's, the dispersion term %s cannot be understood by SSM. Please define it.",t));
		}
	      });
	    }
	  });
	});
	
	// Variables on which SDE's are defined must belong to parameters, but cannot be state variables.
	model.sde.drift.forEach(function(d){
	  if ( Parameters.indexOf(d.name) == -1 ){
	    return fail(util.format("The variable %s on which you define an SDE is not defined in 'inputs'. Please define it.",d.name));
	  }
	  if ( (StateVariables.indexOf(d.name) != -1) || (PopSizes.indexOf(d.name) != -1) ){
	    return fail(util.format("SDE's cannot be defined on a state variable or population size (%s).",d.name));
	  }
	});
      }
	
      return dpkg;
	
    } else {
      return fail(tv4.error);
    }

  } catch (e){
    return fail(e.message);
    return fail("Please check package.json, it does not follow JSON standards.");
  }

}




