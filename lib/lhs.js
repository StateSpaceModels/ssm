var lhs = require('lhs')
  , util = require('util') 
  , links = require('./links')
  , clone = require('clone')
  , rmath = require('rmath.js');

function makeLhs(dpkgRoot, dpkg, options, callback){

  links.resolve(dpkgRoot, dpkg, dpkg.model.inputs, function(err, inputs){

    var t0 = new Date(Math.min.apply(Math, dpkg.model.observations.map(function(x){return new Date(x.start);})));

    try{
      var populations = _patchPopulations(dpkg.model, inputs, t0);
    } catch(e){
      return callback(e);
    }

    var finputs = makeFinputs(inputs, t0, options);
    var pfit = inputs
      .filter(function(x){
        return ('require' in x) && !('fields' in x) && ('data' in x) && ('distribution' in x.data) && (x.data.distribution !== 'fixed');       
      })
      .map(function(x){return x.name;});
    
    var rlhs;
    var cnt = 0;
    var ok = true;
    var jsonrow;


    do{
      rlhs = lhs.random(options.samples, pfit.length);
      ok = true;

      for(var i=0; i<rlhs.length; i++){
        jsonrow = {};
        pfit.forEach(function(pname, j){
          jsonrow[finputs[pname].priorName] = lhs.rescale(rlhs[i][j], finputs[pname].qmin, finputs[pname].qmax);
          finputs[pname].value = jsonrow[finputs[pname].priorName];                 
        });
        rlhs[i] = jsonrow;

        ok = checkIc(finputs, populations);
        if(!ok) break;
      }
      
    } while(!ok && cnt++ < options.trials);

    if(ok){
      callback(null, rlhs);
    } else {
      callback(new Error("can't statisfy constraints"));
    }    

  });

};


/**
 * add _toSum and resolve pop_size to its value (_pop_size_value) in
 * model.populations
 */
function _patchPopulations(model, inputs, t0){

  var populations = clone(model.populations);

  if(populations){
    for(var i=0; i<populations.length; i++){
      if('remainder' in populations[i]){
        var pop = inputs.filter(function(x){return x.name === populations[i].remainder.pop_size})[0];
        if(!pop){
          throw new Error('invalid pop_size');
        }
        populations[i]._pop_size_value = _getValueT0(pop, t0);       
        if(!populations[i]._pop_size_value){
          throw new Error('could not find a value for t0');          
        }
        populations[i]._toSum = populations[i].composition.filter(function(x){return x !== populations[i].remainder.name;});
      }
    }
  }

  return populations;
};


function _getValueT0(input, t0){
  var value;
  var dateName = input.require.fields[0];
  var valueName = input.require.fields[1];
  var j = 0;

  while( (t0 >= input.data[j][dateName]) && (j< input.data.length) ){
    value = input.data[j++][valueName];
  }

  return value; 
};



/**
 * input a resolved input from model.inputs
 * p: a percentile
 */
function _qify(input, p){

  
  if ( !( ('require' in input) && !('fields' in input) && ('data' in input) && ('distribution' in input.data) ) ){
    return undefined;
  }

  var prior = input.data;
  var x;

  switch(prior.distribution){
  case 'uniform':
    return rmath.qunif(p, prior.lower, prior.upper, true, false);
    
  case 'normal':
    x =  rmath.qnorm(p, prior.mean, prior.sd, true, false);
    if('lower' in prior){
      x = Math.max(x, prior.lower);
    }
    if('upper' in prior){
      x = Math.min(x, prior.upper);
    }
    return x;

  case 'fixed':
    return undefined;

  };
};


function makeFinputs(inputs, t0, options){
  var finputs = {};  

  inputs.forEach(function(input){
    finputs[input.name] = {      
      'qmin': _qify(input, options.pmin),
      'qmax': _qify(input, options.pmax)
    };

    if(('require' in input) && !('fields' in input) && ('data' in input)){
      finputs[input.name]['value'] = (input.data.distribution === 'fixed') ? input.data.value : undefined;
      finputs[input.name]['priorName'] = input.require.resource;
    } else if(('require' in input) && ('fields' in input)) {
      finputs[input.name]['value'] = _getValueT0(input, t0);
      finputs[input.name]['priorName'] = input.require.name || input.name;
    } else {
      finputs[input.name]['value'] = undefined;
      finputs[input.name]['priorName'] = input.name;
    }

    if('transformation' in input){
      finputs[input.name]['f'] = function(valuesObj){
        with(Math){
          with(valuesObj){
            return eval(input.transformation);
          }
        }
      };
    } else {
      finputs[input.name]['f'] = function(valuesObj){
        with(Math){
          with(valuesObj){
            return eval(input.name);
          }
        }
      };
    }
  });


  return finputs;
};


/**
 * 2 new properties have been added to populations:
 * - _toSum
 * - _pop_size_value
 * we check than _toSum < _size;
 */
function checkIc(finputs, populations){
  if(!populations) return true;

  var valuesObj = {};
  for(var key in finputs){
    valuesObj[finputs[key].priorName] = finputs[key].value;
  }

  for(var i=0; i< populations.length; i++){

    if('remainder' in populations[i]){
      var sum = 0;
      populations[i]._toSum.forEach(function(pname){
        sum += finputs[pname].f(valuesObj);
      });      
      if (sum >= populations[i]._pop_size_value) {
        return false;  
      }
    }

  }

  return true;
};

module.exports = makeLhs;
