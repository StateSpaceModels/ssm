var lhs = require('lhs')
  , util = require('util') 
  , links = require('./links')
  , clone = require('clone')
  , rmath = require('rmath.js');

function makeLhs(dpkgRoot, dpkg, model, options, callback){

  links.resolve(dpkgRoot, dpkg, model.inputs, function(err, inputs){

    var t0 = new Date(Math.min.apply(Math, model.observations.map(function(x){return new Date(x.start);})));

    try{
      var populations = _patchPopulations(model, inputs, t0);
    } catch(e){
      return callback(e);
    }

    var finputs = makeFinputs(inputs, t0, options);
    var pfit = inputs
      .filter(function(x){
        return ('data' in x) && !Array.isArray(x.data) && ('data' in x.data) && ('distribution' in x.data.data) && (x.data.data.distribution !== 'fixed');       
      })
      .map(function(x){return x.name;});
    
    var rlhs;
    var cnt = 0;
    var ok = true;

    do{
      rlhs = lhs.random(options.samples, pfit.length);
      ok = true;

      for(var i=0; i<rlhs.length; i++){
        pfit.forEach(function(pname, j){
          rlhs[i][j] = lhs.rescale(rlhs[i][j], finputs[pname].qmin, finputs[pname].qmax);
          finputs[pname].value = rlhs[i][j];       
        });

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
        populations[i]._toSum = populations[i].composition.filter(function(x){return x !== populations[i].remainder.name;});
      }
    }
  }

  return populations;
};


function _getValueT0(input, t0){
  var value;
  var dateName = input.schema.fields[0].foreignkey.field;
  var valueName = input.schema.fields[1].foreignkey.field;
  var j = 0;
  while( (t0 >= input.data[j][dateName]) && (j< input.data.length)){
    value = input.data[j++][valueName];
  }

  return value; 
};



/**
 * input a resolved input from model.inputs
 * p: a percentile
 */
function _qify(input, p){

  if (! (('data' in input) && !Array.isArray(input.data) && ('data' in input.data) && ('distribution' in input.data.data))){
    return undefined;
  }

  var prior = input.data.data;
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

    if(('data' in input) && !Array.isArray(input.data) && ('data' in input.data)){
      finputs[input.name]['value'] = (input.data.data.distribution === 'fixed') ? input.data.data.value : undefined;
      finputs[input.name]['priorName'] = input.data.name;
    } else if(('data' in input) && Array.isArray(input.data)) {
      finputs[input.name]['value'] = _getValueT0(input, t0);
      finputs[input.name]['priorName'] = input.schema.fields[1].foreignkey.name || input.name;
    } else {
      finputs[input.name]['value'] = undefined;
      finputs[input.name]['priorName'] = input.name;
    }

    //TODO change grammar so that can be used with with(Math){}...
    if('transformation' in input){
      finputs[input.name]['f'] = function(valuesObj){
        with(valuesObj){
          return eval(input.transformation);
        }
      };
    } else {
      finputs[input.name]['f'] = function(valuesObj){
        with(valuesObj){
          return eval(input.name);
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
  console.log(valuesObj);

  for(var i=0; i< populations.length; i++){

    if('remainder' in populations[i]){
      var sum = 0;
      populations[i]._toSum.forEach(function(pname){
        console.log('aaa', pname, finputs[pname].f(valuesObj));
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

var dpkgRoot = '/Users/team/ssm/examples/foo';
var dpkg = require('/Users/team/ssm/examples/foo/package.json');

makeLhs(dpkgRoot, dpkg, dpkg.models[0], {pmin: 0.05, pmax: 0.95, samples: 10, trials:3}, function(err, mylhs){
  if(err) console.log(err);
  console.log(mylhs);
});
