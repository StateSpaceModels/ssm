var lhs = require('lhs')
  , util = require('util') 
  , links = require('./links')
  , rmath = require('rmath.js');

function makeLhs(dpkgRoot, dpkg, model, options, callback){

  links.resolve(dpkgRoot, dpkg, model.inputs, function(err, inputs){

    var finputs = makeFinputs(inputs, options);
    console.log(util.inspect(finputs, { depth: 3 }));

    var pfit = inputs
      .filter(function(x){
        return ('data' in x) && !Array.isArray(x.data) && ('data' in x.data) && ('distribution' in x.data.data) && (x.data.data.distribution !== 'fixed');       
      })
      .map(function(x){return x.name;});


//    
//    var rlhs;
//    var cnt = 0;
//    var ok = true;
//
//    do{
//      rlhs = lhs.random(options.samples, pfit.length);
//      ok = true;
//
//      for(var i=0; i<rlhs.length; i++){
//        pfit.forEach(function(pname, j){
//          rlhs[i][j] = lhs.rescale(rlhs[i][j], inputs[pname].qmin, inputs[pname].qmax);
//          inputs[pname].value = rlhs[i][j];       
//        });
//
//        ok = checkIc(inputs, populations);
//        if(!ok) break;
//      }
//      
//    } while(!ok && cnt++ < options.trials);
//
//
//    if(ok){
//      callback(null, rlhs);
//    } else {
//      callback(new Error("can't statisfy constraints"));
//    }    

  });


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


function makeFinputs(inputs, options){
  var finputs = {};  


  inputs.forEach(function(input){
    finputs[input.name] = {
      'value': ( ('data' in input) && !Array.isArray(input.data) && ('data' in input.data) && (input.data.data.distribution === 'fixed') ) ? input.data.data.value : undefined,
      'qmin': _qify(input, options.pmin),
      'qmax': _qify(input, options.pmax),
      //'f': ('transformation' in p) ? transf[p.name] = Function( (p.prior && p.prior.name) || p.name, p.transformation ) : fid
    };
  });

  return finputs;
};


/**
 * 2 new properties have been added to populations:
 * - toSum
 * - size
 * we check than toSum < size;
 */

//var populations = clone(dpkg.resources.filter(function(x){return x.name === 'populations';})[0]);
//TODO: if remainders: add toSum and resolve N to it's value


function checkIc(inputs, populations){
  if(!populations) return true;

  for(var i=0; i< populations.length; i++){

    if('remainder' in populations[i]){
      var sum = 0;
      populations[i].toSum.forEach(function(pname){
        sum += inputs[pname].f(inputs);
      });

      if (sum >= populations[i].size) {
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
  console.log(mylhs);
});


