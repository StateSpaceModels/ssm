var lhs = require('lhs')
  , rmath = require('rmath.js');

function makeInputs(parameters, options){
  var inputs = {};  

  function fid(x) {
    return x;
  }

  function qify(prior, p){
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
      return x;

    };
  };

  parameters.forEach(function(p){
    inputs[p.name] = {
      'value': ( ('prior' in x) && (x.prior.distribution === 'fixed') ) ? x.prior.value : undefined,
      'qmin': ('prior' in  x) ? qify(prior, options.pmin): undefined,
      'qmax': ('prior' in  x) ? qify(prior, options.pmax): undefined,
      'f': ('transformation' in p) ? transf[p.name] = Function( (p.prior && p.prior.name) || p.name, p.transformation ) : fid
    };
  });

  return inputs;
};


/**
 * 2 new properties have been added to populations:
 * - toSum
 * - size
 * we check than toSum < size;
 */
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


module.exports = function(dpkg, options, callback){

  var populations = clone(dpkg.resources.filter(function(x){return x.name === 'populations';})[0]);
  //TODO: if remainders: add toSum and resolve N to it's value

  //TODO getStuff!!!
  getStuff(dpkg, function(err, dpkg){

    var parameters = dpkg.resources.filter(function(x){return x.name === 'parameters';})[0];

    var inputs = makeInputs(parameters, options);

    var pfit = parameters
      .filter(function(x){return ('prior' in x) && (x.prior.distribution !== 'fixed');})
      .map(function(x){return x.id;});
    
    var rlhs;
    var cnt = 0;
    var ok = true;

    do{
      rlhs = lhs.random(options.samples, pfit.length);
      ok = true;

      for(var i=0; i<rlhs.length; i++){
        pfit.forEach(function(pname, j){
          rlhs[i][j] = lhs.rescale(rlhs[i][j], inputs[pname].qmin, inputs[pname].qmax);
          inputs[pname].value = rlhs[i][j];       
        });

        ok = checkIc(inputs, populations);
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
