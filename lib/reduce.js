var async = require('async')
  , fs = require('fs')
  , path = require('path');

/**
 * return the best dpkg and its path
 */
function bestDpkg(root, options, callback){
  fs.readdir(root, function(err, filenames){
    if(err) return callback(err);

    var jsonPaths = filenames
      .filter(function(x){
        console.log(x);
        var ext = path.extname(x);
        return ext === '.json' && x[0] !== '.';
      })
      .map(function(x){
        return path.join(root, x);
      });

    if(!jsonPaths.length){
      return callback(new Error('could not find any files containing a summary resource'));
    }

    async.map(jsonPaths, function(jsonPath, cb){
      fs.readFile(jsonPath, function(err, jsonString){
        if(err) return cb(err);
        
        try{
          var obj = JSON.parse(jsonString);
        }catch (e){
          return cb(err);
        }

        return cb(null, {obj:obj, path: jsonPath});
      });
    }, function(err, objsAndPaths){
      if(err) return callback(err);

      //filter objects to valide datapackages with summary resources      
      var dpkgsAndPaths = objsAndPaths.filter(function(x){
        return ('resources' in x.obj) && Array.isArray(x.obj.resources) && x.obj.resources.filter(function(x){return x.name === 'summary';})[0];
      });

      if(!dpkgsAndPaths.length){
        return callback(new Error('could not find any files containing a summary resource'));
      }

      //always sort from big to small
      //-1 smaller is better, 1 bigger is better
      var multipliers = {
        DIC: -1,
        AICc: -1, 
        AIC: -1,
        log_ltp: 1,
        log_likelihood: 1,
        sum_squares: -1 
      };

      var bestDpkgAndPath = dpkgsAndPaths.sort(function(a, b){
        var m = multipliers[options.by]
          , ca = a.obj.resources.filter(function(x){return x.name === 'summary';})[0]['data'][options.by]*m
          , cb = b.obj.resources.filter(function(x){return x.name === 'summary';})[0]['data'][options.by]*m;
        
        if(ca > cb) return -1;
        if(ca < cb) return 1;
        return 0
      })[0];
           
      callback(null, bestDpkgAndPath.obj, bestDpkgAndPath.path);
    });

  });

};

module.exports = bestDpkg;
