var Ctnr = require('data-streams').Ctnr
  , fs = require('fs')
  , once = require('once')
  , path = require('path')
  , clone = require('clone')
  , async = require('async');

/**
 * inputs is a list of object comming from model.data and/or model.inputs
 */
function resolve(dpkgRoot, dpkg, inputs, callback){
  
  var links = clone(inputs);

  async.eachSeries(links, function(l, cb){

    cb = once(cb);

    var mydpkg, s;

    if(!('require' in l)){
      return cb(null);
    }

    if('fields' in l.require){ //CSV (data and covariate)

      //create a container and use Ctnr to coerce
      var ctnr = {
        dataset: [
          {
            name: l.name,        
            about: [
              {name: l.require.fields[0], valueType: 'xsd:date'},
              {name: l.require.fields[1], valueType: 'xsd:double'}
            ],
            distribution: { 
              encodingFormat: 'text/csv', 
              contentPath: path.join(dpkgRoot, l.require.path) 
            }
          }
        ]
      };

      l.data = [];
      myctnr = new Ctnr(ctnr, dpkgRoot);
      s = myctnr.createReadStream(l.name, {coerce: true, filter: l.require.fields});
      s.on('error', cb);
      s.on('data', function(row){ l.data.push(row); });
      s.on('end', function(){ cb(null); });        

    } else { //JSON (prior)

      fs.readFile(path.join(dpkgRoot, l.require.path), function(err, file){
        if(err) return cb(err);
        try {
          l.data = JSON.parse(file)
        } catch(e){
          return cb(e);
        }
        cb(null);
      });

    }

  }, function(err){

    callback(err, links);

  });

};

exports.resolve = resolve;
