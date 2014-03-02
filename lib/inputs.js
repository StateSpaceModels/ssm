var Dpkg = require('data-streams').Dpkg
  , fs = require('fs')
  , once = require('once')
  , path = require('path')
  , clone = require('clone')
  , async = require('async');

/**
 * resolve (add data) links from a datapackage metadata
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

    if('fields' in l.require){ //SDF (data and covariate)

      //create a virtual data package and use Dpkg to coerce
      var resources = [{
        name: l.name, 
        path: path.join(dpkgRoot, 'data_modules', l.require.datapackage, 'data', l.require.resource +'.csv'),
        schema: {
          fields: [
            {name: l.require.fields[0], type: 'date'},
            {name: l.require.fields[1], type: 'number'}
          ]
        }, 
        fields: l.require.fields
      }];

      mydpkg = new Dpkg({resources: resources}, dpkgRoot);
      s = mydpkg.createReadStream(l.name, {coerce: true});
      l.data = [];
      s.on('error', cb);
      s.on('data', function(row){ l.data.push(row); });
      s.on('end', function(){ cb(null); });        

    } else { //JSON (prior)

      fs.readFile(path.resolve(dpkgRoot, 'data_modules', l.require.datapackage, 'data', l.require.resource +'.json'), function(err, file){
        if(err) return cb(err);
        try{
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
