var Dpkg = require('data-streams')
  , fs = require('fs')
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

    var hasCallbacked = false;
    var mydpkg, s;

    if('datapackage' in l.require){ //cached

      if('fields' in l.require){ //SDF (data and covariate)

        //create a virtual data package and use Dpkg to coerce
        var resources = [{
          name: l.name, 
          path: path.resolve(dpkgRoot, 'data', l.name +'.csv'),
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
        s.on('error', function(err){ hasCallbacked = true; cb(err); });
        s.on('data', function(row){ l.data.push(row); });
        s.on('end', function(){ if(!hasCallbacked){ cb(null); } });        

      } else { //JSON (prior)

        fs.readFile(path.resolve(dpkgRoot, 'data', l.name +'.json'), function(err, file){
          if(err) return cb(err);
          try{
            l.data = JSON.parse(file)
          } catch(e){
            return cb(e);
          }
          cb(null);
        });

      }

    } else { //inside dpkg

      var r = clone(dpkg.resources.filter(function(x){return x.name === l.require.resource})[0]);
      if(!r) return cb(new Error('invalid link for '+ l.name));
      if( 'fields' in l.require) r.fields = l.require.fields;
      mydpkg = new Dpkg({resources: [r]}, dpkgRoot);
      s = mydpkg.createReadStream(r.name, {coerce: true});
      s.on('error', function(err){ hasCallbacked = true; cb(err); });
      if((r.format === 'json') || (('data' in r) && (typeof r.data !== 'string'))){
        s.on('data', function(row){ l.data = row; });
      } else {       
        l.data = [];
        s.on('data', function(row){ l.data.push(row); });
      }

      s.on('end', function(){ if(!hasCallbacked){ cb(null); } });        

    }

  }, function(err){
    callback(err, links);
  });
  
};

exports.resolve = resolve;
