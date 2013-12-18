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

    if('datapackage' in l.require){ //cached

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

    } else { //inside dpkg

      var r = clone(dpkg.resources.filter(function(x){return x.name === l.require.resource})[0]);

      if(!r) return cb(new Error('invalid link for '+ l.name));
      if('require' in r) return cb(new Error('resource of a model data pacakge cannot use require ('+ r.name + ')'));

      if( 'fields' in l.require) r.fields = l.require.fields;

      mydpkg = new Dpkg({resources: [r]}, dpkgRoot);

      var opts = {};
      if((r.format === 'json') || (('data' in r) && (typeof r.data !== 'string'))){
        r.format = 'json';
      } else {
        opts.coerce = true;
      }

      s = mydpkg.createReadStream(r.name, opts);
      s.on('error', cb);

      if(r.format === 'json'){
        s.on('data', function(row){ 
          try{
            l.data = JSON.parse(row); 
          } catch(e){
            cb(e);
          }
        });
      } else {       
        l.data = [];
        s.on('data', function(row){ l.data.push(row); });
      }

      s.on('end', function(){ cb(null) });        
      
    }

  }, function(err){
    callback(err, links);
  });

};

exports.resolve = resolve;
