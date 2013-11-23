var datapackage = require('datapackage')
  , fs = require('fs')
  , path = require('path')
  , clone = require('clone')
  , async = require('async');

/**
 * resolve (add data) links from a datapackage metadata
 * inputs is a list of object comming from model.data and/or model.inputs
 */
function resolve(dpkgRoot, dpkg, inputs, callback){
  
  var links = _linkify(inputs);

  async.eachSeries(links, function(l, cb){
    if( !(('schema' in l) || ('_link' in l)) ) {
      return cb(null);
    }

    //non official links _link (i.e links from non SDF resources)
    var cachedDeps = {};
    if('_link' in l){
      resolve_link(dpkgRoot, dpkg, cachedDeps, l._link, function(err, data){
        if(err) return cb(err);
        l.data = data;
        return cb(null);
      });
    } else { //SDF link

      var add2Resources = {};
      var fields = l.schema.fields;
      var f, resource;

      for(var i=0; i<fields.length; i++){
        f = fields[i];
        if( !('datapackage' in f.foreignkey) && !(f.foreignkey.resource in add2Resources) ){
          resource = dpkg.resources.filter(function(x){return x.name === f.foreignkey.resource;})[0];
          if(!resource){
            return cb(new Error('invalid link'));
          }
          add2Resources[f.foreignkey.resource] = resource;
        }
      }

      //create a virtual data package.
      var resources = [{"name": "__target", "schema": l.schema}];
      for(var key in add2Resources){
        resources.push(add2Resources[key]);
      }

      l.data = [];
      var s = datapackage.createReadStream(dpkgRoot, {resources: resources}, '__target', {coerce:true, foreignkeys:true});
      var hasCallbacked = false;
      s.on('error', function(err){
        hasCallbacked = true;
        cb(err);
      });
      s.on('data', function(row){
        l.data.push(row);
      });
      s.on('end', function(){
        if(!hasCallbacked){
          cb(null);
        }
      });

    }


  }, function(err){
    callback(err, links);
  });
  
};


/**
 * get links from inputs (model.data or model.inputs)
 */
function _linkify(inputs){
  var links = clone(inputs);

  links.forEach(function(d){
    if( ('data' in d) && Array.isArray(d.data) ){
      d.schema = { fields: d.data.map(function(x){return {foreignkey:x}}) };
    } else if( ('data' in d) && !Array.isArray(d.data) ){
      d._link = d.data;
      delete d.data;
    }
  });

  return links;
};


/**
 * resolve non official (i.e links from non SDF resources) links (_link) 
 * if name is present in link object, overwrite name of the fetched resource
 */
function resolve_link (dpkgRoot, dpkg, cachedDeps, _link, callback){

  if( !('datapackage' in _link) ){
    return _getResource(dpkg, _link, callback);
  } else {
    
    if (_link.datapackage in cachedDeps){
      return _getResource(cachedDeps[_link.datapackage], _link, callback);
    } else {

      fs.readFile(path.join(dpkgRoot, 'node_modules', _link.datapackage, 'package.json'), function(err, data){
        if(err) return callback(err);
        try {
          cachedDeps[_link.datapackage] = JSON.parse(data);
        } catch(e){
          return callback(e);
        }

        return _getResource(cachedDeps[_link.datapackage], _link, callback);
      });
      
    }
  }
};


function _getResource(dpkg, _link, callback){
  var r = clone(dpkg.resources.filter(function(x){return x.name === _link.resource})[0]);
  if(!r){ return callback(new Error('invalid link')); }

  if('name' in _link){
    r.name = _link.name;
  }
  return callback(null, r);
};


exports.resolve = resolve;
