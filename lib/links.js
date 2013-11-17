var datapackage = require('datapackage')
  , async = require('async');

/**
 * resolve (add data) SDF links from a datapackage metadata ([{name: ,
 * schema}])
 */
function resolve(dpkgRoot, links, callback){
  
  async.eachSeries(links, function(l, cb){
    var add2Resources = {};
    var fields = l.schema.fields;
    var f, resource;

    for(var i=0; i<fields.length; i++){
      f = fields[i];
      if( !('datapackage' in f.foreignkey) && !(f.foreignkey.resource in add2Resources) ){
        resource = dpgk.resources.filter(function(x){return x.name === f.foreignkey.resource;})[0];
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

  }, function(err){
    callback(err, links);
  });
  
};


/**
 * get links from a model
 */
function get(model){
  var links = [];
  model.data.concat(model.inputs).forEach(function(d){
    if( ('data' in d) && Array.isArray(d.data) && (d.data.length === 2)){
      links.push({
        name: d.name,
        schema: {
          fields: d.data.map(function(x){return {foreignkey:x}})
        }
      });
    }
  });

  return links;
};


exports.resolve = resolve;
exports.get = get;
