var fs = require('fs')
  , mkdirp = require('mkdirp')
  , rimraf = require('rimraf')
  , path = require('path')
  , spawn = require('child_process').spawn
  , links = require('../lib/links');


module.exports = function(dpkgRoot, dpkg, pathModel, keepSources, emitter, callback){

  function fail(err){
    if(err){
      emitter.emit('error', err.message);
      process.exit(1);
    }
  };

  var model = dpkg.model;

  mkdirp(pathModel, function (err) {
    if(err) return fail(err);

    //get data
    links.resolve(dpkgRoot, dpkg, model.data.concat(model.inputs).filter(function(x){return ('require' in x) && ('fields' in x.require);}), function(err, rlinks){
      if(err) return fail(err);

      fs.writeFile(path.join(pathModel, '.data.json'), JSON.stringify(rlinks), function(err){
        if(err) return fail(err);
        
        var tplter = [
          "import os",
          "import sys",
          "import json",
          "sys.path.append('" + path.resolve(__dirname, '..', 'src') + "')",
          "from Builder import Builder",
          "path_model_coded = '" + pathModel + "'",
          "dpkg = json.load(open('" + path.join(dpkgRoot, 'package.json') + "'))",
          "try:",
          "\tb = Builder(path_model_coded, '"+ dpkgRoot + "', dpkg)",
          "\tb.prepare()",
          "\tb.code()",
          "\tb.write_data()",
          "except:",
          "\tsys.exit(1)"
        ].join('\n');

        var templater = spawn('python', ['-c', tplter]);
        templater.stdout.setEncoding('utf8');
        templater.stderr.setEncoding('utf8');
        templater.stdout.on('data', function(data){
          emitter.emit('logEol', data);
        });
        templater.stderr.on('data', function(data){
          emitter.emit('errorEol', data);
        });

        templater.on('exit', function (code) {

          if(code === 0) {
            var make = spawn('make', ['install'], {cwd: path.join(pathModel, 'C', 'templates')});
            make.stdout.setEncoding('utf8');
            make.stderr.setEncoding('utf8');
            make.stdout.on('data', function(data){
              emitter.emit('logEol', data);
            });
            make.stderr.on('data', function(data){
              emitter.emit('errorEol', data);
            });

            make.on('exit', function (code) {
              if(keepSources){

                if(code === 0){
                  emitter.emit('success','model has been created in ' + pathModel);
                  callback(null);
                } else {
                  callback(new Error('could not build the model ('+ code +').'));
                }

              } else {

                rimraf(path.join(pathModel, 'C'), function(err){
                  if(err) emitter.emit('error', err.message);

                  if(code === 0){
                    emitter.emit('success','model has been created in ' + pathModel);
                    callback(null);
                  } else {
                    callback(new Error('could not build the model ('+ code +').'));
                  }

                });

              }
              
            });
          } else {
            callback(new Error('could not template the model ('+ code +').'));
          }

        });

      });

    });    

  });

};
