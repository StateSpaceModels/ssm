var path = require('path')
  , assert = require('assert')
  , predict = require('../lib/predict');

describe('predict', function(){
  
  it('should create the data object of a prediction resource', function(done){
    predict("2013-01-02", path.resolve(__dirname, 'data','x.csv'), path.resolve(__dirname, 'data','trace.csv'), function(err, data){
      if(err) throw err;

      var expected =  [
        { resources:
          [ { name: 'values', data: { c: 10, d: 60 } },
            { name: 'states', data: { a: 3, b: 7 } } ] },
        { resources:
          [ { name: 'values', data: { c: 40, d: 90 } },
            { name: 'states', data: { a: 4, b: 8 } } ] }
      ];

      assert.deepEqual(data, expected);

      done()      
    });
  });

});
