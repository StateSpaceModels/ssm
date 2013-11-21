var tv4 = require("tv4");
var fs  = require("fs");

var noise = require('./noise_example.json');
var foo = require('./foo_example.json');
var schema = require('./schema.json');

var valid = tv4.validate(noise, schema);
if (valid){
  console.log(valid);
} else {
  console.log(tv4.error);
}

var valid = tv4.validate(foo, schema);
if (valid){
  console.log(valid);
} else {
  console.log(tv4.error);
}
