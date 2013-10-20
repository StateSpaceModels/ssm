var path = require('path');

exports.resolvePath = function (string) {
  if (string.substr(0,1) === '~')
    string = process.env.HOME + string.substr(1)
  return path.resolve(string)
}
