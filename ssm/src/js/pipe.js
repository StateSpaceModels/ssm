process.stdin.setEncoding('utf8');

process.stdin.on('data', function(chunk) {
  console.log('js: ' + chunk);

  setTimeout(function(){
    console.log("js: hello");
  }, 2000);
  

});

