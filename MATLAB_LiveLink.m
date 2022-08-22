% This code is to connect to the COMSOL livelink with MATLAB server

% add the livelink path to the server for comsol5.4
addpath( '/appl/comsol/6.0/mli/' );

% connect to the server, the default port number is 2036 and then 2037 ...
mphstart(2036)

% the MATLAB file you want to run
MAIN;

% exit the program
exit;
