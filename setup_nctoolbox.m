function setup_nctoolbox

home = fileparts(which(mfilename));
addpath(fullfile(home, 'java'));
addpath(fullfile(home, 'cdm'));
addpath(fullfile(home, 'cdm', 'utilities'));


warning off
try 
    setup_nctoolbox_java;
catch me
    ex = MException('MBARI:NCTOOLBOX', 'Failed to setup the Java classpath');
    ex.throw
end
warning on
