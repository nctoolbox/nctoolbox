function setup_ncdataset

home = fileparts(which(mfilename));
addpath(fullfile(home, 'java'));
addpath(fullfile(home, 'cdm'));

warning off
try 
    setup_nctoolbox_java;
catch me
    ex = MException('MBARI:NCTOOLBOX', 'Failed to setup the Java classpath');
    ex.throw
end
warning on
