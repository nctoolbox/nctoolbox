% nctoolbox_info  Provide information on the NCTOOLBOX version and dependencies
%

% David Forrest <drf@vims.edu>
% 2015-03-12

ver('nctoolbox')

fpath = fileparts(which('setup_nctoolbox_java'));
disp(sprintf('Java files in %s',fpath));
ls(fpath)

try
   system(sprintf('cd %s ; git describe --tags; git status',fpath));
end


