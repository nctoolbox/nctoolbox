% nctoolbox_info  Provide information on the NCTOOLBOX version and dependencies
%

% David Forrest <drf@vims.edu>
% 2015-03-12

ver('nctoolbox')

disp(sprintf('Java files in %s',fileparts(which('setup_nctoolbox_java'))));
ls(fileparts(which('setup_nctoolbox_java')))

