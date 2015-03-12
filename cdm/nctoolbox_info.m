% nctoolbox_info  Provide information on the NCTOOLBOX version and dependencies
%

% David Forrest <drf@vims.edu>
% 2015-03-12

ver('nctoolbox')

fpath = fileparts(which('setup_nctoolbox_java'));
disp(sprintf('Java files in %s: ',fpath));
ls(fpath)

fpath=fileparts(which('setup_nctoolbox'));

if (exist(strcat(fpath,'/.git'),'dir'))
  % git information
  try
    disp 'NCtoolbox git install information:'
    [retval,gitdescribe]= system(sprintf('cd %s ; git describe --tags',fpath));
    if ~retval
	disp (sprintf('git describe info (tag-commits-commitID): %s',gitdescribe));
    end
    [retval,gitstatus]= system(sprintf('cd %s ;git status',fpath));
    if ~retval
        disp (gitstatus)
    end
  end
end


