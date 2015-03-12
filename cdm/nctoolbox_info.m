% nctoolbox_info  Provide information on the NCTOOLBOX version and dependencies
%

% David Forrest <drf@vims.edu>
% 2015-03-12

ver('nctoolbox')

fpath = fileparts(which('setup_nctoolbox_java'));
disp(sprintf('Java files in %s: ',fpath));
ls(fpath)

fpath = fileparts(mfilename('fullpath'));
rootpath = fileparts(fpath);
gitpath = fullfile(rootpath, '.git');
[rv, r] = system('git --version');
%useGit = exist(gitpath,'dir') && ~rv;
useGit = false; % for testing without git

if useGit
  % git information
  try
    disp 'NCtoolbox git install information:'
    [rv, ncDescribe] = system(sprintf('cd %s; git describe --tags', rootpath));
    [rv, ncCommit] = system(sprintf('cd %s; git rev-parse --short HEAD', rootpath));
    % git log behave weirdly on old version, like the one shipped by Apple.
    % Disabled date for now
    %[rv, ncDate] = system(sprintf('cd %s;git --no-pager log -1 --format=%cd --date=short', rootpath));
    [rv, ncBranch] = system(sprintf('cd %s;git rev-parse --abbrev-ref HEAD', rootpath));
    [rv, ncDate] = system(sprintf('cd %s;git --no-pager log -1 --format=%cd --date=short', rootpath));
    fprintf(1, 'GIT COMMIT:   %s', ncCommit);
    fprintf(1, 'GIT BRANCH:   %s', ncBranch);
    %fprintf(1, 'GIT DATE:     %s', ncDate);
    fprintf(1, 'GIT DESCRIBE: %s', ncDescribe);
  catch
      useGit = false;
  end
end

if ~useGit
    infofile = fullfile(fpath, 'nctoolbox-version.txt');
    if exist(infofile, 'file')
        type(infofile);
    end
end

fprintf(1, '\n')
