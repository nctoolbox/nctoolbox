% Ramadda Browser

function [nc]= ramadda_browse()
% Usage: >> output = edc_matlab;

mypath=which('ramadda_browse');
myDir=fileparts(mypath); % Find path
% myJava=strcat(myDir,'/Matlab/toolsUI-4.2.jar');
%
% addpath(strcat(myDir,'/Matlab'), '-end'); % Temp add Matlab folder to matlab path
% javaaddpath(myJava);                    % Temp add bundled (newer) Jar to javaclasspath
try
  if ispc
    [~, result] = dos(['cd ', myDir, ' & python ramaddacatalogbrowser.pyc']);
  elseif ismac
    try
      [~, result] = unix(['cd ', myDir, ' ; python ramaddacatalogbrowser.pyc']);
    catch
      [~, result] = system(['cd ', myDir, ' ; python ramaddacatalogbrowser.pyc']);
    end
  elseif isunix
    [~, result] = unix(['cd ', myDir, ' ; python ramaddacatalogbrowser.pyc']);
  end
catch
  error('Problem locating Python interpreter make sure that it is installed and located on the system path of your operating system.')
end

if length(result) > 0
  
  spaces = find(isspace(result));
  index = cat(2, 1, spaces);
  for i = 1:length(index)-1
    if i ~=length(index)-1
      links{i} = result(index(i)+2:index(i+1)-3);
    else
      links{i} = result(index(i)+2:index(i+1)-3);
    end
  end
  
  
  for i = 1:length(links)
    try
      disp('Connecting to dataset...')
        nc.(['dataset_',num2str(i)]) = ncgeodataset(links{i});
      disp('...Done')
    catch
      try
        nc.(['dataset_',num2str(i)]) = cfdataset(links{i});
        disp('...Done')
      catch
          warning(['Dataset ', num2str(i), ' cannot be opened on the server.'])
      end
    end
  end
else
  error('No datasets were selected or Ramadda returned empty resource links.');
end




