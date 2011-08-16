% ncpanelcomparison - Function to compare 2 planer 2d horzitonal netcdf model/sat datasets.
function [model1 model2 section1 section2 profile1 profile2] = ncpanelcomparison(dap1, dap2, varargin)
if nargin > 2
    inputs = arg2hash(varargin);
    colors = value4key(inputs, 'colormap');
    style = value4key(inputs, 'style');
    date = value4key(inputs, 'date');
    scale = value4key(inputs, 'scale'); % user could input 'log' for chl data, 'linear' by default
    depth = value4key(inputs, 'depth');
    pathcoords = value4key(inputs, 'path');
    profilecoords = value4key(inputs, 'profile');
    if isempty(colors)
        colors = 'default';
    end
    if isempty(style)
        style = 'default';
    end
    if isempty(scale)
        style = 'linear';
    end
    if isempty(depth)
        depth = 'surface';
    end
    if isempty(pathcoords)
        automatepath = 0;
    else
        automatepath = 1;
    end
    if isempty(profilecoords)
        automateprofile = 0;
    else
        automateprofile = 1;
    end
else
    style = 'default';
    colors = 'default';
    scale = 'linear';
    depth = 'surface';
    automateprofile = 0;
    automatepath = 0;
end


        
switch nargin
    case 0
        datasets = ramadda_browse; 
    case 1
        error('Less than two arguments. Please enter two dataset opendap links. Or try using modellook.m.')
    otherwise
        switch isempty(dap1) || isempty(dap2)
            case 1
                datasets = ramadda_browse;
            case 0
                datasets.dataset_1 = ncgeodataset(dap1);
                datasets.dataset_2 = ncgeodataset(dap2);
        end
end

date = inputdlg('Input Time As YYYY-MM-DD-hh:mm:ss or YYYY/MM/DD');



% Select var form dataset1
vars = datasets.dataset_1.variables;
for i = 1:length(vars)
    sizes{i} = [vars{i}, mat2str(datasets.dataset_1.size(vars{i}))];
   
end
[select, ok] = listdlg('ListString', sizes,'SelectionMode','single', 'Name', 'ModelLook: Quick nc vis',...
    'PromptString', 'Please Select Variable:' );
var1 = vars{select};

% Select var form dataset1
vars = datasets.dataset_2.variables;
for i = 1:length(vars)
    sizes{i} = [vars{i}, mat2str(datasets.dataset_2.size(vars{i}))];
   
end
[select, ok] = listdlg('ListString', sizes,'SelectionMode','single', 'Name', 'ModelLook: Quick nc vis',...
    'PromptString', 'Please Select Variable:' );
var2 = vars{select};

% Units dance:
units1 = datasets.dataset_1.attribute(var1, 'units');
units2 = datasets.dataset_2.attribute(var2, 'units');
[a b c match1] = regexpi(units1, {'cell', 'mol', 'chl'});
[a b c match2] = regexpi(units2, {'cell', 'mol', 'chl'});
if ~isempty(match1{1})
    error('no conversions available at the moment for cell')
elseif ~isempty(match1{2})
    units1 = 'mg/L';
elseif ~isempty(match1{3})  
    units1 = 'mg/L';
else
    % pass
end
if ~isempty(match2{1})
    error('no conversions available at the moment for cell')
elseif ~isempty(match2{2})
    units2 = 'mg/L';
elseif ~isempty(match2{3})  
    units2 = 'mg/L';
else
    units2 = units1; % get units from first dataset to convert second dataset (usually for temp or salt, or velocities)
end


% Start data access and plotting
figure;
b = subplot(1,2,1);
model1 = panellook(datasets.dataset_1, var1, date, colors, style, units1, scale, depth);
amerc % psudo-project data

% ar = get(b, 'PlotBoxAspectRatio');
% axis square
v = subplot(1,2,2);
model2 = panellook(datasets.dataset_2, var2, date, colors, style, units2, scale, depth);
amerc % psudo-project data

% colorbar('PlotBoxAspectRatio', [1 40 1])
% set(v, 'PlotBoxAspectRatio', ar);
% axis square
linkaxes([b v]);
linkprop([b v], 'CLim');



end