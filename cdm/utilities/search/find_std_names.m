function vars = find_std_names(nc,std_name)
% FIND_STD_NAMES finds variable that have the specified standard_name attribute
% vars = find_std_name(nc,std_name)
% Inputs: nc = ncgeodataset object or url
%         std_name = a CF standard_name  (e.g. 'sea_water_salinity')
% Output: vars = cell array of variable names that have the std_name attribute
    vars={};
    if ~isa(nc,'ncdataset'),
        nc=ncdataset(nc);
    end
    avars=nc.variables;
    k=1;
    for j=1:length(avars); %loop through variables to find standard_names
        vart=nc.variable(avars{j});
        svar=vart.attribute('standard_name');
        if strcmp(svar,std_name),
            vars{k}=avars{j};
        end
    end
