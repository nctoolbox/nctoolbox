function [desired_values] = hypoxia_units(current_values, current_units, desired_units)
% HYPOXIA_UNITS - In this case desired units must be stanard densities mass/vol
%
% Usage:
%             >> newvals = hypoxia_units(current_vals, current_units, desired_units)
%  Expects that the current units if molar are (millimole_x)
%
%
path = which('hypoxia_units.m');
table = load([path,'at']);

[units{1} units{2}] = strtok(current_units);

[scale compound] = strtok(units{1},'_');
compound = regexprep(compound, '_', '');

match = strcmpi(compound, table.hypoxia_units.periodic_table(:,1));
if sum(match) > 0
    amu = table.hypoxia_units.periodic_table{match, 4};
    makeup = '';
else
    clear match 
    match = strcmpi(compound, table.hypoxia_units.compounds(:,2));
    makeup = table.hypoxia_units.compounds{match, 3};
    amu = 0;
    for i = 1:length(makeup);
        if strcmp(makeup, 'chlorophyll')
            % pass
        else
            tempmatch = strcmpi(makeup{i}, table.hypoxia_units.periodic_table(:,1));
            amu = amu + table.hypoxia_units.compounds{tempmatch, 3};
        end
    end
end

if strcmp(makeup, 'chlorophyll')
    desired_values = ncunits(current_values, ['mg', units{2}], desired_units);
else
    
    % 1amu = 1.66053873e-24 g
    % N = 6.02214199e23 indv/mol
    g = (amu*(6.02214199e23))*(1.66053873e-24)*(current_values/1000); % mass in g    (amu/mol)*(g/amu)*(mmol)*(mol/mmol)
    
    desired_values = ncunits(g, ['g', units{2}], desired_units);
end




end