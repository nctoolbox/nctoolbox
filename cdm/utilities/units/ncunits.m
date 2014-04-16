% NCUNITS - Function to do unit conversions using the ucar.units library in the Unidata netcdf-java
% package.
%
% NCTOOLBOX (https://github.com/nctoolbox/nctoolbox)
%
% Usage: converted_data = ncunits(data, 'm/s', 'cm/s');
%
function converted = ncunits(unconverted, originalunits, newunits)
% NCUNITS - Function to convert data in one unit set to another set of units.
%
% Usage: converted_data = ncunits(data, 'm/s', 'cm/s');
%
%
import ucar.nc2.units.*

% Convert uppercase chars to lowercase (required)
% originalunits = lower(originalunits);
% newunits = lower(newunits);
switch originalunits
    case {'MM/S', 'MM/SEC','MMS/SEC', 'MMS/SEC', 'MILLIMETERS/S', 'MILLIMETERS/SEC'};
        originalunits = 'mm/s' 
    case {'CM/S', 'CM/SEC','CMS/S', 'CMS/SEC', 'CENTIMETERS/S', 'CENTIMETERS/SEC'};
        originalunits = 'cm/s';
    case {'M/S', 'M/SEC', 'METERS/S', 'METERS/SEC', 'MS/S', 'MS/SEC'};
        originalunits = 'm/s';
    case {'deg C', 'degree C', 'degrees C'};
        originalunits = 'degC';
    case {'ppt', 'psu', 'PPT', 'PSU'};
        originalunits = 'psu';
    otherwise
        % pass
end

switch newunits
    case {'MM/S', 'MM/SEC','MMS/SEC', 'MMS/SEC', 'MILLIMETERS/S', 'MILLIMETERS/SEC'};
        newunits = 'mm/s'
    case {'CM/S', 'CM/SEC','CMS/S', 'CMS/SEC', 'CENTIMETERS/S', 'CENTIMETERS/SEC'};
        newunits = 'cm/s';
    case {'M/S', 'M/SEC', 'METERS/S', 'METERS/SEC', 'MS/S', 'MS/SEC'};
        newunits = 'm/s';
    case {'deg C', 'degree C', 'degrees C'};
        newunits = 'degC';
    case {'ppt', 'psu', 'PPT', 'PSU'};
        newunits = 'psu';
    otherwise
        % pass
end


% Check if conversion is possible, if not throw error
if SimpleUnit.isCompatible(originalunits, newunits)
%     converter = SimpleUnit.getConversionFactor(originalunits, newunits); % Get conversion
%     converted = unconverted .* converter; % Apply conversion
    converter = SimpleUnit.factory(originalunits);
    newunits = SimpleUnit.factory(newunits);
    unconverted = double(unconverted);
    lengths = size(unconverted);
    lengths = [lengths ones([1 4-length(lengths)])];
    for i = 1:lengths(1)
        for j = 1:lengths(2);
            for k = 1:lengths(3);
                for l = 1:lengths(4);
                    converted(i,j,k,l) = converter.convertTo(unconverted(i,j,k,l), newunits);
                end
            end
        end
    end
    
else
%     try
        
%         originalunits = strrep(originalunits, 's/', '/');
%         newunits = strrep(newunits, 's/', '/');
% %         SimpleUnit.isCompatible(originalunits, newunits);
%         converter = SimpleUnit.getConversionFactor(originalunits, newunits); % Get conversion
%     catch
        error('NCTOOLBOX:NCUNITS', 'Input or output units are illegal, or input units are unable to be converted into the output units.');
%     end
end
converted = squeeze(converted);

end