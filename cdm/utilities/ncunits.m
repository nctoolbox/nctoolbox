% NCUNITS - Function to do unit conversions using the ucar.units library in the Unidata netcdf-java
% package.
%
% NCTOOLBOX (http://code.google.com/p/nctoolbox)
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
originalunits = lower(originalunits);
newunits = lower(newunits);

% Check if conversion is possible, if not throw error
if SimpleUnit.isCompatible(originalunits, newunits)
    converter = SimpleUnit.getConversionFactor(originalunits, newunits); % Get conversion
    
    
else
    try
        originalunits = strrep(originalunits, 's/', '/');
        newunits = strrep(newunits, 's/', '/');
%         SimpleUnit.isCompatible(originalunits, newunits);
        converter = SimpleUnit.getConversionFactor(originalunits, newunits); % Get conversion
    catch
        error('NCTOOLBOX:NCUNITS', 'Input or output units are illegal, or input units are unable to be converted into the output units.');
    end
end

converted = unconverted .* converter; % Apply conversion
end