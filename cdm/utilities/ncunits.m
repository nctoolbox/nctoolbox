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

if SimpleUnit.isCompatible(originalunits, newunits)
    converter = SimpleUnit.getConversionFactor(originalunits, newunits);
    
    converted = unconverted .* converter;
else
    error('NCTOOLBOX:NCUNITS', 'Input or output units are illegal, or input units are unable to be converted into the output units.');
end

% That was easier than I thought!

end