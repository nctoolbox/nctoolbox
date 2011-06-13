% NCUNITS - Function to do unit conversions using the ucar.units library in the Unidata netcdf-java 
% package. 
% 
% NCTOOLBOX (http://code.google.com/p/nctoolbox)
%
%
function converted = ncunits(unconverted, originalunits, newunits)
    % NCUNITS - Function to convert data in one unit set to another set of units.
    %
    % Usage: converted_data = ncunits(data, 'm/s', 'cm/s');
    %
    %
    import ucar.nc2.units.*
    

converter = SimpleUnit.getConversionFactor(originalunits, newunits);
    
converted = unconverted .* converter;
    
    
% That was easier than I thought!
    
end