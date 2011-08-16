function [temp] = timeseriesfrom3d(temp, out)
%
%
%
%
%
[temp.temp temp.temp.temp temp.profile] = nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, out(1), out(2));
try
    temp.data = double(squeeze(temp.var.data((temp.ind-100):(temp.ind+100), temp.profile(1), temp.profile(2))));
    temp.depths = temp.var.grid_interop((temp.ind-100):(temp.ind+100), temp.profile(1), temp.profile(2));
catch
    try
        temp.data = double(squeeze(temp.var.data((temp.ind-100):end, temp.profile(1), temp.profile(2))));
        temp.depths = temp.var.grid_interop((temp.ind-100):end, temp.profile(1), temp.profile(2));
    catch
        try
            temp.data = double(squeeze(temp.var.data(1:(temp.ind+100), temp.profile(1), temp.profile(2))));
            temp.depths = temp.var.grid_interop(1:(temp.ind+100), temp.profile(1), temp.profile(2));
        catch
            temp.data = double(squeeze(temp.var.data(1:end, temp.profile(1), temp.profile(2))));
            temp.depths = temp.var.grid_interop(1:end, temp.profile(1), temp.profile(2));
        end
    end
end

end
