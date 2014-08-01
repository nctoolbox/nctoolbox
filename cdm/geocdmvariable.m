classdef geocdmvariable < cdmvariable

    properties (SetAccess = private, GetAccess = private)
        axes_info % list of axes names and dinensions in one cell dict
    end

    methods

        function obj = geocdmvariable(cdmInstance, variableName)
            obj = obj@cdmvariable(cdmInstance, variableName)
        end

        function v = grid_interop(obj, varargin)
            ignore = {}; % TODO - implement ignore
            try 
                griddataset = ucar.nc2.dt.grid.GridDataset.open(obj.dataset.location);
            catch me

            end

            g = obj.grid(varargin{:})
            names = fieldnames(g);

            for i = 1:length(names) % Loop through fields returned by grid

                tempname = names{i};
                javaaxisvar = obj.dataset.netcdf.findVariable(tempname);
                try 
                    % Expecting a uvar.nc2.dataset.CoordinateAxis here
                    type = char(javaaxisvar.getAxisType());

                catch me 
                    % If not a CorrdinateAxis use empty string
                    type = '';
                end

                if isempty(type)
                    v.(tempname) = g.(tempname);
                else
                    switch type
                    case 'GeoX'
                        % Do Nothing. Handled in GeoY as GeoX and GeoY alwasy travel together
                    case 'GeoY'
                        [x y] = obj.handleGeoXY(tempname, g, griddataset);
                        if ~ismember('lon', ignore)
                            v.lon = x;
                        end
                        if ~ismember('lat', ignore)
                            v.lat = y;
                        end
                    case 'GeoZ'
                        if ~ismember('z', ignore)
                            [first, last, stride] = obj.toncindex(varargin);
                            v.z = obj.handleGeoZ(tempname, g, griddataset, first, last, stride);
                        end
                    case 'Height'
                        if ~ismember('z', ignore)
                            v.z = obj.handleHeight(tempname, g);
                        end
                    case 'Lat'
                        if ~ismember('lat', ignore)
                            v.lat = double(g.(tempname));
                        end
                    case 'Lon'
                        if ~ismember('lon', ignore)
                            v.lon = obj.handleLon(tempname, g);
                        end
                    case 'Time'
                        if ~ismember('time', ignore)
                            v.time = obj.handleTime(tempname, g);
                        end
                    otherwise
                        if ~ismember(tempname, ignore)
                            v.(tempname) = g.(tempname);
                        end
                    end
                end

            end
        end


    end

    methods (Access = protected)

        function v = handleHeight(obj, variableName, gridData) 
            javaaxisvar = obj.dataset.netcdf.findVariable(variableName);
            pos_z = char(javaaxisvar.getPositive());
            if strcmpi(pos_z, 'POSITIVE_DOWN') || strcmpi(pos_z, 'down')
                tmp = g.(variableName);
                v = temp * -1 % adjust for positive direction
            else
                v = g.(variableName);
            end
        end

        function [vx, vy] = handleGeoXY(obj, variableName, gridData, griddataset)
            g = griddataset.findGridByName(obj.name)
            cs = g.getCoordinateSystem();
            try 
                [x, y] = meshgrid(gridData.x, gridData.y);
                s = size(x);
                x = reshape(x, [1 numel(x)]);
                y = reshape(y, [1 numel(y)]);
                txy = [x; y];
                projection = g.getProjection();
                tmpLatLon = projection.projToLatLon(txy);
                vx = reshape(tmpLatLon(1, :), s);
                vy = reshape(tmpLatLon(2, :), s);
            catch me
                warning('NCTOOLBOX:geocdmvariable:handleGeoXY', 'Unable to process GeoX and GeoY');
                vx = [];
                vy = [];
            end
        end

        function v = handleGeoZ(obj, variableName, gridData, griddataset, first, last, stride)
            javaaxisvar = obj.dataset.netcdf.findVariable(variableName);
            pos_z = char(javaaxisvar.getPositive());
            z_sn = obj.dataset.attribute(variableName, 'standard_name');
            k = strfind(z_sn, 'ocean_s');
            if isempty(k)
                trange = ucar.ma2.Range(first(1) - 1, last(1) - 1, stride(1));
                zrange = ucar.ma2.Range(first(2) - 1, last(2) - 1, stride(2));
                xrange = ucar.ma2.Range(first(4) - 1, last(4) - 1, stride(4));
                yrange = ucar.ma2.Range(first(3) - 1, last(3) - 1, stride(3));
                grid = griddataset.findGridByName(obj.name);
                grid = grid.getCoordinateSystem();
                subgrid = grid.getVerticalTransform();
                subgrid = subgrid.subset(trange, zrange, yrange, xrange);
                % It looks like this dataset (http://geoport.whoi.edu/thredds/dodsC/usgs/vault0/models/examples/bora_feb.nc)
                % works when the vertical transform isnt subset. NJTBX uses the same methodology but subsets the resulting
                % z coordinate field in matlab after calling the getCoordinateArray method, not before in the java.
                % Should I implement a try and then default to subsetting afterwards on catch? Or should I just always do it
                % after? I am choosing the former. Anyone have any thoughts on the matter? -acrosby
                try
                    try
                        try
                            for q = first(1):stride(1):last(1)
                                array = subgrid.getCoordinateArray(q - 1);
                                v(q, :, :, :) = array.copyToNDJavaArray();
                            end
                        catch me
                            c = 1;
                            for q = first(1):stride(1):last(1)
                                subgrid = grid.getVerticalTransform();
                                array = subgrid.getCoordinateArray(q - 1); % Issue 27 is failing here...
                                v(c, :, :, :) = array.copyToNDJavaArray();
                                c = c + 1;
                            end
                            v = v(:, first(2):stride(2):last(2), ...
                                first(3):stride(3):last(3), ...
                                first(4):stride(4):last(4));
                        end
                    catch me
                        array = subgrid.getCoordinateArray(0);
                        v = array.copyToNDJavaArray();
                    end
                catch me
                    disp('Could you please add the code you are trying to run to Issue 27 at the nctoolbox issue tracking site.');
                    web http://code.google.com/p/nctoolbox/issues/detail?id=27
                    me.throw()
                    % me.error('There is a problem applying the vertical coordinate tranform and subsetting the resuting values.');
                end
            else
                v = g.(variableName);
                if strcmp(pos_z, 'POSITIVE_DOWN')
                    v = v * -1;
                end
            end
        end

        function v = handleTime(obj, variableName, gridData)
            tmp = g.(variableName);
            v = obj.dataset(time, variableName, tmp);
        end

        function v = handleLon(obj, variableName, gridData)
            tmp = g.(variableName);
            if max(max(tmp)) > 360
                if min(min(tmp)) > 0
                    tmp(tmp > 360) = tmp(tmp > 360) - 360;
                    tmp(tmp > 360) = tmp(tmp > 360) - 360;
                    tmp(tmp > 360) = tmp(tmp > 360) - 360;
                    tmp(tmp > 180) = tmp(tmp > 180) - 360;
                else
                    error('GEOCDMVARIABLE:grid_interop',...
                        'Longitude contains values that follow both -180/180 and 0/360+ conventions; can not subset.');
                end
            elseif min(tmp) <~ 0
                % Do nothing. Assume input is -180/+180 convention
            elseif max(tmp) > 180
                tmp = tmp - 360;
            end
            v = double(tmp);
        end

    end

end