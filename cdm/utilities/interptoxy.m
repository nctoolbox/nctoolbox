        
        function interpolatedvalues = interptoxy(data, lon, lat, x, y, method)
          % NCGEOVARIABLE.INTERPTOXY - Interpolate gridded data at an x/y (lon/lat) location or set of locat-
          % ions when inputs are vectors of the same length.
          % Useage: >> data = var.interptoxy(data, lat, lon, -77, 45, 'linear');
          s = size(data);
          if length(x) ~= length(y)
            error('NCTOOLBOX:interptoxy','Input coordinate vectors do not agree.')
          end
          
          % Need good method for dealing with vector ans matrix lat/lon then to vector for scattered 
          % interpretation.
          if isvector(lat)
            lat = repmat(lat, 1, length(lon));
            lon = repmat(lon, length(lat), 1);
          end
          lat = reshape(lat, numel(lat), 1);
          lon = reshape(lon, numel(lon), 1);
          
          % Convert to double
          data = double(data);
         
          if length(s) > 2
            for j = 1:s(1)
              
              if length(s) > 3
                for i=1:s(2)
                  redata = reshape(data(j,i,:,:), numel(data(j,i,:,:)), 1);
                  int = TriScatteredInterp(lon, lat, redata, method);
                  clear redata
                  for k = 1:length(x)
                    interpolatedvalues(j,i,k) = int(x(k),y(k));
                  end
                end
                
              else
                redata = reshape(data(j,:,:), numel(data(j,:,:)), 1);
                int = TriScatteredInterp(lon, lat, redata, method);
                clear redata
                for k = 1:length(x)
                  interpolatedvalues(j,k) = int(x(k),y(k)); 
                end
              end

            end
            
          else
            redata = reshape(data, numel(data), 1);
            int = TriScatteredInterp(lon, lat, redata, method);
            clear redata
            for k = 1:length(x)
              interpolatedvalues(k) = int(x(k),y(k));
            end
          end
          
        end