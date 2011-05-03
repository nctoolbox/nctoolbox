        
        function int = getinterp(data, lon, lat, method)
          % NCGEOVARIABLE.INTERPTOXY - Interpolate gridded data at an x/y (lon/lat) location or set of locat-
          % ions when inputs are vectors of the same length.
          % Useage: >> data = interptoxy(data, datalon, datalat, -77, 45, 'linear');
          s = size(data);
%           if length(x) ~= length(y)
%             error('NCTOOLBOX:interptoxy','Input coordinate vectors do not agree.')
%           end
          
          
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
          lat = double(lat);
          lon = double(lon);
%           x = double(x);
%           y = double(y);
          
          if length(s) > 2
            for j = 1:s(1)
              
              if length(s) > 3
                for i=1:s(2)
                  redata = reshape(data(j,i,:,:), numel(data(j,i,:,:)), 1);
                  relat = lat;
                  relon = lon;
                  flags1 = isnan(redata) | isnan(relat) | isnan(relon);
                  redata(flags1) = [];
                  relat(flags1) = [];
                  relon(flags1) = [];
                  int = TriScatteredInterp(relon, relat, redata, method);
                  clear redata
%                   if isvector(x)
%                     for k = 1:length(x)
%                       interpolatedvalues(j,i,k) = int(x(k),y(k));
%                     end
%                   else
%                     num = size(x);
%                     for k = 1:num(1)
%                       for l = 1:num(2)
%                         interpolatedvalues(j, i, k, l) = int(x(k, l), y(k, l));
%                       end
%                     end
%                   end
                end
                
              else
                redata = reshape(data(j,:,:), numel(data(j,:,:)), 1);
                relat = lat;
                relon = lon;
                flags1 = isnan(redata) | isnan(relat) | isnan(relon);
                redata(flags1) = [];
                relat(flags1) = [];
                relon(flags1) = [];
                int = TriScatteredInterp(relon, relat, redata, method);
                clear redata
%                 if isvector(x)
%                   for k = 1:length(x)
% %                     interpolatedvalues(j,k) = int(x(k),y(k));
%                   end
%                 else
%                   num = size(x);
%                   for k = 1:num(1)
%                     for l = 1:num(2)
% %                       interpolatedvalues(j, k, l) = int(x(k, l), y(k, l));
%                     end
%                   end
%                 end
              end

            end
            
          else
            redata = reshape(data, numel(data), 1);
            relat = lat;
            relon = lon;
            flags1 = isnan(redata) | isnan(relat) | isnan(relon);
            redata(flags1) = [];
            relat(flags1) = [];
            relon(flags1) = [];
            int = TriScatteredInterp(relon, relat, redata, method);
            clear redata
%             if isvector(x)
%               for k = 1:length(x)
% %                 interpolatedvalues(k) = int(x(k),y(k));
%               end
%               num = size(x);
%               for k = 1:num(1)
%                 for l = 1:num(2)
% %                   interpolatedvalues(k, l) = int(x(k, l), y(k, l));
%                 end
%               end
%             end
          end
          
        end