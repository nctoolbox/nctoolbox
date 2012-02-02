classdef complex_ongrid < handle
    
    properties (SetAccess = private)
        uvar
        vvar
        rhovar
        anglevar
    end % end properties
    
    
    methods
        function obj = complex_ongrid(rhovar, ncgeovariable_u, ...
                ncgeovariable_v, ncgeovariable_angle)
            obj.uvar = ncgeovariable_u;
            obj.vvar = ncgeovariable_v;
            if nargin < 4
                obj.anglevar = -1;
            else
                obj.anglevar = ncgeovariable_angle;
            end
            obj.rhovar = rhovar;
        end % end constructor
        
        function mag = magnitude(obj, first, last, stride)
            vec = obj.vectors(first, last, stride);
            mag = abs(vec);
        end % end mag
        
        function vec = vectors(obj, first, last, stride)
            % get lon,lat size from hvar
            if nargin < 2
                hsize = obj.rhovar.size;
                horder = obj.rhovar.getaxesorder;
                if horder.lat == horder.lon
                    jj = 1:hsize(horder.lat);
                    ii = 1:hsize(horder.lat  + 1);
                else
                    %         ii = hsize(order.lat);
                    %         jj = hsize(order.lon);
                    error('order.lat != order.lon');
                end
            else
                jj = first(end-1):last(end-1);
                ii = first(end):last(end);
            end
            ujj = (jj(1) + 1):(jj(end) - 1);
            uii = ii(1):ii(end) - 1;
            vjj = jj(1):jj(end) - 1;
            vii = (ii(1) + 1):(ii(end) - 1);
            
            u = double(squeeze(obj.uvar.data(first(1), first(2), ujj, uii))); % get u
            v = double(squeeze(obj.vvar.data(first(1), first(2), vjj, vii))); % get v
            
            g = obj.rhovar.grid_interop(jj, ii); % get lon,lat from rho-point variable (like 'h' in ROMS)
            
            U = ones(size(g.lon)) * nan; % template for U at rho points
            U(2:end-1, 2:end-1) = complex(av2(u.').', av2(v)); %average u,v to rho
            
            if obj.anglevar ~= -1
                ang = obj.anglevar.data(jj, ii); % get angle
                U = U .* exp(sqrt(-1) * ang); % rotate
            end
            
            if nargin < 4
                stride = ones(length(first));
            end
            
            vec = U(1:stride(end-1):end, 1:stride(end):end); % output
            
        end % end vec
        
        function g = grid(obj, first, last, stride)
            if nargin < 2
                hsize = obj.rhovar.size;
                horder = obj.rhovar.getaxesorder;
                if horder.lat == horder.lon
                    jj = 1:hsize(horder.lat);
                    ii = 1:hsize(horder.lat  + 1);
                else
                    %         ii = hsize(order.lat);
                    %         jj = hsize(order.lon);
                    error('order.lat != order.lon');
                end
            else
                if nargin < 4
                    stride = ones(length(first));
                end
                jj = first(end-1):stride(end-1):last(end-1);
                ii = first(end):stride(end-1):last(end);
            end
            g = obj.rhovar.grid_interop(jj, ii); % get lon,lat from rho-point variable (like 'h' in ROMS)
            
            % output
            tim = obj.uvar.timewindowij( double(first(1)) );
            g.time = tim.time;
            g.klevel = first(2);
            g.itime = first(1);
            
        end % end grid
        
        function e = end(obj, k, n)
            usize = size(obj.uvar);
            hsize = size(obj.rhovar);
            n = [usize(1:2) hsize(end-1:end)];
            e = n(k);
        end % Added to deal with end indexing functionality,
        % otherwise the indexing arugment is ignored.
        
        function sref = subsref(obj,s)
            switch s(1).type
                % Use the built-in subsref for dot notation
                case '.'
                    usize = size(obj.uvar);
                    hsize = size(obj.rhovar);
                    nums = [usize(1:2) hsize(end-1:end)];
                    switch s(1).subs
                        case 'grid'
                            
                            if ~isempty(nums)
                                switch length(s)
                                    case 1
                                        sref = obj;
                                    case 2
                                        [first last stride] = indexing(s(2).subs, double(nums));
                                        sref = obj.grid(first, last, stride);
                                end
                                
                            else
                                sref = obj.data;
                                warning('NCTOOLBOX:ncvariable:subsref', ...
                                    ['Variable "' obj.name '" has no netcdf dimension associated with it. Errors may result from non CF compliant files.'])
                            end
                        case 'vectors'
                            %                             nums = size(obj.rhovar);
                            if ~isempty(nums)
                                switch length(s)
                                    case 1
                                        sref = obj;
                                    case 2
                                        [first last stride] = indexing(s(2).subs, double(nums));
                                        sref = obj.vectors(first, last, stride);
                                end
                                
                            else
                                warning('NCTOOLBOX:ncvariable:subsref', ...
                                    ['Variable "' name '" has no netcdf dimension associated with it. Errors may result from non CF compliant files.'])
                            end
                        case 'magnitude'
                            %                             nums = size(obj.rhovar);
                            if ~isempty(nums)
                                switch length(s)
                                    case 1
                                        sref = obj;
                                    case 2
                                        [first last stride] = indexing(s(2).subs, double(nums));
                                        sref = obj.magnitude(first, last, stride);
                                end
                                
                            else
                                warning('NCTOOLBOX:ncvariable:subsref', ...
                                    ['Variable "' name '" has no netcdf dimension associated with it. Errors may result from non CF compliant files.'])
                            end
                        otherwise
                            sref = builtin('subsref',obj,s);
                    end
                case '()'
                    if length(s)<2
                        % Note that obj.Data is passed to subsref
                        sref = builtin('subsref',obj.data,s);
                    else
                        sref = builtin('subsref',obj,s);
                    end
                    % No support for indexing using '{}'
                case '{}'
                    error('NCTOOLBOX:ncvariable:subsref', ...
                        'Not a supported subscripted reference, "{}" are not permitted to call variable object methods');
            end
        end
    end % end methods
    
end % end classdef