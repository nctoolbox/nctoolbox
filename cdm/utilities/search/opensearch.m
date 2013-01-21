function [ links, params, requesturl, box] = opensearch( q )
%OPENSEARCH Queries GI-cat using the OpenSearch interface
% Returns a list of available links from the matching resources
%
% Usage:  [ links, params, requesturl ] = opensearch( q )
%
%   q.endpoint = URL to opensearch server
%   q.string_text = text string to match 
%   q.bbox= bounding box, vector [lon_min lon_max lat_min lat_max]
%       or string 'lon_min,lat_min,lon_max,lat_max' 
%   q.time_start = start time, datenum/datestr or ISO string
%   q.time_end  = ending time, datenum/datestr or ISO string
% 
% For full options and more explanation, see the GI-CAT OpenSearch doc
%    http://essi-lab.eu/cgi-bin/twiki/view/GIcat/OpenSearchGuide
%
% Example 1:
% q.endpoint='http://geoport.whoi.edu/gi-cat/services/opensearch';
% q.string_text='sea_water_temperature';
% q.bbox=[-74.0 -66.0 39.0 45.];
% q.time_start='1999-01-20 04:33:45';
% q.time_end='1999-08-01 00:00:00';
% [links,params]=opensearch(q);   % make the query
% dap=links2dap(links); % find only the OPeNDAP links
%
% Example 2:
% q.endpoint='http://testbedapps.sura.org/gi-cat/services/opensearch';
% q.string_text='sea_water_salinity';
% q.bbox=[-82 -73 36 40];
% q.time_start=[2000 06 20 04 33 45];
% q.time_end=[2006 08 01 00 00 00];
% [links,params]=opensearch(q);   % make the query
% dap=links2dap(links); % find only the OPeNDAP links
% 
% Example 3:
% q.endpoint='http://geoport.whoi.edu/gi-cat/services/opensearch';
% q.bbox=[-180 180 -30 70];  % whole globe from -30 to 70 north
% q.time_start=now_utc-28; % start 28 days ago
% q.time_end=now_utc;    % stop now
% q.string_text='relative_humidity'; %search for datasets containing 'relative_humidity'
% [links,params]=opensearch(q);   % make the query
% dap=links2dap(links); % find only the OPeNDAP links

if(nargin<1);help opensearch; return;end
if ~isfield(q,'endpoint'); disp('must define endpoint (URL)'); return; end
if ~isfield(q,'string_text'); q.string_text='';end
if ~isfield(q,'bbox')
  q.bbox='';
else
  if length(q.bbox)==4, % matlab style bbox [lon_min lon_max lat_min lat_max]
    q.bbox=sprintf('%g, %g, %g, %g',q.bbox([1 3 2 4]));
  end
end
if ~isfield(q,'time_start')
  q.time_start='';
else
  try 
    dn=datenum(q.time_start); % handle datenum/datestr input
    q.time_start=datestr(dn,'yyyy-mm-ddTHH:MM:SSZ');
  catch
  end 
end
if ~isfield(q,'time_end')
  q.time_end='';
else
  try 
    dn=datenum(q.time_end);
    q.time_end=datestr(dn,'yyyy-mm-ddTHH:MM:SSZ');
  catch
  end
end
if ~isfield(q,'loc'); q.loc='';end
if ~isfield(q,'si'); q.si='1';end

% Default "&ct=&" returns only 10 records so specify default to be very 
% large if not supplied.  
if ~isfield(q,'ct'); q.ct='1000';end

% In GI-CAT 8.4, default "&rel=&" for rel:{geo:relation} returns 
% records outside bounding box, so here we specify the default 
% to be 'overlaps' if not supplied.  Options are
% {disjoint|overlaps|contains}
if ~isfield(q,'rel'); q.rel='overlaps';end  

params = {
'si';   q.si;
'ct';   q.ct; 
'rel';  q.rel;   
'st';   q.string_text;
'bbox'; q.bbox;
'loc';  q.loc;
'ts';   q.time_start;
'te';   q.time_end; 
'outputFormat';'application/atom+xml';
'enableCrawling'; 'false';
};

[results] = urlread(q.endpoint,'get', params);

% find all URL links
[~,~,~,~,links1]=regexp(results,'<gmd:URL>(.*?)</gmd:URL>','match'); % early versions of nciso/gicat
[~,~,~,~,links2]=regexp(results,'<dm:srvOp xmlns:dm="http://floraresearch.eu/sdi/services/7.0/dataModel/schema" name="(.*?)</dm:srvOp>','match'); % modern versions of nciso/gicat
for ind = 1:length(links2)
   link = links2(ind);
   st = regexp(link{1}, 'http');
   links2{ind}{1}(1:st{1}-1) = [];
end

% find all box values (bounding box)
[~,~,~,~,box]=regexp(results,'<box xmlns="http://www.georss.org/georss">(.*?)</box>','match');  

% Optionally return the request URL used for the search
requesturl = [q.endpoint, '?&', params{1},'=',params{2},'&',...
    params{3},'=',params{4},'&',...    
    params{5},'=',params{6},'&',...
    params{7},'=',params{8},'&',...
    params{9},'=',params{10},'&',...
    params{11},'=',params{12},'&',...
    params{13},'=',params{14},'&',...
    params{15},'=',params{16},'&',...
    params{17},'=',params{18},'&',...
    params{19},'=',params{20}];

if ~isempty(links2)
    links = links2;
else
    links = links1;
end

end




