function [ links, params, requesturl] = opensearch( q )
%OPENSEARCH Queries GI-cat using the opensearch interface
% returns a list of available links from the matching resources
% Usage:  [ links, params ] = opensearch(q )
%
%   q.endpoint = URL to opensearch server
%   q.string_text = text string to match 
%   q.bbox= 'xmin,ymin,xmax,ymax'
%   q.time_start = start time string
%   q.time_end  = ending time string
% 
% For full options and more explanation, see the GI-CAT OpenSearch doc
%    http://essi-lab.eu/cgi-bin/twiki/view/GIcat/OpenSearchGuide
%
% Example 1:
% q.endpoint='http://geoport.whoi.edu/gi-cat/services/opensearch';
% q.string_text='sea_water_temperature';
% q.bbox='-74,39,-66,45';
% q.time_start='1999-01-20T04:33:45Z';
% q.time_end='1999-08-01T00:00:00Z';
% [links,params]=opensearch(q);   % make the query
% dap=links2dap(links); % find only the OPeNDAP links
%
% Example 2:
% q.endpoint='http://testbedapps.sura.org/gi-cat/services/opensearch';
% q.string_text='sea_water_salinity';
% q.bbox='-82,36,-73,40';
% q.time_start='2000-06-20T04:33:45Z';
% q.time_end='2006-08-01T00:00:00Z';
% [links,params]=opensearch(q);   % make the query
% dap=links2dap(links); % find only the OPeNDAP links
% 
% Example 3:
% q.endpoint='http://geoport.whoi.edu/gi-cat/services/opensearch';
% q.bbox='-180,0,180,90'; % northern hemisphere
% start=now_utc-7;  %start one week ago
% stop=now_utc;     %stop now
% q.time_start=datestr(start,'yyyy-mm-ddTHH:MM:SSZ');% convert to ISO
% q.time_end=datestr(stop,'yyyy-mm-ddTHH:MM:SSZ');% convert to ISO
% q.string_text='relative_humidity'; %search for datasets containing 'relative_humidity'
% [links,params]=opensearch(q);   % make the query
% dap=links2dap(links); % find only the OPeNDAP links

if(nargin<1);help opensearch; return;end
if ~isfield(q,'endpoint'); disp('must define endpoint (URL)'); return; end
if ~isfield(q,'string_text'); q.string_text='';end
if ~isfield(q,'bbox'); q.bbox='';end
if ~isfield(q,'time_start'); q.time_start='';end
if ~isfield(q,'time_end'); q.time_end='';end
if ~isfield(q,'loc'); q.loc='';end
if ~isfield(q,'si'); q.si='1';end

% Default "&ct=&" returns only 10 records so specify default to be very 
% large if not supplied.  
if ~isfield(q,'ct'); q.ct='400';end

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

[startIndex, endIndex, tokIndex, matchStr, links, exprNames,splitStr] = ...
    regexp(results, '<gmd:URL>(.*?)</gmd:URL>','match');

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

end



