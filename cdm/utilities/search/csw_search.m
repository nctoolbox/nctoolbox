function [records,paramString,s_input]=csw_search(s_input)
% CSW_SEARCH finds data web services for datasets in specified query
% query is limited here to bbox,start,stop,any_text
% and all these must be specified.
% enter bbox=[] to omit bounding box
%
% Usage:
% [records,paramString,s]=csw_search(s)
% Inputs: s.endpoint = csw endpoint url (string)
%         s.bbox = bounding box [lon_min lon_max lat_min lat_max]
%         s.start = start time (string or datenum)
%         s.stop = stop time (string or datenum)
%         s.any_text = text string to look for (e.g. 'sea_water_salinity')
%
% Example:
%s.endpoint = 'http://www.ngdc.noaa.gov/geoportal/csw';
%s.bbox = [-75.0 -71.0 39.0 41.0];
%s.start = '2014-11-12 18:00';
%s.stop = '2014-11-18 18:00';
%s.any_text = 'sea_water_salinity';
%records = csw_search(s);
% see also: CSW_SCHEME, TEST_CSW_SEARCH

% this script loads a mat file containing a template of the XML request
% that we previously saved thusly:
% >> paramString=fileread('simple_csw_request.xml');
% >> save paramString simple_csw_request.mat

csw_endpoint=s_input.endpoint;
bbox = s_input.bbox;
if ~exist('s.max_records'),
    s_input.max_records=1000;
end
max_records = s_input.max_records;
any_text = s_input.any_text;
start=datestr(s_input.start,'yyyy-mm-ddTHH:MM:SSZ');
stop=datestr(s_input.stop,'yyyy-mm-ddTHH:MM:SSZ');

if ~isempty(bbox),
    load('simple_csw_request.mat','paramString');
    paramString = strrep(paramString,'LON_MIN',num2str(bbox(1),'%f'));
    paramString = strrep(paramString,'LON_MAX',num2str(bbox(2),'%f'));
    paramString = strrep(paramString,'LAT_MIN',num2str(bbox(3),'%f'));
    paramString = strrep(paramString,'LAT_MAX',num2str(bbox(4),'%f'));
else
    load('simpler_csw_request.mat','paramString');
end
paramString = strrep(paramString,'ANY_TEXT_STRING',any_text);
paramString = strrep(paramString,'TIME_START',start);
paramString = strrep(paramString,'TIME_STOP',stop);
paramString = strrep(paramString,'MAX_RECORDS',num2str(max_records));


[output,extras] = urlread2(csw_endpoint,'POST',paramString);

% split on records</csw:Record>
[~,~,~,~,records]=regexp(output,'<csw:Record>(.*?)</csw:Record>','match');



