%% TEST 1 (NGDC Geoportal)
s.bbox=[-75.0 -71.0 39.0 41.0]; % new york region
%s.start='2014-11-12 18:00';
%s.stop='2014-11-18 18:00';
s.max_records=10
s.start= now-7
s.stop = now
s.any_text='sea_water_salinity';
s.endpoint = 'http://www.ngdc.noaa.gov/geoportal/csw';
[records]=csw_search(s);
scheme='urn:x-esri:specification:ServiceType:odp:url';
s2 = csw_scheme(records,scheme);
length(s2)
s2{1}
%% SURA Testbed CSW test (pycsw)
% find water level results for hurricane ike
s.endpoint = 'http://comt.sura.org:8000';
s.bbox=[-100.0 -97.0 27.0 30.0]; %near the mississippi delta
s.start='2008-09-01 00:00';
s.stop='2008-09-15 00:00';
s.max_records=100
%s.start= now-3
%s.stop = now
s.any_text='surface_height';
[records]=csw_search(s);
scheme='OPeNDAP:OPeNDAP';
s2 = csw_scheme(records,scheme);
length(s2)
s2{1}


