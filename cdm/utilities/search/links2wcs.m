function wcs=links2wcs(links);
% LINKS2WCS extract WCS Data URLs from OpenSearch links
% Usage: wcs=links2wcs(links);
%    [links,struc]=opensearch(q);
%    wcs=links2wcs(links);
% Example wcs GetCapabilities request:
%    url=[wcs{i} '?service=wcs&version=1.1.0&request=GetCapabilities'];
% Rich Signell (rsignell@usgs.gov)
k=0;
wcs=[];
for i=1:length(links)
    s=char(links{i});
    if ~isempty(regexpi(s,'wcs')),
        k=k+1;
        wcs{k}=s;
    end
end