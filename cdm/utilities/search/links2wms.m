function wms=links2wms(links);
% LINKS2WMS extract WMS Data URLs from OpenSearch links
% Usage: wms=links2wms(links);
%    [links,struc]=opensearch(q);
%    wms=links2wms(links);
% Example WMS GetCapabilities request:
%    url=[wms{i} '?service=WMS&version=1.3.0&request=GetCapabilities'];
% Rich Signell (rsignell@usgs.gov)
k=0;
wms=[];
for i=1:length(links)
    s=char(links{i});
    if ~isempty(regexpi(s,'wms')),
        k=k+1;
        wms{k}=s;
    end
end