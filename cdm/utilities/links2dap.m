function dap=links2dap(links);
% LINKS2DAP extract OpenDAP Data URLs from OpenSearch links
% Usage: dap=links2dap(links);
%    [links,struc]=opensearch(q);
%    dap=links2dap(links);

% Rich Signell (rsignell@usgs.gov)
k=0;
dap=[];
for i=1:length(links)
    s=char(links{i});
    if ~isempty(regexpi(s,'dodsC')),
        k=k+1;
        e=regexpi(s,'\.html');
        dap{k}=s(1:e-1);
    end
end