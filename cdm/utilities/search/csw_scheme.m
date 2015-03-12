function [s]=csw_scheme(records,scheme)
% CSW_SCHEME parse records returned by csw_search for scheme
% Usage: [s]=CSW_SCHEME(records,scheme)
% Inputs: records = records returned by csw_search
%         scheme = data endpoint scheme to find 
% Example:
% scheme = 'urn:x-esri:specification:ServiceType:odp:url';
%% scheme = 'OPeNDAP:OPeNDAP';
% [s]=csw_scheme(records,scheme)
% s{1}
% see also: CSW_SEARCH,TEST_CSW_SEARCH
s={};
k=0;
for i=1:length(records)
    record_xml=char(records{i});
    service_template='<dct:references scheme="DATA_SCHEME">(.*?)</dct:references>';
    scheme_expr=strrep(service_template,'DATA_SCHEME',scheme);
    [~,~,~,~,data_url]=regexp(record_xml,scheme_expr,'match');
    if ~isempty(data_url),
        [~,~,~,~,data_title]=regexp(record_xml,'<dc:title>(.*?)</dc:title>','match');
        k=k+1;
        s{k}.title=char(data_title{1});
        s{k}.scheme=scheme;
        s{k}.url=char(data_url{1});
        
    end
end