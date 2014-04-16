function setproxy(proxyname, port)
% SETPROXY - Set the proxy to use when making all remote calls with
% NCTOOLBOX (https://github.com/nctoolbox/nctoolbox)
%
% Use as: setproxy(proxyname, port)
%
% Inputs:
%  proxyname = the proxy name that all connections will be routed to
%  port      = The port on the proxy to connect through
%

% Brian Schlining
% 2010-02-16
%
% httpClient = org.apache.commons.httpclient.HttpClient();
% httpClient.getHostConfiguration().setProxy(proxyname, port);
% ucar.nc2.dataset.NetcdfDataset.setHttpClient(httpClient);

%ucar.nc2.util.net.HTTPSession.setGlobalSimpleProxy(proxyname, port)
ucar.nc2.util.net.HTTPSession.setGlobalProxy(proxyname, port)
