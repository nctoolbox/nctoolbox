function setpassword(username, password)
% SETPASSWORD - Set credentials needed to access basic password protected URLS
%
% Sets the GLOBAL credentials used for authentication against servers. 
%
% Use As:
%    setpassword(username, password)
%
% Inputs:
%
%    username - The username for accessing a server
%    password = The password for the given username needed to access the server
%

% Brian Schlining
% 2012-01-11

credentialProvider = BasicCredentialsProvider(username, password);
ucar.nc2.util.net.HTTPSession.setGlobalCredentialsProvider(credentialProvider);

