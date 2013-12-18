% ACCESS_TEST1 - comparison of different styles of access
echo on

url = 'http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/OS_M1_20081008_TS.nc'

% Object hierarchy: nc (root) <- cf <- geo
nc = ncdataset(url);     % Very basic access. Makes no assumptions about data
cf = cfdataset(url);     % Adds functionality. Assumes data is CF-compliant
geo = ncgeodataset(url); % Even more goodness. Assumes geographic oriented or CF-compliant data

% Basic data access (first, last)
nc.data('TEMP', [1 1 1 1], [10 2 1 1]) % works, 10x2
cf.data('TEMP', [1 1 1 1], [10 2 1 1])
geo.data('TEMP', [1 1 1 1], [10 2 1 1]) % works, 10x2

% variable access
% v1 = nc.variable('TEMP')  % Does not exist on ncdataset
v2 = cf.variable('TEMP')
v3 = geo.variable('TEMP')  % works, returns ncvariable

% Using ncvariable constructor. Should work for all ncdataset classes
n1 = ncvariable(nc, 'TEMP')% works: returns ncvariable
n2 = ncvariable(cf, 'TEMP') % works, returns ncvariable
n3 = ncvariable(geo, 'TEMP') % works, returns ncvariable

% Creating an ncgeovariable
% g1 = ncgeovariable(geo,'TEMP') % fails: not enough input arguments.
% g1 should never work!!
% g2 = geo.ncgeovariable('TEMP') % fails, no method
g3 = geo.geovariable('TEMP') % works, returns ncgeovariable  
% g4 = geo.ncvariable('TEMP') % fails, no method
% Use: g4 = geo.variable('TEMP') instead

g3(1:10,1:2) % nothing
g3.data(1:10,1:2) % 10x2 array
% g3.data([1 1 1 1],[10 2 1 1]) % ERROR: 1 value, strides must be positive and constant

n3(1:10,1:2)  % 10x2 array
% n3.data(1:10,1:2)  % ERROR: Index exceeds matrix dimensions.
n3.data([1 1 1 1],[10 2 1 1]) % 10x2 array

data(n3,[1 1 1 1],[10 2 1 1])  % 10x2
data(g3,[1 1 1 1],[10 2 1 1])  % 10x2
% data(g3,1:10,1:2)  % ERROR: data exceeds matrix dimensions
% data(n3,1:10,1:2)  % ERROR: data exceeds matrix dimensions

% ncgeodataset{varname} subsref access
g5 = geo{'TEMP'} % works, returns ncgeovariable
geo{'TEMP'}(1:10,1:2)       % 10x2
geo{'TEMP'}(1:10,1:2).data  % 10x2 (works with patch)
geo{'TEMP'}(1:10,1:2).grid  % structure with subsetted grid

echo off