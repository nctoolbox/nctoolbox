% DEMO6 - Gratuitis demo showing AUV path

echo('on') 
% Starting DEMO6 ----------------------------------------------------------
% Gratuitis demo showing AUV path

nav = ncdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/ssds/generated/netcdf/files/ssds.shore.mbari.org/auvbi/missionlogs/2009/2009173/2009.173.03/navigation.nc');
y = nav.data('mPos_y');
x = nav.data('mPos_x');
z = nav.data('mDepth');
plot3(x, y, z, '.')
set(gca, 'ZDir', 'reverse')
axis('equal')
grid('on')
xlabel('meters')
ylabel('meters')
zlabel('meters')
echo('off') % Ending DEMO6 ------------------------------------------------
