
%variable = 'salt';
variable = 'temp';


% will use dynamic field names with structures to process variable named
% above


% -----------------------------------------------------------------------
% Obs: Single glider transect
%obs.url = ['http://tds.marine.rutgers.edu/thredds/dodsC/' ...
%    'cool/glider/mab/Gridded/20101025T000000_20101117T000000_marcoos_ru22.nc']
%obs.url=['http://tds.marine.rutgers.edu/thredds/dodsC/' ...
%    'cool/glider/mab/Gridded/20100402T000000_20100511T000000_cara_ru26d.nc']
obs.url = ['http://tds.marine.rutgers.edu:8080/thredds/dodsC/' ...
    'cool/glider/mab/Gridded/20130911T000000_20130920T000000_gp2013_modena.nc']
%obs.url = ['http://tds.marine.rutgers.edu/thredds/dodsC/' ...
%    'cool/glider/mab/Gridded/20130813T000000_20130826T000000_njdep_ru28.nc']
obs.file = obs.url;
obs.temp.name = 'temperature';
obs.salt.name = 'salinity';
obs.lonname  = 'longitude';
obs.latname  = 'latitude';
obs.zname    = 'depth';

nc = ncgeodataset(obs.url)
% Load the observations
disp(['Loading obs from ' obs.url])
disp(['  Variable is ' obs.(variable).name])
data     = nc{obs.(variable).name}(:);
obs.data = double(data);
obs.lon  = nc{obs.lonname}(:);
obs.lat  = nc{obs.latname}(:);
obs.z    = nc{obs.zname}(:);
obs.dist = cumsum([0; sw_dist(obs.lat,obs.lon,'km')]);
obs.(variable).data = data;
obs.(variable).dist = obs.dist;
obs.(variable).z    = obs.z;
obs.time = nj_time(nc,'temperature');
tstart = min(obs.time);
tend = max(obs.time);
disp('  Time interval of obs:')
disp(['    ' datestr(tstart) ' to ' datestr(tend)])

% -----------------------------------------------------------------------


% -----------------------------------------------------------------------
% Model: Global NCOM CF-compliant aggregation
ncom.name = 'global_ncom';
ncom.url = ['http://ecowatch.ncddc.noaa.gov/thredds/dodsC/' ...
    'ncom/ncom_reg1_agg/NCOM_Region_1_Aggregation_best.ncd'];
ncom.file = 'ncom.nc';
ncom.temp.name = 'water_temp';
ncom.salt.name = 'salinity';

% -----------------------------------------------------------------------
% Model: US-EAST (NCOM) CF-compliant aggregation
useast.name = 'useast';
useast.url = ['http://ecowatch.ncddc.noaa.gov/thredds/dodsC/' ...
    'ncom_us_east_agg/US_East_Apr_05_2013_to_Current_best.ncd']
useast.file = 'useast.nc';
useast.temp.name = 'water_temp';
useast.salt.name = 'salinity';
% -----------------------------------------------------------------------
% Model: MERCATOR CF-compliant nc file extracted at myocean.eu
% Registration and user/password authentication required
mercator.name = 'mercator';
mercator.url = ['./dataset-psy2v4-pgs-nat-myocean-bestestimate_' ...
    '1295878772263.nc'];
mercator.file = mercator.url;
mercator.temp.name = 'temperature'; % <<< in KELVIN !!!!!!!!!!!!!!!!!!!!!
mercator.salt.name = 'salinity';

% -----------------------------------------------------------------------
% Model: COAWST CF-compliant ROMS aggregation
coawst.name = 'coawst';
% before june 25, 2012
if tend <= datenum(2012,6,25),
    coawst.url = ['http://geoport.whoi.edu/thredds/dodsC/' ...
    'coawst_2_2/fmrc/coawst_2_2_best.ncd'];
else% after june 25, 2012
    coawst.url = ['http://geoport.whoi.edu/thredds/dodsC/' ...
     'coawst_4/use/fmrc/coawst_4_use_best.ncd'];
end
coawst.file = 'coawst.nc';
coawst.temp.name = 'temp';
coawst.salt.name = 'salt';

% -----------------------------------------------------------------------
% Model: ESPreSSO CF-compliant ROMS aggregation
espresso.name = 'espresso';
espresso.url = ['http://tds.marine.rutgers.edu:8080/thredds/' ...
    'dodsC/roms/espresso/2009_da/his'];
espresso.file = 'espresso.nc';
espresso.temp.name = 'temp';
espresso.salt.name = 'salt';

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% Model: SABGOM CF-compliant ROMS aggregation
sabgom.name = 'sabgom';
sabgom.url = ['http://omgsrv1.meas.ncsu.edu:8080/thredds/dodsC/' ...
    'fmrc/sabgom/SABGOM_Forecast_Model_Run_Collection_best.ncd'];
sabgom.file = 'sabgom.nc';
sabgom.temp.name = 'temp';
sabgom.salt.name = 'salt';
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% Model: AMSEAS CF-compliant NCOM aggregation
amseas.name = 'amseas';
amseas.url = ['http://edac-dap3.northerngulfinstitute.org/thredds/dodsC/' ...
    'ncom_amseas_agg/AmSeas_Aggregation_best.ncd'];
amseas.file = 'amseas.nc';
amseas.temp.name = 'water_temp';
amseas.salt.name = 'salinity';
% -----------------------------------------------------------------------
% Model: Global HYCOM RTOFS (HYCOM) Region 1
hycom.name = 'hycom';
hycom.url = ['http://ecowatch.ncddc.noaa.gov/thredds/dodsC/' ...
    'hycom/hycom_reg1_agg/HYCOM_Region_1_Aggregation_best.ncd'];
hycom.file = 'hycom.nc';
hycom.temp.name = 'water_temp';
hycom.salt.name = 'salinity';

%% models to compare with data
%model_list = {'USEAST','ESPreSSO','HYCOM'};  %MARACOOS
model_list = {'USEAST','SABGOM','HYCOM'};     %SECOORA


ncks = 0;

for m = 1:length(model_list)
    
    tic
    
    mname = char(model_list{m});
    
    % work with a temporary structure named 'model'
    
    eval(['model = ' lower(mname)])
    
    if ncks
        str = nc_genslice(model.url,model.(variable).name,...,
            obs.lon,obs.lat,obs.time,'ncks');
        disp([str ' ' model.name '.nc'])
        return
    end
    
    [Tvar,Tdis,Tzed] = nc_genslice(model.url,model.(variable).name,...
        obs.lon,obs.lat,obs.time,'verbose');
    
    if ~isempty(findstr(model.url,'myocean')) && strcmp(variable,'temp')
        Tvar = Tvar - 272.15;
    end
    model.(variable).data = Tvar;
    model.(variable).dist = Tdis;
    model.(variable).z = Tzed;
    
    % copy 'model' back to the oroginal named strucutre for this model
    eval([model.name ' = model;'])
    tocs(m)=toc;
    disp('----------------------------------------------------------------')
    disp(['  Elapsed time processing ' mname])
    disp(['  was ' num2str(toc,3) ' seconds'])
    disp('----------------------------------------------------------------')
    
end
disp('----------------------------------------------------------------')
disp(['  Total Elapsed Time' ])
disp(['  was ' num2str(sum(tocs),3) ' seconds'])
disp('----------------------------------------------------------------')

clear nc
save secoora_models.mat
