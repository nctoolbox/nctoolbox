
% Plot sections distance versus time of glider and model temp and salt
%
% load previously saved data file, or just
% call this script after calling "load_secoora_models.m"
% load_secoora_models
% load secoora_models.mat

[obs.server,obs.name] = fileparts(obs.url);
%obsname = strrep_(obs.name);
obsname = obs.name

switch variable
  case 'temp'
    bins = 22:0.25:30.;
    c = range(bins);
  case 'salt'
 %   bins = 29:0.2:36.5;  %MARACOOS
     bins = 35.5:.1:36.3; % SECOORA
    c = range(bins);
  otherwise
    c = 'auto';
end

%model_list = {'OBS','ESPreSSo','USEAST','HYCOM'}; %MARACOOS
model_list = {'OBS','SABGOM','USEAST','HYCOM'}; %SECOORA

clf
hax = nfigaxes([length(model_list) 1],[0 0.02],[0.1 0.95],[0.1 0.95]);

for m = 1:length(model_list)
  
  mname = char(model_list{m});
  eval(['model = ' lower(mname)])
  
  axes(hax(m))
  
  pcolorjw(model.(variable).dist,model.(variable).z,model.(variable).data)
  axis([0 200 -80 0]);
  %axis([0 300 -50 0]);
  caxis(c)
  colorbar
  
  title_pos=-60
  if m==1
    title(upper(variable))
    text(10,title_pos,{mname,(model.name)},'interpreter','none');
  else
    text(10,title_pos,mname);
  end
  
end

set(hax(1:(end-1)),'xticklabel',[])
xlabel('distance (km)')


