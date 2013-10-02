
% Plot sections distance versus time of glider and model temp and salt

[obs.server,obs.name] = fileparts(obs.url);
%obsname = strrep_(obs.name);
obsname = obs.name

switch variable
  case 'temp'
    bins = 20:0.25:30.5;
    c = range(bins);
  case 'salt'
 %   bins = 29:0.2:36.5;
     bins = 35.5:.1:36.3;
    c = range(bins);
  otherwise
    c = 'auto';
end

%model_list = {'OBS','NCOM','Mercator','COAWST','HyCOM','ESPreSSO'};
%model_list = {'OBS','COAWST','ESPreSSO'};
%model_list = {'OBS','ESPreSSO','NCOM','HyCOM'};
%model_list = {'OBS','ESPreSSO','NCOM'}
model_list = {'OBS','USEAST','SABGOM','RTOFS'}

clf
hax = nfigaxes([length(model_list) 1],[0 0.02],[0.1 0.95],[0.1 0.95]);

for m = 1:length(model_list)
  
  mname = char(model_list{m});
  eval(['model = ' lower(mname)])
  
  axes(hax(m))
  
  pcolorjw(model.(variable).dist,model.(variable).z,model.(variable).data)
  axis([0 200 -100 0])
  caxis(c)
  colorbar
  
  if m==1
    title(upper(variable))
    text(10,-90,{mname,(model.name)});
  else
    text(10,-80,mname);
  end
  
end

set(hax(1:(end-1)),'xticklabel',[])
xlabel('distance (km)')


