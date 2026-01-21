function [] = step_22_summarize_deployment()
% function [] = step_22_summarize_deployment()
% 
% GEOMAR SVN $Id: step_22_summarize_deployment.m 557 2019-01-25 10:53:26Z gkrahmann@geomar.de $
%
% summarize the deployment
%
% requires :   processing_parameters.m in the current working folder
%
% version 5  last change 16.01.2019

% G.Krahmann, GEOMAR, Feb 2014

% catch non-GEOMAR case                            GK, 02.03.2017  1-->2
% add processing diary                             GK, 30.01.2018  2-->3
% extended diary                                   GK, 15.01.2019  3-->4
% adapt to new cond parameters                     GK, 16.01.2019  4-->5


%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 22  at  ',datestr(now)])
disp('summarize the deployment')
diary off



%
% look for and load processing parameters which contain the location of
% the flash card copy of the glider
%
if exist('./processing_parameters.m')
  processing_parameters
else
  error(['could not find file ''processing_parameters.m'' in current folder ',pwd])
end

%
% load data
%
load([op.deplname,'_yos.mat']);
load([op.deplname,'_dynamics.mat']);

diary deployment_report.txt
disp(' ')
disp(['Deployment : ',op.deplname])
disp(['Yos        : ',int2str(length(m_present_time))])

coulomb1 = nmin([m_coulomb_amphr{:}]);
coulomb2 = nmax([m_coulomb_amphr{:}]);
coul = coulomb2-coulomb1;
if isnan(coul)
  coul = '-';
else
  coul = int2str(round(coul));
end

tim1 = func_webbtime2mattime(nmin(m_present_time{1}));
tim2 = func_webbtime2mattime(nmax(m_present_time{end}));
dv1 = datestr(datevec(tim1),'DD.mm.YYYY');
dv2 = datestr(datevec(tim2),'DD.mm.YYYY');
dv1b = datestr(datevec(tim1),'YYYY-mm-DD');
dv2b = datestr(datevec(tim2),'YYYY-mm-DD');
 
dx = 0;
ca2 = 0;
ca1 = 1000000;
di2 = 0;
di1 = 1000000;
for n=1:length(m_present_time)
  la = n_lat{n};
  lo = n_lon{n};
  good = find(~isnan(la));
  if ~isempty(good)
    la1 = la(good(1));
    lo1 = lo(good(1));
    la2 = la(good(end));
    lo2 = lo(good(end));
    dx = dx + sqrt((la2-la1)^2 + cosd(la1)^2*(lo2-lo1)^2);
  end
end
dx = dx*1.852*60;
dx2 = 0;
for n=3:2:length(m_present_time)
  la = n_lat{n-2};
  lo = n_lon{n-2};
  good = find(~isnan(la));
  if ~isempty(good)
    la1 = la(good(1));
    lo1 = lo(good(1));
  else
    la1 = nan;
    lo1 = nan;
  end
  la = n_lat{n};
  lo = n_lon{n};
  good = find(~isnan(la));
  if ~isempty(good)
    la2 = la(good(1));
    lo2 = lo(good(1));
  else
    la2 = nan;
    lo2 = nan;
  end
  dx2 = dx2 + sqrt((la2-la1)^2 + cosd(la1)^2*(lo2-lo1)^2);
end
dx2 = dx2*1.852*60;
if dx2>dx
  dx = dx2;
end
dx = int2str(dx);

if ~exist('m_iridium_call_num')
  ca = 'x';
  di = 'x';
else
  ca2 = nmax([m_iridium_call_num{:}]);
  ca1 = nmin([m_iridium_call_num{:}]);
  di2 = nmax([m_iridium_dialed_num{:}]);
  di1 = nmin([m_iridium_dialed_num{:}]);
  ca = int2str(ca2-ca1);
  di = int2str(di2-di1);
end

disp(' ')
disp(' ')
disp('Deployment information: ')
disp(' ')
disp(['Days       : ',int2str(round(tim2-tim1))])
disp(['Energy     : ',coul,' Ahr'])
disp(['Distance   : ',dx,' km'])
disp(['Dialed     : ',di])
disp(['Call       : ',ca])
disp(' ')

disp(['| ',op.deplname(end+[-1,0]),' | ',dv1,' - ',dv2,' | ',int2str(round(tim2-tim1)),...
	' | ',int2str(length(m_present_time)),' | ',dx,' km | ',coul,' Ah | ',ca,' / ',di,' | ? |'])

disp(' ')
disp(' ')
disp('Flight model: ')
disp(' ')
if isempty(fixed_params)
  disp('no parameters were determined')
else
  disp(['Drag coeff : ',num2str(fixed_params.C_D_0(1)),' - ',num2str(fixed_params.C_D_0(end))])
  disp(['epsilon    : ',num2str(fixed_params.epsilon(1))])
  disp(['alpha_T    : ',num2str(fixed_params.alpha_T(1))])
  disp(['temp_filter: ',num2str(fixed_params.temperature_filter(1))])
  disp(['mass       : ',num2str(fixed_params.m_g(1)),' - ',num2str(fixed_params.m_g(end)) ])
  disp(['volume     : ',num2str(op.glider_volume(1)*1000)])
  disp(['density    : ',num2str(fixed_params.m_g(1)/op.glider_volume),' - ',...
    num2str(fixed_params.m_g(end)/op.glider_volume)])
  disp(['temperature: ',num2str(nmean([sci_water_temp{:}]))])
  disp(['front area : ',num2str(op.frontal_area)])
  disp(['wing area  : ',num2str(op.effective_wing_area)])
  disp(['air in oil : ',num2str(op.air_V)])
end


%
% get the pressure correction information
%
load([op.deplname,'_yos'],'scaling_factor','nav_pressure_offset','sci_pressure_offset',...
  'pressure_jump_correction');
disp(' ')
disp('Pressure correction:')
disp(' ')
disp(['scaling factor multiplied onto nav_pressure : ',num2str(scaling_factor)])
disp(['nav_pressure_offset added to nav_pressure   : ',num2str(nmean(nav_pressure_offset))])
disp(['sci_pressure_offset added to sci_pressure   : ',num2str(nmean(sci_pressure_offset))])
disp(['pressure_jump_correction                    : ',num2str(pressure_jump_correction)])


%
% get optode calibration information
%
geo_o_delay = 0;
geo_o_t_delay = 0;
if isfield(op,'force_optode_delay_o')
  disp('optode delays were set, not determined')
  disp(['applied o-delay : ',int2str(op.force_optode_delay_o),' secs'])
  disp(['applied t-delay : ',int2str(op.force_optode_delay_t),' secs'])
  aa_o_delay = nan;
  aa_o_t_delay = nan;
  geo_o_delay = op.force_optode_delay_o;
  geo_o_t_delay = op.force_optode_delay_t;
else
  load(['aanderaa_optode_delays']);
  disp(' ')
  disp('Aanderaa Optode calibration:')
  disp(' ')
  disp(['applied o-delay : ',int2str(best_delay_o),' secs'])
  disp(['applied t-delay : ',int2str(best_delay_t),' secs'])
  aa_o_delay = best_delay_o;
  aa_o_t_delay = best_delay_t;
  if exist(['geomar_optode_delays.mat']);
    load(['geomar_optode_delays']);
    disp(' ')
    disp('Geomar Optode calibration:')
    disp(' ')
    disp(['applied o-delay : ',int2str(best_delay_o),' secs'])
    disp(['applied t-delay : ',int2str(best_delay_t),' secs'])
    geo_o_delay = best_delay_o;
    geo_o_t_delay = best_delay_t;
  end
end
if isempty(op.file_with_geomar_optode_calibration_parameters)
  disp('using Aanderaa coefficients')
else
  disp('using Geomar coefficients')
end


%
% find the best set of optimization parameters
%
load([op.deplname,'_best_optimizer'])

disp(' ')
disp('Conductivity optimization:')
disp(' ')
disp(['best deviation : ',num2str(bestdevi)])
disp(['parameters     : ',num2str(res(1:5))])
disp(['flow type      : ',flow_type])
res1 = res;




disp(' ')
disp('Calibration offsets:')
disp(' ')
disp(['T  : ',num2str(op.t_offset)])
disp(['S  : ',num2str(op.s_offset)])
disp(['OA : ',num2str(op.oa_offset)])
disp(['OG : ',num2str(op.og_offset)])

diary off


%
% submit information to database
% this is GEOMAR stuff only
%
if ~exist('po_connect_db')
  return
end
po_connect_db;
if po_connection_flag==1
  res = po_query('glider_deployment',op.deplname);
  po_insert('glider',res.id,'number_of_yos',length(m_present_time));
  po_insert('glider',res.id,'days',round(tim2-tim1));
  po_insert('glider',res.id,'distance_in_km',dx);
  po_insert('glider',res.id,'date_start',dv1b);
  po_insert('glider',res.id,'date_end',dv2b);
  po_insert('glider',res.id,'drag_coeff_start',fixed_params.C_D_0(1));
  po_insert('glider',res.id,'drag_coeff_end',fixed_params.C_D_0(end));
  po_insert('glider',res.id,'compressibility_epsilon',fixed_params.epsilon(1));
  po_insert('glider',res.id,'thermal_expansion_alpha_t',fixed_params.alpha_T(1));
  po_insert('glider',res.id,'thermal_expansion_t_filter',fixed_params.temperature_filter(1));
  po_insert('glider',res.id,'mass_start',fixed_params.m_g(1));
  po_insert('glider',res.id,'mass_end',fixed_params.m_g(end));
  po_insert('glider',res.id,'volume',op.glider_volume(1)*1000);
  po_insert('glider',res.id,'density_start',fixed_params.m_g(1)/op.glider_volume*1000);
  po_insert('glider',res.id,'density_end',fixed_params.m_g(end)/op.glider_volume*1000);
  po_insert('glider',res.id,'at_temp',nmean([sci_water_temp{:}]));
  po_insert('glider',res.id,'front_area_in_m2',op.frontal_area);
  po_insert('glider',res.id,'wing_area_in_m2',op.effective_wing_area);
  po_insert('glider',res.id,'cc_air_oil',op.air_V);
  po_insert('glider',res.id,'nav_pressure_scale_factor',scaling_factor);
  po_insert('glider',res.id,'nav_pressure_bogus_step',pressure_jump_correction);
  po_insert('glider',res.id,'nav_pressure_offset',nmean(nav_pressure_offset));
  po_insert('glider',res.id,'sci_pressure_offset',nmean(sci_pressure_offset));
  po_insert('glider',res.id,'optode_delay_aanderaa',aa_o_delay);
  po_insert('glider',res.id,'optode_delay_geomar',geo_o_delay);
  po_insert('glider',res.id,'optode_temp_delay_aanderaa',aa_o_t_delay);
  po_insert('glider',res.id,'optode_temp_delay_geomar',geo_o_t_delay);
  po_insert('glider',res.id,'best_flow_type',flow_type);
  if isnan(res1(4))
    po_insert('glider',res.id,'cond_param_alpha',res1(1));
    po_insert('glider',res.id,'cond_param_tau',res1(2));
    po_insert('glider',res.id,'cond_param_time_offset',res1(3));
  else
    po_insert('glider',res.id,'cond_param_alpha',res1(1));
    po_insert('glider',res.id,'cond_param_slope_alpha',res1(2));
    po_insert('glider',res.id,'cond_param_tau',res1(3));
    po_insert('glider',res.id,'cond_param_slope_tau',res1(4));
    po_insert('glider',res.id,'cond_param_time_offset',res1(5));
  end
  po_insert('glider',res.id,'temp_offset',op.t_offset);
  po_insert('glider',res.id,'salt_offset',op.s_offset);
  po_insert('glider',res.id,'oxygen_offset_aanderaa',op.oa_offset);
  po_insert('glider',res.id,'oxygen_offset_geomar',op.og_offset);
else
  disp('Could not establish connection to PO database. No info submitted.')
end


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 22  at  ',datestr(now)])
disp(' ')
diary off

