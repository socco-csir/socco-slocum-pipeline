function [] = func_derive_temp_for_cond_and_get_salinity(op)
% function [] = func_derive_temp_for_cond_and_get_salinity(op)
% 
% GEOMAR SVN $Id: func_derive_temp_for_cond_and_get_salinity.m 714 2020-09-02 11:25:52Z gkrahmann@geomar.de $
%
% derive a temperature for the conductivity cell
% and afterwards calculate salinity
%
% version 4.3.0  last change 25.01.2023

% G.Krahmann, GEOMAR  Aug 2012

% G.Krahmann:
% This is modified from freely available code downloaded from http://www.socib.es/;glider/doco/gliderToolbox/ctdTools/thermalLagTools
% in 2011. See 'Availability of the code' statement in DOI:10.1175/JTECH-D-10-05030.1

% use CTD pressure instead of filtered navigational          GK, 12.12.2012  2-->3
% handle unpumped/pumped res difference                      GK, 04.01.2019  3-->4
% added possibility for multiple force_tomeu_params          GK, 02.09.2020  4-->4.1.0
% fix bug with multiple                                      GK, 26.11.2020  4.1.0-->4.1.1
% added code source statement                                GK, 23.08.2022  4.1.1-->4.2.0
% changed tomeu to garau                                     
% changed from CSIRO seawater to TEOS-10 library             GK, 25.01.2023  4.2.0-->4.3.0 

%
% find the best set of optimization parameters 
% or a set preset parameters
%
if isfield(op,'force_garau_params')
  if size(op.force_garau_params,1)==1
    res = op.force_garau_params;
    multiple_res = [];
  else
    res = [];
    multiple_res = op.force_garau_params;
  end
  flow_type = 'dpdt';
else
  load([op.deplname,'_',op.use_flow_speed])
  res = repmat(nan,[1,6]);
  bestdevi = 1e6;
  for n=1:length(allres)
    disp(' ')
    disp(allres{n}.res)
    disp(allres{n}.devi)
    disp(allres{n}.devi(2)+abs(allres{n}.devi(1))+allres{n}.devi(4))
    %if allres{n}.devi(2)<bestdevi
    if allres{n}.devi(2)+abs(allres{n}.devi(1))+allres{n}.devi(4)<bestdevi
      res = allres{n}.res;
      %bestdevi = allres{n}.devi(2);
      bestdevi = allres{n}.devi(2)+abs(allres{n}.devi(1))+allres{n}.devi(4);
      flow_type = allres{n}.flow;
      op = allres{n}.op;
      if ~isfield(op,'area_function')
        op.area_function = 1;
      end
    end
  end
  multiple_res = [];
end
disp(' ')
disp(['parameters     : ',num2str(res)])
disp(['flow type      : ',flow_type])


%
% load 1 second interpolated data
%
disp(['loading 1sec data for deployment ',op.deplname])
data = load([op.deplname,'_1sec']);
data2 = load([op.deplname,'_1sec_derived']);
if strcmp(op.use_flow_speed,'dpdt')
  data.flow_speed = func_cell_flow_speed_from_glider_speed(...
    data2.glider_speed_in_cell_direction_dpdt,op);
  flow = '_dpdt';
else
  data.flow_speed = func_cell_flow_speed_from_glider_speed(...
    data2.glider_speed_in_cell_direction_model,op);
  flow = '_model';
end
disp(['using cell-flow-speed : ',flow])
disp(['found ',int2str(length(data.main_datenum)),' yos'])


%
% apply the best coefficients to calculate the temperature in the conductivity cell
%
for m=1:length(data.main_datenum)
  good = find(~isnan(data.ctd_temperature{m}+data.flow_speed{m}));
  data.tempInCell{m} = nan*data.ctd_temperature{m};
  data.ctd_salinity{m} = nan*data.ctd_temperature{m};
  if op.is_pumped_ctd==1
    if isempty(multiple_res)
      data.tempInCell{m}(good) = func_correct_thermal_lag_garau(data.main_datenum{m}(good),...
        data.ctd_temperature{m}(good),data.flow_speed{m}(good),res(1:3),op.area_function);
    else
      ind = find(multiple_res(:,1)<=m & multiple_res(:,2)>=m);
      data.tempInCell{m}(good) = func_correct_thermal_lag_garau(data.main_datenum{m}(good),...
        data.ctd_temperature{m}(good),data.flow_speed{m}(good),multiple_res(ind,3:5),op.area_function);
    end
  else
    data.tempInCell{m}(good) = func_correct_thermal_lag_garau(data.main_datenum{m}(good),...
      data.ctd_temperature{m}(good),data.flow_speed{m}(good),res(1:5),op.area_function);
  end
  data.ctd_salinity{m}(good) = gsw_SP_from_R(data.ctd_conductivity{m}(good)/gsw_C3515,data.tempInCell{m}(good),...
    data.ctd_pressure{m}(good));
end

eval(['data2.tempInCell',flow,'=data.tempInCell;'])
eval(['data2.ctd_salinity',flow,'=data.ctd_salinity;'])
save([op.deplname,'_1sec_derived'],'-struct','data2')
