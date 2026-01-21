function [upcasts,downcasts,meanflow,data] = func_prep4calc(data,op);
% function [upcasts,downcasts,meanflow,data] = func_prep4calc(data,op);
% 
% GEOMAR SVN $Id: func_prep4calc.m 703 2020-07-28 14:50:02Z gkrahmann@geomar.de $
%
% rearrange data into up and downcasts
%
% version 3.1.0  last change 25.01.2023

% G.Krahmann, GEOMAR 2012

% catch extra short deployment                           GK, 29.05.2019  1-->2
% bug in s_yo_number handling                            GK, 10.03.2020  2-->3.0.0
% changed from CSIRO seawater to TEOS-10 library         GK< 25.01.2023  3.0.0-->3.1.0

%data = fix_problem_case(op.deplname,data);

disp(['Found ',int2str(length(data.main_datenum)/2),' profiles.'])

if isfield(op,'min_depth4cond_temp')
  if isempty(op.min_depth4cond_temp)
    min_depth4cond_temp = -10;
  else
    min_depth4cond_temp = op.min_depth4cond_temp;
  end
else
  min_depth4cond_temp = -10;
end

count = 1;
res = repmat(nan,[length(data),6]);
flows = [];
maxflows = [];
downcasts = struct([]);
upcasts = struct([]);

s_yo_numbers = [1+2:op.subsample_profile*2:length(data.main_datenum)];
if ~isempty(op.s_yo_numbers)
  for n=1:length(s_yo_numbers)
    if ~any(op.s_yo_numbers==s_yo_numbers(n))
      s_yo_numbers(n) = nan;
    end
  end
  s_yo_numbers = s_yo_numbers( find( ~isnan( s_yo_numbers ) ) );
end
if ~isempty(s_yo_numbers)
  if s_yo_numbers(1)==1
    s_yo_numbers = s_yo_numbers(3:end);
  end
end
for ind = s_yo_numbers

  clear dtd dtu

  %down.datenum = data.datenum{ind};
  %up = data(ind-1);

  downcasts(count).datenum = data.main_datenum{ind};
  downcasts(count).ctd_pressure = data.ctd_pressure{ind};
  downcasts(count).ctd_temperature = data.ctd_temperature{ind};
  downcasts(count).ctd_conductivity_ratio = data.ctd_conductivity{ind}/gsw_C3515;
  if strcmp(op.use_flow_speed,'dpdt')
    downcasts(count).flow_speed = func_cell_flow_speed_from_glider_speed(...
       data.glider_speed_in_cell_direction_dpdt{ind},op);
  elseif all(isnan([data.glider_speed_in_cell_direction_model{:}]))
    error('you have to run the dynamical model first')
  else
    downcasts(count).flow_speed = func_cell_flow_speed_from_glider_speed(...
       data.glider_speed_in_cell_direction_model{ind},op);
  end
  downcasts(count).pitch = data.pitch{ind};
  yonumbers(1,count) = ind;

  bad = find(downcasts(count).flow_speed<op.minimum_flow_speed);
  if ~isempty(bad)
    downcasts(count).flow_speed(bad) = op.minimum_flow_speed;
  end

  upcasts(count).datenum = data.main_datenum{ind-1};
  upcasts(count).ctd_pressure = data.ctd_pressure{ind-1};
  upcasts(count).ctd_temperature = data.ctd_temperature{ind-1};
  upcasts(count).ctd_conductivity_ratio = data.ctd_conductivity{ind-1}/gsw_C3515;
  if strcmp(op.use_flow_speed,'dpdt')
    upcasts(count).flow_speed = func_cell_flow_speed_from_glider_speed(...
       data.glider_speed_in_cell_direction_dpdt{ind-1},op);
  elseif all(isnan([data.glider_speed_in_cell_direction_model{:}]))
    error('you have to run the dynamical model first')
  else
    upcasts(count).flow_speed = func_cell_flow_speed_from_glider_speed(...
       data.glider_speed_in_cell_direction_model{ind-1},op);
  end
  upcasts(count).pitch = data.pitch{ind-1};
  yonumbers(2,count) = ind-1;

  bad = find(upcasts(count).flow_speed<op.minimum_flow_speed);
  if ~isempty(bad)
    upcasts(count).flow_speed(bad) = op.minimum_flow_speed;
  end

  % blank parts of profile that exist only in up or down parts
  depth_range = [nmax([1.5;nmin(upcasts(count).ctd_pressure);...
	nmin(downcasts(count).ctd_pressure)]),...
	nmin([nmax(upcasts(count).ctd_pressure);...
	nmax(downcasts(count).ctd_pressure)])];

  depth_range(1) = nmax([depth_range(1),min_depth4cond_temp]);

  bad_down1 = find(downcasts(count).ctd_pressure<depth_range(1));
  bad_down2 = find(downcasts(count).ctd_pressure>depth_range(2));
  if ~isempty(bad_down1)
    bad_down1 = [1:bad_down1(end)];
    downcasts(count).ctd_conductivity_ratio(bad_down1) = nan;
    downcasts(count).ctd_temperature(bad_down1) = nan;
    %downcasts(count).flow_speed(bad_down1) = nan;
    downcasts(count).ctd_pressure(bad_down1) = nan;
  end
  if ~isempty(bad_down2)
    bad_down2 = [bad_down2(1):length(downcasts(count).ctd_pressure)];
    downcasts(count).ctd_conductivity_ratio(bad_down2) = nan;
    downcasts(count).ctd_temperature(bad_down2) = nan;
    %downcasts(count).flow_speed(bad_down2) = nan;
    downcasts(count).ctd_pressure(bad_down2) = nan;
  end

  bad_up1 = find(upcasts(count).ctd_pressure<depth_range(1));
  bad_up2 = find(upcasts(count).ctd_pressure>depth_range(2));
  if ~isempty(bad_up1)
    bad_up1 = [bad_up1(1):length(upcasts(count).ctd_pressure)];
    upcasts(count).ctd_conductivity_ratio(bad_up1) = nan;
    upcasts(count).ctd_temperature(bad_up1) = nan;
    %upcasts(count).flow_speed(bad_up1) = nan;
    upcasts(count).ctd_pressure(bad_up1) = nan;
  end
  if ~isempty(bad_up2)
    bad_up2 = [1:bad_up2(end)];
    upcasts(count).ctd_conductivity_ratio(bad_up2) = nan;
    upcasts(count).ctd_temperature(bad_up2) = nan;
    %upcasts(count).flow_speed(bad_up2) = nan;
    upcasts(count).ctd_pressure(bad_up2) = nan;
  end

  flows = [flows;nmean(downcasts(count).flow_speed)];
  flows = [flows;nmean(upcasts(count).flow_speed)];
  maxflows = [maxflows;nmax(downcasts(count).flow_speed)];
  maxflows = [maxflows;nmax(upcasts(count).flow_speed)];

  count = count+1;
end

meanflow = nmean(flows);
maxflow = nmax(maxflows);

disp(['mean flow : ', num2str(meanflow)])
disp(['max flow  : ', num2str(maxflow)])

good = 0*[1:length(downcasts)];
for m=1:length(downcasts)
  if sum(~isnan(upcasts(m).ctd_conductivity_ratio))>100 &...
	 sum(~isnan(downcasts(m).ctd_conductivity_ratio))>100
    good(m) = 1;
  end
  if sum(~isnan(upcasts(m).ctd_conductivity_ratio)) ~= sum(~isnan(upcasts(m).ctd_temperature))
    good(m) = 0;
  end
  if sum(~isnan(downcasts(m).ctd_conductivity_ratio)) ~= sum(~isnan(downcasts(m).ctd_temperature))
    good(m) = 0;
  end
end
clear downcast2 upcast2
if any(good==0)
  count = 1;
  for m=1:length(good)
    if good(m)==1
      downcast2(count) = downcasts(m);
      upcast2(count) = upcasts(m);
      %disp(['updown yos ',int2str(yonumbers(:,m)')])
      count = count+1;
    end
  end
  if exist('downcast2')
    downcasts = downcast2;
    upcasts = upcast2;
  else
    downcasts = struct([]);
    upcasts = struct([]);
  end
end

disp(' ')
disp(['Optimization preparation found ',int2str(length(downcasts)),' profiles with up AND down data'])

