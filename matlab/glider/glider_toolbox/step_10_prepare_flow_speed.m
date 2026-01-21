function [] = step_10_prepare_flow_speed()
% function [] = step_10_prepare_flow_speed()
% 
% GEOMAR SVN $Id: step_10_prepare_flow_speed.m 557 2019-01-25 10:53:26Z gkrahmann@geomar.de $
%
% prepare flow speed data for salinity calculation and filter pressure time series
%
% version 6  last change 15.01.2019

% G.Krahmann, GEOMAR  Aug 2012

% properly handle variable pitch                        GK, 19.02.2013  2-->3
% output measured_w   and removed again                 GK, 15.02.2016  3-->4
% add processing diary, remove v6 mat comp              GK, 30.01.2017  4-->5
% fix file header and extend diary                      GK, 15.01.2019  5-->6

%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 10  at  ',datestr(now)])
disp('derive some additional 1sec data')
diary off




%
% look for and load processing parameters
%
if exist('./processing_parameters.m') | exist('.\processing_parameters.m')
  processing_parameters
else
  error(['could not find file ''processing_parameters.m'' in current folder ',pwd])
end


%
% load 1 sec data
%
load([op.deplname,'_1sec'],'pressure','pitch','main_datenum','heading','latitude');
  

% despike and filter navigational pressure 
% this is done be running a lowpass filter, look for deviations
% and linearly interpolate a [-15:15] second interval around the
% spikes
% the result is called .p , the original .p_orig and
% the at the moment best pressure series .p_filtered
disp('filtering pressure')
[b,a] = butter(3,0.03,'low');
for n=1:length(pressure)
  fprintf(1,'.',[]);
  nav_p = pressure{n};
  good = find(~isnan(nav_p));
  nav_p_filtered = nav_p;
  nav_p_filtered(good) = filtfilt(b,a,nav_p(good));
  noise_p = nstd(nav_p-nav_p_filtered);
  for k=1:2
    bad = find(abs(nav_p-nav_p_filtered)>noise_p*3);
    dx = 15;
    if ~isempty(bad)
      for m=1:length(bad)
        if bad(m)>dx & bad(m)<length(nav_p)-dx
          indi = bad(m)+[-dx:dx];
          good = find(~isnan(nav_p(indi)));
          nav_p(indi) = interp1(indi(good),nav_p(indi(good)),indi);
        end
      end
    end
  end
  pressure_orig{n} = pressure{n};
  pressure{n} = nav_p;

  good = find(~isnan(nav_p));
  nav_p(good) = filtfilt(b,a,nav_p(good));
  pressure_filtered{n} = nav_p;

if 0
  figure(1)
  clf
  plot(pressure{n})
  hold on
  plot(pressure_filtered{n},'r')
  pause
end

end
fprintf(1,'\n',[]);

%
% calculate and store the flow speed through the conductivity cell
%
disp('estimating dpdt glider speed')
[glider_speed_in_glider_direction_dpdt,glider_speed_in_movement_direction_dpdt,...
  movement_direction,glider_direction,glider_speed_in_cell_direction_dpdt,...
  angle_of_attack_dpdt] =...
  func_ctd_cell_flow_speed_estimate(pressure_filtered,pitch,op,latitude);
for n=1:length(glider_speed_in_movement_direction_dpdt)
  glider_speed_in_movement_direction_model{n} = nan*glider_speed_in_movement_direction_dpdt{n};
  glider_speed_in_glider_direction_model{n} = nan*glider_speed_in_movement_direction_dpdt{n};
  glider_speed_in_cell_direction_model{n} = nan*glider_speed_in_movement_direction_dpdt{n};
end


%
% save the data
%
save([op.deplname,'_1sec_derived'],'pressure_filtered','glider_speed_in_movement_direction_dpdt',...
  'glider_speed_in_movement_direction_model','glider_speed_in_glider_direction_dpdt',...
  'glider_speed_in_glider_direction_model','glider_speed_in_cell_direction_dpdt',...
  'glider_speed_in_cell_direction_model','angle_of_attack_dpdt');


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 10  at  ',datestr(now)])
disp(' ')
diary off

