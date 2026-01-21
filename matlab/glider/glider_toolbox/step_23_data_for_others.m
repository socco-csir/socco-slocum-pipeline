function [] = step_23_data_for_others()
% function [] = step_23_data_for_others()
% 
% GEOMAR SVN $Id: step_23_data_for_others.m 557 2019-01-25 10:53:26Z gkrahmann@geomar.de $
%
% extract some data sets for other purposes (e.g. Microrider)
%
% version 5  last change 15.01.2019

% G.Krahmann, GEOMAR,  Sep 2014

% create T,S  data for SUNA                        GK, 22.07.2015  1-->2
% cleanup SUNA part                                GK, 31.03.2016  2-->3
% add processing diary                             GK, 30.01.2018  3-->4
% extended diary                                   GK, 15.01.2019  4-->5


%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 23  at  ',datestr(now)])
disp('create files for others')
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
% create Microrider extra data
%
if isfield(op,'microrider_installed')
  if op.microrider_installed == 1

    %
    % load data
    %
    load([op.deplname,'_final_1sec'])
    load([op.deplname,'_1sec'],'pitch','roll','fin')
    load([op.deplname,'_1sec_derived'],'best_salinity_was_derived_with',...
      'glider_speed_in_movement_direction_dpdt',...
      'glider_speed_in_movement_direction_model','angle_of_attack_model','angle_of_attack')
    angle_of_attack_dpdt = angle_of_attack;
    if strcmp(best_salinity_was_derived_with,'dpdt')
      best_salinity_was_derived_with = 'dpdt';
    else
      best_salinity_was_derived_with = 'model';
    end
    load([op.deplname,'_dynamics.mat']);
    if isempty(fixed_params)
      angle_of_attack_model = [];
    end


    %
    % store data for MR
    %

    % ISABELLE EDIT 4/2/2025 - adding wetlabs variables
    save([op.deplname,'_for_mr_processing'],'temperature','time_datenum','roll','pitch','pressure',...
      'salinity','turbidity','oxygen','chlorophyll','bb470','bb700','longitude','latitude','fin',...
      'angle_of_attack_dpdt','best_salinity_was_derived_with','glider_speed_in_movement_direction_model',...
      'glider_speed_in_movement_direction_dpdt','angle_of_attack_model');

  end
end


%
% create SUNA T,S data
%
if isfield(op,'suna_installed')
  if op.suna_installed == 1
    fid = fopen([op.deplname,'_suna_ts.csv'],'wt');
    load([op.deplname,'_final_1sec'])
    for n=1:length(temperature)
      disp(n)
      d = time_datenum{n};
      s = salinity{n};
      t = temperature{n};
      if strcmp(op.deplname,'ifm13_depl01')
        d = d-365+1/24;
      end
      for m=1:length(s)
        if ~isnan(s(m)+t(m))
          fprintf(fid,[datestr(d(m),'yyyy-mm-dd HH:MM:SS'),',',...
            num2str(t(m),'%7.4f'),',',num2str(s(m),'%7.4f'),'\n'],[]);
        end
      end
    end
    fclose(fid);
  end
end


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 23  at  ',datestr(now)])
disp(' ')
diary off

