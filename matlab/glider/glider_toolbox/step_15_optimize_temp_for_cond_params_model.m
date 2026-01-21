function [] = step_15_optimize_temp_for_cond_params_model(no_optim)
% function [] = step_15_optimize_temp_for_cond_params_model(no_optim)
% 
% GEOMAR SVN $Id: step_15_optimize_temp_for_cond_params_model.m 557 2019-01-25 10:53:26Z gkrahmann@geomar.de $
%
% optimize the parameters to estimate the temperature of the conductivity cell
% to accomplish that run several different attempts of the optimization
%
% version 7  last change 15.01.2019

% G.Krahmann, GEOMAR  Aug 2012

% modify minimum subsample_profile numbers             GK, 21.12.2012  3-->4
% allow to change ts_data_reduction through processing_parameter.m file
%                                                      GK, 18.03.2013  4-->5
% add processing diary                                 GK, 30.01.2018  5-->6
% extended diary                                       GK, 15.01.2019  6-->7



%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 15  at  ',datestr(now)])
disp('optimize parameters for salinity calculation based on the flight model')
diary off


%
% allow for an argument to suppress the optimization
%
if nargin==1
  no_optim = 1;
else
  no_optim = 0;
end



%
% look for and load processing parameters
%
if exist('./processing_parameters.m') | exist('.\processing_parameters.m')
  processing_parameters
else
  error(['could not find file ''processing_parameters.m'' in current folder ',pwd])
end


%
% check whether the parameters have been preset
%
if ~isfield(op,'force_garau_params') & op.is_pumped_ctd~=1

  %
  % set that we are doing the model optimization
  %
  op.use_flow_speed = 'model';


  %
  % call the rest of the script
  %
  if no_optim==0
    func_optimization_dummy_script;
  end

elseif op.is_pumped_ctd==1
  
  disp('skipping flow speed optimization because of pumped CTD')
  op.use_flow_speed = 'dpdt';

else

  disp('using preset parameters for conductivity cell')
  op.use_flow_speed = 'model';

end


%
% calculate the optimal salinity
%
func_derive_temp_for_cond_and_get_salinity(op)


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 15  at  ',datestr(now)])
disp(' ')
diary off

