function [] = step_11_optimize_temp_for_cond_params_dpdt(no_optim)
% function [] = step_11_optimize_temp_for_cond_params_dpdt(no_optim)
% 
% GEOMAR SVN $Id: step_11_optimize_temp_for_cond_params_dpdt.m 557 2019-01-25 10:53:26Z gkrahmann@geomar.de $
%
% optimize the parameters to estimate the temperature of the conductivity cell
% to accomplish that run several different attempts of the optimization
%
% version 7.1.0  last change 25.01.2023

% G.Krahmann, GEOMAR  Aug 2012

% modify minimum subsample_profile numbers             GK, 21.12.2012  3-->4
% allow to change ts_data_reduction through processing_parameter.m file
%                                                      GK, 18.03.2013  4-->5
% add processing diary, remove v6 mat comp             GK, 30.01.2018  5-->6
% extended diary                                       GK, 15.01.2019  6-->7
% changed path to Garau's functions                    GK, 25.01.2023  7-->7.1.0

%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 11  at  ',datestr(now)])
disp('optimize parameters for salinity calculation based on dp/dt')
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
% set some paths
%
external_functions = [fileparts(which('step_11_optimize_temp_for_cond_params_dpdt')),filesep,'external'];
addpath(external_functions)


%
% check whether the parameters are set, or whether we need to optimize them
%
if ~isfield(op,'force_garau_params')

  %
  % set that we are doing the dpdt optimization
  %
  op.use_flow_speed = 'dpdt';


  %
  % call the rest of the script
  %
  if no_optim==0
    func_optimization_dummy_script;
  end

else
 
  disp('using preset parameters for conductivity cell') 
  op.use_flow_speed = 'dpdt';

end


%
% calculate the optimal salinity
%
func_derive_temp_for_cond_and_get_salinity(op)


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 11  at  ',datestr(now)])
disp(' ')
diary off

