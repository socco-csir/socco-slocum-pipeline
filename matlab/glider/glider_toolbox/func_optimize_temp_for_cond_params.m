function [] = func_optimize_temp_for_cond_params(op)
% function [] = func_optimize_temp_for_cond_params(op)
% 
% GEOMAR SVN $Id: func_optimize_temp_for_cond_params.m 557 2019-01-25 10:53:26Z gkrahmann@geomar.de $
%
% optimize the parameters to estimate the temperature of the conductivity cell
%
% input  :  op                     - structure with options controlling the optimization
%
% version 5  last change 28.03.2018

% G.Krahmann, GEOMAR  Sep 2012

% add automatic selection of optimization toolbox            GK, 28.02.2017
% add automatic selection of area_function                   GK, 06.03.2017
% handle PC/UNIX and folder plots                            GK, 22.03.2017
% bug in mkdir call                                          GK, 28.03.2018

%
% load 1 second interpolated data
%
disp(['loading 1sec data for deployment ',op.deplname])
data = load([op.deplname,'_1sec']);
load([op.deplname,'_1sec_derived'],'glider_speed_in_cell_direction_dpdt');
data.glider_speed_in_cell_direction_dpdt = glider_speed_in_cell_direction_dpdt;
load([op.deplname,'_1sec_derived'],'glider_speed_in_cell_direction_model');
data.glider_speed_in_cell_direction_model = glider_speed_in_cell_direction_model;


%
% prepare data structures for optimization
%
[upcasts,downcasts,meanflow] = func_prep4calc(data,op);


%
% pick toolboxes
%
tv = ver('optim');
if isempty(tv)
  optimizer_mode = 2;
else
  optimizer_mode = 1;
end
mv = ver('map');
if isempty(mv)
  op.area_function = 2;
else
  op.area_function = 1;
end
if op.area_function==2
  disp(' ')
  disp('Using GEOMAR''s area function to calculate the minimization.')
else
  disp(' ')
  disp('Using Garau et al.''s area function to calculate the minimization.')
end


%
% Select the number of parameters to optimize and call the relevant optimization
% routines. To speed up the process and get closer to the final optimal values
% we start optimizing only 2 parameters and leave the other 2 a typical values.
% Next we add a third free parameter and last we optimize all 4 free parameters.
% This last optimization is repeated once as sometimes this gives an
% improved result.
%
% 2 free parameters
disp(' ')
disp('First guess for optimization #1 is:')
if op.is_pumped_ctd==1
  firstGuess = [0.015,30,0.3];
  disp(['Alpha-offset: ',num2str(firstGuess(1)),'  Tau-offset: ',num2str(firstGuess(2)),...
    '  Time-offset: ',num2str(firstGuess(3))])
else
  firstGuess = [0.05,0.05,8,12,0.01];
  disp(['Alpha-offset: ',num2str(firstGuess(1)),'  Alpha-slope: ',num2str(firstGuess(2)),...
    '  Tau-offset: ',num2str(firstGuess(3)),'  Tau-slope: ',num2str(firstGuess(4)),...
    '  Time-offset: ',num2str(firstGuess(5))])
end
disp(' ')
[res,devi] = func_optimizer_garau(downcasts,upcasts,firstGuess,optimizer_mode,op);
disp(' ')
disp('Optimization result #1 is:')
if op.is_pumped_ctd==1
  disp(['Alpha-offset: ',num2str(res(1)),'  Tau-offset: ',num2str(res(2)),...
    '  Time-offset: ',num2str(res(3))])
else
  disp(['Alpha-offset: ',num2str(res(1)),'  Alpha-slope: ',num2str(res(2)),...
    '  Tau-offset: ',num2str(res(3)),'  Tau-slope: ',num2str(res(4)),...
    '  Time-offset: ',num2str(res(5))])
end
disp(' ')

secondGuess = res;
[res,devi] = func_optimizer_garau(downcasts,upcasts,secondGuess,optimizer_mode,op);

disp(' ')
disp('Final result :')
if op.is_pumped_ctd==1
  disp(['Alpha-offset: ',num2str(res(1)),'  Tau-offset: ',num2str(res(2)),...
    '  Time-offset: ',num2str(res(3))])
else
  disp(['Alpha-offset: ',num2str(res(1)),'  Alpha-slope: ',num2str(res(2)),...
    '  Tau-offset: ',num2str(res(3)),'  Tau-slope: ',num2str(res(4)),...
    '  Time-offset: ',num2str(res(5))])
end
disp(' ')

fifthGuess = [res];
[res,devi] = func_optimizer_garau(downcasts,upcasts,fifthGuess,optimizer_mode,op);
res

%
% append results to cell array allres and save in a 
% mat file
%
op.res = res;
if exist([op.deplname,'_',op.use_flow_speed,'.mat'])
  load([op.deplname,'_',op.use_flow_speed,'.mat'])
  if ~exist('allres')
    allres{1} = op;
  else
    allres{end+1} = op;
  end
else
  allres{1} = op;
end
[v1,v2,v3,v4,s_field] = func_check_result(op);
allres{end}.devi = [v1,v2,v3,v4];
allres{end}.s_field = s_field;
allres{end}.ts_data_reduction = op.ts_data_reduction;
allres{end}.s_yo_numbers = op.s_yo_numbers;
allres{end}.subsample_profile = op.subsample_profile;
allres{end}.flow_speed_filter = op.flow_speed_filter;
allres{end}.flow = op.use_flow_speed;
allres{end}.area_function = op.area_function;
allres{end}.op = op;
op.devi = [v1,v2,v3,v4];
save([op.deplname,'_',op.use_flow_speed],'allres')

%
% display result
%
figure
func_plot_delay_function(op.res)
if ~exist('plots','dir')
  mkdir('plots')
end
print('-djpeg',['plots',filesep,'last_conductivity_delay_function_',op.use_flow_speed,'.jpg'])

