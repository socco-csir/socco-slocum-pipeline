function [] = step_16_pick_model_or_dpdt()
% function [] = step_16_pick_model_or_dpdt()
% 
% GEOMAR SVN $Id: step_16_pick_model_or_dpdt.m 714 2020-09-02 11:25:52Z gkrahmann@geomar.de $
%
% select the best salinity result, model or dpdt derived
%
% version 5.1.0  last change 25.01.2023

% G.Krahmann, GEOMAR,  May 2014

% correct handling of forced parameters                      GK, 06.03.2017  1-->2
% add processing diary, remove v6 mat comp                   GK, 30.01.2018  2-->3
% for a pumped CTD always use the dpdt optimization          GK, 14.01.2019  3-->4
% extended diary                                             GK, 15.01.2019  4-->5
% added a condition to handle force_tomeu_params             GK, 02.09.2020  5-->5.0.1
% changed tomeu to garau                                     GK, 25.01.2023  5.0.1-->5.1.0

%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 16  at  ',datestr(now)])
disp('pick the best salinity calculation')
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
% check whether the conductivity cell parameters are preset
% in which case we always use dpdt
%
if ~isfield(op,'force_garau_params') & op.is_pumped_ctd~=1

  %
  % load results for model and dpdt flowspeeds and pick the best one
  %
  load([op.deplname,'_model']);
  res_m = repmat(nan,[1,6]);
  bestdevi_m = 1e6;
  disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  disp('model based salinity optimizations:') 
  for n=1:length(allres)
    disp(' ')
    disp(allres{n}.res)
    disp(allres{n}.devi)
    disp(allres{n}.devi(2)+abs(allres{n}.devi(1))+allres{n}.devi(4))
    if allres{n}.devi(2)+abs(allres{n}.devi(1))+allres{n}.devi(4)<bestdevi_m
      res_m(1:length(allres{n}.res)) = allres{n}.res;
      bestdevi_m = allres{n}.devi(2)+abs(allres{n}.devi(1))+allres{n}.devi(4);
    end
  end
  disp('best model based one is:')
  disp(' ')
  disp(['parameters     : ',num2str(res_m(1:4))])
  disp(' ')
  
  res_p = repmat(nan,[1,6]);
  bestdevi_p = 1e6;
  disp('-------------------------------------------------------------------')
  disp('dpdt based salinity optimizations:') 
  load([op.deplname,'_dpdt']);
  for n=1:length(allres)
    disp(' ')
    disp(allres{n}.res)
    disp(allres{n}.devi)
    disp(allres{n}.devi(2)+abs(allres{n}.devi(1))+allres{n}.devi(4))
    if allres{n}.devi(2)+abs(allres{n}.devi(1))+allres{n}.devi(4)<bestdevi_p
      res_p(1:length(allres{n}.res)) = allres{n}.res;
      bestdevi_p = allres{n}.devi(2)+abs(allres{n}.devi(1))+allres{n}.devi(4);
    end
  end
  disp(' ')
  disp('best dpdt based one is:')
  disp(' ')
  disp(['parameters     : ',num2str(res_p(1:4))])
  disp(' ')
  disp('-------------------------------------------------------------------')
  disp(' ')
 
  if bestdevi_p>bestdevi_m 
    disp('model based salinity is better than dpdt based one')
    disp(' ')
    disp(['parameters     : ',num2str(res_m(1:4))])
    flow_type = '_model';
    bestdevi = bestdevi_m;
    res = res_m;
  else
    disp('dpdt based salinity is better than model based one')
    disp(' ')
    disp(['parameters     : ',num2str(res_p(1:4))])
    flow_type = '_dpdt';
    bestdevi = bestdevi_p;
    res = res_p;
  end

elseif ~isfield(op,'force_garau_params') & op.is_pumped_ctd==1

  flow_type = '_dpdt';
  res = repmat(nan,[1,6]);
  bestdevi = 1e6;
  disp('-------------------------------------------------------------------')
  disp('dpdt based salinity optimizations:') 
  load([op.deplname,'_dpdt']);
  for n=1:length(allres)
    disp(' ')
    disp(allres{n}.res)
    disp(allres{n}.devi)
    disp(allres{n}.devi(2)+abs(allres{n}.devi(1))+allres{n}.devi(4))
    if allres{n}.devi(2)+abs(allres{n}.devi(1))+allres{n}.devi(4)<bestdevi
      res(1:length(allres{n}.res)) = allres{n}.res;
      bestdevi = allres{n}.devi(2)+abs(allres{n}.devi(1))+allres{n}.devi(4);
    end
  end
  disp(' ')
  disp('best dpdt based one is:')
  disp(' ')
  disp(['parameters     : ',num2str(res(1:4))])
  disp(' ')
  disp('-------------------------------------------------------------------')
  disp(' ')

else

  if isfield(op,'force_garau_dpdt_or_model')
    flow_type = ['_',op.force_garau_dpdt_or_model];
  else
    flow_type = '_dpdt';
  end
  res = op.force_garau_params;
  bestdevi = nan;

end

%
% select the better one
%
disp(['using flow_type : ',flow_type])
data = load([op.deplname,'_1sec_derived']);
if strcmp(flow_type,'_model') & ~isfield(op,'force_garau_params')
  data.ctd_salinity = data.ctd_salinity_model;
  data.angle_of_attack = data.angle_of_attack_model;
else
  data.ctd_salinity = data.ctd_salinity_dpdt;
  data.angle_of_attack = data.angle_of_attack_dpdt;
end
data.best_salinity_was_derived_with = flow_type(2:end);
save([op.deplname,'_1sec_derived'],'-struct','data');
save([op.deplname,'_best_optimizer'],'res','flow_type','bestdevi');


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 16  at  ',datestr(now)])
disp(' ')
diary off

