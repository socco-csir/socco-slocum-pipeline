% this is a script that is being used twice
% by the steps    11 and 15
% 
% GEOMAR SVN $Id: func_optimization_dummy_script.m 303 2017-03-06 08:52:37Z gkrahmann@geomar.de $
%
% version 1  last change 19.05.2014

% G.Krahmann, GEOMAR  Aug 2012

%
% load the data and figure out how many up/down pairs we have
%
load([op.deplname,'_1sec'],'ctd_temperature');
for n=1:length(ctd_temperature)
  mean_temp(n) = nmean(ctd_temperature{n});
end
updown_pairs = mean_temp(1:2:end)+mean_temp(2:2:end);
updown_pairs = length(find(~isnan(updown_pairs)));


ts_data_reduction = [20,10];

if ~isfield(op,'subsample_profile')
  subsample_profile = [ceil(updown_pairs/30),ceil(updown_pairs/60),ceil(updown_pairs/100)];
  if subsample_profile(2)==1
    subsample_profile(2) = 2;
    ts_data_reduction = [10,5];
  end
  if subsample_profile(3)==1
    subsample_profile(3) = 3;
    ts_data_reduction = [5,3];
  end
elseif isempty(op.subsample_profile)
  subsample_profile = [ceil(updown_pairs/30),ceil(updown_pairs/60),ceil(updown_pairs/100)];
  if subsample_profile(2)==1
    subsample_profile(2) = 2;
    ts_data_reduction = [10,5];
  end
  if subsample_profile(3)==1
    subsample_profile(3) = 3;
    ts_data_reduction = [5,3];
  end
else
  subsample_profile = op.subsample_profile;
end

if isfield(op,'ts_data_reduction')
  ts_data_reduction = op.ts_data_reduction;
end

if ~isfield(op,'flow_speed_filter')
  op.flow_speed_filter = 30;
end

for m=1:length(ts_data_reduction)
  for k=1:length(subsample_profile)
    disp(' ')
    disp(['Optimizing with ',int2str(ts_data_reduction(m)),' data point and ',...
      int2str(subsample_profile(k)),' profile reduction'])
    disp(' ')
    op.ts_data_reduction = ts_data_reduction(m);
    op.subsample_profile = subsample_profile(k);
     
    func_optimize_temp_for_cond_params(op);

  end
end
