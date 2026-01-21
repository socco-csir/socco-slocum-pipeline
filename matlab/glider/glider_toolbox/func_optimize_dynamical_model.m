function [all_res,yo_number,result] = func_optimize_dynamical_model(data,op)
% function [] = func_optimize_dynamical_model()
% 
% GEOMAR SVN $Id: func_optimize_dynamical_model.m 703 2020-07-28 14:50:02Z gkrahmann@geomar.de $
%
% optimize parameters that model the gliders dynamical flight behaviour
%
% requires :   processing_parameters.m in the current working folder
%
% version 4.1.0  last change 24.03.2020

% G.Krahmann, GEOMAR, Aug 2012

% references: Merckelbach, Lucas, David Smeed, Gwyn Griffiths, 2010: Vertical Water 
% Velocities from Underwater Gliders. J. Atmos. Oceanic Technol., 27, 547–563. 

% only use chunks which have a sufficient amount of down AND up data
%                                                      G.Krahmann, 20.2.2014  1-->2
% output final difference result                       GK, 03.01.2019  2-->3
% optimize_dynamical_model_n_yos_together              GK, 23.01.2019  3-->4
% catch no oil_vol data when pitch value exists,
% never happened before.                               GK, 24.03.2020  4-->4.1.0

global val plotnew handle_model_graph w w_glider_model nnn facs fixed_params param_list no_plot

no_plot = op.no_plot;


param_list = {'C_D_0','epsilon','m_g','alpha_T','temperature_filter'};
start_par = [];
facs = [];
fixed_params
for n=1:length(param_list)
  if isempty(getfield(fixed_params,param_list{n}))
    if strcmp('C_D_0',param_list{n})
      start_par = [start_par,0.11];
      facs = [facs,1];
    elseif strcmp('epsilon',param_list{n})
      start_par = [start_par,0.5e-9];
      facs = [facs,1e9];
    elseif strcmp('m_g',param_list{n})
      start_par = [start_par,data.V_g*1024];
      facs = [facs,0.01];
    elseif strcmp('alpha_T',param_list{n})
      start_par = [start_par,1e-6];
      facs = [facs,100000];
    elseif strcmp('temperature_filter',param_list{n})
      start_par = [start_par,0.005];
      facs = [facs,100];
    end
  end
end

 
% this parameter sets how many yos should be optimized in one step
% this should be more than 2 to even out distortions by internal waves 
% anything from 4 to 10 appears to be reasonable.
% If too large, the optimization will slow down and can not account for time changes
% in drag coefficient or glider mass.
if ~isfield(op,'optimize_dynamical_model_n_yos_together')
  nn = 18;
else
  nn = op.optimize_dynamical_model_n_yos_together;
end
plotnew = 1;
count = 1;
all_res = [];
for n=[1:nn:floor(length(data.pressure)/nn)*nn]

  nnn = n;
  
  val.pitch = [data.pitch{n+[0:nn-1]}]/180*pi; 
  val.w = [data.w{n+[0:nn-1]}];
  val.rho = [data.rho{n+[0:nn-1]}];
  val.P = [data.pressure{n+[0:nn-1]}];
  val.T = [data.filtered_ctd_temperature{n+[0:nn-1]}];
  val.salinity = [data.ctd_salinity{n+[0:nn-1]}];
  val.delta_V_bp = [data.delta_V_bp{n+[0:nn-1]}];
  val.battery_position = [data.battery_position{n+[0:nn-1]}];
  val.time = [data.time{n+[0:nn-1]}];
  val.fin = [data.fin{n+[0:nn-1]}];
  val.shear_lift = [data.shear_lift{n+[0:nn-1]}];
  val.V_g = data.V_g;
  val.S = data.S;
  val.A_h = data.A_h;
  val.latitude = nmean([data.latitude{n+[0:nn-1]}]);
  val.yo_index = [];
  val.thruster_power = [data.thruster_power{n+[0:nn-1]}];
  for m=1:nn
    val.yo_index = [val.yo_index,m*ones(1,length(data.w{n+m-1}))];
  end
  val.max_yo_index = nn;

  % define different excluding criteria
  bad = zeros(1,length(val.time));

  % pitch smaller 15 deg
  bad(find(abs(val.pitch)<15/180*pi)) = 1;

  % slow vertical speed
  bad(find(abs(val.w)<0.05)) = 1;

  % no good science data
  bad(find(isnan([val.rho]))) = 1;

  % no good pitch angle
  bad(find(isnan([val.pitch]))) = 1;

  % no good w
  bad(find(isnan([val.w]))) = 1;

  % no good w
  bad(find(isnan([val.delta_V_bp]))) = 1;

  % thruster is on
  bad(find(val.thruster_power>0)) = 1;

  % 60 seconds after a moving battery
  dummy = find(abs(diff([val.battery_position]))>0.01);
  for m=1:length(dummy)
    bad(find(val.time-val.time(dummy(m))>0 & val.time-val.time(dummy(m))<60/86400)) = 1;
  end

  % 60 seconds after changing volume
  dummy = find(abs(diff([val.delta_V_bp]))>1);
  for m=1:length(dummy)
    bad(find(val.time-val.time(dummy(m))>0 & val.time-val.time(dummy(m))<60/86400)) = 1;
  end

  % upper seven dbar
  dummy = find(val.P<7);
  for m=1:length(dummy)
    bad(dummy) = 1;
  end

  % lower ten dbar
  pmax = nmax(val.P);
  dummy = find(pmax-val.P<10);
  for m=1:length(dummy)
    bad(dummy) = 1;
  end 

  % determine remaining good values
  bad = find(bad==1);
  bad = bad(find(bad<=length(val.P)));
  good = [1:length(val.P)];
  good(bad) = nan;
  good = good(find(~isnan(good)));

  updown = 0;
  if any(val.w(good)>0) & any(val.w(good)<0)
    length_up = length(find(val.w(good)>0));
    length_dn = length(find(val.w(good)<0));
    if length_up>0.25*length(good) & length_dn>0.25*length(good)
      updown = 1;
    end
  end

  if length(good)>1 & updown==1
    val.pitch = val.pitch(good); 
    val.w = val.w(good);
    val.salinity = val.salinity(good);

    val.P = val.P(good)*10^4; 
    val.delta_V_bp = val.delta_V_bp(good)/(1*10^6);
    val.rho = val.rho(good);
    val.T = val.T(good);
    val.good = good;
    val.time = val.time(good);
    val.fin = val.fin(good);
    val.battery_position = val.battery_position(good);
    val.yo_index = val.yo_index(good);
    val.shear_lift = val.shear_lift(good);

    % use data at the beginning and end of a yo multiple times
if 0
    yo = unique(val.yo_index);
    all_ind = [];
    for m=yo
      ind = find(m==val.yo_index);
      if length(ind)>300
        all_ind = [all_ind,ind(1:300),ind,ind(end-300:end)];
      end
   end
   val.pitch = val.pitch(all_ind);
   val.w = val.w(all_ind);
   val.salinity = val.salinity(all_ind);
   val.P = val.P(all_ind);
   val.delta_V_bp = val.delta_V_bp(all_ind);
   val.rho = val.rho(all_ind);
   val.good = val.good(all_ind);
   val.time = val.time(all_ind);
   val.fin = val.fin(all_ind);
   val.battery_position = val.battery_position(all_ind);
   val.yo_index = val.yo_index(all_ind);
   val.shear_lift = val.shear_lift(all_ind);
end

    if no_plot~=1
      func_sfigure(1);
      clf
      plot(val(1).w,'g');
      hold on
      handle_model_graph = plot(val(1).w);
      xlabel('Time [sec]')
      ylabel('Vertical velocity [m/s]')
      title(['Dive #',int2str((n-1)/2+1),' to #',int2str((n-1)/2+1+nn/2-1),'       g: dp/dt  b: model'])
    end

    opt = optimset('tolfun',1e-6,'tolx',1e-6);
    res = fminsearch('func_dynamic_model_deviation',start_par.*facs,opt);
    res = fminsearch('func_dynamic_model_deviation',res,opt);

    result = func_dynamic_model_deviation(res);
    if n~=0
      all_res(count,:) = res./facs;
      yo_number(count) = n;
      count = count+1;
    else
      disp('all data:')
    end
    disp(' ')
    disp(['optimization result for dive #',int2str((n-1)/2+1),' to #',int2str((n-1)/2+1+nn/2-1),' :'])
    disp(res./facs)
    drawnow update

  end
end
