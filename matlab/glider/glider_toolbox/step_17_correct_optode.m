function [] = step_17_correct_optode(no_optim)
% function [] = step_17_correct_optode(no_optim)
% 
% GEOMAR SVN $Id: step_17_correct_optode.m 703 2020-07-28 14:50:02Z gkrahmann@geomar.de $
%
% First run a reverse filter on the oxygen data to undo the time shift of the
% slow sensor. Then do a 0-phase filter to smooth the now noisy data but have it
% at the correct times.
%
% version 7.1.0  last change 21.07.2020

% G.Krahmann, GEOMAR  sep 2012

% fixed bug in handling different calculations                   GK, 15.03.2017  2-->3
% add processing diary, remove v6 mat comp                       GK, 30.01.2018  3-->4
% run calculation with previous results when choosing no_optim   GK, 31.01.2018  4-->5
% extended diary                                                 GK, 15.01.2019  5-->6
% better error message for forcing optode delays                 GK, 29.05.2019  6-->7
% handle NaN in optode calibration coefficients                  GK, 21.07.2020  7-->7.1.0

%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 17  at  ',datestr(now)])
disp('find optimal optode parameters')
diary off


%
% input arguments
%
if nargin<1
  no_optim = 0;
else
  no_optim = 1;
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
optode_functions = [fileparts(which('step_17_correct_optode')),filesep,'optode_functions'];
addpath(optode_functions)
optode_calibrations = [fileparts(which('step_17_correct_optode')),filesep,'optode_calibrations'];
addpath(optode_calibrations)


%
% check which optimizations have already been done
%
runs = [];
if ~exist('simple_aanderaa_optode_delays.mat')
  runs = 1;
else
  disp('Skipping simple Aanderaa delay optimization as result already exists. To rerun delete simple_aanderaa_optode_delays.mat')
end
if ~exist('aanderaa_optode_delays.mat')
  runs = [runs,2];
else
  disp('Skipping Aanderaa delay optimization as result already exists. To rerun delete aanderaa_optode_delays.mat')
end
if ~exist('geomar_optode_delays.mat')
  runs = [runs,3];
else
  disp('Skipping GEOMAR delay optimization as result already exists. To rerun delete geomar_optode_delays.mat')
end
if no_optim==1
  runs = [1,2,3];
end

%
% loop over the three variants of calibration parameters
% 1: simple Aanderaa   2: Aanderaa with corrections  3: GEOMAR
%
for nn=runs

  clear res2

  %
  % load the data and reduce to the pairs when we have up and down oxygen data
  %
  data = load([op.deplname,'_1sec']);
  load([op.deplname,'_1sec_derived'],'ctd_salinity_dpdt');
  data.ctd_salinity_dpdt = ctd_salinity_dpdt;

  ind = [1,length(data.oxygen)];
  for n=2:length(data.oxygen)-1
    if any(~isnan(data.oxygen{n})) & any(~isnan(data.oxygen{n+1}))
      if ~isempty(op.o_yo_numbers)
        if any(n==op.o_yo_numbers)
          ind = [ind,n,n+1];
        end
      else
        ind = [ind,n,n+1];
      end
    end
  end
  ind = unique(ind);
  if length(ind)==2 & ~isfield(op,'force_optode_delay_o')
    disp('No up-down oxygen pairs found. Please set optode delays in processing_params.m .')
    disp('e.g.:')
    disp('op.force_optode_delay_o = 30;')
    disp('op.force_optode_delay_t = 35;')
    return
  end
  fn = fieldnames(data);
  for n=1:length(fn)
    dummy = getfield(data,fn{n});
    dummy = {dummy{ind}};
    data = setfield(data,fn{n},dummy);
  end


  % 
  % look whether we can find a set of coefficients for this particular optode
  %
  cal_coeff.calib_param.tau_O2 = 7;
  cal_coeff.calib_param.tau_fact = -0.4;
  cal_coeff.calib_param.tempref_tauO2 = 20;
  cal_coeff.calib_param.cal_method = '';
  if nn==1
      found_coeff_type = 100;
  elseif nn==2
    if ~isempty(op.file_with_aanderaa_optode_calibration_parameters)
      disp(' ')
      disp('Found file with Aanderaa optode calibration parameters.')
      eval(['C = ',op.file_with_aanderaa_optode_calibration_parameters,';']);
      foil_coeff = convert_optfoilcoefstruc2mat(C);
      found_coeff_type = 1;
      cal_coeff.foil_coeff = foil_coeff;
    else
      disp(' ')
      found_coeff_type = -1;
      cal_coeff.foil_coeff = [];
      disp('Found no file with Aanderaa optode calibration parameters.')
    end
  elseif nn==3
    if ~isempty(op.file_with_geomar_optode_calibration_parameters)
      disp(' ')
      disp('Found file with GEOMAR optode calibration parameters.')
      cal_coeff = load([optode_calibrations,filesep,op.file_with_geomar_optode_calibration_parameters]);
      if isnan(cal_coeff.calib_param.tau_fact)
        cal_coeff.calib_param.tau_fact = -0.4;
      end
      if isnan(cal_coeff.calib_param.tempref_tauO2)
        cal_coeff.calib_param.tempref_tauO2 = 20;
      end
      cal_coeff.foil_coeff = [];
      found_coeff_type = 2;
    else
      disp(' ')
      disp('Found no file with GEOMAR optode calibration parameters.')
      found_coeff_type = -1;
    end
  end

  if found_coeff_type>0
    disp(' ')
    disp('Using:')
    disp(['tau_fact      = ',num2str(cal_coeff.calib_param.tau_fact),' sec/degC'])
    disp(['tempref_tauO2 = ',num2str(cal_coeff.calib_param.tempref_tauO2),' degC'])

    if isfield(op,'force_optode_delay_o')

      best_delay_o = op.force_optode_delay_o;
      best_delay_t = op.force_optode_delay_t;
      disp('using set optode delays')

    elseif no_optim==1
      
      if nn==1
        load simple_aanderaa_optode_delays
        best_delay_t = 0;
      elseif nn==2
        load aanderaa_optode_delays
      elseif nn==3
        load geomar_optode_delays
      end

    else

      delays_o = [10:4:50];
      if nn>1
        delays_t = [0:40:200];
      else
        delays_t = [0,0];
      end
      for k=[1:length(delays_o)]

        cal_coeff.calib_param.tau_O2 = delays_o(k);

        disp(['tau_O2        = ',num2str(cal_coeff.calib_param.tau_O2),' sec'])
        [res] = func_optimize_optode_delay(data,cal_coeff,delays_t,found_coeff_type);

        if length(res)>1
          for n=1:length(res)
            res2(n,k) = nsum(nsum(abs(res{n}.oxygen_calculated_undelayed(:,2:2:end-1)-...
                              res{n}.oxygen_calculated_undelayed(:,3:2:end))));
          end
        else
          res2(1,k) = nsum(nsum(abs(res{1}.oxygen_undelayed(:,2:2:end-1)-...
                              res{1}.oxygen_undelayed(:,3:2:end))));
        end
      end
      x = [min(delays_o):max(delays_o)];
      y = [min(delays_t):max(delays_t)];
      if length(y)==1
        res3 = interp1(delays_o,res2,x,'pchip');
        [minres3,ind] = nmin(res3);
        best_delay_o = x(ind);
        best_delay_t = 0;
      else
        res3 = interp2(delays_o,delays_t',res2,x,y','cubic');
        [x2,y2] = meshgrid(x,y);
        [minres3,ind] = nmin(res3(:));
        best_delay_o = x2(ind);
        best_delay_t = y2(ind);
      end
  
      if op.no_plot~=1
        figure(nn)
        clf
        if length(y)>1
          pcolor(x2,y2,res3)
          xlabel('Oxygen delay [sec]')
          ylabel('Temperature delay [sec]')
          colorbar
        else
          plot(x,res3)
          xlabel('Oxygen delay [sec]')
        end
        if nn==1
          title('Simple Aanderaa oxygen delay optimization')
        elseif nn==1
          title('Aanderaa oxygen delay optimization')
        elseif nn==1
          title('GEOMAR oxygen delay optimization')
        end
        print('-djpeg',['plots',filesep,'step_17_optode_delays_',int2str(nn),'.jpg'])
      end

      if nn==1
        save simple_aanderaa_optode_delays best_delay_o res3 x y res2 delays_o
      elseif nn==2
        save aanderaa_optode_delays best_delay_o best_delay_t res3 x y res2 delays_o delays_t
      elseif nn==3
        save geomar_optode_delays best_delay_o best_delay_t res3 x y res2 delays_o delays_t
      end

    end

    disp(' ')
    disp('Minimal down-up oxygen difference were found for:')
    disp(['o-delay of ',int2str(best_delay_o),' secs'])
    disp(['t-delay of ',int2str(best_delay_t),' secs'])

    cal_coeff.calib_param.tau_O2 = best_delay_o;
    data = load([op.deplname,'_1sec']);
    load([op.deplname,'_1sec_derived'],'ctd_salinity_dpdt');
    data.ctd_salinity_dpdt = ctd_salinity_dpdt;

    for mm=1:length(data.oxygen)
      data.oxygen_mumol_l_orig{mm} = data.oxygen{mm};
      d = sw_dens(data.ctd_salinity_dpdt{mm},data.ctd_temperature{mm},data.pressure{mm});
      data.oxygen{mm} = data.oxygen{mm}./d*1000;
    end

    [res1] = func_optimize_optode_delay(data,cal_coeff,best_delay_t,found_coeff_type);

    %
    % apply a 0-phase filter to the undelayed oxygen data to get rid of the introduced noise
    %
    data = load([op.deplname,'_1sec_derived']);
    if nn==1
      data.aanderaa_oxygen_undelayed = res1{1}.oxygen_undelayed;
      [b,a] = butter(3,0.05);
      for n=1:length(data.aanderaa_oxygen_undelayed)
        good = find(~isnan(data.aanderaa_oxygen_undelayed{n}));
        data.aanderaa_oxygen_undelayed_filtered{n} = data.aanderaa_oxygen_undelayed{n};
        if length(good)>10
          data.aanderaa_oxygen_undelayed_filtered{n}(good) =...
            filtfilt(b,a,data.aanderaa_oxygen_undelayed{n}(good));
        end
      end
    elseif nn==2
      data.aanderaa_oxygen_calculated = res1.oxygen_calculated;
      data.aanderaa_oxygen_calculated_undelayed = res1.oxygen_calculated_undelayed;
      [b,a] = butter(3,0.05);
      for n=1:length(data.aanderaa_oxygen_calculated_undelayed)
        good = find(~isnan(data.aanderaa_oxygen_calculated_undelayed{n}));
        data.aanderaa_oxygen_calculated_undelayed_filtered{n} = data.aanderaa_oxygen_calculated_undelayed{n};
        if length(good)>10
          data.aanderaa_oxygen_calculated_undelayed_filtered{n}(good) = ...
            filtfilt(b,a,data.aanderaa_oxygen_calculated_undelayed{n}(good));
        end
      end
    elseif nn==3
      data.geomar_oxygen_calculated_undelayed = res1.oxygen_calculated_undelayed;
      [b,a] = butter(3,0.05);
      for n=1:length(data.geomar_oxygen_calculated_undelayed)
        good = find(~isnan(data.geomar_oxygen_calculated_undelayed{n}));
        data.geomar_oxygen_calculated_undelayed_filtered{n} = data.geomar_oxygen_calculated_undelayed{n};
        if length(good)>10
          data.geomar_oxygen_calculated_undelayed_filtered{n}(good) = ...
            filtfilt(b,a,data.geomar_oxygen_calculated_undelayed{n}(good));
         end
      end
    end

    %
    % store the result
    %
    save([op.deplname,'_1sec_derived'],'-struct','data');

  end
end


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 17  at  ',datestr(now)])
disp(' ')
diary off

