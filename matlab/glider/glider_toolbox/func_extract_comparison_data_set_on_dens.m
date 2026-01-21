function [] = func_extract_comparison_data_set_on_dens(op)
% function [] = func_extract_comparison_data_set_on_dens(op)
% 
% GEOMAR SVN $Id: func_extract_comparison_data_set_on_dens.m 703 2020-07-28 14:50:02Z gkrahmann@geomar.de $
%
% function to load the gridded glider data and extract a reduced data set of T,S,O on certain
% density levels
% this is done both for data with and without offsets applied
%
% input  :  op                    - processing option structure
%
% version 5.2.0  last change 25.01.2023

% G.Krahmann, GEOMAR, Oct 2014

% handle third oxygen variant (simple_aanderaa)                     GK, 03.03.2017  1-->2
% restrict comparison range to good salinities via op.s_yo_range    GK, 21.02.2018  2-->3
% added op.oa_factor                                                GK, 06.03.2020  3-->4.0.0
% added op.og_factor and op.og_pfactor                              GK, 13.03.2020  4.0.0-->5.0.0
% convert s_offset into a c_factor correction                       GK, 22.07.2022  5.0.0-->5.1.0
% change from CSIRO seawater to TEOS-10 library                     GK, 25.01.2023  5.1.0-->5.2.0


%
% load gridded data
%
load([op.deplname,'_gridded'])


%
% define reduced pressure set
%
% check_dens = [[1022:0.05:1025],[1025.01:0.01:1028]];

check_dens = [[1022:0.05:1025],[1025.01:0.01:1033]];

%
% loop over the two variants with and without offsets applied
% once without offsets applied and once with
%
for nn=1:2

  if nn==2
    ctd_temperature = ctd_temperature - op.t_offset;
    if op.use_c_factor==0
      good_salinity = good_salinity - op.s_offset;
    else
      % determine a conductivity factor that matches the given s_offset
      dummy_cond1 = gsw_R_from_SP(good_salinity,ctd_temperature,ctd_pressure);
      dummy_cond2 = gsw_R_from_SP(good_salinity-op.s_offset,ctd_temperature,ctd_pressure);
      cond_factor = nmean(dummy_cond2(:)./dummy_cond1(:));
      disp(['conductivity factor corr/uncorr : ',num2str(cond_factor)])
      good_salinity = gsw_SP_from_R(dummy_cond1*cond_factor,ctd_temperature,ctd_pressure);
      save cond_factor cond_factor
    end

%     oa = aanderaa_oxygen_undelayed_filtered * op.oa_factor - op.oa_offset;
%     if exist('aanderaa_oxygen_calculated_undelayed_filtered')
%       oa = aanderaa_oxygen_calculated_undelayed_filtered * op.oa_factor - op.oa_offset;
%     end
%     if exist('geomar_oxygen_calculated_undelayed_filtered')
%       og = geomar_oxygen_calculated_undelayed_filtered * op.og_factor -...
%         ctd_pressure * op.og_pfactor - op.og_offset;
%     else
%       og = nan*aanderaa_oxygen_undelayed_filtered;
%     end
%   else
%     oa = aanderaa_oxygen_undelayed_filtered;
%     if exist('aanderaa_oxygen_calculated_undelayed_filtered')
%       oa = aanderaa_oxygen_calculated_undelayed_filtered;
%     end
%     if exist('geomar_oxygen_calculated_undelayed_filtered')
%       og = geomar_oxygen_calculated_undelayed_filtered;
%     else
%       og = nan*aanderaa_oxygen_undelayed_filtered;
%     end
  end

  SA = gsw_SA_from_SP(good_salinity,pressure*ones(1,size(ctd_temperature,2)),0,0);
  CT = gsw_CT_from_t(SA,ctd_temperature,pressure*ones(1,size(ctd_temperature,2)));
  density = gsw_rho(SA,CT,pressure*ones(1,size(ctd_salinity,2)));

%  density = sw_pden(good_salinity,ctd_temperature,pressure*ones(1,size(ctd_temperature,2)),0);
  if 1
  fdensity = meanfilt(density,11);
  else
    fdensity = density;
    for n=1:size(density,2)
      dd = diff(fdensity(:,n));
      ind = find(dd<0);
      if ~isempty(ind)
        fdensity(1:ind(end)+1,n) = nan;
      end
    end
  end

  deep_salinity = good_salinity;
  deep_salinity(1:150,:) = nan;
  deep_temperature = ctd_temperature;
  deep_temperature(1:150,:) = nan;
%   deep_oa = oa;
%   deep_oa(1:150,:) = nan;
%   deep_og = og;
%   deep_og(1:150,:) = nan;


  %
  % prepare result arrays
  %
  if ~isempty(op.s_yo_numbers)
    indi = op.s_yo_numbers;
  else
    indi = [1:size(main_datenum,2)];
  end
  np = repmat(nan,length(check_dens),length(indi));
  nt = repmat(nan,length(check_dens),length(indi));
  ntd = repmat(nan,length(check_dens),length(indi));
  ns = repmat(nan,length(check_dens),length(indi));
  nsd = repmat(nan,length(check_dens),length(indi));
  noa = repmat(nan,length(check_dens),length(indi));
  noad = repmat(nan,length(check_dens),length(indi));
  nog = repmat(nan,length(check_dens),length(indi));
  nogd = repmat(nan,length(check_dens),length(indi));
  nlat = nmean(latitude);
  nlon = nmean(longitude);
  ntim = nmean(main_datenum);
  nlat = nlat(indi);
  nlon = nlon(indi);
  ntim = ntim(indi);


  %
  % loop over yos
  %
  for n=indi
    fprintf(1,'.');
    dd = diff(fdensity(:,n));
    ind = find(dd<=0);
    if ~isempty(ind)
      if ind(end)<length(dd)/2
        ind = [ind(end)+1:size(density,1)];
        ind = ind(find(~isnan(fdensity(ind,n))));
      else
        ind = [];
      end
    else
      ind = find(~isnan(fdensity(:,n)));
    end
    if length(ind)>1
      np(:,n) = interp1(fdensity(ind,n),pressure(ind,1),check_dens');
      nt(:,n) = interp1(fdensity(ind,n),ctd_temperature(ind,n),check_dens');
      ntd(:,n) = interp1(fdensity(ind,n),deep_temperature(ind,n),check_dens');
      if any(n==op.s_yo_numbers) | isempty(op.s_yo_numbers)
        ns(:,n) = interp1(fdensity(ind,n),good_salinity(ind,n),check_dens');
        nsd(:,n) = interp1(fdensity(ind,n),deep_salinity(ind,n),check_dens');
      end
%       if any(n==op.o_yo_numbers) | isempty(op.o_yo_numbers)
%         noa(:,n) = interp1(fdensity(ind,n),oa(ind,n),check_dens');
%         noad(:,n) = interp1(fdensity(ind,n),deep_oa(ind,n),check_dens');
%         nog(:,n) = interp1(fdensity(ind,n),og(ind,n),check_dens');
%         nogd(:,n) = interp1(fdensity(ind,n),deep_og(ind,n),check_dens');
%       end
    end
  end
  fprintf(1,'\n');

  %
  % save result
  %
  if nn==1
    save([op.deplname,'_comparison_data_on_dens'],'nt','ntd','ns','nsd','noa','noad',...
      'nog','nogd','np','nlat','nlon','ntim')
  else
    save([op.deplname,'_comparison_data_on_dens_with_offsets_applied'],...
      'nt','ntd','ns','nsd','noa','noad','nog','nogd','np','nlat','nlon','ntim')
  end

end
