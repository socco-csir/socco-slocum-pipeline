function [] = func_extract_comparison_data_set_on_press(op)
% function [] = func_extract_comparison_data_set_on_press(op)
% 
% GEOMAR SVN $Id: func_extract_comparison_data_set_on_press.m 303 2017-03-06 08:52:37Z gkrahmann@geomar.de $
%
% function to load the gridded glider data and extract a reduced depth data set of T,S,O
%
% input  :  op                    - processing option structure
%
% version 4  last change 16.10.2014

% G.Krahmann, GEOMAR, Mar 2014

% use only op.yo_numbers for S and op.o_yo_numbers for O comparison data
% bug in yo_numbers                                                        GK, 17.07.2014  2-->3
% rename function and store as 'on_press' file                             GK, 16.10.2014  3-->4

%
% load gridded data
%
load([op.deplname,'_gridded'])


%
% define reduced pressure set
%
check_press = [10:10:1000];


%
% loop over applied offsets
%
for nn=1:2


  %
  % prepare result arrays
  %
  np = repmat(nan,length(check_press),size(main_datenum,2));
  nt = repmat(nan,length(check_press),size(main_datenum,2));
  ns = repmat(nan,length(check_press),size(main_datenum,2));
  noa = repmat(nan,length(check_press),size(main_datenum,2));
  nog = repmat(nan,length(check_press),size(main_datenum,2));
  nlat = nmean(latitude);
  nlon = nmean(longitude);
  ntim = nmean(main_datenum);


  %
  % loop over yos
  %
  for n=1:size(main_datenum,2)
    fprintf(1,'.')
    for m=1:length(check_press)
      [dummy,ind] = nmin(abs(pressure-check_press(m)));
      if dummy<1 
        if nn==1
          nt(m,n) = ctd_temperature(ind,n);
        else
          nt(m,n) = ctd_temperature(ind,n) - op.t_offset;
        end
        np(m,n) = pressure(ind);
        if any(n==op.s_yo_numbers) | isempty(op.s_yo_numbers)
          if nn==1
            ns(m,n) = good_salinity(ind,n);
          else
            ns(m,n) = good_salinity(ind,n) - op.s_offset;
          end
        end
%         if any(n==op.o_yo_numbers) | isempty(op.o_yo_numbers)
%           if nn==1
%             noa(m,n) = aanderaa_oxygen_calculated_undelayed_filtered(ind,n);
%           else
%             noa(m,n) = aanderaa_oxygen_calculated_undelayed_filtered(ind,n) - op.oa_offset;
%           end
%           if ~exist('geomar_oxygen_calculated_undelayed_filtered')
%             nog(m,n) = nan;
%           else
%             if nn==1
%               nog(m,n) = geomar_oxygen_calculated_undelayed_filtered(ind,n);
%             else
%               nog(m,n) = geomar_oxygen_calculated_undelayed_filtered(ind,n) - op.og_offset;
%             end
%           end
%         end
      end
    end
  end
  fprintf(1,'\n')

  ntd = nt;
  ntd(1:15,:) = nan;
  nsd = ns;
  nsd(1:15,:) = nan;
  nogd = nog;
  nogd(1:15,:) = nan;
  noad = noa;
  noad(1:15,:) = nan;


  %
  % save result
  %
  if nn==1
    save([op.deplname,'_comparison_data_on_press'],'nt','ntd','ns','nsd','nog','nogd','noa','noad',...
      'np','nlat','nlon','ntim')
  else
    save([op.deplname,'_comparison_data_on_press_with_offsets_applied'],...
      'nt','ntd','ns','nsd','nog','nogd','noa','noad',...
      'np','nlat','nlon','ntim')
  end
end
