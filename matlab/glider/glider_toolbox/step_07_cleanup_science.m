function [] = step_07_cleanup_science()
% function [] = step_07_cleanup_science()
% 
% GEOMAR SVN $Id: step_07_cleanup_science.m 981 2022-12-16 16:52:53Z gkrahmann@geomar.de $
%
% SLOCUM science data has some spurious values and we try to get rid of them here
%
% also set to NaN, data that has been found bad
%
% version 10.1.0  last change 16.12.2022

% G.Krahmann, GEOMAR  Aug 2012

% read bad_data_list.txt file for bad data and blank them out     GK, 20.12.2012  2-->3
% small bug in <20 data points check                              GK, 21.03.2013  3-->4
% allow for yo-range in bad_data_list instead of single yos only  GK, 10.10.2013  4-->5
% catch 0 science pressure when CTD is being turned on            GK, 25.06.2014  5-->6
% remove exactly 0 timestamps                                     GK, 07.07.2014  6-->7
% add processing diary, remove v6 mat comp                        GK, 30.01.2018  7-->8
% extended diary                                                  GK, 15.01.2019  8-->9
% added variable and made code less variable name explicit,
% changed allowed CTD timestamps                                  GK, 09.12.2020  9-->10.0.0
% copy new SUNA variable into old SUNA variable                   GK, 16.12.2022  10.0.0-->10.1.0

%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 07  at  ',datestr(now)])
disp('clean up science data and set bad data to NaN')
diary off


%
% the flash card copy of the glider
%
if exist('./processing_parameters.m') | exist('.\processing_parameters.m')
  processing_parameters
else
  error(['could not find file ''processing_parameters.m'' in current folder ',pwd])
end


%
% load list of variables
%
slocum_variable_list = func_load_variable_list;

%
% load data
%
load([op.deplname,'_yos'])

%
% loop over yos and apply several corrections and cleanups to the data
%
for n=1:length(m_present_time)

  fprintf(1,'.',[]);
  
  %    
  % setting negative conductivity and temperatures (<-3) to nan
  %
  bad=find(sci_water_cond{n} < 0);
  if ~isempty(bad)
    sci_water_cond{n}(bad) = nan;
  end
  bad=find(sci_water_temp{n} <= -3);
  if ~isempty(bad)
    sci_water_temp{n}(bad) = nan;
  end

  %
  % remove 0 timestamps
  %
  for m=1:length(slocum_variable_list)
    if ~isempty(findstr(slocum_variable_list{m},'sci_')) &...
       ~isempty(findstr(slocum_variable_list{m},'_timestamp'))
      eval([slocum_variable_list{m},'{n}=nans(',slocum_variable_list{m},'{n},nan,0,''=='');'])
    end
  end


  %
  % there are flntu and optode data sets with just a single value of 0 as the first value. These
  % then fail the all-NaN tests. Set these single values to 0
  % 
  for m=1:length(slocum_variable_list)
    if ~isempty(findstr(slocum_variable_list{m},'sci_')) &...
       ~isempty(findstr(slocum_variable_list{m},'_units'))
      eval(['ind=find(~isnan(',slocum_variable_list{m},'{n}));'])
      if length(ind)==1
        eval([slocum_variable_list{m},'{n}(ind) = nan;'])
      end
    end
    if ~isempty(findstr(slocum_variable_list{m},'sci_')) &...
       ~isempty(findstr(slocum_variable_list{m},'_oxygen'))
      eval(['ind=find(~isnan(',slocum_variable_list{m},'{n}));'])
      if length(ind)==1
        eval([slocum_variable_list{m},'{n}(ind) = nan;'])
      end
    end
  end
if 0
  ind = find(~isnan(sci_flntu_chlor_units{n}));
  if length(ind)==1
    sci_flntu_chlor_units{n}(ind) = nan;
  end
  ind = find(~isnan(sci_flntu_turb_units{n}));
  if length(ind)==1
    sci_flntu_turb_units{n}(ind) = nan;
  end
  ind = find(~isnan(sci_flbbcd_bb_units{n}));
  if length(ind)==1
    sci_flbbcd_bb_units{n}(ind) = nan;
  end
  ind = find(~isnan(sci_flbbcd_chlor_units{n}));
  if length(ind)==1
    sci_flbbcd_chlor_units{n}(ind) = nan;
  end
  ind = find(~isnan(sci_flbbcd_cdom_units{n}));
  if length(ind)==1
    sci_flbbcd_cdom_units{n}(ind) = nan;
  end
  ind = find(~isnan(sci_oxy4_oxygen{n}));
  if length(ind)==1
    sci_oxy4_oxygen{n}(ind) = nan;
  end
  ind = find(~isnan(sci_oxy3835_oxygen{n}));
  if length(ind)==1
    sci_oxy3835_oxygen{n}(ind) = nan;
  end
  ind = find(~isnan(sci_oxy3835_wphase_oxygen{n}));
  if length(ind)==1
    sci_oxy3835_wphase_oxygen{n}(ind) = nan;
  end
  ind = find(~isnan(sci_fl2urrh_rhod_units{n}));
  if length(ind)==1
    sci_fl2urrh_rhod_units{n}(ind) = nan;
  end
  ind = find(~isnan(sci_bb2flsV2_units{n}));
  if length(ind)==1
    sci_bb2flsV2_units{n}(ind) = nan;
  end
  ind = find(~isnan(sci_bb2flsV8_units{n}));
  if length(ind)==1
    sci_bb2flsV8_units{n}(ind) = nan;
  end
end



  % 
  % Some science pressures are identical 0 but the glider is definitely not at the
  % surface. Possibly these occur when the CTD is just turned on during an inflection.
  %
  ind = find(~isnan(sci_water_pressure{n}));
  ind2 = find(sci_water_pressure{n}(ind)==0);
  if ~isempty(ind2)
    if ind2(1)==1
      sci_water_pressure{n}(ind(ind2(1))) = nan;
      sci_water_temp{n}(ind(ind2(1))) = nan;
      sci_water_cond{n}(ind(ind2(1))) = nan;
    end
  end


  %
  % fill in missing science sensor timestamps with science timestamp
  %
  for m=1:length(slocum_variable_list)
    if ~isempty(findstr(slocum_variable_list{m},'sci_')) &...
       ~isempty(findstr(slocum_variable_list{m},'_timestamp'))
      eval(['dummy=',slocum_variable_list{m},'{n};'])
      if all(isnan(dummy))
        eval([slocum_variable_list{m},'{n}=sci_m_present_time{n};'])
      end
    end
  end
if 0
  if all(isnan(sci_ctd41cp_timestamp{n}))
    sci_ctd41cp_timestamp{n} = sci_m_present_time{n};
  end
  if all(isnan(sci_oxy3835_wphase_timestamp{n}))
    sci_oxy3835_wphase_timestamp{n} = sci_m_present_time{n};
  end
  if all(isnan(sci_oxy3835_timestamp{n}))
    sci_oxy3835_timestamp{n} = sci_m_present_time{n};
  end
  if all(isnan(sci_oxy4_timestamp{n}))
    sci_oxy4_timestamp{n} = sci_m_present_time{n};
  end
  if all(isnan(sci_flntu_timestamp{n}))
    sci_flntu_timestamp{n} = sci_m_present_time{n};
  end
  if all(isnan(sci_flbbcd_timestamp{n}))
    sci_flbbcd_timestamp{n} = sci_m_present_time{n};
  end
end


  %
  % remove CTD data for early CTD timestamps
  %
  if ~isempty(sci_ctd41cp_timestamp{n})
    bad = find(sci_ctd41cp_timestamp{n}<0.9e9);
    if ~isempty(bad)
      sci_water_pressure{n}(bad) = nan;
      sci_water_temp{n}(bad) = nan;
      sci_water_cond{n}(bad) = nan;
      sci_ctd41cp_timestamp{n}(bad) = nan;
    end  
  end


  %
  % Discard all science data, if the first data set contains
  % science data, but the following 9 do not.
  % This seems to be the case of a leftover data point.
  % Good profiles that contain data in the first data set, will
  % also contain data in the following 9 data sets.
  %
  if length(sci_water_pressure{n})>=10
    if ~isnan(sci_water_pressure{n}(1)) & ...
	all( isnan(sci_water_pressure{n}(2:10)) )
      for m=1:length(slocum_variable_list)
        if ~isempty(findstr('sci_',slocum_variable_list{m})) &... 
           ~strcmp('sci_m_present_time',slocum_variable_list{m})
          eval([slocum_variable_list{m},'{n}(1) = nan;'])
        end
      end
    end
  end


  %
  % Remove all science data for which the CTD timestamp is more
  % than 60 seconds before the first main data set.
  % This takes care of cases when the CTD timestamp is exactly 0.
  %
  if length(m_present_time{n})>1 % there was one case with an empty yo
    bad = find( sci_ctd41cp_timestamp{n} < m_present_time{n}(1)-60 );
    if ~isempty(bad)
      for m=1:length(slocum_variable_list)
        if ~isempty(findstr('sci_',slocum_variable_list{m})) &... 
           ~strcmp('sci_m_present_time',slocum_variable_list{m})
          eval([slocum_variable_list{m},'{n}(bad) = nan;'])
        end
      end
    end
  end


  %
  % Remove data from the profiles before and after
  % these points are included because of the 'extension' of the
  % profile into the surrounding one.
  % One can recognize them by going to the central point of
  % the main profile and then look forward and backward for
  % gaps > 200sec in the science data set.
  % data before or after such gaps is deleted as it likely belongs
  % to another profile.
  %
  sp = sci_water_pressure{n};
  st = sci_ctd41cp_timestamp{n};
  if any(~isnan(sp))
    ind = find(~isnan(sp));
    mind = floor(nmedian(ind'));
    dst = diff(st(ind));

    longind = find(dst>200);
    if ~isempty(longind)
      longind1 = find( ind(longind)<mind );
      longind2 = find( ind(longind)>mind );
      if ~isempty(longind1)
        longind1 = ind(longind(longind1(end))+1);
        bad = [1:longind1-1];
        for m=1:length(slocum_variable_list)
          if ~isempty(findstr('sci_',slocum_variable_list{m})) &... 
             ~strcmp('sci_m_present_time',slocum_variable_list{m})
            eval([slocum_variable_list{m},'{n}(bad) = nan;'])
          end
        end
      end
      if ~isempty(longind2)
        longind2 = ind(longind(longind2(1)));
        bad = [longind2+1:length(st)];
        for m=1:length(slocum_variable_list)
          if ~isempty(findstr('sci_',slocum_variable_list{m})) &... 
             ~strcmp('sci_m_present_time',slocum_variable_list{m})
            eval([slocum_variable_list{m},'{n}(bad) = nan;'])
          end
        end
      end
    else
      longind1 = [];
      longind2 = [];
    end  
  end

  %
  % If science data is shorter than 20, remove all science data.
  %
  sp = sci_water_pressure{n};
  ind = find(~isnan(sp));
  if length(ind)<20
    for m=1:length(slocum_variable_list)
      if ~isempty(findstr('sci_',slocum_variable_list{m})) &... 
         ~strcmp('sci_m_present_time',slocum_variable_list{m})
        eval([slocum_variable_list{m},'{n} = nan*',slocum_variable_list{m},'{n};'])
      end
    end
  end 

  
  % 
  % there exists data where the time stamps of the whole 
  % science computer or the CTD are doubtful since they
  % do not change any more or are even reversed
  % here these should be captured and set to NaN
  %
  good = find(~isnan(sci_m_present_time{n}));
  dt = diff(sci_m_present_time{n}(good));
  ind = find(dt<=0);
  while ~isempty(ind)
    sci_m_present_time{n}(good(ind(end)+1)) = nan;
    good = find(~isnan(sci_m_present_time{n}));
    dt = diff(sci_m_present_time{n}(good));
    ind = find(dt<=0);
  end
  good = find(~isnan(sci_ctd41cp_timestamp{n}));
  dt = diff(sci_ctd41cp_timestamp{n}(good));
  ind = find(dt<=0);
  while ~isempty(ind)
    sci_ctd41cp_timestamp{n}(good(ind(end)+1)) = nan;
    good = find(~isnan(sci_ctd41cp_timestamp{n}));
    dt = diff(sci_ctd41cp_timestamp{n}(good));
    ind = find(dt<=0);
  end

  % copy new SUNA variable into old variable
  if any( ~isnan(sci_suna_nitrate_concentration{n}) )
    sci_suna_nitrate_um{n} = sci_suna_nitrate_concentration{n};
  end

end
fprintf(1,'\n',[]);

%
% look for bad data list
%
diary diary_processing.txt
if exist('bad_data_list.txt')
  disp('found bad_data_list.txt')
  fid = fopen('bad_data_list.txt','rt');
  while ~feof(fid)
    li = fgetl(fid);
    [vname,yo1,yo2,ind1,ind2] = strread(li,'%s %d %d %d %d');
    for m=[yo1:yo2]
      ind2_local = min([length(m_present_time{m}),ind2]);
      disp([vname{1},'{',int2str(m),'}(',int2str(ind1),':',int2str(ind2_local),')=nan;'])
      eval([vname{1},'{',int2str(m),'}(',int2str(ind1),':',int2str(ind2_local),')=nan;'])
    end
  end
end
diary off


%
% save the cleaned data
%
str = ['save(''',op.deplname,'_yos'''];
for m=1:length(slocum_variable_list)
  str = [str,',''',slocum_variable_list{m},''''];
end
str = [str,');'];
eval(str)


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 07  at  ',datestr(now)])
disp(' ')
diary off

