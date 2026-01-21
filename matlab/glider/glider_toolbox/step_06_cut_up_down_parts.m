function [] = step_06_cut_up_down_parts()
% function [] = step_06_cut_up_down_parts()
% 
% GEOMAR SVN $Id: step_06_cut_up_down_parts.m 950 2022-10-14 15:42:13Z gkrahmann@geomar.de $
%
% take the merged data from the whole deployment and separate the 
% different up/down segments
%
% version 6.2.0  last change 14.10.2022

% G.Krahmann, GEOMAR  Aug 2012

% removed n_pressure                                    GK, 10.12.2012  2-->3
% add processing diary, remove v6 mat comp              GK, 30.01.2018  3-->4
% extended diary, save plots, glider clock correction   GK, 15.01.2019  4-->5
% fixed bugs in clock correction variables              GK, 23.01.2019  5-->6
% fixed bug in storing uncorrected clock                GK, 09.09.2020  6-->6.1.0
% warn against time reversals                           GK, 14.10.2022  6.1.0-->6.2.0

%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 06  at  ',datestr(now)])
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
% load list of variables
%
slocum_variable_list = func_load_variable_list;


%
% recut the merged data so that every single data set contains
% one Yo (down+up are separate elements)
%
load([op.deplname,'_1vector'])
if length(find(~isnan(m_depth)))>length(find(~isnan(m_pressure)))+100
  disp('This should not occur. You might have to hardwire replacing')
  disp('m_pressure with m_depth in module_6')
  warning('found more than 100 more usable m_depth values than m_pressure values')
end
if strcmp(op.deplname,'ifm11_depl02') | strcmp(op.deplname,'ifm10_depl01') |... 
    strcmp(op.deplname,'ifm07_depl02')
  l_pressure = m_depth;
  m_pressure = m_depth/10;
else
  l_pressure = m_pressure*10;
end
if length(~isnan(l_pressure))==0
  l_pressure = m_depth;
end
good = find(~isnan(l_pressure));
bad = find(diff(m_present_time(good))==0);
if ~isempty(bad)
  disp('found not changing time vector')
  good(bad+1) = [];
end
l_pressure(good(1):good(end)) = interp1(m_present_time(good),l_pressure(good),...
 	[m_present_time(good(1):good(end))]); 
l_pressure(3:end-2) = nmean([l_pressure(1:end-4);l_pressure(2:end-3);l_pressure(3:end-2);...
          l_pressure(4:end-1);l_pressure(5:end)]);
[b,a] = butter(3,0.05,'low');
good = find(~isnan(l_pressure));
l_pressure(good) = filtfilt(b,a,l_pressure(good));
starts = [];
bottoms = [];
ends = [];
diveids = [];

ids = unique(dos_id);
for n=1:length(ids)
  good = find(dos_id==ids(n));

  if nmax(l_pressure(good))>op.minimum_deepest_depth_of_dive &...
	( nmin(m_de_oil_vol(good))<0 | nmin(m_ballast_pumped(good))<0 )
    [pmax,imax] = nmax(l_pressure(good));
    dn = find(l_pressure(good(1:end-1))<pmax*op.top_turnaround_limit &... 
        l_pressure(good(2:end))>pmax*op.top_turnaround_limit);
    up = find(l_pressure(good(1:end-1))>pmax*op.top_turnaround_limit &... 
            l_pressure(good(2:end))<pmax*op.top_turnaround_limit);
    if 0
            disp(ids(n))
            disp(length(good))
            clf
            plot(good,l_pressure(good))
            hold on
            plot(good(dn),l_pressure(good(dn)),'xr')
            plot(good(up),l_pressure(good(up)),'xr')
             pause
    end
    if length(dn)==1 && length(up)==1
      bottoms(end+1) = good(imax);
      starts(end+1) = good(1);
      ends(end+1) = good(end);
      diveids(end+1) = ids(n);
    elseif length(dn)==length(up) && length(dn)~=0
      starts(end+1) = good(1);
      diveids(end+1) = ids(n);
      for m=1:length(dn)-1
        [pmin,imin] = nmin(l_pressure(good(up(m):dn(m+1))));
        ends(end+1) = good(up(m)+imin-1);
        starts(end+1) = ends(end)+1;
        diveids(end+1) = ids(n);
      end
      ends(end+1) = good(end);
      for m=1:length(dn)
        [pmax,imax] = nmax(l_pressure(good(dn(m):up(m))));
        bottoms(end+1) = good(dn(m)+imax-1);
      end
    elseif length(dn)==length(up) && length(dn)==0
      disp('no real down and up branches, skipping these')
    else
      disp('unequal number of down and up branches, skipping these')
      disp(ids(n))
      disp([n,length(dn),length(up)])
    end
  end
end

if op.no_plot~=1
  clf
  plot(m_present_time,l_pressure)
  hold on
  for n=1:length(starts)
    plot(m_present_time([starts(n),bottoms(n),ends(n)]),...
           l_pressure([starts(n),bottoms(n),ends(n)]),'r')
  end
end


yo_yos = [starts',bottoms',ends',m_present_time(starts)',...
          m_present_time(bottoms)',m_present_time(ends)',diveids'];
 
lt = length(m_present_time);
for m=1:length(slocum_variable_list)
  disp(slocum_variable_list{m})
  eval(['dummy = ',slocum_variable_list{m},';'])
  clear(slocum_variable_list{m})
  count = 1;
  for n=1:size(yo_yos,1)
    if length(dummy)==lt
      eval([slocum_variable_list{m},'{count}=dummy(yo_yos(n,1):yo_yos(n,2));'])
      eval([slocum_variable_list{m},'{count+1}=dummy(yo_yos(n,2)+1:yo_yos(n,3));'])
    else
      eval([slocum_variable_list{m},'{count}=[];'])
      eval([slocum_variable_list{m},'{count+1}=[];'])
    end
    count = count+2;
  end
end
diary diary_processing.txt
disp(['separated glider data into ',int2str(length(m_present_time)),' yos'])
diary off


%
% look whether there is an m-file that applies corrections to the raw yo data
% if so, execute this file
%
if exist('correct_bad_data_yos.m')
  disp('found m-file correcting bad data')
  correct_bad_data_yos;
end


%
% determine the yo-median time lag between glider and GPS and correct for
% each yo
%
last_clock_lag = 0;
if ~isfield(op,'do_not_correct_gps_glider_time_lags')
  correct_time = 1;
elseif op.do_not_correct_gps_glider_time_lags==1
  correct_time = 1;
else
  correct_time = 0;
end
for n=1:length(m_present_time)
  n_present_time_before_clock_correction{n} = m_present_time{n};
end
if correct_time==1
  disp(' ')
  disp('Glider clock correction')
  disp(' ')
  for n=1:length(m_system_clock_lags_gps)
    clock_lag(n) = nmean(m_system_clock_lags_gps{n});
    clock_lag_std(n) = nstd(m_system_clock_lags_gps{n});
    clock(n) = nmean(m_present_time{n});
    if isnan(clock_lag(n))
      clock_lag(n) = last_clock_lag;
    end
    last_clock_lag = clock_lag(n);
if 0
figure(1)
clf
plot(m_system_clock_lags_gps{n},'x')
hold on
plot(x_system_clock_adjusted{n},'r+')
n
pause
end
  end
  ind_jump = find(abs(diff(clock_lag))>=7);
  diary diary_processing.txt
  if length(ind_jump)==0
    disp('found no time jump because of a glider internal clock correction')
  else
    disp(['found time jumps because of a glider internal clock correction at yos ',int2str(ind_jump)])
  end
  diary off
  disp(' ')
  ind = [0,ind_jump];
  for n=1:length(ind)-1
    parts(n,1) = ind(n)+1;
    parts(n,2) = ind(n+1);
  end
  parts(length(ind),1) = ind(end)+1;
  parts(length(ind),2) = length(clock_lag);
  range = [-10:10];
  for n=1:size(parts,1)
    for m=parts(n,1):parts(n,2)
      ind = m+range;
      ind = ind( find(ind>=parts(n,1) & ind<=parts(n,2)) );  
      filtered_clock_lag(m) = median(clock_lag(ind));
    end
  end
  figure(2)
  clf
  plot(clock_lag)
  hold on
  plot(filtered_clock_lag,'g')
  plot([ind_jump;ind_jump+1],[clock_lag(ind_jump);clock_lag(ind_jump+1)],'r')
  xlabel('Yo')
  ylabel('System Clock lags GPS [seconds]')
  grid on
  if op.correct_glider_times==1
    disp('correcting glider times to better match GPS times')
    w = whos;
    for n=1:length(w)
      if ( ~isempty(strfind(w(n).name,'_timestamp')) | ~isempty(strfind(w(n).name,'present_time')) ) &...
        isempty(strfind(w(n).name,'before_clock_correction'))
        for m=1:length(m_present_time)
          eval([w(n).name,'{m}=',w(n).name,'{m}+filtered_clock_lag(m);'])
        end
      end
    end
    title('Glider clock correction (the filtered lags (green) are corrected)')
  else
    title('Glider clock correction (the filtered lags (green) are NOT corrected)')
    disp('NOT correcting glider times to better match GPS times')
  end
  figure(2)
  clf
  plot(clock_lag)
  hold on
  plot(filtered_clock_lag,'g')
  plot([ind_jump;ind_jump+1],[clock_lag(ind_jump);clock_lag(ind_jump+1)],'r')
  xlabel('Yo')
  ylabel('System Clock lags GPS [seconds]')
  title('Glider clock correction (the filtered lags (green) are corrected)')
  grid on
  print('-djpeg',['plots',filesep,'step_06_clock_lag.jpg'])
else
  diary diary_processing.txt
  disp('NOT correcting glider times to better match GPS times')
end
diary off
disp(' ')


%
% check for time reversals in the corrected time vector
%
bad = find(diff([m_present_time{:}])<=0);
if ~isempty(bad)
  warning('found time reversals in corrected time vector')
end


%
% create the save command and then execute it
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
disp(['Ending step 06  at  ',datestr(now)])
disp(' ')
diary off

