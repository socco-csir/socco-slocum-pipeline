function [] = step_04_merge_mat_files()
% function [] = step_04_merge_mat_files()
% 
% GEOMAR SVN $Id: step_04_merge_mat_files.m 768 2021-01-04 14:26:19Z gkrahmann@geomar.de $
%
% merge the data from one or two mat files per call to a single mat file
% variable names remain the same and all data is concatenated
%
% the mat-filename of the merged data is by default the name of the current folder
%
% uses : func_load_variable_list.m  func_slocum_sparse.m
%
% version 6.2.0  last change 11.12.2020

% G.Krahmann, GEOMAR  Aug 2012

% add processing diary, remove v6 mat comp                           GK, 30.01.2018  2-->3
% add filling of gappy state variables                               GK, 09.01.2019  3-->4
% extended diary                                                     GK, 15.01.2019  4-->5
% added '_' to mat file name for R2020 compatibility                 GK, 23.10.2020  5-->6.0.0
% changed handling of '_'  so that old and new works                 GK, 08.12.2020  6.0.0-->6.0.1
% catch and fix GPS overruns and too long data sets, fix typo        GK, 08.12.2020  6.0.1-->6.1.0
% fix GPS overruns                                                   GK, 11.12.2020  6.1.0-->6.2.0


%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 04  at  ',datestr(now)])
disp('merge the single mat files to one mat file')
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
% load variable list
%
[slocum_variable_list] = func_load_variable_list;


%
% the single mat files are now loaded and the variables are appended 
%
if ispc
  d = dir('raw_mat\*bd_.mat');
  if isempty(d)
    d = dir('raw_mat\*bd.mat');
  end
else
  d = dir('raw_mat/*bd_.mat');
  if isempty(d)
    d = dir('raw_mat/*bd.mat');
  end
end
id_all = [];
for m=1:length(slocum_variable_list)
  disp(slocum_variable_list{m})
  eval([slocum_variable_list{m},'=[];'])
  for n=1:length(d)
    if ispc
      dummy = load(['raw_mat\',d(n).name],slocum_variable_list{m});
    else
      dummy = load(['raw_mat/',d(n).name],slocum_variable_list{m});
    end
    dummy = getfield(dummy,slocum_variable_list{m})';
    dummy = func_slocum_sparse(dummy);
    if ~isempty(dummy)
      eval([slocum_variable_list{m},'=[',slocum_variable_list{m},',dummy(:)''];'])
    else
      eval([slocum_variable_list{m},'=[',slocum_variable_list{m},',nan*m_present_time{m}(:)''];'])
    end
  end
  if strcmp(slocum_variable_list{m},'m_present_time')
    if nmax(m_present_time)-nmin(m_present_time)>365*86400 & ~isfield(op,'gps_overrun')
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
      disp('WARNING')
      disp(['First data point   : ',datestr(webbtime2mattime(nmin(m_present_time)))])
      disp(['Last data point    : ',datestr(webbtime2mattime(nmax(m_present_time)))])
      disp('This is longer than one year !')
      disp('Check whether you have included old bd-files or whether you had a GPS overrun.')
      disp('To fix the problem set   op.gps_overrun=YYYY  to the correct year of your deployment') 
      disp('Paused.   Hit ctrl-C to stop  or any other key to continue')
      pause
    end
  end
  if ~isempty(findstr(slocum_variable_list{m},'_time'))
    eval(['dummy=',slocum_variable_list{m},';'])
    ind = find( dummy>0.9e9 & dummy<1.1e9 );
    if ~isempty(ind)
      dummy(ind) = dummy(ind)+1024*7*86400;
      eval([slocum_variable_list{m},'=dummy;'])
    end
    disp(['added GPS overrun to  ',slocum_variable_list{m}])
  end
end


%
% here the m_present_time vector is being sorted 
% and this same sorting is applied to all other variables
%
% at the same time, turn the variables back to the sparse state to save space
%
[dummy,ind] = sort(m_present_time);
for m=1:length(slocum_variable_list)
  eval([slocum_variable_list{m},' = ',slocum_variable_list{m},'(ind);'])
end


%
% fill gaps in some variables
% Slocum variables are not stored again unless they change. This leads to some
% problems later on. Here we fill the NaN after a known state of the variable
% with the last known state,
%
disp(' ')
disp('filling gappy state variables')
complete_variable_list = {'m_air_pump'};
for n=1:length(complete_variable_list)
  disp(complete_variable_list{n})
  eval(['dummy = ',complete_variable_list{n},';'])
  for m=2:length(dummy)
    if isnan(dummy(m))
      dummy(m) = dummy(m-1);
    end
  end
  eval([complete_variable_list{n},' = dummy;'])
end


%
% the merged data is saved
%
str = ['save(''',op.deplname,'_1vector'''];
for m=1:length(slocum_variable_list)
  str = [str,',''',slocum_variable_list{m},''''];
end
str = [str,');'];
eval(str)


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 04  at  ',datestr(now)])
disp(' ')
diary off

