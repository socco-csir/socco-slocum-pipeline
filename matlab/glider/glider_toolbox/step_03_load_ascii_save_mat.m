function [] = step_03_load_ascii_save_mat()
% function [] = step_03_load_ascii_save_mat()
% 
% GEOMAR SVN $Id: step_03_load_ascii_save_mat.m 893 2022-01-10 15:07:34Z gkrahmann@geomar.de $
%
% load the merged and ascii/matlab converted files into matlab and resave as mat files
%
% uses : func_load_variable_list.m  func_slocum_sparse.m
%
% version 5.2.0  last change 10.01.2022

% G.Krahmann GEOMAR  Aug 2012

% add processing diary                                              GK, 30.01.2018  2-->3
% extended diary                                                    GK, 15.01.2019  3-->4
% added '_' to mat file name for R2020 compatibility                GK, 23.10.2020  4-->5.0.0
% there is a problem with the file names                            GK, 07.12.2020  5.0.0-->5.0.1
% use matlab internal gzip/gunzip                                   GK, 11.11.2021  5.0.1-->5.1.0
% add a variable existence check                                    GK, 10.01.2022  5.1.0-->5.2.0

%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 03  at  ',datestr(now)])
disp('converting ASCII to mat file')
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
% load list of variables that will be stored
%
slocum_variable_list = func_load_variable_list;


%
% figure out where the supplied executables are
%
%exe_dir = fileparts(which('step_03_load_ascii_save_mat'));
%if ~exist([exe_dir,filesep,'win_gzip.exe']) 
%  error('Could not find Webb''s executables')
%end



%
% since the loading of the ASCII data is slow we load it once and
% save it as mat file
% in this step this happens for each of the original files
%
cd raw_mat
d = dir('*bd_.m');
for n=1:length(d)

  clear dummy

  fname2 = d(n).name		% old file name

  [~,fname2] = fileparts(fname2);	% new file name
%   fname2 = [fname2,'_'];	% new file name   %isabelle edit from Gerd 2019 version


  % uncompress ascii files (*.dat), in case they are compressed
  if ~exist([fname2,'.mat']) & exist([fname2(1:end-1),'.dat.gz'])
    gunzip([fname2(1:end-1),'.dat.gz'])
%    if ispc
%      system([exe_dir,'\gzip -d ',fname2,'.dat.gz']);
%    else
%      system(['gunzip ',fname2,'.dat.gz']);
%    end
  end

  % catch file names and variables
  if exist(fname2,'var')
    whos
    error(['found variable  ',fname2,'  . It should not exist!!!'])
  end

  % only do the following, if the mat file has not yet been created
  if ~exist([fname2,'.mat'])

    disp(fname2)

    % execute the m-file, which loads the ascii data
    eval(fname2)

    % At this data all data is stored in the array 'data'
    % and all variables contain the column of 'data' in which the
    % data is stored.
    % Since this is not very handy we stick the data itself into the
    % variables.
    % Should the desired variable not exist in the data file, we fill
    % it with NaN's.
    for m=1:length(slocum_variable_list)
      if exist(slocum_variable_list{m})==1
        eval([slocum_variable_list{m},'=data(:,',slocum_variable_list{m},');'])
      else
        eval([slocum_variable_list{m},'=nan*m_present_time;'])
      end
    end  

    % run_name is a variable from the header of the m-file that loads the
    % original data
    % it contains the 8.3 and the full file name
    % here we extract the '8' of the 8.3 file name and store it as 'dos_id' variable
    rr1 = findstr('(',run_name);
    rr2 = findstr(')',run_name);
    dos_id = str2num(run_name(rr1+1:rr2-1))+m_present_time'*0;

    % look for cases when ALL science time information appears to be lacking
    % in these cases we try to interpolate  sci_m_present_time
    %
    % first we build a cell structure with all vectors that contain science time info
    % that is all time stamps and sci_m_present_time
    w = whos('sci*timestamp');
    all_time_infos = {};
    for m=1:length(w)
      if prod(w(m).size)>1
        all_time_infos{end+1} = w(m).name;
      end
    end

    % then we check whether any of these contain usable info
    dummy = sci_m_present_time;
    for m=1:length(all_time_infos)
      eval(['dummy = [dummy,',all_time_infos{m},'];'])
    end
    dummy = nsum(dummy,2);
    bad = find( dummy==0 );

    % If there are lines without time info we try to interpolate.
    % We simply replace  sci_m_present_time  with  m_present_time  
    % in case of all differences between them are exactly 0.
    if ~isempty(bad)
      warning('Found science data lines with no usable science information. Creating time stamps.')
      disp(run_name)
      disp('')
      dt = sci_m_present_time-m_present_time;
      good = find(~isnan(dt));
      if any(dt(good)~=0)
        tgood = find(~isnan(sci_m_present_time));
        if length(tgood)>1
          sci_m_present_time(bad) = interp1(tgood,sci_m_present_time(tgood),bad,'linear','extrap');
        elseif length(tgood)==1
          sci_m_present_time(bad) = interp1(tgood+[0,0.1],sci_m_present_time(tgood)+[0,0.1],bad,'linear','extrap');
        else
disp('needs to be programmed')
keyboard
        end
      else
        sci_m_present_time = m_present_time; 
      end
      sci_m_present_time(bad) = m_present_time(bad);
    end

    % save all data
    m_present_time = func_slocum_sparse( m_present_time );
    save([fname2,'.mat'],'m_present_time')
    clear m_present_time
    for m=1:length(slocum_variable_list)
      if ~strcmp(slocum_variable_list{m},'m_present_time')
        eval([slocum_variable_list{m},'=func_slocum_sparse(',slocum_variable_list{m},');'])
        save([fname2,'.mat'],'-append',slocum_variable_list{m})
        clear(slocum_variable_list{m})
      end
    end

    % recompress data file to save space
    gzip([fname2(1:end-1),'.dat'])
%    if ispc
%      system([exe_dir,'\gzip -1 ',fname2,'.dat']);
%    else
%      system(['gzip -1 ',fname2,'.dat']);
%    end

  end
end

func_slocum_irregular_geomar_names;

%
% switch back to main folder
%
cd ..


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 03  at  ',datestr(now)])
disp(' ')
diary off

