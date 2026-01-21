function [] = step_01_copy_bd_to_processing_folder()
% function [] = step_01_copy_bd_to_processing_folder()
% 
% GEOMAR SVN $Id: step_01_copy_bd_to_processing_folder.m 893 2022-01-10 15:07:34Z gkrahmann@geomar.de $
%
% copy all sorts of bd and log files from the states folders to the processing folder
% also copy all compression information that is stored in the cache folders on the gliders
%
% requires :   processing_parameters.m in the current working folder
%
% version 9.1.1  last change 10.01.2022

% G.Krahmann, GEOMAR, Aug 2012

% better PC/Linux handling                GK, 27.10.2012  2-->3
% fixed bug in removing file 0000*        GK, 19.12.2012  3-->4
% upper case CAC file handling on unix    GK, 13.02.2013  4-->5
% more case problems solved               GK, 18.03.2013  5-->6
% create folder plots                     GK, 14.03.2017  6-->7
% add diary                               GK, 30.01.2018  7-->8
% extended diary                          GK, 15.01.2019  8-->9
% use matlab internal gunzip              GK, 11.11.2021  9-->9.1.0
% fixed bugs with internal gunzip (it crashes when there are no input files)
% added another path to CAC files         GK, 10.01.2022  9.1.0-->9.1.1

%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 01  at  ',datestr(now)])
disp('copying *bd files to processing folder')
diary off


%
% look for and load processing parameters which contain the location of
% the flash card copy of the glider
%
if exist('./processing_parameters.m')
  processing_parameters
else
  error(['could not find file ''processing_parameters.m'' in current folder ',pwd])
end


%
% check whether we are in the correct folder and whether all required folders exist
%
[pth,nam] = fileparts(pwd);
if ~exist('./plots')
  mkdir plots
end
if ~strcmp(op.deplname,nam)
  error(['current working folder does agree with entry in processing_parameters.m'])
elseif exist('./raw_bd')
  cd raw_bd
else
  mkdir raw_bd
  cd raw_bd
end
if ~exist('./cache')
  mkdir cache
end


%
% get paths
%
md = op.path_to_main_state_copy;
sd = op.path_to_science_state_copy;
disp(['copying files from  ',md])

if ispc
  eval(['!copy ',md,'\LOGS\*.*bd* .'])
  eval(['!copy ',md,'\LOGS\*.*BD* .'])
  eval(['!copy ',md,'\logs\*.*bd* .'])
  eval(['!copy ',md,'\logs\*.*BD* .'])
  eval(['!copy ',md,'\SENTLOGS\*.*bd* .'])
  eval(['!copy ',md,'\SENTLOGS\*.*BD* .'])
  eval(['!copy ',md,'\sentlogs\*.*bd* .'])
  eval(['!copy ',md,'\sentlogs\*.*BD* .'])
  eval(['!copy ',md,'\STATE\CACHE\*.CAC* cache'])
  eval(['!copy ',md,'\STATE\CACHE\*.cac* cache'])
  eval(['!copy ',md,'\State\CACHE\*.CAC* cache'])
  eval(['!copy ',md,'\State\CACHE\*.cac* cache'])
  eval(['!copy ',md,'\state\cache\*.CAC* cache'])
  eval(['!copy ',md,'\state\cache\*.cac* cache'])
  eval(['!copy ',md,'\State\cache\*.CAC* cache'])
  eval(['!copy ',md,'\State\cache\*.cac* cache'])
  eval(['!copy ',md,'\state\CACHE\*.CAC* cache'])
  eval(['!copy ',md,'\state\CACHE\*.cac* cache'])

  eval(['!copy ',sd,'\LOGS\*.*bd* .'])
  eval(['!copy ',sd,'\LOGS\*.*BD* .'])
  eval(['!copy ',sd,'\logs\*.*bd* .'])
  eval(['!copy ',sd,'\logs\*.*BD* .'])
  eval(['!copy ',sd,'\SENTLOGS\*.*bd* .'])
  eval(['!copy ',sd,'\SENTLOGS\*.*BD* .'])
  eval(['!copy ',sd,'\sentlogs\*.*bd* .'])
  eval(['!copy ',sd,'\sentlogs\*.*BD* .'])
  eval(['!copy ',sd,'\STATE\CACHE\*.CAC* cache'])
  eval(['!copy ',sd,'\STATE\CACHE\*.cac* cache'])
  eval(['!copy ',sd,'\state\cache\*.CAC* cache'])
  eval(['!copy ',sd,'\state\cache\*.cac* cache'])
  eval(['!copy ',sd,'\State\cache\*.CAC* cache'])
  eval(['!copy ',sd,'\State\cache\*.cac* cache'])
  eval(['!copy ',sd,'\state\CACHE\*.CAC* cache'])
  eval(['!copy ',sd,'\state\CACHE\*.cac* cache'])

else

  eval(['!cp ',md,'/LOGS/*.*bd* .'])
  eval(['!cp ',md,'/LOGS/*.*BD* .'])
  eval(['!cp ',md,'/logs/*.*bd* .'])
  eval(['!cp ',md,'/logs/*.*BD* .'])
  eval(['!cp ',md,'/SENTLOGS/*.*bd* .'])
  eval(['!cp ',md,'/SENTLOGS/*.*BD* .'])
  eval(['!cp ',md,'/sentlogs/*.*bd* .'])
  eval(['!cp ',md,'/sentlogs/*.*BD* .'])
  eval(['!cp ',md,'/STATE/CACHE/*.CAC* cache'])
  eval(['!cp ',md,'/STATE/CACHE/*.cac* cache'])
  eval(['!cp ',md,'/State/CACHE/*.CAC* cache'])
  eval(['!cp ',md,'/State/CACHE/*.cac* cache'])
  eval(['!cp ',md,'/state/cache/*.CAC* cache'])
  eval(['!cp ',md,'/state/cache/*.cac* cache'])
  eval(['!cp ',md,'/state/CACHE/*.CAC* cache'])
  eval(['!cp ',md,'/state/CACHE/*.cac* cache'])

  eval(['!cp ',sd,'/LOGS/*.*bd* .'])
  eval(['!cp ',sd,'/LOGS/*.*BD* .'])
  eval(['!cp ',sd,'/logs/*.*bd* .'])
  eval(['!cp ',sd,'/logs/*.*BD* .'])
  eval(['!cp ',sd,'/SENTLOGS/*.*bd* .'])
  eval(['!cp ',sd,'/SENTLOGS/*.*BD* .'])
  eval(['!cp ',sd,'/sentlogs/*.*bd* .'])
  eval(['!cp ',sd,'/sentlogs/*.*BD* .'])
  eval(['!cp ',sd,'/STATE/CACHE/*.CAC* cache'])
  eval(['!cp ',sd,'/STATE/CACHE/*.cac* cache'])
  eval(['!cp ',sd,'/State/CACHE/*.CAC* cache'])
  eval(['!cp ',sd,'/State/CACHE/*.cac* cache'])
  eval(['!cp ',sd,'/state/cache/*.CAC* cache'])
  eval(['!cp ',sd,'/state/cache/*.cac* cache'])
  eval(['!cp ',sd,'/state/CACHE/*.CAC* cache'])
  eval(['!cp ',sd,'/state/CACHE/*.cac* cache'])
  eval(['!cp ',sd,'/State/cache/*.CAC* cache'])
  eval(['!cp ',sd,'/State/cache/*.cac* cache'])

  ddd = dir('*.gz');
  if length(ddd)>0 
    gunzip('*.gz')
  end
  ddd = dir('*.GZ');
  if length(ddd)>0 
    gunzip('*.GZ')
  end
  ddd = dir('cache/*.gz');
  if length(ddd)>0 
    cd cache
    gunzip('*.gz')
    cd ..
  end
  ddd = dir('cache/*.GZ');
  if length(ddd)>0 
    cd cache
    gunzip('*.GZ')
    cd ..
  end

end


%
% make sure that we have the cache files in lower case
% as the unix version wants this
%
cd cache
d = dir('*.CAC');
for n=1:length(d)
  eval(['!cp ',d(n).name,' ',lower(d(n).name)])
end
cd ..

%
% in case there is a minimum mission number larger than 0 
% or the maximum mission number smaller than 9999 we have a look at
% all files and delete the ones that have lower mission numbers
%
if op.minimum_mission_number>0 | op.maximum_mission_number<9999
  diary diary_processing.txt
  disp('removing some *bd files according to processing_parameters.m')
  diary off
  d = dir('*');
  lastfname = 'xxxx';
  for n=1:length(d)
    fname = d(n).name;
    if length(fname)==12
      if (str2num(fname(1:4))<op.minimum_mission_number | ...
	  str2num(fname(1:4))>op.maximum_mission_number) & ~strcmp(lastfname,fname(1:4)) 
        disp(['!rm ',fname(1:4),'*'])
        eval(['!rm ',fname(1:4),'*'])
        lastfname = fname(1:4);
      end
    end
  end
end


%
% rename from DOS 8.3 format file names to full bd file names
%
func_slocum_rename_all_bd_files;


%
% go back to deployment folder
%
cd ..


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 01  at  ',datestr(now)])
disp(' ')
diary off
