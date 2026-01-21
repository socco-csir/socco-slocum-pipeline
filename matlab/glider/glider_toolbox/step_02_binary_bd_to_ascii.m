function [] = step_02_binary_bd_to_ascii()
% function [] = step_02_binary_bd_to_ascii()
% 
% GEOMAR SVN $Id: step_02_binary_bd_to_ascii.m 999 2023-01-11 11:28:54Z gkrahmann@geomar.de $
%
% convert binary bd files to ASCII
%
% At the moment this happens only for DBD and MBD file versions.
% Should we need to also convert SBD file versions, the code needs to
% be adapted!
%
% version 7.3.0  last change 11.01.2023

% G.Krahmann, GEOMAR  Aug 2012

% bug in mbd handling                                         GK, 04.12.2012  2-->3
% do not repeat calculations when they already had been done  GK, 01.03.2013  3-->4
% renamed executables and move to svn server                  GK, 04.07.2016  4-->5
% add processing diary                                        GK, 30.01.2018  5-->6
% extended diary                                              GK, 15.01.2019  6-->7
% use matlab internal gzip                                    GK, 11.11.2021  7-->7.1.0
% changed gzip call                                           GK, 10.01.2022  7.1.0-->7.1.1
% revert back to system gzip for Linux compression            GK, 05.08.2022  7.1.1-->7.2.0
% changed exe path                                            GK, 11.01.2023  7.2.0-->7.3.0

%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 02  at  ',datestr(now)])
disp('converting *bd to ASCII')
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
% check whether we are in the correct folder and whether all required folders exist
%
[pth,nam] = fileparts(pwd);
if ~strcmp(op.deplname,nam)
  error(['current working folder does agree with entry in processing_parameters.m'])
elseif exist('./raw_bd') | exist('.\raw_bd')
  if ~exist('./raw_mat') & ~exist('.\raw_mat')
    mkdir raw_mat
  end
  cd raw_bd
else
  error(['current working folder does not contain raw_bd folder'])
end


%
% figure out where Webb's executables are
%
exe_dir = fileparts(which('step_02_binary_bd_to_ascii'));
exe_dir = [exe_dir,filesep,'proprietary'];
if ~exist([exe_dir,filesep,'win_dbd2asc.exe'])
  error('Could not find Webb''s executables')
end


%
% first we create a list of dbd and mbd files and merge the lists
% so that we have one full list
%
% files where we do not have main-data (d/m/sbd files) are not used !
%
d1 = dir(['*.dbd']);
d2 = dir(['*.mbd']);
for n=1:length(d2)
  is_already_in_list = 0;
  for m=1:length(d1)
    if strcmp(d2(n).name(1:end-4),d1(m).name(1:end-4))
      is_already_in_list = 1;
    end
  end
  if is_already_in_list==0
    d1(end+1).name = d2(n).name;
    d1(end).bytes = d2(n).bytes;
    d1(end).date = d2(n).date;
    d1(end).isdir = d2(n).isdir;
    d1(end).datenum = d2(n).datenum;
  end
end
d = d1; 
all_bytes = sum([d.bytes]);
done_bytes = 0;

d = d(end:-1:1);

%
% now we have the complete list in structure d and figure out the
% matching science files (e/n/tbd files)
%
for n=1:length(d)
  m_name = d(n).name;
  s_name = m_name;
  s_name(end-2) = 'e';	
  if ~exist(s_name)		% Try existence of 'ebd' file. If not, switch to 'ndb' file.
    s_name(end-2) = 'n';
  end
  if  ~exist(s_name)		% Try existence of 'ebd' or 'nbd' file. If not, complain and empty.
    disp(['found no science file for ',d(n).name])
    s_name = '';
  end
  d(n).s_name = s_name;
end


%
% create the matlab-external system commands
%
if ispc
  comm1 = [exe_dir,'\win_dbd2asc '];
  comm2 = [exe_dir,'\win_dba_merge '];
  comm3 = ['move '];
  comm4 = [exe_dir,'\win_dba2_orig_matlab < dummy'];
else
  comm1 = [exe_dir,'/linux_dbd2asc '];
  comm2 = [exe_dir,'/linux_dba_merge '];
  comm3 = ['mv '];
  comm4 = [exe_dir,'/linux_dba2_orig_matlab < dummy'];
  % check whether they are executable
  [~,att] = fileattrib(deblank(comm1));
  if att.UserExecute==0
    fileattrib(deblank(comm1),'+x','u');
    [~,att] = fileattrib(deblank(comm1));
    if att.UserExecute==0
      error(['Could not make ',comm1(1:end-1),' executable. You have to do this by hand.']);
    end
  end
  [~,att] = fileattrib(deblank(comm2));
  if att.UserExecute==0
    fileattrib(deblank(comm2),'+x','u');
    [~,att] = fileattrib(deblank(comm2));
    if att.UserExecute==0
      error(['Could not make ',deblank(comm2),' executable. You have to do this by hand.']);
    end
  end
  [~,att] = fileattrib(comm4(1:end-8));
  if att.UserExecute==0
    fileattrib(comm4(1:end-8),'+x','u');
    [~,att] = fileattrib(comm4(1:end-8));
    if att.UserExecute==0
      error(['Could not make ',comm4(1:end-8),' executable. You have to do this by hand.']);
    end
  end
end


%
% execute the external commands
%
for n=1:length(d)

  % display file names
  disp(' ')
  done_bytes = done_bytes + d(n).bytes;
  disp(['file ',int2str(n),' of ',int2str(length(d)),'   ',...
	num2str(done_bytes/all_bytes*100,'%5.2f'),...
	'% of data converted'])

  % Look for existence of already converted files. 
  % Only if not already converted, convert the merged binary files to matlab-ascci files
  % Also compress the ascii files as they are large
  clear fname2
  fname2 = d(n).name;
  fname2 = fname2(1:end-4);
  if strcmp(d(n).name(end-2),'d')
    fname2 = [fname2,'_dbd.m'];
  elseif strcmp(d(n).name(end-2),'m')
    fname2 = [fname2,'_mbd.m'];
  else
    fname2 = [fname2,'_ebd.m'];
  end
  ind = findstr('-',fname2);
  fname2(ind) = '_';
  if ~exist(['../raw_mat/',fname2]) & ~exist(['..\raw_mat\',fname2])
    % merge or rename main/science files
    if ~isempty(d(n).s_name)		
      system([comm1,d(n).name,' > dummy1']);	% convert main file to ascii
      system([comm1,d(n).s_name,' > dummy2']);	% convert science file to ascii
      system([comm2,'dummy1 dummy2 > dummy']);	% merge the two ascii files into 'dummy'
    else
      system([comm1,d(n).name,' > dummy1']);	% convert science file to ascii
      system([comm3,' dummy1 dummy']);		% rename file to 'dummy'
    end
    system(comm4);
    if ispc
      ppd = pwd;
      ppd = fileparts(ppd);
      system([comm3,' ',fname2(1:end-6),'_dbd.dat ',ppd,'\raw_mat']);
      system([comm3,' ',fname2(1:end-6),'_dbd.m ',ppd,'\raw_mat']);
      ddd = dir('..\raw_mat\*.dat');
      if length(ddd)>0
        gzip(['..\raw_mat\*.dat'])
        delete(['..\raw_mat\*.dat'])
      end
    else
      system([comm3,' ',fname2(1:end-1),'* ../raw_mat']);
%      system(['gzip -1 -f ../raw_mat/*.dat']);
      ddd = dir('../raw_mat/*.dat');
      if length(ddd)>0
        system(['gzip -f ../raw_mat/*.dat']);
        %gzip(['../raw_mat/*.dat'])
      end
    end
  end  

  % if this pause statement is missing, it is difficult to ctrl-c the processing !
  pause(0.1)

end

%
% switch back to main folder
%
cd ..


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 02  at  ',datestr(now)])
disp(' ')
diary off

