function [] = func_slocum_rename_all_bd_files()
% function [] = func_slocum_rename_all_bd_files()
% 
% GEOMAR SVN $Id: func_slocum_rename_all_bd_files.m 999 2023-01-11 11:28:54Z gkrahmann@geomar.de $
%
% rename all DOS-named bd files to full file names
%
% uses : func_append_struct.m
%
% version 3.1.0  last change 11.01.2023

% G.Krahmann, GEOMAR  Aug 2012

% now called from within copy_               GK, 11.12.2012  2-->3
% changed exe path                           GK, 11.01.2023  3-->3.1.0

%
% figure out where Webb's executables are
%
exe_dir = fileparts(which('step_01_copy_bd_to_processing_folder'));
exe_dir = [exe_dir,filesep,'proprietary'];
if ~exist([exe_dir,filesep,'win_dbd2asc.exe'])
  error('Could not find Webb''s executables')
end


%
% rename 8.3 format dbd files to full name dbd files
%
% first we need to make a list-structure of ALL bd files
%
disp('renaming files')
d1 = dir('*.DBD');
d2 = dir('*.dbd');
d = func_append_struct(d1,d2);
d2 = dir('*.EBD');
d = func_append_struct(d,d2);
d2 = dir('*.ebd');
d = func_append_struct(d,d2);
d2 = dir('*.MBD');
d = func_append_struct(d,d2);
d2 = dir('*.mbd');
d = func_append_struct(d,d2);
d2 = dir('*.NBD');
d = func_append_struct(d,d2);
d2 = dir('*.nbd');
d = func_append_struct(d,d2);
d2 = dir('*.SBD');
d = func_append_struct(d,d2);
d2 = dir('*.sbd');
d = func_append_struct(d,d2);
d2 = dir('*.TBD');
d = func_append_struct(d,d2);
d2 = dir('*.tbd');
d = func_append_struct(d,d2);
for n=1:length(d)
  if length(d(n).name)==12	% only these are proper 8.3 length DOS file names
    if ispc
      system([exe_dir,filesep,'win_rename_dbd_files.exe ',d(n).name]);
    else
      [~,att] = fileattrib([exe_dir,filesep,'linux_rename_dbd_files']);
      if att.UserExecute==0
        fileattrib([exe_dir,filesep,'linux_rename_dbd_files'],'+x','u');
        [~,att] = fileattrib([exe_dir,filesep,'linux_rename_dbd_files']);
        if att.UserExecute==0
          error(['Could not make ',exe_dir,filesep,'linux_rename_dbd_files',...
            ' executable. You have to do this by hand.']);
        end
      end
      system([exe_dir,filesep,'linux_rename_dbd_files ',d(n).name]);
    end
  end
end
