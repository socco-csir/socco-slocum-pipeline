function [] = step_00_prepare_processing()
% function [] = step_00_prepare_processing()
% 
% GEOMAR SVN $Id: step_01_copy_bd_to_processing_folder.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% Ask for a deployment name, create a folder and copy the template file.
%
% requires :   -
%
% version 1  last change 22.02.2017

% G.Krahmann, GEOMAR, Feb 2017

%
% first make sure we are not the in the toolbox folder
%
pd = pwd;
td = which('step_00_prepare_processing.m');
td = fileparts(td);
if strcmp(td,pd)
  error('You should not mix the toolbox with the actual processing. Please cd to another folder and restart.')
end


% 
% ask for folder name
%
fname = input('Please enter a name for the deployment. A folder with this name will be created. ','s');
suc = mkdir(fname);
if suc==1
  if ispc
    eval(['!copy ',td,filesep,'processing_parameters_template.m ',fname,filesep,'processing_parameters.m']);
  else
    eval(['!cp ',td,filesep,'processing_parameters_template.m ',fname,filesep,'processing_parameters.m']);
  end
  disp('  ')
  disp(['Before continuing you have to edit  ',fname,filesep,'processing_parameters.m .'])
else
  error('Folder could not be created.')
end


%
% check whether there is any reference data folder
%
if exist(['ctd_reference_data'])
  warning('Could not find a folder for calibrated CTD reference data. If you have such data, create a folder ''ctd_reference_data'' and put the data there.')
end
