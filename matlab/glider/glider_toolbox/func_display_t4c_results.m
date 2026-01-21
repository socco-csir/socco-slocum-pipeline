function [] = func_display_t4c_results()
% function [] = func_display_t4c_results()
% 
% GEOMAR SVN $Id: func_display_t4c_results.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% display the results for the temperature for conductivity optimization
%
% version 1  last change 20.02.2013

% G.Krahmann, GEOMAR Feb 2013

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
% check whether we are in the correct folder and whether the results mat-file exists
%
[pth,nam] = fileparts(pwd);
if ~strcmp(op.deplname,nam)
  error(['current working folder does agree with entry in processing_parameters.m'])
end
if ~exist([op.deplname,'.mat'])
  error(['could not find optimization results'])
end


%
% load results file
%
load(op.deplname)


%
% loop over results and display them
%
for n=1:length(allres)
  allres{n}.devi
end

