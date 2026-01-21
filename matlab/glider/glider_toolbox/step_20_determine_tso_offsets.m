function [] = step_20_determine_tso_offsets()
% function [] = step_20_determine_tso_offsets()
% 
% GEOMAR SVN $Id: step_20_determine_tso_offsets.m 557 2019-01-25 10:53:26Z gkrahmann@geomar.de $
%
% compare the glider T,S,O data against a database of other well calibrated data
%
% version 3  last change 15.01.2019

% G.Krahmann, GEOMAR,  Feb 2014

% add processing diary                                                     GK, 30.01.2018  1-->2
% extended diary                                                           GK, 15.01.2019  2-->3



%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 20  at  ',datestr(now)])
disp('calculate offsets for T, S, O against CTD and others gliders')
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
% loop over possible comparisons
%
!rm tso_offsets.txt
diary tso_offsets.txt
for n_ref = 1:length(op.glider_reference_data)

%  func_determine_tso_glider_offsets_on_press(op,op.glider_reference_data{n_ref})

  func_determine_tso_glider_offsets_on_dens(op,op.glider_reference_data{n_ref})

end

if ~isempty(op.ctd_dens_reference_data)
  func_determine_tso_ctd_offsets_on_dens(op,op.ctd_dens_reference_data);
end

%if ~isempty(op.ctd_press_reference_data)
%  func_determine_tso_ctd_offsets_on_press(op,op.ctd_press_reference_data);
%end

if ~isempty(op.mooring_reference_data)
  func_determine_tso_mooring_offsets(op.mooring_reference_data)
end

diary off


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 20  at  ',datestr(now)])
disp(' ')
diary off

