function [] = func_determine_tso_glider_offsets(op,reference_file)
% function [] = func_determine_tso_glider_offsets(op,reference_file)
% 
% GEOMAR SVN $Id: func_determine_tso_glider_offsets_on_press.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% compare the glider T,S,O data against a database of other well calibrated data
%
% version 1  last change 20.02.2014

% G.Krahmann, GEOMAR,  Feb 2014


%
% load reference data
%
if ~exist(['../',reference_file,'/',reference_file,'_comparison_data_on_press.mat']);
  return
end
ref = load(['../',reference_file,'/',reference_file,'_comparison_data_on_press']);
ref_off = load(['../',reference_file,'/',reference_file,'_comparison_data_on_press_with_offsets_applied']);
check_press = [10:10:1000];


% 
% screen display
%
disp('--------------------------------------------------------------------------------- ')
disp(' ')
disp('Differences glider data minus reference data')
disp([op.deplname,'  minus  ',reference_file,'   on pressure levels'])
disp(' ')

%
% load glider data
%
gli = load([op.deplname,'_comparison_data_on_press']);
comp_type = ['glider_on_press_',reference_file];
func_analyze_tso_data(gli,ref,comp_type,op,reference_file);

gli = load([op.deplname,'_comparison_data_on_press_with_offsets_applied']);
comp_type = ['glider_on_press_with_offsets_applied_',reference_file];
func_analyze_tso_data(gli,ref_off,comp_type,op,reference_file);

