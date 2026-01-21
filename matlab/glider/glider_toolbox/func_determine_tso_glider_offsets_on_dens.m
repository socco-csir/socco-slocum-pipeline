function [] = func_determine_tso_glider_offsets(op,reference_file)
% function [] = func_determine_tso_glider_offsets(op,reference_file)
% 
% GEOMAR SVN $Id: func_determine_tso_glider_offsets_on_dens.m 468 2018-01-31 16:52:24Z gkrahmann@geomar.de $
%
% compare the glider T,S,O data against a database of other well calibrated data
%
% version 1  last change 20.02.2014

% G.Krahmann, GEOMAR,  Feb 2014


%
% load reference data
%
if ~exist(['../',reference_file,'/',reference_file,'_comparison_data_on_dens.mat']);
  return
end
ref = load(['../',reference_file,'/',reference_file,'_comparison_data_on_dens']);
ref_off = load(['../',reference_file,'/',reference_file,'_comparison_data_on_dens_with_offsets_applied']);


% 
% screen display
%
disp('--------------------------------------------------------------------------------- ')
disp(' ')
disp('Differences glider data minus reference data')
disp([op.deplname,'  minus  ',reference_file,'   on density levels'])
disp(' ')

%
% load glider data
%
if 0
gli = load([op.deplname,'_comparison_data_on_dens']);
comp_type = ['glider_on_dens_',reference_file];
func_analyze_tso_data(gli,ref,comp_type,op,reference_file);
end

gli = load([op.deplname,'_comparison_data_on_dens_with_offsets_applied']);
comp_type = ['glider_on_dens_with_offsets_applied_',reference_file];
func_analyze_tso_data(gli,ref_off,comp_type,op,reference_file);

