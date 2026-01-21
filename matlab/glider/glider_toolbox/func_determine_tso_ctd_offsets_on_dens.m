function [] = func_determine_tso_ctd_offsets_on_dens(op,reference_file)
% function [] = func_determine_tso_ctd_offsets_on_dens(op,reference_file)
% 
% GEOMAR SVN $Id: func_determine_tso_ctd_offsets_on_dens.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% compare the glider T,S,O data against a database of other well calibrated data on density levels
%
% version 2  last change 16.10.2014

% G.Krahmann, GEOMAR,  Feb 2014

% rewrite from original version to use prepared comparison data      GK, 16.10.2014  1-->2


%
% load reference data
%
ref = load(reference_file);
% if ~isfield(ref,'noa')
%   ref.noa = ref.no;
%   ref.nog = ref.no;
%   ref.noad = ref.nod;
%   ref.nogd = ref.nod;
% end


%
% load glider data and ananlyze
%
gli = load([op.deplname,'_comparison_data_on_dens']);
comp_type = 'ctd_on_dens';
func_analyze_tso_data(gli,ref,comp_type,op,reference_file);

gli = load([op.deplname,'_comparison_data_on_dens_with_offsets_applied']);
comp_type = 'ctd_on_dens_with_offsets_applied';
func_analyze_tso_data(gli,ref,comp_type,op,reference_file);
