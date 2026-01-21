function [] = step_19_create_comparison_data()
% function [] = step_19_create_comparison_data()
% 
% GEOMAR SVN $Id: step_19_create_comparison_data.m 703 2020-07-28 14:50:02Z gkrahmann@geomar.de $
%
% load gridded data and create reduced comparison data sets
%
% version 5.0.0  last change 13.03.2020

% G.Krahmann, GEOMAR,  Oct 2012

% add processing diary                                              GK, 30.01.2018  1-->2
% extended diary                                                    GK, 15.01.2019  2-->3
% added op.oa_factor                                                GK, 06.03.2020  3-->4.0.0
% added op.og_factor and op.og_pfactor                              GK, 13.03.2020  4.0.0-->5.0.0



%
% look for and load processing parameters
%
if exist('./processing_parameters.m') | exist('.\processing_parameters.m')
  processing_parameters
else
  error(['could not find file ''processing_parameters.m'' in current folder ',pwd])
end


%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 19  at  ',datestr(now)])
disp('apply offsets to T, S, O and create comparison data')
disp(['t-off  = ',num2str(op.t_offset)])
disp(['s-off  = ',num2str(op.s_offset)])
disp(['use_c_factor  = ',num2str(op.use_c_factor)])
disp(['oa-off  = ',num2str(op.oa_offset)])
disp(['oa-fac  = ',num2str(op.oa_factor)])
disp(['og-off  = ',num2str(op.og_offset)])
disp(['og-fac  = ',num2str(op.og_factor)])
disp(['og-pfac = ',num2str(op.og_pfactor)])
diary off

%
% extract comparison data set on pressure
%
func_extract_comparison_data_set_on_press(op);


%
% extract comparison data set on density
%
func_extract_comparison_data_set_on_dens(op);


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 19  at  ',datestr(now)])
disp(' ')
diary off

