function [] = step_05_best_position_vector()
% function [] = step_05_best_position_vector()
% 
% GEOMAR SVN $Id: step_05_best_position_vector.m 768 2021-01-04 14:26:19Z gkrahmann@geomar.de $
%
% load merged data and derive the best estimate of the glider's position 
% throughout the whole deployment
%
% this will use GPS positions at the surface and the glider's internal
% dead reckoned position under water
% the dead reckoned position will be corrected for difference between 
% the expected and the experienced surfacing positions
%
% version 5  last change 15.01.2019

% G.Krahmann, GEOMAR  2012

% apply magnetic deviation to m_water_vx/y and create a variable in *_1vector that this
% has been done.                                                  GK, 03.01.2018  2-->3
% removed magnetic correction part and added it to the finalization step
% this keeps variable names consistent                        
% add processing diary, remove v6 mat comp                        GK, 30.01.2018  3-->4
% extended diary and save plots                                   GK, 15.01.2019  4-->5

%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 05  at  ',datestr(now)])
disp('derive the best position vector by checking the glider''s GPS data interpretation')
disp('  and correcting the dead reckoned data to match GPS fixes')
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
% a 'best' guess of the position of the glider is developed.
% for deepy on deployment 9 this was a real mess
% not really sure whether positions are now really good
%
data = load([op.deplname,'_1vector']);


%
% if there is a hand-fix file, call it now
%
if exist('fix_gps_problems.m')
  disp('executing fix_gps_problems.m')
  fix_gps_problems;
end


%
% call the external routine for the best positions
%
% they will be stored under the new fieldnames n_lat and n_lon
%
[data] = func_slocum_best_position(data,op);


%
% fill in gaps in the position time series
%
bad = find( isnan(data.n_lat) | (data.n_lat==0 & data.n_lon==0) );
good = find( ~isnan(data.n_lat) & (data.n_lat~=0 | data.n_lon~=0) );
if length(good)>1
  data.n_lat(bad) = interp1(data.m_present_time(good),data.n_lat(good),...
    data.m_present_time(bad),'linear','extrap');
  data.n_lon(bad) = interp1(data.m_present_time(good),data.n_lon(good),...
	data.m_present_time(bad),'linear','extrap');
end


%
% resave the data
%
save([op.deplname,'_1vector'],'-struct','data')


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 05  at  ',datestr(now)])
disp(' ')
diary off

