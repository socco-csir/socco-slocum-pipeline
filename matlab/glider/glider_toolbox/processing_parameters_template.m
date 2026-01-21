op.no_plot = 1;

%
% this m-file contains all required parameters to process the following glider deployment
%
op.deplname = 'ifm13_depl07';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin Slocum relevant paramters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% location of a copy of the two Slocum flash cards
%
op.path_to_main_state_copy = '/data2/states/ifm13/main/ifm13_main_2019_12_19';
op.path_to_science_state_copy = '/data2/states/ifm13/science/ifm13_science_2019_12_19';

%
% minimum and maximum mission number
% here you can restrict the bd files being copied, in case there are missions from
% other deployments or simulations on the copies of the flash cards
%
% default (copy all bd files) is 0 and 9999
%
op.minimum_mission_number = 0;
op.maximum_mission_number = 187;

%
% To separate up and down casts we look for times at which a certain fraction of
% the maximum dive depth is passed either on the way down or up.
% Using these indices we go forward and backward to identify the inflection points.
% This is based solely on the pressure record and is independent of internal Slocum
% variables which in cases of aborts might be missing.
% Usually a fraction of 1/3 works well. In cases of W-shaped dives with a very deep
% center inflection point, one might need to use something different.
%
op.top_turnaround_limit = 1/3;

%
% to get rid of very shallow test dives one can set a minimum deepest depth for
% a dive to be included as a separate dive
% For deep gliders 16m works well, for very shallow deployments one needs to lower
% this value
%
op.minimum_deepest_depth_of_dive = 16;

%
% in the module to estimate the best position during the whole dive, 
% a special treatment is required for data from deployment ifm02_depl09
% 
% there a median is used to
%
% usually this should be set to 0
op.use_median = 0;

% 
% the best position estimate procedure is only used when there is a gap between
% observed GPS values for longer than 10 minutes
% If the gap is shorter we do not apply the procedure.
%
% usually this should be set to 250
op.min_dive_length = 250;

op.correct_glider_times = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End Slocum relevant paramters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin parameters for the optimization of the conductivity cell temperature estimate
% for unpumped Seabird CTDs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% this estimate is derived following
% ...
% Therein a glider speed ('surgespeed') dependent delay is applied to the measured temperature
% and this temperature then used to derive a salinity.
% A number of parameters control the details of the delay and are optimized through
% a non-linear optimization. The criterion is a minimum between the salinities of
% subsequent dives. Subsequent here means pairs of 'up' and subsequent 'down' data.
% This order of up and then down is important as it minimizes the time between
% data collected near the surface.

% An important parameter when using a glider speed dependent delay, is a minimum
% speed. The original code allowed for 0 speed which results in division by 0.
% Here a minimum speed (or flow through the cell, 'flowspeed' in the paper) 
% is used making the optimization more stable. 
% 0.05 m/s appears to be a good value.
op.minimum_flow_speed = 0.15;

op.flow_speed_degree = 0;

% We use both dp/dt and a dynamical flight model to derive the speed of the
% glider through the water that is then used to estimate the flow through
% the conductivity cell. In the model two parameters can be set. One to
% account for air in the oil and one for air bubbles that are usually
% present in the tail cone after a surfacing.
%op.first_dive_bubbles = 45;
%op.first_dive_bubbles = 350;
%op.sum_of_negative_angles_to_release_bubble = -30;

op.start_dive_parameter = [2,0.5,10];

op.is_pumped_ctd = 1;


% if there is air in the oil system, one can try to enter this here
% this is a number in cubic centimeters
% so if you think you 20cc air in the system, add a 20 here
op.air_V = 0;

% The optimization makes only sense when the measured conductivity data is good.
% During some deployment we had strong bio-fouling in the second half of the 
% deployment.
op.s_yo_numbers = [];
op.o_yo_numbers = [];

% The optimization over all profiles of a long deployment is likely to be too
% much for most computers. Here one can reduce the number of profiles included
% in the optimization. A number of 40 to 100 deep profiles is managable by
% a fast computer.
% Please note that only profiles for which down AND up cast exist are included in
% the optimization.
op.subsample_profile = [3,5,10];

% In the optimization procedure the area between the up and the down trace in a
% theta-S diagram is minimized. It appears that the apparent up and down shifts
% in the S data relative to the T data are quite large. This means that the temperature
% of the conductivity cell is lagging quite a bit. To speed up the calculation we
% here can subsample the data used in the theta-S diagram.
% It appears that reduction values up to 10 still deliver consistent results.
% If you want to be on the safe side or have only few profiles with up AND down data,
% set this to 1.
op.ts_data_reduction = [5];

% Aanderaa coefficients
op.file_with_aanderaa_optode_calibration_parameters = 'optode_00196_1206E_coeff';

% GEOMAR coefficients
op.file_with_geomar_optode_calibration_parameters = 'opt_calcoeff_M4831_00196_01206E_2017033103_M135_CTD113_uch_CALID03';

% file with CTD reference data
op.glider_reference_data = {};
op.mooring_reference_data = '';
op.ctd_press_reference_data = '../ctd_reference_data/glider_check_data_press.mat';
op.ctd_dens_reference_data = '../ctd_reference_data/glider_check_data_dens.mat';



% offsets
op.t_offset = 0;
op.s_offset = 0;
op.use_c_factor = 0;
op.oa_offset = 0;
op.oa_factor = 1;
op.og_offset = 0;
op.og_factor = 1;
op.og_pfactor = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End parameters for the optimization of the conductivity cell temperature estimate
% for unpumped Seabird CTDs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

op.microrider_installed = 0;
op.suna_installed = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings for the glider dynamical model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Slocum Volume list:
% G1 short shallow ?
% G1 short deep      55.2
% G1 long deep       63.8
% G1 MR              61.0
% G2 short deep      57.8
% G2 long deep       72.4
% G2 SUNA            59.3
op.glider_volume = (59.3)/1000;   % volume in m^3 

% frontal area of the glider
% number taken from Merckelbach et al.
% 0.038 m^2 for a regular glider, 1.5 * 0.038 m^2 for a microrider glider, 1.3 * 0.038 m^2 for a SUNA glider
op.frontal_area = 0.038; % in m^2 

% effective wing area 
% number taken from Merckelbach et al.
op.effective_wing_area = 0.1;  % in m^3
