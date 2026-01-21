function [jump_ind,jump_correction,down_or_up] = func_find_pressure_jumps(nav_pressure_old,m_present_time,sci_pressure,ctd_time)
% function [jump_ind,jump_correction,down_or_up] = func_find_pressure_jumps(nav_pressure_old,m_present_time,sci_pressure,ctd_time)
% 
% GEOMAR SVN $Id: func_find_pressure_jumps.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% The navigational pressure sensor of Slocum gliders in some cases has a regular jump to
% larger pressures at about 5m depth. This jump reverses at the end of dives.
% Here we try to identify and correct the problem cases.
%
% input  : nav_pressure_old             - navigational pressure time series in dbar
%          m_present_time               - time stamp for the navigational pressure measurement
%          sci_pressure                 - science pressure time series in dbar
%          ctd_time                     - time stamp for the science pressure measurement
%
% output : jump_ind                     - index of the value before the jump
%          jump_correction              - how much to subtract
%          down_or_up                   - is this a down or up yo
%
% version 1, last change 21.03.2013

% G.Krahmann, GEOMAR, Mar 2013

%
% first we interpolate the science pressure onto the time stamps of the navigational pressure sensor
%
good_sci = find(~isnan(ctd_time+sci_pressure));
good_nav = find(~isnan(nav_pressure_old));
sci_pressure_int = nan*nav_pressure_old;
if length(good_sci)>1
  sci_pressure_int = interp1(ctd_time(good_sci),sci_pressure(good_sci),m_present_time);
else
  sci_pressure_int = nan*m_present_time;
end

%
% prepare data
%
delta_ti = diff(m_present_time(good_nav));
delta_p_nav = diff(nav_pressure_old(good_nav));
delta_p_sci = diff(sci_pressure_int(good_nav));
jump_nav = delta_p_nav./delta_ti;
jump_sci = delta_p_sci./delta_ti;
down_or_up = 0;

%
% try to determine jumps to higher pressure
%
if nav_pressure_old(good_nav(1))<nav_pressure_old(good_nav(end))
  ind = find(jump_nav>0.05 & nav_pressure_old(good_nav(1:end-1))<1.5);
  good_ind = [];
  for n=1:length(ind)
    ind1 = ind(n)+[-20:0];
    ind1 = ind1(find(ind1>0));
    ind2 = ind(n)+[0:20];
    ind2 = ind2(find(ind2<=length(good_nav)));
    if nmean(nav_pressure_old(good_nav(ind1)))<0.5 & nmean(nav_pressure_old(good_nav(ind2)))>0.5 &...
        nmean(nav_pressure_old(good_nav(ind2))) - nmean(nav_pressure_old(good_nav(ind1))) > 0.2
      good_ind = [good_ind,ind(n)];
    end
  end
  down_or_up = -1;
end

%
% try to determine jumps to lower pressure
%
if nav_pressure_old(good_nav(1))>nav_pressure_old(good_nav(end))
  ind = find(jump_nav<-0.05 & nav_pressure_old(good_nav(1:end-1))<1.5);
  good_ind = [];
  for n=1:length(ind)
    ind1 = ind(n)+[-20:0];
    ind1 = ind1(find(ind1>0));
    ind2 = ind(n)+[0:20];
    ind2 = ind2(find(ind2<=length(good_nav)));
    if nmean(nav_pressure_old(good_nav(ind1)))>0.5 & nmean(nav_pressure_old(good_nav(ind2)))<0.5 &...
        nmean(nav_pressure_old(good_nav(ind1))) - nmean(nav_pressure_old(good_nav(ind2))) > 0.2
      good_ind = [good_ind,ind(n)];
    end
  end
  down_or_up = 1;
end
%keyboard

if 0
clf
plot(m_present_time,nav_pressure_old,'+',ctd_time,sci_pressure,'x')
hold on
plot(m_present_time(good_nav(good_ind)),nav_pressure_old(good_nav(good_ind)),'og')
pause
end

%
% determine the correction value
%
if ~isempty(good_ind)
  bad_diff = nav_pressure_old(good_nav(good_ind(1)+1)) - nav_pressure_old(good_nav(good_ind(1)));
  correct_diff = nmedian(jump_nav) * (m_present_time(good_nav(good_ind(1)+1)) - m_present_time(good_nav(good_ind(1))));
  jump_correction = bad_diff - correct_diff;
  if 0
    disp([good_ind(1),bad_diff,jump_correction])
    pause
  end
  jump_ind = good_nav(good_ind(1));
else
  jump_ind = 1;
  jump_correction = 0;
end

