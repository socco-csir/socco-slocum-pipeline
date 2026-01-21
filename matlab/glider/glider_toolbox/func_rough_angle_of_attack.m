function [angle_of_attack,flight_angle] = func_rough_angle_of_attack(pitch)
% function [angle_of_attack,flight_angle] = func_rough_angle_of_attack(pitch)
% 
% GEOMAR SVN $Id: func_rough_angle_of_attack.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% Very crude estimate of the angle of attack simply from
% the pitch angle. This does not take the descent or ascent
% rate into account and should be used only to limit the
% forward speed calculated from pitch angle and vertical speed.
%
% input  :  pitch               - pitch angle in degrees
%
% output :  angle_of_attack     - angle of attack in degrees to be
%                                 added to the pitch angle to get
%                                 the total movement angle
%           flight_angle        - pitch + angle_of_attack limited to 90 degrees
%
% version 1  last change 25.02.2013

% G.Krahmann, GEOMAR, Feb 2013

% references:  Cooney, L.: Angle of Attack and Slocum Vehicle Model, 7pp, 2013.

%
% establish a simple fit to the data from the referenced manuscript
% but importantly include an AoA of 90 degrees (vertical fall) at 0 degree
% pitch angle, that is not included in the manuscript.
%
if 0
  x = [0,5,10,15,20,25,30];
  y = [90,5.2,4,3.3,2.7,2.4,2];
  figure(1)
  clf
  plot(x,1./(y.^1.5))
end



if ~iscell(pitch)
  si = sign(pitch);
  angle_of_attack = pitch/30*0.35;
  zero_s = find(pitch==0);
  if ~isempty(zero_s)
    angle_of_attack(zero_s) = 0.0000000001;
    si(zero_s) = 1;
  end
  angle_of_attack = 1./abs(angle_of_attack);
  angle_of_attack = angle_of_attack.^(1/1.5);
  bad = find(angle_of_attack>90);
  if ~isempty(bad)
    angle_of_attack(bad) = 90;
  end
  angle_of_attack = si.*angle_of_attack;
  flight_angle = pitch + angle_of_attack;
  bad = find(flight_angle>90);
  if ~isempty(bad)
    flight_angle(bad) = 90;
  end
  bad = find(flight_angle<-90);
  if ~isempty(bad)
    flight_angle(bad) = -90;
  end
else
  for n=1:length(pitch)
    [angle_of_attack{n},flight_angle{n}] = func_rough_angle_of_ttack(pitch{n});
  end
end

