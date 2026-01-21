function [ur,vr] = rotate_uv(u,v,ang)
% function [ur,vr] = rotate_uv(u,v,ang)
% 
% GEOMAR SVN $Id: func_rotate_uv.m 453 2018-01-03 12:54:59Z gkrahmann@geomar.de $
%
% Calculate the new components (ur,vr) of a velocity vector (u,v) 
% in a new coordinate system that is rotated relative to the 
% coordinate system in which (u,v) was measured.
%
% A positive angle 'ang' will give the numbers for a new
% coordinate system that is rotated counter-clockwise.
%
% The typical application is (u,v) measured by an instrument
% with a magnetic compass. When entering the local magnetic deviation
% given by the function 'magdev.m' as 'ang' you will receive
% the velocity vector (ur,vr) in true coordinates (ur is eastward, 
% vr is northward).   E.g.  Kiel/Germany has  magdev(54,10,0,2017)=2.72 degrees
%
% input  : u                        - u velocity component in instrument(magnetic) coordinate system
%          v                        - v velocity component in instrument(magnetic) coordinate system
%          ang                      - local magnetic deviation in degrees at time and place of measurement
%
% output : ur                       - rotated u velocity component
%          vr                       - rotated v velocity component
%
% 'ang' can fewer dimensions than u/v and will be inflated to match the size of u/v
%
% version 1  last change 02.01.2018

% G.Krahmann, GEOMAR Jan 2018    based on older functions with no clear explanations

% automatic help
if nargin==0
  help rotate_uv
  return
end

% inflate angle argument to match the size of the u/v matrices
if size(ang,1)==1
  ang = ones(size(u,1),1)*ang;
end
if size(ang,2)==1
  ang = ang*ones(1,size(u,2),1);
end

% convert to radian
ang = -ang*pi/180;

% calculate factors
cr = cos(ang);
sr = sin(ang);

% calculate rotated components
ur = u.*cr-v.*sr;
vr = u.*sr+v.*cr;

