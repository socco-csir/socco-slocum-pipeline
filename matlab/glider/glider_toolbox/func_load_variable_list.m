function [slocum_variable_list,allglider_variable_list] = func_load_variable_lists()
% function [slocum_variable_list,allglider_variable_list] = func_load_variable_lists()
% 
% GEOMAR SVN $Id: func_load_variable_list.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% load the list with all variables contained in the file
% 'variables_list.txt'
%
% version 1  last change 16.08.2012

% G.Krahmann GEOMAR  2012

%
% load the variable_list file for slocum gliders
%
fn = which('slocum_variable_list.txt');
fid = fopen(fn,'rt');
slocum_variable_list = {};
li = fgetl(fid);	% read header line
while ~feof(fid)
  li = fgetl(fid);
  slocum_variable_list{end+1} = deblank(fliplr(deblank(fliplr(li))));  
end
fclose(fid);

%
% load the processed variable_list file for all (?) gliders
%
fn = which('allglider_variable_list.txt');
fid = fopen(fn,'rt');
allglider_variable_list = {};
li = fgetl(fid);	% read header line
while ~feof(fid)
  li = fgetl(fid);
  allglider_variable_list{end+1} = deblank(fliplr(deblank(fliplr(li))));  
end
fclose(fid);
