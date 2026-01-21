function [] = func_slocum_irregular_geomar_names()
% function [] = func_slocum_irregular_geomar_names()
% 
% GEOMAR SVN $Id: func_slocum_irregular_geomar_names.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% rename all resaved mat-files
%
% uses : func_append_struct.m
%
% version 2  last change 11.12.2012

% G.Krahmann, GEOMAR  Oct 2012

% now called from load_ascii_save_mat         GK, 11.12.2012  1-->2


%
% handle the special case of irregular names
%
% the following is only required for GEOMAR gliders
% it should however not significantly slow down processing of others
% so we leave this code part in here
%
ir_names = {'deepy','bonpland'};		% these are the 'wrong' names
re_names = {'ifm02','ifm99'};                   % which will be replaced by these 'correct' names
for n_names = 1:length(ir_names)
  d1 = dir([upper(ir_names{n_names}),'*.MAT']);
  d2 = dir([lower(ir_names{n_names}),'*.mat']);
  d = func_append_struct(d1,d2);
  for n=1:length(d)
    dummy = [re_names{n_names},d(n).name(length(ir_names{n_names})+1:end)];
    if ispc
      system(['rename ',d(n).name,' ',dummy]);
    else
      system(['mv ',d(n).name,' ',dummy]);
    end
  end
end
