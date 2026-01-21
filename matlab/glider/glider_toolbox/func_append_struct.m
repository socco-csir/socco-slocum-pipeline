function [newstruct] = append_struct(struct1,struct2)
% append two structures
% 
% GEOMAR SVN $Id: func_append_struct.m 890 2021-12-15 15:38:48Z gkrahmann@geomar.de $
%
% function [newstruct] = append_struct(struct1,struct2)

% append two (vector) structures along their elements
% fieldnames existing in only one structure will be created empty
% for the other
%
% input  :  struct1         - first structure
%           struct2         - second structure
%
% output :  newstruct       - appended structure
%
% version 1.0.0   last change 17.11.2021

% G.Krahmann, IFM-GEOMAR, Dec 2006

% fixed poss problem with empty struct       GK, 25.06.2008  0.1-->0.2
% something else might now be broken !!!
% changed header                             GK, 17.11.2021  0.2-->1.0.0


% start with first structure
newstruct = struct1;

% add missing fieldnames
fnames = fieldnames(struct2);
for n=1:length(fnames)
  if ~isfield(newstruct,fnames{n})
    newstruct = setfield(newstruct,{1},fnames{n},[]);
  end
end

% add elements of second structure
count = length(struct1);
for n=1:length(struct2)
  count = count+1;
  for m=1:length(fnames)
    newstruct = setfield(newstruct,{count},fnames{m},...
	getfield(struct2,{n},fnames{m}));
  end
end
