function [dataout] = func_slocum_sparse(datain)
% function [dataout] = func_slocum_sparse(datain)
% 
% GEOMAR SVN $Id: func_slocum_sparse.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% 'sparse' the data and reverse
% in difference to standard matlab, NaN are treated as non-elements
%
% input  :	datain		- structure or cell or array of data
%			
% output :	dataout		- sparsed or filled structure or array of data
%
% uses :	nans.m
%
% version 2	last change 20.08.2012

% G.Krahmann, IFM-GEOMAR, July 2006

% handle cell objects                           GK, 20.08.2012  0.1-->2

if isstruct(datain)
  fnames = fieldnames(datain);
  if ~issparse(getfield(datain,fnames{1}));
    dataout = datain;
    for n=1:length(fnames)
      if ~strcmp(fnames{n},'fname') & ~strcmp(fnames{n},'sname') &...
		~strcmp(fnames{n},'start')
        d = getfield(datain,fnames{n});
        d = nans(d,-987654321,0,'==');
        d = nans(d,0,nan,'==');  
        d = sparse(d);
        dataout = setfield(dataout,fnames{n},d);
      end
    end
  else  
    dataout = datain;
    for n=1:length(fnames)
      if ~strcmp(fnames{n},'fname') & ~strcmp(fnames{n},'sname') &...
		~strcmp(fnames{n},'start')
        d = getfield(datain,fnames{n});
        d = full(d);
        d = nans(d,nan,0,'==');
        d = nans(d,0,-987654321,'==');  
        dataout = setfield(dataout,fnames{n},d);
      end
    end
  end  
elseif iscell(datain)
  for n=1:length(datain)
    if ~issparse(datain{n})
      d = datain{n};
      d = nans(d,-987654321,0,'==');
      d = nans(d,0,nan,'==');  
      dataout{n} = sparse(d);
    else
      d = full(datain{n});
      d = nans(d,nan,0,'==');
      dataout{n} = nans(d,0,-987654321,'==');  
    end
  end
else
  if ~issparse(datain);
    d = datain;
    d = nans(d,-987654321,0,'==');
    d = nans(d,0,nan,'==');  
    dataout = sparse(d);
  else  
    d = datain;
    d = full(d);
    d = nans(d,nan,0,'==');
    dataout = nans(d,0,-987654321,'==');  
    dataout = nans(dataout,0,-9876543210,'==');  
  end  
end
