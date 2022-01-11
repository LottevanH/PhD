function [variable, ind] = search_variable(varname)
%SEARCH_VARIABLE Summary of this function goes here
%   Detailed explanation goes here
  global datasetinfo dataset
  variable=[];
  ix=regexp(datasetinfo.paralist, varname);
  ix=~cellfun('isempty',ix);
  ind = find(ix, 1);
  if isempty(ind)
      return 
  end
  variable = dataset{ind};
end

