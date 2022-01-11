function [NaN] = savedata(varargin)
%SAVEDATA Summary of this function goes here
%   Detailed explanation goes here
  global dataset datasetinfo
  t_st=datasetinfo.dateran(1);
  t_ed=datasetinfo.dateran(2);
    
  fp_res=[pwd '/results/'];
  prefn=strrep(datasetinfo.projname,'/','-');
  if floor(t_st)==floor(t_ed)
    fn=[prefn '_' datestr(t_st,'yyyymmdd_HHMM') 	...
      '-' datestr(t_ed,'HHMM')];
  else
    fn=[prefn '_' datestr(t_st,'yyyymmdd_HHMM') 	...
      '-' datestr(t_ed,'yyyymmdd_HHMM')];
  end
  
  save([fp_res fn '.mat'], 'dataset', 'datasetinfo');
end

