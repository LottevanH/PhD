function [label]=axes_axis_label(varargin)
  nvar=length(varargin);
  if nvar==1
    name=varargin{1};
    unit=[];
  elseif nvar==2
    name=varargin{1};
    unit=varargin{2};
  end
  if isempty(unit)
    label=name;
  else
    label=[name ' (' unit ')'];
  end
end