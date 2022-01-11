function [label]=axes_axis_label(varargin)
  nvar=length(varargin);
  mode=1;
  if nvar==1
    name=varargin{1};
    unit=[];
  elseif nvar==2
    name=varargin{1};
    unit=varargin{2};
  elseif nvar==3
    name=varargin{1};
    unit=varargin{2};
    mode=varargin{3};
  end
  if isempty(unit)
    label=name;
    if mode==2
      label={label,[]};
    end
  else
    if mode==1
      label=[name ' (' unit ')'];
    elseif mode==2
      label={name, ['(' unit ')']};
    end
  end
end