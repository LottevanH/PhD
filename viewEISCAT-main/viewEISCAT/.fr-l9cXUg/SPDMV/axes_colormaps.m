function [cmap]=vis_colormaps(go)
  addpath([pwd '/colormaps/']);
  if ~go
    cmap{1}='jet';
  else
    cmap{1}=vis_colorpalette(1,'whole',80,1);
    cmap{2}=myb(7,'whole',120,1); % wind magnitude
    cmap{3}=myb(4,'whole',100,0); % neutral winds
    cmap{4}=myb(6,'whole',180,0); % wind clock angles?
  end
end