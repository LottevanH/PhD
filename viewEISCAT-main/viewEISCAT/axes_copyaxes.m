function [ ha_new ] = axes_copyaxes( hf,ha )
%AXES_COPYAXES Summary of this function goes here
%   Detailed explanation goes here
  set(0,'currentfigure',hf)
  pos=get(ha,'position');
  xlim=get(ha,'xlim');
  ylim=get(ha,'ylim');
  ha_new=axes('position',pos,'xlim',xlim,'ylim',ylim);
end

