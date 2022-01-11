function [lim]=axes_axis_limits(lim_in,lim_setting)

lim=[];
if isempty(lim_setting)
  lim = lim_in;
  return
end
if lim_setting(1)==-inf && lim_setting(2)==inf
  limmax=max(abs(lim_in));
  lim=[-limmax limmax];
elseif lim_setting(1)==0 && lim_setting(2)==inf && lim_in(2)>0
  lim=[0 lim_in(2)];
elseif  lim_setting(1)==-inf && lim_setting(2)==0 && lim_in(1)<0
  lim=[lim_in(1) 0];
elseif isfinite(lim_setting(1)) && isfinite(lim_setting(1))
  lim=lim_setting;
else
  lim=lim_in;
end
end