function [ cmp ] = ColorWind( cs )
%KBGYRM Summary of this function goes here
%   Detailed explanation goes here
  f=[ 80   0   0   0   0 255 255 255 105 255 255
       0   0 255  80 255 255 255 77   77   0  69
      80 255 255   0   0 255   0   0  20   0  84 ]'/255;
if strcmp(cs,'positive')
  cmp=f(6:11,:);
elseif strcmp(cs,'negative')
  cmp=f(1:6,:);
else
  cmp=f;
end   

end
