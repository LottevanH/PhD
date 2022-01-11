function [ cmp ] = ColorWheel( cs )
%KBGYRM Summary of this function goes here
%   Detailed explanation goes here
%   f=[ 80   0   0   0   0 255 255 255 105 255 255
%        0   0 255  80 255 255 255 77   77   0  69
%       80 255 255   0   0 255   0   0  20   0  84 ]'/255;
f=[  0 156   0
     0 174 104
     0 147 154 
    11  95 165
    17  52 171
    28  26 178
    55  20 176
    84  15 172
   115   9 170
   170   0 162
   205   0 116
   228   0  69
   255   0   0
   255  74   0
   255 112   0
   255 144   0
   255 170   0
   255 193   0
   255 213   0
   255 255   0
   207 247   0
   155 237   0
   103 227   0
     0 156   0
  ]/255;
if strcmp(cs,'positive')
  cmp=f(6:11,:);
elseif strcmp(cs,'negative')
  cmp=f(1:6,:);
else
  cmp=f;
end   

end
