function [ cmp ] = KBGYRM( cs )
%KBGYRM Summary of this function goes here
%   Detailed explanation goes here
%   f=[0 0 0 0 0 1 2 2 2 2 2 2 2
%      0 0 0 1 2 2 2 1 0 0 0 1 2
%      0 1 2 1 0 0 0 0 0 1 2 2 2]'/2;
f=[  0   0   0
     0   0 153
     0   0 255
     0 102 255
     0 255 255
     0 255 102
     0 255   0
   153 255   0
   255 255   0
   255 102   0
   255   0   0
   255   0 102
   255   0 153
   255   0 255
   255 153 255
   255 255 255
   ]/255;
   
if strcmp(cs,'positive')
  cmp=f(6:13,:);
elseif strcmp(cs,'negative')
  cmp=f(1:6,:);
else
  cmp=f;
end   

end

