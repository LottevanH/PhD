function [ cmp ] = DarkKBGYRM( opt )
%DARKKBGYRM Summary of this function goes here
%   Detailed explanation goes here

%   cmp=flipud([0 48  16  39  0   84  85  139 139 135 85  125 155 255 
%               0 48  78  64  128 139 107 90  71  38  26  38  48  255
%               0 48  139 139 128 84  47  0   93  87  139 205 255 255 ]'/255);
%   cmp=flipud([0 0 0 0 0 1 2 2 2 2 2 2 
%               0 0 0 1 2 2 2 1 0 0 0 1 
%               0 1 2 1 0 0 0 0 0 1 2 2 ]'/5);
cmp=flipud([   0   0   0
       25  25  51
       25  25 102
        0  51 102
        0 102 102
        0 102  51
        0 102   0
       51 102   0
      102 102   0
      102  51   0
      102   0   0
      102   0  51
      102   0 102
      102  51 102
      102  51 153
     ]/255/1.5);
end

