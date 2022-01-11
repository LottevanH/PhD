function [ output_args ] = vis_init_EISCAT_lv0( input_args )
%VIS_INIT_EISCAT_LV0 Summary of this function goes here
%   Detailed explanation goes here

global database datasetinfo

drawopttemplate1=vis_drawopttemplates(1);
drawopttemplate2=vis_drawopttemplates(2);
% NOTE: Add the parameter information whenever a new parameters is added!

%% EISCAT level 0 data
parainfo={'EF_E_NWIND',            ...
          'EF_N_NWIND'};
        
cmap=axes_colormaps(0);       
%% Item 1: Eastward electric field
itemname = 'EF_E_NW';
database.(itemname).drawopt=drawopttemplate1;
item={                                  ...
      'xlabel',         'UT';       	...
      'ylim',           [];     ...
      'yscale',         'linear'};
for i=1:size(item,1)
  database.(itemname).drawopt.(item{i,1})=item{i,2};
end

%% Item 2: Northward electric field
itemname = 'EF_N_NW';
database.(itemname).drawopt=drawopttemplate1;
item={                                  ...
      'xlabel',         'UT';       	...
      'ylim',           [];     ...
      'yscale',         'linear'};
for i=1:size(item,1)
  database.(itemname).drawopt.(item{i,1})=item{i,2};
end

%%
end

