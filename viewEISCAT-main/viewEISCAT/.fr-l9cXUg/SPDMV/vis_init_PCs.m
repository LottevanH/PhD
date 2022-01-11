global database

drawopttemplate=vis_drawopttemplates(1);

% NOTE: Add the parameter information whenever a new parameters is added!

parainfo={'PCN',  ...
          'PCS'};
        
%% Item 1: PCN
itemname='PCN';
database.(itemname).drawopt=drawopttemplate;
item={                                  ...
      'xlabel',         'UT';       	...
      'ylim',           [];     ...
      'yscale',         'linear'};
for i=1:size(item,1)
  database.(itemname).drawopt.(item{i,1})=item{i,2};
end

%% Item 2: PCS
itemname='PCS';
database.(itemname).drawopt=drawopttemplate;
item={                                  ...
      'xlabel',         'UT';       	...
      'ylim',           [];     ...
      'yscale',         'linear'};
for i=1:size(item,1)
  database.(itemname).drawopt.(item{i,1})=item{i,2};
end
