global database

drawopttemplate=vis_drawopttemplates(1);

% NOTE: Add the parameter information whenever a new parameters is added!

parainfo={'IL',  ...
          'IU',  ...
          'IE'};
        
%% Item 1: IL
itemname='IL';
database.(itemname).drawopt=drawopttemplate;
item={                                  ...
      'xlabel',         'UT';       	...
      'ylim',           [];     ...
      'yscale',         'linear'};
for i=1:size(item,1)
  database.(itemname).drawopt.(item{i,1})=item{i,2};
end

%% Item 2: IU
itemname='IU';
database.(itemname).drawopt=drawopttemplate;
item={                                  ...
      'xlabel',         'UT';       	...
      'ylim',           [];     ...
      'yscale',         'linear'};
for i=1:size(item,1)
  database.(itemname).drawopt.(item{i,1})=item{i,2};
end

%% Item 2: IE
itemname='IE';
database.(itemname).drawopt=drawopttemplate;
item={                                  ...
      'xlabel',         'UT';       	...
      'ylim',           [];       ...
      'yscale',         'linear'};  
for i=1:size(item,1)
  database.(itemname).drawopt.(item{i,1})=item{i,2};
end