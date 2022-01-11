%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% database for IL, IU, IE index from FMI/IMAGE    %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global database datasetinfo drawopt

% NOTE: Add the parameter information whenever a new parameters is added!


databasetemplate=database_variabletemplates(1);

parainfo={'PCN',  ...
          'PCS'};
ix=ismember(parainfo,datasetinfo.paralist);
if ~isempty(parainfo(ix))
  % linked file for reading FPI data
  filereader='assign_PCs.m';
  %% Item 1: PCN
  itemname='PCN';
  database.(itemname)=databasetemplate;
  item={                            ...
    'name',         'PCN';        ...
    'label',        'PCN';          ...
    'unit',         'mV/m';          ...
    'ndim',          1;              ...
    'group',        'PCs';    ...
    'dscp',         'PC index (Northern hemisphere)';      ...
    'filereader',   filereader};
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end
  %% Item 2: PCS
  itemname='PCS';
  database.(itemname)=databasetemplate;
  item={                            ...
    'name',         'PCS';        ...
    'label',        'PCS';          ...
    'unit',         'mV/m';          ...
    'ndim',          1;              ...
    'group',        'PCs';    ...
    'dscp',         'PC index (Southern hemisphere)';      ...
    'filereader',   filereader};
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end 
  
  %% Initialize plot settings
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  vis_init_PCs;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end