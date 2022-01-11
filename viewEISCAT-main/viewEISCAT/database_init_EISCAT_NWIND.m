global database datasetinfo

% NOTE: Add the parameter information whenever a new parameters is added!

databasetemplate1=database_variabletemplates(1);
databasetemplate2=database_variabletemplates(2);

%% Eastward electric field
parainfo={'EF_E_NW',            ...
          'EF_N_NW'};
ix=ismember(parainfo,datasetinfo.paralist);
if ~isempty(parainfo(ix))
  % linked file for reading EISCAT data
  filereader='assign_EISCAT_NWIND.m';
  
  %% Item 1: Electron density
  itemname='EF_E_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'E_E';    ...
    'unit',         'mV/m';          ...
    'ndim',          1;              ... 
    'group',        'EF';    ...
    'dscp',         'Electric field (Eastward)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end
  
  %% Item 2: Northward electric field
  itemname='EF_N_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'E_N';    ...
    'unit',         'mV/m';          ...
    'ndim',          1;              ... 
    'group',        'EF';    ...
    'dscp',         'Electric field (Northward)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end
  
  %% Item 3: uEh_E
  itemname='uEh_E_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'u_E';    ...
    'unit',         'm/s';          ...
    'ndim',          2;              ... 
    'group',        'u';    ...
    'dscp',         'nuetral wind (Eastward)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end
  
    %% Item 4: uEh_N
  itemname='uEh_N_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'u_N';    ...
    'unit',         'm/s';          ...
    'ndim',          2;              ... 
    'group',        'u';    ...
    'dscp',         'nuetral wind (N)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end
  
  %% Item 5: uEh_h
  itemname='uEh_h_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'u_h';    ...
    'unit',         'm/s';          ...
    'ndim',          2;              ... 
    'group',        'u';    ...
    'dscp',         'nuetral wind (h)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end
  
    %% Item 6: uEz_E
  itemname='uEz_E_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'u_E';    ...
    'unit',         'm/s';          ...
    'ndim',          2;              ... 
    'group',        'u';    ...
    'dscp',         'nuetral wind (Eastward)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end
  
    %% Item 7: uEz_N
  itemname='uEz_N_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'u_N';    ...
    'unit',         'm/s';          ...
    'ndim',          2;              ... 
    'group',        'u';    ...
    'dscp',         'nuetral wind (N)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end
  
  %% Item 8: uEz_z
  itemname='uEz_z_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'u_h';    ...
    'unit',         'm/s';          ...
    'ndim',          2;              ... 
    'group',        'u';    ...
    'dscp',         'nuetral wind (h)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end
  
    %% Item 9: vEh_E
  itemname='vEh_E_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'vi_E';    ...
    'unit',         'm/s';          ...
    'ndim',          2;              ... 
    'group',        'v';    ...
    'dscp',         'ion velcity (Eastward)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end
  
    %% Item 10: vEh_N
  itemname='vEh_N_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'v_N';    ...
    'unit',         'm/s';          ...
    'ndim',          2;              ... 
    'group',        'v';    ...
    'dscp',         'ion velocity (N)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end
  
  %% Item 11: vEh_h
  itemname='vEh_h_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'v_h';    ...
    'unit',         'm/s';          ...
    'ndim',          2;              ... 
    'group',        'v';    ...
    'dscp',         'ion velocity (h)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end
  
  %% Item 12: vEz_E
  itemname='vEz_E_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'vi_E';    ...
    'unit',         'm/s';          ...
    'ndim',          2;              ... 
    'group',        'v';    ...
    'dscp',         'ion velcity (Eastward)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end
  
    %% Item 13: vEz_N
  itemname='vEz_N_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'v_N';    ...
    'unit',         'm/s';          ...
    'ndim',          2;              ... 
    'group',        'v';    ...
    'dscp',         'ion velocity (N)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end
  
  %% Item 14: vEz_z
  itemname='vEz_z_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'v_h';    ...
    'unit',         'm/s';          ...
    'ndim',          2;              ... 
    'group',        'v';    ...
    'dscp',         'ion velocity (h)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end  
  
    %% Item 15: vFh_E
  itemname='vFh_E_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'vi_E';    ...
    'unit',         'm/s';          ...
    'ndim',          1;              ... 
    'group',        'v';    ...
    'dscp',         'ion velcity (Eastward)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end
  
    %% Item 16: vFh_N
  itemname='vFh_N_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'v_N';    ...
    'unit',         'm/s';          ...
    'ndim',          1;              ... 
    'group',        'v';    ...
    'dscp',         'ion velocity (N)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end
  
  %% Item 17: vFh_h
  itemname='vFh_h_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'v_h';    ...
    'unit',         'm/s';          ...
    'ndim',          1;              ... 
    'group',        'v';    ...
    'dscp',         'ion velocity (h)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end
  
  %% Item 18: vFz_E
  itemname='vFz_E_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'vi_E';    ...
    'unit',         'm/s';          ...
    'ndim',          1;              ... 
    'group',        'v';    ...
    'dscp',         'ion velcity (Eastward)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end
  
    %% Item 19: vFz_N
  itemname='vFz_N_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'v_N';    ...
    'unit',         'm/s';          ...
    'ndim',          1;              ... 
    'group',        'v';    ...
    'dscp',         'ion velocity (N)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end
  
  %% Item 20: vFz_z
  itemname='vFz_z_NW';
  database.(itemname)=databasetemplate2;
  item={                            ...
    'name',         itemname;        ...
    'label',        'v_h';    ...
    'unit',         'm/s';          ...
    'ndim',          1;              ... 
    'group',        'v';    ...
    'dscp',         'ion velocity (h)';      ...
    'filereader',   filereader      ...
    };
  for i=1:size(item,1)
    database.(itemname).(item{i,1})=item{i,2};
  end  
  %% Initialize plot settings
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  vis_init_EISCAT_NWIND;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end