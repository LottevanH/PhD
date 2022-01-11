function assign_OMNI()
global datasetinfo dataset

%fp_data=['/home/leicai/Projects/Analysis/J_eq/'];

dateran=datasetinfo.dateran;

pnlist={
          'Bx_gsm',     'OMNI.Bx_gsm',      'OMNI.tl';        ...
          'By_gsm',     'OMNI.By_gsm',      'OMNI.tl';        ...
          'Bz_gsm',     'OMNI.Bz_gsm',      'OMNI.tl';        ...
          'Btot_gsm',   'OMNI.Btot_gsm',    'OMNI.tl';        ...
          'v_sw',       'OMNI.v_sw',        'OMNI.tl';        ...
          'E_sw',       'OMNI.E_sw',        'OMNI.tl';        ...
          'Temp',       'OMNI.Temp',        'OMNI.tl';        ...
          'N_p',        'OMNI.N_p',         'OMNI.tl';        ...
          'Pdyn',       'OMNI.Pdyn',        'OMNI.tl';        ...
          'AE',         'OMNI.AE',          'OMNI.tl';        ...
          'AL',         'OMNI.AL',          'OMNI.tl';        ...
          'AU',         'OMNI.AU',          'OMNI.tl';        ...
          'PC',         'OMNI.PC',          'OMNI.tl'};

pararec=datasetinfo.filereader(datasetinfo.filereaderrec).pararec;
npara=length(pararec);
para=cell(npara,1);
for i=1:npara
  para{i}.val=[];
  para{i}.tl=[];
end

%% load data
for dn=floor(dateran(1)):floor(dateran(2))
  if dn==floor(dateran(2)) && dateran(2)==floor(dateran(2)) 
    continue
  end
  
  OMNI=loaddata_OMNI(dn,datasetinfo.paralist(pararec));
  
  for i=1:length(pararec)
    ix=regexp(pnlist(:,1),datasetinfo.paralist(pararec(i)));
    ix=~cellfun('isempty',ix);
    val=eval(pnlist{ix,2});
    tl=eval(pnlist{ix,3}); 
    
    para{i}.val=[para{i}.val val];
    para{i}.tl=[para{i}.tl tl];
  end
end
%% restrict time range
for i=1:length(pararec)
  [tl,indtst,indted]=ts_confinetimeline(para{i}.tl,dateran(1),dateran(2));
  para{i}.val=para{i}.val(:,indtst:indted);
  para{i}.tl=tl;
end
%% Assign values to dataset
names=fieldnames(para{i});
for i=1:length(pararec)
  for j=1:length(names)
    dataset{pararec(i)}.(names{j})=para{i}.(names{j});
  end
end
end