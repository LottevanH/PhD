function [NaN]=assign_ElSpecGLOW()
global datasetinfo dataset

dateran=datasetinfo.dateran;

pnlist={'nu_ni', 'CDT.nu_ni.val', 'CDT.nu_ni.err', 'CDT.tl', 'CDT.alt';
        'tau_ni', 'CDT.tau_ni.val', 'CDT.tau_ni.err', 'CDT.tl', 'CDT.alt'};

pararec=datasetinfo.filereader(datasetinfo.filereaderrec).pararec;
npara=length(pararec);
para=cell(npara,1);
for i=1:npara
  para{i}.val=[];
  para{i}.tl=[];
  para{i}.err=[];
  para{i}.alt = [];
end

%% load data
for dn=floor(dateran(1)-0.5):floor(dateran(2)-0.5)
  if dn==floor(dateran(2)) && dateran(2)==floor(dateran(2)) 
    continue
  end
  
  [CDT]=loaddata_EISCAT_cdt(dn);
  
  if isempty(CDT)
    continue
  end
  for i=1:length(pararec)
    ix=regexp(pnlist(:,1),datasetinfo.paralist(pararec(i)));
    ix=~cellfun('isempty',ix);
    val=eval(pnlist{ix,2});
    tl=eval(pnlist{ix,4});
    err=eval(pnlist{ix,3});
    if isempty(pnlist{ix, 5})
      alt = [];
    else
      alt=eval(pnlist{ix,5});
    end
    
    para{i}.val=[para{i}.val val];
    para{i}.tl=[para{i}.tl tl];
    para{i}.err=[para{i}.err err];
    para{i}.alt=alt;

  end
end
%% restrict time range
for i=1:length(pararec)
  if isempty(para{i}.tl); continue; end
  [tl,indtst,indted]=ts_confinetimeline(para{i}.tl,dateran(1),dateran(2));
  para{i}.val=para{i}.val(:,indtst:indted);
  if ~isempty(para{i}.err)
      para{i}.err=para{i}.err(:,indtst:indted);
  end
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