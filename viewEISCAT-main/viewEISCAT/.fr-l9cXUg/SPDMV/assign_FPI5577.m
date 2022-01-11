%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% read FPI data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NaN]=assign_FPI5577()
global dataset datasetinfo

emline='5577';
dateran=datasetinfo.dateran;


pnlist={'u55_N', 'FPI.uN.val','FPI.uN.err','FPI.uN.tl';     ...
        'u55_E', 'FPI.uE.val','FPI.uE.err','FPI.uE.tl';     ...
        'u55_h', 'FPI.uh.val','FPI.uh.err','FPI.uh.tl';     ...
        'u55_bm1','FPI.beam1.val','FPI.beam1.err','FPI.beam1.tl';   ...
        'u55_bm2','FPI.beam2.val','FPI.beam2.err','FPI.beam2.tl';   ...
        'u55_bm3','FPI.beam3.val','FPI.beam3.err','FPI.beam3.tl';   ...
        'u55_bm4','FPI.beam4.val','FPI.beam4.err','FPI.beam4.tl';   ...
        'u55_bm5','FPI.beam5.val','FPI.beam5.err','FPI.beam5.tl';   ...
        'a55_N','FPI.aN.val','FPI.aN.err','FPI.aN.tl';      ...
        'a55_E','FPI.aE.val','FPI.aE.err','FPI.aE.tl'};
      
pararec=datasetinfo.filereader(datasetinfo.filereaderrec).pararec;
npara=length(pararec);
para=cell(npara,1);
for i=1:npara
  para{i}.val=[];
  para{i}.err=[];
  para{i}.tl=[];
end

%% check if los velocity need to be calculated
ind_los=ismember(datasetinfo.paralist(pararec(:)),pnlist(4:8,1));
if ~isempty(find(ind_los)>0)
  los=1;
else
  los=0;
end

%% check if acceleration rate need to be calculated
ind_los=ismember(datasetinfo.paralist(pararec(:)),pnlist(9:10,1));
if ~isempty(find(ind_los)>0)
  au=1;
else
  au=0;
end

%% load data
for dn=floor(dateran(1)-0.5):floor(dateran(2)-0.5)
  if los
    FPI=loaddata_FPIwithLOS(dn,emline);
  else
    FPI=loaddata_FPI(dn,emline);
  end
  
  if au
    FPI.aE=cal_FPIa(FPI.uE.val, FPI.uE.tl*3600);
    FPI.aE.tl=FPI.aE.tl/3600;
    FPI.aN=cal_FPIa(FPI.uN.val, FPI.uN.tl*3600);
    FPI.aN.tl=FPI.aN.tl/3600;
  end
    
  for i=1:length(pararec)
      ix=regexp(pnlist(:,1),datasetinfo.paralist(pararec(i)));
      ix=~cellfun('isempty',ix);
      val=eval(pnlist{ix,2});
      err=eval(pnlist{ix,3});
      tl=dn+eval(pnlist{ix,4})/24;
      if size(val,2)==1; val=val';end
      if size(err,2)==1; err=err';end 
      if size(tl,2)==1; tl=tl';end 
      
      para{i}.val=[para{i}.val val];
      para{i}.err=[para{i}.err err];
      para{i}.tl=[para{i}.tl tl];
  end
  
end
%% restrict time range
for i=1:length(pararec)
  [tl,indtst,indted]=ts_confinetimeline(para{i}.tl,dateran(1),dateran(2));
  if ~isempty(para{i}.val)
    para{i}.val=para{i}.val(:,indtst:indted);
  end
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

