function assign_PCs()
global datasetinfo dataset

%fp_data=['/home/leicai/Projects/Analysis/J_eq/'];

dateran=datasetinfo.dateran;

pnlist={'PCN','PCs.N','PCs.tl';       ...
        'PCS','PCs.S','PCs.tl'};

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
  
  PCs=loaddata_PCs(dn);
%   [yy, mm, dd, HH, MM, SS]=datevec(dn);
%   syy=sprintf('%04d',yy);
%   syy_2=syy(3:4);
%   smm=sprintf('%02d',mm);
%   sdd=sprintf('%02d',dd);
  
  %   fdir=[root_dat];
  %   flist=dir([fdir '*.mat']);
  %   fnlist={flist.name};
  %   ix=regexp(fnlist,['wind_' syy_2 smm sdd]);
  %   ix=~cellfun('isempty',ix);
  %   fn=fnlist(ix);
  %
  %   fp_dat=fdir;
  %   fn_dat=fn{1};
  
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
  if ~isempty(para{i}.tl)
  [tl,indtst,indted]=ts_confinetimeline(para{i}.tl,dateran(1),dateran(2));
  para{i}.val=para{i}.val(:,indtst:indted);
  para{i}.tl=tl;
  end
end
%% Assign values to dataset
names=fieldnames(para{i});
for i=1:length(pararec)
  for j=1:length(names)
    dataset{pararec(i)}.(names{j})=para{i}.(names{j});
  end
end
end