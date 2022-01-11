function [NaN]=assign_IMAGEMAG_EQ2D()
global datasetinfo dataset

fp_data=['/home/leicai/Projects/Analysis/J_eq/'];

dateran=datasetinfo.dateran;

pnlist={'EQ1Eext', 'EQ1D.Elat.val','EQ1D.Elat.tl','EQ1D.Elat.lat','EQ1D.Elat.lon';     ...
  'EQ1Eext_mono', 'EQ1D.Elat.monosite.val','EQ1D.Elat.tl','EQ1D.Elat.monosite.lat','EQ1D.Elat.monosite.lon';    ...
  'EQ1Next', 'EQ1D.Nlat.val','EQ1D.Nlat.tl','EQ1D.Nlat.lat','EQ1D.Nlat.lon';     ...
  'EQ1Next_mono', 'EQ1D.Nlat.monosite.val','EQ1D.Nlat.tl','EQ1D.Nlat.monosite.lat','EQ1D.Nlat.monosite.lon'};

pararec=datasetinfo.filereader(datasetinfo.filereaderrec).pararec;
npara=length(pararec);
para=cell(npara,1);
for i=1:npara
  para{i}.val=[];
  para{i}.tl=[];
  para{i}.lat=[];
  para{i}.lon=[];
end

%% load data
for dn=floor(dateran(1)):floor(dateran(2))
  if dn==floor(dateran(2)) && dateran(2)==floor(dateran(2)) 
    continue
  end
  [yy, mm, dd, HH, MM, SS]=datevec(dn);
  syy=sprintf('%04d',yy);
  syy_2=syy(3:4);
  smm=sprintf('%02d',mm);
  sdd=sprintf('%02d',dd);
  
  %   fdir=[root_dat];
  %   flist=dir([fdir '*.mat']);
  %   fnlist={flist.name};
  %   ix=regexp(fnlist,['wind_' syy_2 smm sdd]);
  %   ix=~cellfun('isempty',ix);
  %   fn=fnlist(ix);
  %
  %   fp_dat=fdir;
  %   fn_dat=fn{1};
  
  fn_data=['Jeq_' datestr(dn,'yyyymmdd') '_Heikki_ExtOnly.mat'];
%   if isempty(dir([fp_data fn_data]))
%     continue
%   end
  
  [EQ1D]=loaddata_IMAGEMAG_EQ2D(fp_data,fn_data);
  
  for i=1:length(pararec)
    ix=regexp(pnlist(:,1),datasetinfo.paralist(pararec(i)));
    ix=~cellfun('isempty',ix);
    val=eval(pnlist{ix,2});
    tl=eval(pnlist{ix,3});
    lat=eval(pnlist{ix,4});
    lon=eval(pnlist{ix,5}); 
    
    para{i}.val=[para{i}.val val];
    para{i}.tl=[para{i}.tl tl];
    para{i}.lat=lat;
    para{i}.lon=lon;
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