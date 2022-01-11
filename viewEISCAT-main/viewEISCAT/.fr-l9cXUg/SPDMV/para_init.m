%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Initialize parameters %%%%%%%%%%%%%%%%
function [NaN]=para_init
global dataset datasetinfo database 
% fn_db='databasebib.m';
% fnl_dbb=dir('*.dbb');
% for i=1:length(fnl_dbb);
%   fdbb_name{i}=fnl_dbb(i).name;
%   fdbb_info=dir(fdbb_name{i});
%   fdbb_dn[]=datenum(fdbb_info.date);
% end
% f_db=dir(fn_dn);
% if isempty(f_db)
%   loaddatabaseinfo;
% else
%   load(f_db);
%   fdbb_date=databaseinfo.dn;
%   f_dbbib=dir'';
%% load database

%databaseinfo=[];  

% Initialize values in the database
flist=dir('*.m');
fnlist={flist.name};
ix=regexp(fnlist,'database_init_');
ix=~cellfun('isempty',ix);
fn=fnlist(ix);
for i=1:length(fn)
  [pathstr, name, ext] = fileparts(fn{i});
  eval(name)
end

%% find parameters in the database
paralist=datasetinfo.paralist;
npara=length(paralist);
readerrec={};
pararec={};

dataset=cell(npara,1);
for i=1:npara
  if ~isfield(database,paralist{i})
    errordlg(['''' paralist{i} '''' ' not found in the database'],  ...
      'Initialization error');
    continue;
  end
  eval(['dataset{i}=database.' paralist{i} ';'])
  
  rec=ismember(readerrec, dataset{i}.filereader);
  
  if isempty(readerrec(rec))
    readerrec=[readerrec {dataset{i}.filereader}];
    pararec=[pararec {i}];
  else
    rloc=find(rec==1);
    pararec{rloc}=[pararec{rloc} i];
  end
end

for i=1:length(readerrec)
  datasetinfo.filereader(i).name=readerrec{i};
  datasetinfo.filereader(i).pararec=pararec{i};
end
        
end