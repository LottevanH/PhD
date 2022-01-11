function PCs=loaddata_PCs(varargin)

  if isempty(varargin)
    dn=datenum([2011 12 18]);
  else
    dn=varargin{1};
  end
  
  global datasetinfo
  pwd1=pwd;
  cd ..
  fp_root=[pwd '/data/PC/'];
  cd(pwd1)
  
  [yy, mm, dd, HH, MM, SS]=datevec(dn);
  syy=sprintf('%04d',yy);
  syy_2=syy(3:4);
  smm=sprintf('%02d',mm);
  sdd=sprintf('%02d',dd);
  
  fdir=fp_root;
  flist=dir([fdir '*.mat']);
  fnlist={flist.name};
  ix=regexp(fnlist,['PCs_' syy]);
  ix=~cellfun('isempty',ix);
  fid=-1;
  if isempty(find(ix==1))
    disp(['Warning: PCs_' syy '.mat is not found!'])
  else
    fn=fnlist(ix);
    fp_dat=fdir;
    fn_dat=fn{1};
  end   
  
  load([fp_dat fn_dat])
  
  ind=find(floor(dnlist)==floor(dn));
  
  PCs.N=PCN(ind)';
  PCs.S=PCS(ind)';
  PCs.tl=dnlist(ind)';
  
end