function [IMGInd]=loaddata_IMAGEMAG_IL_level0(varargin)
  
  if isempty(varargin)
    dn=datenum([2011 12 18]);
  else
    dn=varargin{1};
  end
  
  global datasetinfo
  pwd1=pwd;
  cd ..
  fp_root=[pwd '/data/IL/'];
  cd(pwd1)
  
  [yy, mm, dd, HH, MM, SS]=datevec(dn);
  syy=sprintf('%04d',yy);
  syy_2=syy(3:4);
  smm=sprintf('%02d',mm);
  sdd=sprintf('%02d',dd);
  
  fdir=fp_root;
  flist=dir([fdir '*.dat']);
  fnlist={flist.name};
  ix=regexp(fnlist,['IL' syy smm sdd]);
  ix=~cellfun('isempty',ix);
  fid=-1;
  if isempty(find(ix==1))
    disp(['Warning: IL' syy smm sdd '.dat is not found!'])
  else
    fn=fnlist(ix);
    fp_dat=fdir;
    fn_dat=fn{1};
    fid=fopen([fp_dat fn_dat]);
  end   
  
  if fid==-1
    IL=[];
    IU=[];
    IE=[];
    tl=[];
  else
    C=textscan(fid,'%d %d %d %d %d %d %f %f %f','commentstyle','%');
    fclose(fid);
    dnlist=datenum(double(cell2mat(C(1,1:6))));
    IL=C{7};
    IU=C{8};
    IE=C{9};
    tl=dnlist';
  end
  
  IMGInd.IL=IL';
  IMGInd.IU=IU';
  IMGInd.IE=IE';
  IMGInd.tl=tl;
end