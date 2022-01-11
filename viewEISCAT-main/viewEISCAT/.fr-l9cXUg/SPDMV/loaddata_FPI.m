function FPIdata=loaddata_FPI(dn,emline)
  
  root_dat='/mnt/eiscatshin7/analysis/fp01';
  
  [yy, mm, dd, HH, MM, SS]=datevec(dn);
  syy=sprintf('%04d',yy);
  syy_2=syy(3:4);
  smm=sprintf('%02d',mm);
  sdd=sprintf('%02d',dd);
    
  fdir=[root_dat '/' syy '/' emline '/'];
  flist=dir([fdir '*.mat']);
  fnlist={flist.name};
  ix=regexp(fnlist,['wind_' syy_2 smm sdd]);
  ix=~cellfun('isempty',ix);
  fn=fnlist(ix);
  
  fp_dat=fdir;
  fn_dat=fn{1};
    
  eval(['load ' fp_dat fn_dat ';']) 
  
  %% assign wind components
  uN.val=MeanWindArr_north;
  uN.err=StdWindArr_north;
  uN.tl=TimeArr_north;
  uE.val=MeanWindArr_east;
  uE.err=StdWindArr_east;
  uE.tl=TimeArr_east;
  uh.val=MeanWindArr_up;
  uh.err=StdWindArr_up;
  uh.tl=TimeArr_up;
  
  %% assign to mother structure
  FPIdata.uN=uN;       
  FPIdata.uE=uE;
  FPIdata.uh=uh;
end