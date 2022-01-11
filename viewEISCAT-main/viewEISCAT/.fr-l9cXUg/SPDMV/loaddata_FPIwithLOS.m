function FPIdata=loaddata_FPIwithLOS(dn,emline)
  
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
  
  %% calculate LOS velocities
  mfilepath=[pwd '/FPI/'];
  addpath(mfilepath);
  [ beam1.val, beam1.dev, beam1.err, beam1.tl,  ...
    beam2.val, beam2.dev, beam2.err, beam2.tl,  ...
    beam3.val, beam3.dev, beam3.err, beam3.tl,  ...
    beam4.val, beam4.dev, beam4.err, beam4.tl,  ...
    beam5.val, beam5.dev, beam5.err, beam5.tl]  ...
           = func_AnaFPILOSwind_SAP4FPI([fp_dat fn_dat]);
  %% assign to mother structure
  FPIdata.uN=uN;       
  FPIdata.uE=uE;
  FPIdata.uh=uh;
  FPIdata.beam1=beam1;
  FPIdata.beam2=beam2;
  FPIdata.beam3=beam3;
  FPIdata.beam4=beam4;
  FPIdata.beam5=beam5;
end