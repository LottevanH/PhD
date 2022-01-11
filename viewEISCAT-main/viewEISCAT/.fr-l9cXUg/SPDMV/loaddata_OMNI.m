function OMNI=loaddata_OMNI(dn,paranames)

%% find filenames
pwd1=pwd;
cd ..
fp_root=pwd;
cd(pwd1)

% add path
addpath([fp_root '/viewOMNI/'])
fp_lv0=[fp_root '/viewOMNI/data/OMNI_monthly_1min_Lv0/'];
fp_lv1=[fp_root '/viewOMNI/data/OMNI_1min_Lv1/'];
fp_lv2_20m=[fp_root '/viewOMNI/data/OMNI_avg_20min_Lv2/'];
fp_lv2_60m=[fp_root '/viewOMNI/OMNI_avg_60min_Lv2/'];

[yy, mm, dd]=datevec(dn);
syy=sprintf('%4d',yy);smm=sprintf('%02d',mm);sdd=sprintf('%02d',dd);


fn_dat=['OMNI_1min_' syy smm '_Lv1.mat'];
fp_dat=fp_lv1;
f_stat=dir([fp_lv1 fn_dat]);
if isempty(f_stat)
  disp('Wait a moment: Collecting OMNI data!')
  fn_dat=OMNI_read12(fp_lv0, fp_dat, yy, mm);
end

% remove path
rmpath([fp_root '/viewOMNI/']);



%% load data
eval(['load ' fp_dat fn_dat ';'])

%% assign parameters
%mat2cell(scomp,ones(1,46),5)

%  scomp=['ID MF';'ID PL';'n  MF';'n  PL';'Inter';'Dt   ';'RMSDt';'RMSMi';'DBOT ';'<B>  ';'BxGSE';'ByGSE';'BzGSE';'ByGSM';'BzGSM';'RMSBs';'RMSBv';'vel  ';'vxGSE';'vyGSE';'vzGSE';...
%        'n    ';'Temp ';'Pdyn ';'Ey   ';'beta ';'MA   ';'XGSE ';'YGSE ';'ZGSE ';'XBow ';'YBow ';'ZBow ';'AE   ';'AL   ';'AU   ';'SYM/D';'SYM/H';'ASY/D';...
%        'ASY/H';'PC-N ';'Mms  ';'clock';'New07';'Bt   ';'E_sw '];
pnlist={'Bx_gsm',   11;     ...
        'By_gsm',   14;     ...
        'Bz_gsm',   15;     ...
        'Btot_gsm'   ,   45;     ...
        'vx_gse',   19;     ...
        'vy_gse',   20;     ...
        'vz_gse',   21;     ...
        'v_sw',     18;     ...
        'E_sw',     46;     ...
        'Temp',     23;     ...
        'Pdyn',     24;     ...
        'n_p',      22;     ...
        'cagl',     43;     ...
        'New07',    44;     ...
        'AE',       34;     ...
        'AL',       35;     ...
        'AU',       36;     ...
        'PC',       41};

npara=length(paranames);
para=cell(npara,1);
for i=1:npara
  para{i}.val=[];
end

dn_st=datenum(sdate,'yyyy-mm-dd');
dnlist=edesh'/86400+dn_st;
comp1=comp(floor(dnlist)==dn,:)';
for i=1:npara
  ix=regexp(pnlist(:,1),paranames{i});
  ix=~cellfun('isempty',ix);
  OMNI.(paranames{i})=comp1(pnlist{ix,2},:);
end
OMNI.tl=dnlist(floor(dnlist)==dn);
end