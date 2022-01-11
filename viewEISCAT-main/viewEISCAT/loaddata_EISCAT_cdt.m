function CDT = loaddata_EISCAT_cdt(dn)
  
%% find filenames
pwd1=pwd;
cd ..
fp_root=[pwd '/EISCAT_cdt/res_cdt'];
cd(pwd1)

dstr = datestr(dn, 'yyyymmdd');

dstr1 = datestr(dn, 'yyyymm');

flist = dir([fp_root  '/TRO'  dstr1 '/*.mat']);
fnlist = {flist.name};
ix = regexp(fnlist, ['cdt_' dstr '_TRO_T_ie']);
ix = ~cellfun('isempty', ix);
fn = fnlist(ix);

fpfn = fullfile(fp_root, ['TRO' dstr1], fn{1});
eval(['load ' fpfn]);

tl = (datenum(t1') + datenum(t2'))/2;
CDT.tl = tl';
CDT.alt = height_c';

CDT.nu_ni.val = vni;
CDT.nu_ni.err = [];
CDT.tau_ni.val = tau_ni/3600.;
CDT.tau_ni.err = [];

% EF_tl = (datenum(EF_time1') + datenum(EF_time2'))/2;
% NWIND = [];
% scale = 1e3;
% NWIND.EF_E.val = EF_E'*scale;
% NWIND.EF_E.err = sigma_EF_E'*scale;
% NWIND.EF_E.tl = EF_tl';
% NWIND.EF_N.val = EF_N'*scale;
% NWIND.EF_N.err = sigma_EF_N'*scale;
% NWIND.EF_N.tl = EF_tl';

end
