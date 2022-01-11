function NWIND = loaddata_EISCAT_NWIND(dn)
  
%% find filenames
pwd1=pwd;
cd ..
fp_root=[pwd '/Data/EISCAT/lv2_NWIND'];
cd(pwd1)

dstr = datestr(dn, 'yyyy-mm-dd');

flist = dir([fp_root, '/*.mat']);
fnlist = {flist.name};
ix = regexp(fnlist, ['u_ef_', dstr]);
ix = ~cellfun('isempty', ix);
fn = fnlist(ix);

fpfn = fullfile(fp_root, fn{1});
eval(['load ' fpfn]);

v_E(abs(v_E)>5e3) = nan;
v_N(abs(v_E)>5e3) = nan;
v_h(abs(v_E)>5e3) = nan;

NWIND = [];
scale = 1e3;
NWIND.EF_E.val = EF_E'*scale;
NWIND.EF_E.err = sigma_EF_E'*scale;
NWIND.EF_N.val = EF_N'*scale;
NWIND.EF_N.err = sigma_EF_N'*scale;

NWIND.uEh_N.val = u_N';
NWIND.uEh_N.err = sigma_uN';
NWIND.uEh_E.val = u_E';
NWIND.uEh_E.err = sigma_uE';
NWIND.uEh_h.val = u_h';
NWIND.uEh_h.err = sigma_uh';

EF_tl = (datenum(EF_time1') + datenum(EF_time2'))/2;
u_tl = (datenum(u_time1') + datenum(u_time2'))/2;

NWIND.tl_EF = EF_tl';
NWIND.tl_u = u_tl';
NWIND.alt_EF = hgF;
NWIND.alt_u = hgE';

if isempty(search_variable('vFz_N_NW'))
    return
end
NWIND.vFz_N.val = Fvv(1, :);
NWIND.vFz_N.err = sigma_Fvv(1, :);
NWIND.vFz_E.val = Fvv(2, :);
NWIND.vFz_E.err = sigma_Fvv(2, :);
NWIND.vFz_z.val = Fvv(3, :);
NWIND.vFz_z.err = sigma_Fvv(3, :);
NWIND.vEh_N.val = v_N';
NWIND.vEh_N.err = sigma_vN';
NWIND.vEh_E.val = v_E';
NWIND.vEh_E.err = sigma_vE';
NWIND.vEh_h.val = v_h';
NWIND.vEh_h.err = sigma_vh';

Binc = Binc_F;
[NWIND.vFh_E, NWIND.vFh_N, NWIND.vFh_h] =   ...
    transform_z2h_witherror(NWIND.vFz_N, NWIND.vFz_E, NWIND.vFz_z, Binc);

[NWIND.vEz_N, NWIND.vEz_E, NWIND.vEz_z] =   ...
    transform_h2z_witherror(NWIND.vEh_E, NWIND.vEh_N, NWIND.vEh_h, Binc);

[NWIND.uEz_N, NWIND.uEz_E, NWIND.uEz_z] =   ...
    transform_h2z_witherror(NWIND.uEh_E, NWIND.uEh_N, NWIND.uEh_h, Binc);

end
