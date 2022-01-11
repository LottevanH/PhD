function [NaN]=assign_ElSpecGLOW()
global datasetinfo dataset

dateran=datasetinfo.dateran;

pnlist={'EF_E_NW', 'NWIND.EF_E.val', 'NWIND.EF_E.err', 'NWIND.tl_EF', 'NWIND.alt_EF';
        'EF_N_NW', 'NWIND.EF_N.val', 'NWIND.EF_N.err', 'NWIND.tl_EF', 'NWIND.alt_EF';
        'uEh_E_NW', 'NWIND.uEh_E.val', 'NWIND.uEh_E.err', 'NWIND.tl_u', 'NWIND.alt_u';
        'uEh_N_NW', 'NWIND.uEh_N.val', 'NWIND.uEh_N.err', 'NWIND.tl_u', 'NWIND.alt_u';
        'uEh_h_NW', 'NWIND.uEh_h.val', 'NWIND.uEh_h.err', 'NWIND.tl_u', 'NWIND.alt_u';
};
if ~isempty(search_variable('vFz_N_NW'))
    pnlist = [pnlist;   ...
        {
        'uEz_E_NW', 'NWIND.uEz_E.val', 'NWIND.uEz_E.err', 'NWIND.tl_u', 'NWIND.alt_u';
        'uEz_N_NW', 'NWIND.uEz_N.val', 'NWIND.uEz_N.err', 'NWIND.tl_u', 'NWIND.alt_u';
        'uEz_z_NW', 'NWIND.uEz_h.val', 'NWIND.uEz_h.err', 'NWIND.tl_u', 'NWIND.alt_u';
        'vEh_E_NW', 'NWIND.vEh_E.val', 'NWIND.vEh_E.err', 'NWIND.tl_u', 'NWIND.alt_u';
        'vEh_N_NW', 'NWIND.vEh_N.val', 'NWIND.vEh_N.err', 'NWIND.tl_u', 'NWIND.alt_u';
        'vEh_h_NW', 'NWIND.vEh_h.val', 'NWIND.vEh_h.err', 'NWIND.tl_u', 'NWIND.alt_u';
        'vEz_E_NW', 'NWIND.vEz_E.val', 'NWIND.vEz_E.err', 'NWIND.tl_u', 'NWIND.alt_u';
        'vEz_N_NW', 'NWIND.vEz_N.val', 'NWIND.vEz_N.err', 'NWIND.tl_u', 'NWIND.alt_u';
        'vEz_z_NW', 'NWIND.vEz_h.val', 'NWIND.vEz_h.err', 'NWIND.tl_u', 'NWIND.alt_u';
        'vFh_E_NW', 'NWIND.vFh_E.val', 'NWIND.vFh_E.err', 'NWIND.tl_EF', 'NWIND.alt_EF';
        'vFh_N_NW', 'NWIND.vFh_N.val', 'NWIND.vFh_N.err', 'NWIND.tl_EF', 'NWIND.alt_EF';
        'vFh_h_NW', 'NWIND.vFh_h.val', 'NWIND.vFh_h.err', 'NWIND.tl_EF', 'NWIND.alt_EF';
        'vFz_E_NW', 'NWIND.vFz_E.val', 'NWIND.vFz_E.err', 'NWIND.tl_EF', 'NWIND.alt_EF';
        'vFz_N_NW', 'NWIND.vFz_N.val', 'NWIND.vFz_N.err', 'NWIND.tl_EF', 'NWIND.alt_EF';
        'vFz_z_NW', 'NWIND.vFz_z.val', 'NWIND.vFz_z.err', 'NWIND.tl_EF', 'NWIND.alt_EF';
        }
        ];
end
pararec=datasetinfo.filereader(datasetinfo.filereaderrec).pararec;
npara=length(pararec);
para=cell(npara,1);
for i=1:npara
  para{i}.val=[];
  para{i}.tl=[];
  para{i}.err=[];
  para{i}.alt = [];
end

%% load data
for dn=floor(dateran(1)-0.5):floor(dateran(2)-0.5)
  if dn==floor(dateran(2)) && dateran(2)==floor(dateran(2)) 
    continue
  end
  
  [NWIND]=loaddata_EISCAT_NWIND(dn);
  
  if isempty(NWIND)
    continue
  end
  for i=1:length(pararec)
    ix=regexp(pnlist(:,1),datasetinfo.paralist(pararec(i)));
    ix=~cellfun('isempty',ix);
    val=eval(pnlist{ix,2});
    tl=eval(pnlist{ix,4});
    err=eval(pnlist{ix,3});
    if isempty(pnlist{ix, 5})
      alt = [];
    else
      alt=eval(pnlist{ix,5});
    end
    
    para{i}.val=[para{i}.val val];
    para{i}.tl=[para{i}.tl tl];
    para{i}.err=[para{i}.err err];
    para{i}.alt=[para{i}.alt alt];
  end
end
%% restrict time range
for i=1:length(pararec)
  if isempty(para{i}.tl); continue; end
  [tl,indtst,indted]=ts_confinetimeline(para{i}.tl,dateran(1),dateran(2));
  para{i}.val=para{i}.val(:,indtst:indted);
  para{i}.err=para{i}.err(:,indtst:indted);
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