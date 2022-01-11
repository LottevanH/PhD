function ISR = loaddata_EISCAT_cdt(dn)
  
%% find filenames
global datasetinfo

try
    dir_data = datasetinfo.EISCAT.dir_data;
catch
    dir_data = '';
end

try
    sitename = datasetinfo.EISCAT.sitename;
catch
    sitename = 'TRO';
end

try
    antenna = datasetinfo.EISCAT.sitename;
catch
    antenna = 'UHF';
end

try
    pulsecode = datasetinfo.EISCAT.pulsecode;
catch
    pulsecode = 'beata';
end

try
    beammode = datasetinfo.EISCAT.beammode;
catch
    beammode = '';
end

try
    hostcode = datasetinfo.EISCAT.hostcode;
catch
    hostcode = '';
end

try
    filemode = datasetinfo.EISCAT.filemode;
catch
    filemode = 'EISCAT_cdt';    % autoselect, EISCAT_cdt, manual
end



if strcmp(filemode,  'EISCAT_cdt')
    if isempty(dir_data)
        pwd1=pwd;
        cd ..
        fp_root=[pwd '/EISCAT_cdt/dat_cmb/'];
        cd(pwd1)
    end

    dstr = datestr(dn, 'yyyymmdd');

    dstr1 = datestr(dn, 'yyyymm');

    flist = dir([fp_root  '/' sitename '/'  dstr1 '/*.mat']);
    fnlist = {flist.name};
    ix = regexp(fnlist, [sitename '_' dstr '_cmb']);
    ix = ~cellfun('isempty', ix);
    fn = fnlist(ix);

    fpfn = fullfile(fp_root, [sitename dstr1], fn{1});
    eval(['load ' fpfn]);
elseif strcmp(filemode, 'autoselect')
    if isempty(dir_data)
        fp_root=['/kaappi/EISCAT/RESULTS/'];
    end
    
    yystr = datestr(dn, 'yyyy');
    fs = dir([fp_root '/res_' yystr '_OY']);
    ind_subdir = [fs.isdir];
    fps = {fs(ind_subdir).name};
    ix = regexp(fps, [datestr(dn, 'yyyy-mm-dd') '_' pulsecode]);
    ix = ~cellfun('isempty', ix);
    fps = fps(ix);
    ix = regexp(fps, lower(antenna));
    ix = ~cellfun('isempty', ix);
    fp = fps(ix);
    fp = fp{1};
    combine_mat14;
    datasetinfo.EISCAT.antenna = name_ant;
    datasetinfo.EISCAT.pulsecode = name_expr;
    datasetinfo.EISCAT.sitename = fn_site;
    
    if isfile([fp '/gfd_setup.m'])
        datasetinfo.EISCAT.key = get_EISCAT_key(fp);
    else
        datasetinfo.EISCAT.key = [fn_site '_' name_ant '_' name_expr];
    end
elseif strcmp(filemode, 'manual')
    filterSpec=[pwd '/*.mat'];
    filterSpec=['/kaappi/EISCAT/RESULTS/*.mat'];
    [filenamein,pathname]=uigetfile(filterSpec,         ... 
        [datestr(dn, 'yyyy-mm-dd')  ...
         ',Select any *.mat file in the GUISDAP result folder:']);
    if filenamein==0; return, end
    nfd=1;
    fp=pathname;
    combine_mat14;
    datasetinfo.EISCAT.antenna = name_ant;
    datasetinfo.EISCAT.pulsecode = name_expr;
    datasetinfo.EISCAT.sitename = fn_site;
    if isfile([fp '/gfd_setup.m'])
        datasetinfo.EISCAT.key = get_EISCAT_key(fp);
    else
        datasetinfo.EISCAT.key = [fn_site '_' name_ant '_' name_expr];
    end
end
    

%ind=find(abs(az-186.2)<1 & abs(el-77.5)<1);
  ind=1:size(n_e,2);
  ne=n_e(:,ind);
  ti=ti(:,ind);
  te=te(:,ind);
  vel=vel(:,ind);
  height=height(:,ind);
  stat=stat(:,ind);
  resid=resid(:,ind);
  ti_err=ti_err(:,ind);
  
  ind2=find(stat>1);
  ne(ind2)=nan;
  ti(ind2)=nan;
  te(ind2)=nan;
  vel(ind2) = nan;
  
%   ind3=find(resid>10);
%   ne(ind3)=nan;
%   ti(ind3)=nan;
%   te(ind3)=nan;
%   vel(ind3) = nan;
  
%   ind4=find(ti_err./ti>1);
%   ne(ind4)=nan;
%   ti(ind4)=nan;
%   te(ind4)=nan;
%   vel(ind4) = nan;
%   
%   % T_i < 50 K? velocity is not independent from other fittings, bad ion
%   % velocity may be related to wired temperature
%     ind4=find(ti>2500);
%   ne(ind4)=nan;
%   ti(ind4)=nan;
%   te(ind4)=nan;
%   vel(ind4)=nan;
  
  
  t1=datenum(t1');
  t2=datenum(t2');
  
  
  t1=t1(ind);
  t2=t2(ind);
  midtime=(t1+t2)/2;
  sectime=(midtime-floor(midtime))*86400;
  
  
  ISR.ne.val = ne;
  ISR.ne.err = [];
  ISR.Ti.val = ti;
  ISR.Ti.err = [];
  ISR.Te.val = te;
  ISR.Te.err = [];
  ISR.vi.val = vel;
  ISR.vi.err = [];
  ISR.az = az(1, ind);
  ISR.el = el(1, ind);
  ISR.tl = midtime';
  ISR.alt = height;
  
end

function txt=get_EISCAT_key(fp)
    run([fp '/gfd_setup.m'])
    s = regexp(data_path, '/', 'split');
    
    txt = s{end};
end
    