
%%%%%% combine GUISDAP analyzed data into a single file
%%%%%% Input: pathlist - the path list of the GUISDAP data files
%%%%%% Output: Combined files in folder "dat_cmb"
%%%%%% Modified by Lei Cai on 31.03.2014

res_root_path=[pwd '/dat_cmb/' ];
itgtime_lim=5; % threshold for the integration time in second; 

fl=dir([fp '*.mat']);
nfiles=length(fl);
  
% load all files and read parameters in matrices

%%% Comment by L. Cai: For CP2 experiment the number of range gates in
%%% different beam direction are different. To unify the number of range
%%% gates, we set the smallest number of gates in the measurements as the
%%% final number of range gates.
num_gates=0;
npf=0;
outliers=[];
for i=1:nfiles
  data=load([fp fl(i).name]);
  num_gates=max([num_gates size(data.r_range,1)]);
  if size(data.r_range,1)>50
      i;
  end
  if round((datenum(data.r_time(2,:))-datenum(data.r_time(1,:)))*86400)<itgtime_lim
    outliers=[outliers i];
    continue
  else
    npf=npf+1;
    ind_pf(npf)=i;
  end
%   if i==1
%     % Find number of �real��PP gates
%     ind=find(diff(r_pprange)<0); 			 % r_pp incl. real PP plus MP and LP zero lags in succession
%     if ind ~=0;
%       ng_pp=ind(1); 							 % # of PP gates
%     else
%       ng_pp=length(r_pprange)
%     end
%   end
end
if ~isempty(outliers)
  disp(sprintf('Note: data within a short integration time (<%d s) has been taken out!',itgtime_lim)) 
  disp(outliers)
end


%Disp(['Warning: the upper boundary set to the range gate ' num2str(num_gates) '!'])
gran=1:num_gates;
ran=nan(num_gates,npf);height=nan(num_gates,npf);n_e=nan(num_gates,npf);
ti=nan(num_gates,npf);te=nan(num_gates,npf);coll=nan(num_gates,npf);
vel=nan(num_gates,npf);ne_err=nan(num_gates,npf);ti_err=nan(num_gates,npf);
te_err=nan(num_gates,npf);coll_err=nan(num_gates,npf);vel_err=nan(num_gates,npf);
stat=nan(num_gates,npf);resid=nan(num_gates,npf);
comp=nan(num_gates,npf);

for i=1:npf
  ii=ind_pf(i);
 
  data=load([fp fl(ii).name]);

  ng=size(data.r_range,1);
%   if size(data.r_range,1)<num_gates
%     return
%   end
  az(i)=data.r_az;
  el(i)=data.r_el;
  tx_power(i)=data.r_Pt;
%  t_sys1(i)=r_Tsys(1);
%  t_sys2(i)=r_Tsys(2);
  t1(:,i)=data.r_time(1,:)';
  t2(:,i)=data.r_time(2,:)';
%  ran_pp(:,i)=r_pprange(1:ng_pp,1);
%  ne_pp(:,i)=r_pp(1:ng_pp,1);
  ran(1:ng,i)=data.r_range(:,1);
  height(1:ng,i)=data.r_h(:,1);
  n_e(1:ng,i)=data.r_param(:,1);			% Ne in m^-3
  ti(1:ng,i)=data.r_param(:,2);			% Ti in K
  te(1:ng,i)=data.r_param(:,3).*data.r_param(:,2);  % Te = (Te/Ti)*Ti in K
  coll(1:ng,i)=data.r_param(:,4);		% ion-neutral coll. freq.in Hz
  vel(1:ng,i)=data.r_param(:,5); 		% velocity in m/s
  ne_err(1:ng,i)=data.r_error(:,1);		% errors of parameters
  ti_err(1:ng,i)=data.r_error(:,2);
  te_err(1:ng,i)=data.r_error(:,3).*data.r_param(:,2)+data.r_error(:,2).*data.r_param(:,3); % the error of product!
  coll_err(1:ng,i)=data.r_error(:,4);
  vel_err(1:ng,i)=data.r_error(:,5);
  stat(1:ng,i)=data.r_status(:);			% 0 = fit OK, 1 = max nmb of iterations exceeded, 2 = data too noisy, no fit done
  resid(1:ng,i)=data.r_res(:,1);		% residual (good values close to 1)
  comp(1:ng,i)=data.r_dp(:);				% O+/Ne ratio  
end

r_m0=data.r_m0;
r_XMITloc=data.r_XMITloc;
r_RECloc=data.r_RECloc;
r_SCangle=data.r_SCangle;
name_expr=data.name_expr;
name_ant = data.name_ant;
r_ver=data.r_ver;
switch data.name_site
  case 'T'
    site=2;
    fn_site='TRO';
  case 'V'
    site=2;
    fn_site='TRO'
  case 'K'
    site=1;
    fn_site='KIR';
  case 'S'
    site=3;
    fn_site='SOD';
  case 'L'
    site=4;
    fn_site='ESR';
end            
  

% range, az, el => height lat lon is calculated for each point
% Calculation takes time, so a check is made not to calc[     ulate same values twice
% Let's make calculation for the first profile

[height(:,1),lat(:,1),lon(:,1)]=loc2gg_mat(site,ran(:,1),az(1),el(1));
%[height_pp(:,1),lat_pp(:,1),lon_pp(:,1)]=loc2gg_mat(site,ran_pp(:,1),az(1),el(1));

% Other profiles
for i=2:npf
	ind=0;
	for j=1:i-1	
		if (az(i)==az(j) && el(i)==el(j)); ind=j; end		
	end
	if ind~=0
		height(:,i)=height(:,ind);
		lat(:,i)=lat(:,ind);
		lon(:,i)=lon(:,ind);
%		height_pp(:,i)=height_pp(:,ind);
%		lat_pp(:,i)=lat_pp(:,ind);
%		lon_pp(:,i)=lon_pp(:,ind);
else
		[height(:,i),lat(:,i),lon(:,i)]=loc2gg_mat(site,ran(:,i),az(i),el(i));
%		[height_pp(:,i),lat_pp(:,i),lon_pp(:,i)]=loc2gg_mat(site,ran_pp(:,i),az(i),el(i));
	end
end

% Write result file

dates=[num2str(t1(1,floor(npf/2)),'%4.4d') num2str(t1(2,floor(npf/2)),'%2.2d') num2str(t1(3,floor(npf/2)),'%2.2d')];
starttime=t1(4,1)*10000+t1(5,1)*100+t1(6,1);
stoptime=t1(4,npf)*10000+t1(5,npf)*100+t1(6,npf);
name_site=data.name_site;
disp([fn_site ', ' dates])

%res_path=[res_root_path fn_site dates(1:6) '/'];
%if ~isdir(res_path)
%  mkdir(res_root_path,[fn_site dates(1:6)]);
%end

%filenameout=[fn_site '_' dates '_cmb' '.mat'];
%eval(['save ',res_path,filenameout,...
%' az el tx_power t1 t2 ran height lat lon n_e ti te coll vel ',...
%'ne_err ti_err te_err coll_err vel_err stat resid comp r_m0 ',...
%'name_site r_XMITloc r_RECloc r_SCangle name_expr r_ver dates starttime stoptime']);


