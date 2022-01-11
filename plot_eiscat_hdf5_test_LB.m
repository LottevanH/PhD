clear all
% program to use the function in guisdap load_param_hdf5.m to load up data
% from a hdf5 file which has been downloaded from madrigal

% Instructions:
% 1. Download the hdf5 data file from madrigal
% 2. start guisdap in matlab and then run the following command (e.g):
% [Time,par2D,par1D,rpar2D,err2D,errs2d]=load_param_hdf5('D:\EISCAT_data\analyzed\eiscat_website_hdf5\EISCAT_2019-11-20_folke_64@42mb.hdf5')
% Variable names: 
% par2D: [Ran,Alt,Ne,Te,Ti,Vi,Coll,Comp,Res]
% par1D: [Az,El,Pt,Tsys]
% rpar2D: [ppRan,ppAlt,ppNe]
% err2d_id  = {'var_Ne' 'var_Tr' 'var_Ti' 'var_Vi' 'var_Collf'};
%   Note 1: 'Tr' as input will generate 'Te'. In this case 'Ti' is needed as input as well.
%   Note 2: 'Te' will generate NaN...
% Input:
%   pars1d  = names of parameters wanted in par1D
%   pars2d  = names of parameters wanted in par2D
%   do_rpar = true/false (extract pp-data or not)
%   errs2d  = generated in function: those parameters in pars2d with corresponding variances in the
%     HDF5 file will be
% Default (load_param_hdf5(hdf5file)) sets
%   pars1d  = {'az' 'el' 'Pt' 'Tsys'};
%   par2d_id  = {'range' 'h' 'Ne' 'Tr' 'Ti' 'Vi' 'Collf' 'po+' 'res'};
%   do_rpars = true
%   errs2d_id  = {'var_Ne' 'var_Tr' 'var_Ti' 'var_Vi' 'var_Collf'}

% [Time,par2D,par1D,rpar2D,err2D,errs2d]=load_param_hdf5('C:\Svalbard\PhD\Data\ULF_waves_events_currently_unexamined\2017_12_18\EISCAT_2017-12-18_bella_60@vhf.hdf5');%'D:\EISCAT_data\analyzed\eiscat_website_hdf5\EISCAT_2019-11-20_folke_64@42mb.hdf5');
[Time,par2D,par1D,rpar2D,err2D,errs2d]=load_param_hdf5('C:\Users\lotte\OneDrive - Universitetssenteret på Svalbard AS (1)\Svalbard\PhD\Data\Paper_Lisa_2017\EISCAT_2007-12-27_ipy_60@42m.hdf5');
% [Time,par2D,par1D,rpar2D,err2D,errs2d]=load_param_hdf5('C:\Svalbard\PhD\Data\Test_data\EISCAT_2019-11-20_folke_64@42mb.hdf5');

[y,m,d,h,mn,s]=datevec(Time(1,:));
ut_time=h+mn/60.;
alt=par2D(:,:,2);
ne=par2D(:,:,3);
te=par2D(:,:,4);%.*par2D(:,:,5);
ti=par2D(:,:,5);
alt_min = min(alt,[],'all');
Time_datetime = datetime(Time,'ConvertFrom','datenum'); %does something weird for 
%% Introduce NaN values when radar is off
% for i = 1:size(Time,2)-1
%     Time_diff(:,i) = datetime(Time(:,i+1),'ConvertFrom','datenum') - datetime(Time(:,i),'ConvertFrom','datenum');
% end
% if isempty(find(Time_diff(1,:) > seconds(65))) == 0
%     idx_time_gap = find(Time_diff(1,:) > seconds(65));
%     time_gap = Time_diff(1,idx_time_gap);
%     time_missing = floor(time_gap/seconds(64));
%     ne_new = NaN(size(ne,1),size(ne,2)+sum(time_missing));
%     [alt_new,te_new,ti_new] = deal(ne_new);
%     Time_new = NaN(size(Time,1),size(Time,2)+sum(time_missing));
%     ne_new(:,1:idx_time_gap(1)) = ne(:,1:idx_time_gap(1));    
%     ne_new(:,idx_time_gap(end)+sum(time_missing)+1:end) = ne(:,idx_time_gap(end)+1:end);
%     alt_new(:,1:idx_time_gap(1)) = alt(:,1:idx_time_gap(1));    
%     alt_new(:,idx_time_gap(end)+sum(time_missing)+1:end) = alt(:,idx_time_gap(end)+1:end);
%     Time_new(:,1:idx_time_gap(1)) = Time(:,1:idx_time_gap(1));    
%     Time_new(:,idx_time_gap(end)+sum(time_missing)+1:end) = Time(:,idx_time_gap(end)+1:end);
% %     for j = idx_time_gap(1)+1:idx_time_gap(1)+1+time_missing(1)
% %         Time_new(:,j) = datenum(datetime(Time_new(:,idx_time_gap(1)),'ConvertFrom','datenum')+seconds(64)*j);
% %     end
%     te_new(:,1:idx_time_gap(1)) = te(:,1:idx_time_gap(1));    
%     te_new(:,idx_time_gap(end)+sum(time_missing)+1:end) = te(:,idx_time_gap(end)+1:end);
%     ti_new(:,1:idx_time_gap(1)) = ti(:,1:idx_time_gap(1));    
%     ti_new(:,idx_time_gap(end)+sum(time_missing)+1:end) = ti(:,idx_time_gap(end)+1:end);    
%     if size(idx_time_gap,2) > 1
%         for i = 2:size(idx_time_gap,2)
%             ne_new(:,idx_time_gap(i-1)+time_missing(i-1):idx_time_gap(i)+time_missing(i-1)) = ne(:,idx_time_gap(i-1):idx_time_gap(i));
%             alt_new(:,idx_time_gap(i-1)+time_missing(i-1):idx_time_gap(i)+time_missing(i-1)) = alt(:,idx_time_gap(i-1):idx_time_gap(i));
%             te_new(:,idx_time_gap(i-1)+time_missing(i-1):idx_time_gap(i)+time_missing(i-1)) = te(:,idx_time_gap(i-1):idx_time_gap(i));
%             ti_new(:,idx_time_gap(i-1)+time_missing(i-1):idx_time_gap(i)+time_missing(i-1)) = ti(:,idx_time_gap(i-1):idx_time_gap(i));
%             Time_new(:,idx_time_gap(i-1)+time_missing(i-1):idx_time_gap(i)+time_missing(i-1)) = Time(:,idx_time_gap(i-1):idx_time_gap(i));
% %                 for j = idx_time_gap(i)+1:idx_time_gap(i)+1+time_missing(i)
% %                     Time_new(:,j) = datenum(datetime(Time_new(:,idx_time_gap(i)),'ConvertFrom','datenum')+seconds(64)*j);
% %                 end
%         end
%     else ;
%     end
% else ;
% end
% 
% [Time_new_1,T_Time_new_1] = fillmissing(Time_new,'spline')
% %% smooth the Time_new array
% for i = idx_time_gap(1)+1:idx_time_gap(1)+1+time_missing(1)
% Time_new(i) = datenum(datetime(Time_new(idx_time_gap(1)),'ConvertFrom','datenum')+seconds(64)*i);
% end

% Time_diff1 = Time_datetime(2,:) - Time_datetime(1,:);
%% Plot
n1 = Time(1,1); %datenum([2019 11 20 16 00 00]);
n2 = Time(2,end); %datenum([2019 11 20 20 00 00]);
xLimits=[n1 n2];
yLimits=[alt_min 300];
f1=figure();
colormap jet;

subplot(3,1,1);
% pcolor(datenum([y(:),m(:),d(:),h(:),mn(:),s(:)]),alt,log10(ne)),shading flat;
% pcolor(Time(1,:)',alt,log10(ne)),shading flat;
% imagesc(Time',alt,log10(ne))
surface(Time(1,:)',alt,log10(ne),'EdgeColor','none');%,shading flat;

% give the limits of the colourbar
caxis([10 11.4]);
% % Add colorbar to the right of the plot
colorbar;
set(get(colorbar,'label'),'FontSize', 8,'string','Electron density (m^{-3})')
xlim(xLimits);
% xlim([datenum(2017,12,18,02,00,00) datenum(2017,12,18,12,00,00)]);
ylim(yLimits) ;
datetick('x',13,'keeplimits')
% % Add labels to axes and colorbar
ylabel('Altitude [km]','FontSize', 8,'FontName','Arial')
 xlabel(['Time [UT] ',num2str(y(1)),'-',num2str(m(1)),'-',num2str(d(1))],'FontSize', 8,'FontName','Arial')
 
subplot(3,1,2)
pcolor(Time(1,:)',alt,te),shading flat;
% Limits to colorbar
caxis([0 3000]);
% % Add colorbar to the right of the plot
colorbar;
set(get(colorbar,'label'),'FontSize', 8,'string','Electron temperature (K)')
xlim(xLimits);
% xlim([datenum(2017,12,18,02,00,00) datenum(2017,12,18,12,00,00)]);
ylim([min(alt,[],'all') 300])% ylim(yLimits) ;
datetick('x',13,'keeplimits')
% % Add labels to axes and colorbar
ylabel('Altitude [km]','FontSize', 8,'FontName','Arial')
 xlabel(['Time [UT] ',num2str(y(1))],'FontSize', 8,'FontName','Arial')

subplot(3,1,3)
pcolor(Time(1,:)',alt,ti),shading flat;
% Limits to colorbar
caxis([0 3000]);
% % Add colorbar to the right of the plot
colorbar;
set(get(colorbar,'label'),'FontSize', 8,'string','Ion temperature (K)')
xlim(xLimits);
% xlim([datenum(2017,12,18,02,00,00) datenum(2017,12,18,12,00,00)]);
ylim([100 300])%ylim(yLimits) ;
datetick('x',13,'keeplimits')
% % Add labels to axes and colorbar
ylabel('Altitude [km]','FontSize', 8,'FontName','Arial')
 xlabel(['Time [UT] ',num2str(y(1))],'FontSize', 8,'FontName','Arial')

% %%
% n1 = Time_new(1,1); %datenum([2019 11 20 16 00 00]);
% n2 = Time_new(2,end); %datenum([2019 11 20 20 00 00]);
% xLimits=[n1 n2];
% yLimits=[alt_min 300];
% f1=figure();
% colormap jet;
% 
% subplot(3,1,1);
% % pcolor(datenum([y(:),m(:),d(:),h(:),mn(:),s(:)]),alt,log10(ne)),shading flat;
% % pcolor(Time(1,:)',alt,log10(ne)),shading flat;
% % imagesc(Time',alt,log10(ne))
% surface(Time_new(1,:)',alt_new,log10(ne_new),'EdgeColor','none');%,shading flat;
% 
% % give the limits of the colourbar
% caxis([10 11.4]);
% % % Add colorbar to the right of the plot
% colorbar;
% set(get(colorbar,'label'),'FontSize', 8,'string','Electron density (m^{-3})')
% xlim(xLimits);
% % xlim([datenum(2017,12,18,02,00,00) datenum(2017,12,18,12,00,00)]);
% ylim(yLimits) ;
% datetick('x',13,'keeplimits')
% % % Add labels to axes and colorbar
% ylabel('Altitude [km]','FontSize', 8,'FontName','Arial')
%  xlabel(['Time [UT] ',num2str(y(1)),'-',num2str(m(1)),'-',num2str(d(1))],'FontSize', 8,'FontName','Arial')
%  
% subplot(3,1,2)
% pcolor(Time_new(1,:)',alt_new,te_new),shading flat;
% % Limits to colorbar
% caxis([0 3000]);
% % % Add colorbar to the right of the plot
% colorbar;
% set(get(colorbar,'label'),'FontSize', 8,'string','Electron temperature (K)')
% xlim(xLimits);
% % xlim([datenum(2017,12,18,02,00,00) datenum(2017,12,18,12,00,00)]);
% ylim([min(alt_new,[],'all') 300])% ylim(yLimits) ;
% datetick('x',13,'keeplimits')
% % % Add labels to axes and colorbar
% ylabel('Altitude [km]','FontSize', 8,'FontName','Arial')
%  xlabel(['Time [UT] ',num2str(y(1))],'FontSize', 8,'FontName','Arial')
% 
% subplot(3,1,3)
% pcolor(Time_new(1,:)',alt_new,ti_new),shading flat;
% % Limits to colorbar
% caxis([0 3000]);
% % % Add colorbar to the right of the plot
% colorbar;
% set(get(colorbar,'label'),'FontSize', 8,'string','Ion temperature (K)')
% xlim(xLimits);
% % xlim([datenum(2017,12,18,02,00,00) datenum(2017,12,18,12,00,00)]);
% ylim([100 300])%ylim(yLimits) ;
% datetick('x',13,'keeplimits')
% % % Add labels to axes and colorbar
% ylabel('Altitude [km]','FontSize', 8,'FontName','Arial')
%  xlabel(['Time [UT] ',num2str(y(1))],'FontSize', 8,'FontName','Arial')
 %% test FFT
% %  %find altitude between lower limit (low_lim) and higher limit (high_lim) (in km) in which you want to find ULF waves:
% low_lim = 100;
% high_lim = 300;
% te_lim = te;
% te_lim(alt<low_lim | alt>high_lim) = 0;
% 
% alt_lim = 150:5:300;
% te_lim1 = te;
% 
% %divide per hour?
% delta_T = round(3600/64);
% delta_T_array = 1:delta_T:size(te,2);
% 
% T = 64; %period is 64 seconds
% Fs = 1/T; 
% for i = 1:size(alt_lim,2)-1
%     for j = 1:size(delta_T_array,2)-1
%         test{i,j} = te_lim1(alt>alt_lim(i) & alt<alt_lim(i+1) & Time_datetime(1,:) > Time_datetime(1,delta_T_array(j)) & Time_datetime(1,:) < Time_datetime(1,delta_T_array(j+1)));
%     if size(test{i,j},1) > 45%isempty(test{i}) == 0
%         disp(num2str(j))
% %         FFT1{i,j} = fft(test{i,j});
% %         L(i,j) = size(test{i,j},1);
% %         P2{i,j} = abs(FFT1{i,j}/L(i));
% %         P1{i,j} = P2{i}(1:L(i,j)/2+1);
% %         P1_a = P1{i,j};
% %         P1_a(2:end-1) = 2*P1_a(2:end-1);
% %         f{i,j} = Fs*(0:(L(i,j)/2))/L(i,j);
% %         figure()
% %         plot(f{i,j},P1_a)
% %         xlim([0 6e-3])
% %         clear P1_a
%     else ;
%     end
%     end
% end
% 
% % for i = 1:size(alt_lim,2)-1
% %         test{i} = te_lim1(alt>alt_lim(i) & alt<alt_lim(i+1));
% %     if size(test{i},1) > size(te,2)/2%isempty(test{i}) == 0
% %         FFT1{i} = fft(test{i});
% %         L(i) = size(test{i},1);
% %         P2{i} = abs(FFT1{i}/L(i));
% %         P1{i} = P2{i}(1:L(i)/2+1);
% %         P1_a = P1{i};
% %         P1_a(2:end-1) = 2*P1_a(2:end-1);
% %         f{i} = Fs*(0:(L(i)/2))/L(i);
% %         figure()
% %         plot(f{i},P1_a)
% %         xlim([0 6e-3])
% %         clear P1_a
% %     else ;
% %     end
% % end


 %% test FFT 2007.12.27
%  %find altitude between lower limit (low_lim) and higher limit (high_lim) (in km) in which you want to find ULF waves:



low_lim = 100;
high_lim = 300;
te_lim = te;
te_lim(alt<low_lim | alt>high_lim) = 0;

alt_lim = 150:5:300;
te_lim1 = te;

%divide per hour?
delta_T = round(3600/64);
delta_T_array = 1:delta_T:size(te,2);

T = 64; %period is 64 seconds
Fs = 1/T; 
% for i = 1:size(alt_lim,2)-1
%     for j = 1:size(delta_T_array,2)-1
%         test{i,j} = te_lim1(alt>alt_lim(i) & alt<alt_lim(i+1) & Time_datetime(1,:) > Time_datetime(1,delta_T_array(j)) & Time_datetime(1,:) < Time_datetime(1,delta_T_array(j+1)));
%     if size(test{i,j},1) > 45%isempty(test{i}) == 0
%         disp(num2str(j))
% %         FFT1{i,j} = fft(test{i,j});
% %         L(i,j) = size(test{i,j},1);
% %         P2{i,j} = abs(FFT1{i,j}/L(i));
% %         P1{i,j} = P2{i}(1:L(i,j)/2+1);
% %         P1_a = P1{i,j};
% %         P1_a(2:end-1) = 2*P1_a(2:end-1);
% %         f{i,j} = Fs*(0:(L(i,j)/2))/L(i,j);
% %         figure()
% %         plot(f{i,j},P1_a)
% %         xlim([0 6e-3])
% %         clear P1_a
%     else ;
%     end
%     end
% end
% 
% % for i = 1:size(alt_lim,2)-1
% %         test{i} = te_lim1(alt>alt_lim(i) & alt<alt_lim(i+1));
% %     if size(test{i},1) > size(te,2)/2%isempty(test{i}) == 0
% %         FFT1{i} = fft(test{i});
% %         L(i) = size(test{i},1);
% %         P2{i} = abs(FFT1{i}/L(i));
% %         P1{i} = P2{i}(1:L(i)/2+1);
% %         P1_a = P1{i};
% %         P1_a(2:end-1) = 2*P1_a(2:end-1);
% %         f{i} = Fs*(0:(L(i)/2))/L(i);
% %         figure()
% %         plot(f{i},P1_a)
% %         xlim([0 6e-3])
% %         clear P1_a
% %     else ;
% %     end
% % end


 