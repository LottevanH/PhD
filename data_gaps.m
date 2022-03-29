%% Test effect interpolation
%take a time interval without data gaps --> create data gaps of 5, 10, 15,
%20, 25 and 30 minutes and see the effect on the spectrogram

%data without time gap: 
% - 2007-12-27 00:00:00 - 07:10:00 
% - 2017-12-18 04:03:00 - 10:00:00 (ULF waves between 02-07 UT) --> try
% with this one


%% start script
clear all
close all

addpath('C:\Github\PhD\functions');
% filename = 'C:\Users\lotte\OneDrive - Universitetssenteret på Svalbard AS (1)\Svalbard\PhD\Data\Paper_Lisa_2017\EISCAT_2007-12-27_ipy_60@42m.hdf5';
filename = 'C:\Users\lotte\OneDrive - Universitetssenteret på Svalbard AS (1)\Svalbard\PhD\Data\ULF_waves_events_currently_unexamined\2017_12_18\EISCAT_2017-12-18_bella_60@vhf.hdf5';

%% extract variables and get metadeta info
[Time,par2D,par1D,rpar2D,err2D,errs2d]=load_param_hdf5(filename);
metadata_par2d = deblank(h5read(filename,'/metadata/par2d')); %to see which variables are available. In the load_param_hdf5 script the parameters saved for the par2d parameter are given.
metadata_par1d = deblank(h5read(filename,'/metadata/par1d'));
metadata_par0d = deblank(h5read(filename,'/metadata/par0d'));
par0D = h5read(filename,'/data/par0d');
%%
[y,m,d,h,mn,s]=datevec(Time(1,:));
ut_time=h+mn/60.;
alt=par2D(:,:,2);
ne=par2D(:,:,3);
ne(find(ne<0)) = NaN; %in case of negative electron densities (physically impossible..)
te=par2D(:,:,4);%.*par2D(:,:,5); %according to the load_param_hdf5 script (line 42), the 4th column of par2D gives the Te parameter instead of Tr (Te/Ti) 
ti=par2D(:,:,5);
Time_datetime = datetime(Time,'ConvertFrom','datenum'); 

%% Time period of interest
n11 = datenum(datetime(2017,12,18,02,00,00)); %in case of limited time
n22 = datenum(datetime(2017,12,18,07,00,00)); %in case of limited time
% n11 = Time(1,1); %does not work if there is a time gap
% n22 = Time(1,end);

idx_time_start = find(Time(1,:) > n11,1);
idx_time_end = find(Time(1,:) < n22,1,'last');
 
alt_poi = alt(:,idx_time_start:idx_time_end); %poi stands for: period of interest
ne_poi = ne(:,idx_time_start:idx_time_end);
te_poi = te(:,idx_time_start:idx_time_end);
ti_poi = ti(:,idx_time_start:idx_time_end);
Time_poi = Time(:,idx_time_start:idx_time_end);
Time_datetime_poi = Time_datetime(:,idx_time_start:idx_time_end);

%% Add NaNs in case of time gaps
for i = 1:size(Time_datetime_poi,2)-1
    delta_T(i) = Time_datetime_poi(1,i+1) - Time_datetime_poi(1,i);
end
idx_time_gap = find(delta_T > median(delta_T)+seconds(15)); %find time gaps if there are any
if length(idx_time_gap) == 1
    i = 1;
    number_of_NaN(i) = round((delta_T(idx_time_gap(i))+delta_T(idx_time_gap(i) + 1 ))/median(delta_T));
    Time_datetime_new(:,1:idx_time_gap(i)-1) = Time_datetime_poi(:,1:idx_time_gap(i)-1);
    Time_datetime_new(:,idx_time_gap(i):idx_time_gap(i)+number_of_NaN(i)-1) = NaT;
    Time_datetime_new(:,idx_time_gap(i)+number_of_NaN(i):length(Time_datetime_poi)+number_of_NaN-2) = Time_datetime_poi(:,idx_time_gap(i)+2:end);
    ne_new(:,1:idx_time_gap(i)-1) = ne_poi(:,1:idx_time_gap(i)-1);
    ne_new(:,idx_time_gap(i):idx_time_gap(i)+number_of_NaN(i)-1) = NaN;
    ne_new(:,idx_time_gap(i)+number_of_NaN(i):length(Time_datetime_poi)+number_of_NaN-2) = ne_poi(:,idx_time_gap(i)+2:end);
    te_new(:,1:idx_time_gap(i)-1) = te_poi(:,1:idx_time_gap(i)-1);
    te_new(:,idx_time_gap(i):idx_time_gap(i)+number_of_NaN(i)-1) = NaN;
    te_new(:,idx_time_gap(i)+number_of_NaN(i):length(Time_datetime_poi)+number_of_NaN-2) = te_poi(:,idx_time_gap(i)+2:end);    
    ti_new(:,1:idx_time_gap(i)-1) = ti_poi(:,1:idx_time_gap(i)-1);
    ti_new(:,idx_time_gap(i):idx_time_gap(i)+number_of_NaN(i)-1) = NaN;
    ti_new(:,idx_time_gap(i)+number_of_NaN(i):length(Time_datetime_poi)+number_of_NaN-2) = ti_poi(:,idx_time_gap(i)+2:end);
    Time_new(:,1:idx_time_gap(i)-1) = Time_poi(:,1:idx_time_gap(i)-1);
    Time_new(:,idx_time_gap(i):idx_time_gap(i)+number_of_NaN(i)-1) = NaN;
    Time_new(:,idx_time_gap(i)+number_of_NaN(i):length(Time_datetime_poi)+number_of_NaN-2) = Time_poi(:,idx_time_gap(i)+2:end);
    alt_new(:,1:idx_time_gap(i)-1) = alt_poi(:,1:idx_time_gap(i)-1);
    alt_new(:,idx_time_gap(i):idx_time_gap(i)+number_of_NaN(i)-1) = NaN;
    alt_new(:,idx_time_gap(i)+number_of_NaN(i):length(Time_datetime_poi)+number_of_NaN-2) = alt_poi(:,idx_time_gap(i)+2:end);

elseif length(idx_time_gap) > 1 %more than one time gap
    for i = 1:length(idx_time_gap)
        if i == 1
            number_of_NaN(i) = round((delta_T(idx_time_gap(i))+delta_T(idx_time_gap(i) + 1))/median(delta_T));
            Time_datetime_new(:,1:idx_time_gap(i)-1) = Time_datetime_poi(:,1:idx_time_gap(i)-1);
            Time_datetime_new(:,idx_time_gap(i):idx_time_gap(i)+number_of_NaN(i)-1) = NaT;
            ne_new(:,1:idx_time_gap(i)-1) = ne_poi(:,1:idx_time_gap(i)-1);
            ne_new(:,idx_time_gap(i):idx_time_gap(i)+number_of_NaN(i)-1) = NaN;
            te_new(:,1:idx_time_gap(i)-1) = te_poi(:,1:idx_time_gap(i)-1);
            te_new(:,idx_time_gap(i):idx_time_gap(i)+number_of_NaN(i)-1) = NaN;
            ti_new(:,1:idx_time_gap(i)-1) = ti_poi(:,1:idx_time_gap(i)-1);
            ti_new(:,idx_time_gap(i):idx_time_gap(i)+number_of_NaN(i)-1) = NaN;
            Time_new(:,1:idx_time_gap(i)-1) = Time_poi(:,1:idx_time_gap(i)-1);
            Time_new(:,idx_time_gap(i):idx_time_gap(i)+number_of_NaN(i)-1) = NaN;
            alt_new(:,1:idx_time_gap(i)-1) = alt_poi(:,1:idx_time_gap(i)-1);
            alt_new(:,idx_time_gap(i):idx_time_gap(i)+number_of_NaN(i)-1) = NaN;
            Time_datetime_new(:,idx_time_gap(i)+number_of_NaN(i):idx_time_gap(i+1)-1+number_of_NaN(i)-2) = Time_datetime_poi(:,idx_time_gap(i)+2:idx_time_gap(i+1)-1);
            ne_new(:,idx_time_gap(i)+number_of_NaN(i):idx_time_gap(i+1)-1+number_of_NaN(i)-2) = ne_poi(:,idx_time_gap(i)+2:idx_time_gap(i+1)-1);
            te_new(:,idx_time_gap(i)+number_of_NaN(i):idx_time_gap(i+1)-1+number_of_NaN(i)-2) = te_poi(:,idx_time_gap(i)+2:idx_time_gap(i+1)-1);
            alt_new(:,idx_time_gap(i)+number_of_NaN(i):idx_time_gap(i+1)-1+number_of_NaN(i)-2) = alt_poi(:,idx_time_gap(i)+2:idx_time_gap(i+1)-1);
            Time_new(:,idx_time_gap(i)+number_of_NaN(i):idx_time_gap(i+1)-1+number_of_NaN(i)-2) = Time_poi(:,idx_time_gap(i)+2:idx_time_gap(i+1)-1); 
            ti_new(:,idx_time_gap(i)+number_of_NaN(i):idx_time_gap(i+1)-1+number_of_NaN(i)-2) = ti_poi(:,idx_time_gap(i)+2:idx_time_gap(i+1)-1);
        elseif (i > 1) && (i < length(idx_time_gap)) 
            number_of_NaN(i) = round((delta_T(idx_time_gap(i))+delta_T(idx_time_gap(i) + 1 ))/median(delta_T));
            ne_new(:,idx_time_gap(i)+sum(number_of_NaN(1:i-1)-i):idx_time_gap(i)+sum(number_of_NaN(1:i))-i) = NaN;
            te_new(:,idx_time_gap(i)+sum(number_of_NaN(1:i-1)-i):idx_time_gap(i)+sum(number_of_NaN(1:i))-i) = NaN;
            Time_new(:,idx_time_gap(i)+sum(number_of_NaN(1:i-1)-i):idx_time_gap(i)+sum(number_of_NaN(1:i))-i) = NaN; 
            Time_datetime_new(:,idx_time_gap(i)+sum(number_of_NaN(1:i-1))-i:idx_time_gap(i)+sum(number_of_NaN(1:i))-i) = NaT; 
            alt_new(:,idx_time_gap(i)+sum(number_of_NaN(1:i-1)-i):idx_time_gap(i)+sum(number_of_NaN(1:i))-i) = NaN;
            ti_new(:,idx_time_gap(i)+sum(number_of_NaN(1:i-1)-i):idx_time_gap(i)+sum(number_of_NaN(1:i))-i) = NaN;
            Time_datetime_new(:,idx_time_gap(i)+sum(number_of_NaN(1:i))-i:idx_time_gap(i+1)-1+sum(number_of_NaN(1:i))-2-i) = Time_datetime_poi(:,idx_time_gap(i)+2:idx_time_gap(i+1)-1);
            Time_new(:,idx_time_gap(i)+sum(number_of_NaN(1:i))-i:idx_time_gap(i+1)-1+sum(number_of_NaN(1:i))-2-i) = Time_poi(:,idx_time_gap(i)+2:idx_time_gap(i+1)-1);
            ne_new(:,idx_time_gap(i)+sum(number_of_NaN(1:i))-i:idx_time_gap(i+1)-1+sum(number_of_NaN(1:i))-2-i) = ne_poi(:,idx_time_gap(i)+2:idx_time_gap(i+1)-1);
            te_new(:,idx_time_gap(i)+sum(number_of_NaN(1:i))-i:idx_time_gap(i+1)-1+sum(number_of_NaN(1:i))-2-i) = te_poi(:,idx_time_gap(i)+2:idx_time_gap(i+1)-1);
            alt_new(:,idx_time_gap(i)+sum(number_of_NaN(1:i))-i:idx_time_gap(i+1)-1+sum(number_of_NaN(1:i))-2-i) = alt_poi(:,idx_time_gap(i)+2:idx_time_gap(i+1)-1);
            ti_new(:,idx_time_gap(i)+sum(number_of_NaN(1:i))-i:idx_time_gap(i+1)-1+sum(number_of_NaN(1:i))-2-i) = ti_poi(:,idx_time_gap(i)+2:idx_time_gap(i+1)-1);
        elseif i == length(idx_time_gap)
            number_of_NaN(i) = round((delta_T(idx_time_gap(i))+delta_T(idx_time_gap(i) + 1 ))/median(delta_T));
            Time_datetime_new(:,idx_time_gap(i)+sum(number_of_NaN(1:i-1)-i):idx_time_gap(i)+sum(number_of_NaN(:))-1-i) = NaT;
            ne_new(:,idx_time_gap(i)+sum(number_of_NaN(1:i-1)-i):idx_time_gap(i)+sum(number_of_NaN(:))-1-i) = NaN;
            te_new(:,idx_time_gap(i)+sum(number_of_NaN(1:i-1)-i):idx_time_gap(i)+sum(number_of_NaN(:))-1-i) = NaN;
            ti_new(:,idx_time_gap(i)+sum(number_of_NaN(1:i-1)-i):idx_time_gap(i)+sum(number_of_NaN(:))-1-i) = NaN;
            Time_new(:,idx_time_gap(i)+sum(number_of_NaN(1:i-1)-i):idx_time_gap(i)+sum(number_of_NaN(:))-1-i) = NaN; 
            alt_new(:,idx_time_gap(i)+sum(number_of_NaN(1:i-1)-i):idx_time_gap(i)+sum(number_of_NaN(:))-1-i) = NaN;
            Time_datetime_new(:,idx_time_gap(i)+sum(number_of_NaN)-i:length(Time_datetime_poi)+sum(number_of_NaN)-2-i) = Time_datetime_poi(:,idx_time_gap(i)+2:end);
            ne_new(:,idx_time_gap(i)+sum(number_of_NaN)-i:length(Time_datetime_poi)+sum(number_of_NaN)-2-i) = ne_poi(:,idx_time_gap(i)+2:end);
            te_new(:,idx_time_gap(i)+sum(number_of_NaN)-i:length(Time_datetime_poi)+sum(number_of_NaN)-2-i) = te_poi(:,idx_time_gap(i)+2:end);
            ti_new(:,idx_time_gap(i)+sum(number_of_NaN)-i:length(Time_datetime_poi)+sum(number_of_NaN)-2-i) = te_poi(:,idx_time_gap(i)+2:end);
            alt_new(:,idx_time_gap(i)+sum(number_of_NaN)-i:length(Time_datetime_poi)+sum(number_of_NaN)-2-i) = alt_poi(:,idx_time_gap(i)+2:end);
            Time_new(:,idx_time_gap(i)+sum(number_of_NaN)-i:length(Time_datetime_poi)+sum(number_of_NaN)-2-i) = Time_poi(:,idx_time_gap(i)+2:end);
        else
            disp('Something is wrong')
        end            
    end
    
else %so no time gaps
    ne_new = ne_poi;
    te_new = te_poi;
    alt_new = alt_poi;
    Time_new = Time_poi; 
    ti_new = ti_poi;
    Time_datetime_new(:,:) = Time_datetime_poi(:,:);
end

addpath('C:\Github\PhD\functions\Inpaint_nans')
% alt_new = inpaint_nans(alt_new); %interpolate over the NaN values so that these data points are taken into account as well later on
%% Plot period of interest
xLimits=[n11 n22];
yLimits=[110 400];
f1=figure();
colormap jet;

subplot(3,1,1);
pcolor(Time_new(1,:)',alt_new,log10(ne_new)),shading flat;
% imagesc(Time',alt,log10(ne)) %does not work, because needs a vector for Y instead of a matrix (not sure how to do this with a changing altitude) --> this function should make it easier to have empty parts for missing data.
% surface(Time(1,:)',alt,log10(ne),'EdgeColor','none');%,shading flat; %different way to plot --> works

caxis([10 11.4]); % give the limits of the colourbar
colorbar; % Add colorbar to the right of the plot
set(get(colorbar,'label'),'FontSize', 8,'string','Electron density (m^{-3})')
xlim(xLimits); 
ylim(yLimits);
datetick('x',13,'keeplimits')
ylabel('Altitude [km]','FontSize', 8,'FontName','Arial') %Add labels to axes and colorbar
xlabel(['Time [UT] ',num2str(y(1)),'-',num2str(m(1)),'-',num2str(d(1))],'FontSize', 8,'FontName','Arial')
 
subplot(3,1,2)
pcolor(Time_new(1,:)',alt_new,te_new),shading flat;
caxis([0 3000]); % Limits to colorbar
colorbar; % Add colorbar to the right of the plot
set(get(colorbar,'label'),'FontSize', 8,'string','Electron temperature (K)') % Add labels to axes and colorbar
xlim(xLimits); 
ylim(yLimits); %ylim([min(alt,[],'all') 300])
datetick('x',13,'keeplimits')
ylabel('Altitude [km]','FontSize', 8,'FontName','Arial') % Add labels to axes and colorbar
xlabel(['Time [UT] ',num2str(y(1))],'FontSize', 8,'FontName','Arial') % Add labels to axes and colorbar

subplot(3,1,3)
pcolor(Time_new(1,:)',alt_new,ti_new),shading flat;
caxis([0 3000]); % Limits to colorbar
colorbar; % Add colorbar to the right of the plot
set(get(colorbar,'label'),'FontSize', 8,'string','Ion temperature (K)')
xlim(xLimits);
ylim(yLimits) ;
datetick('x',13,'keeplimits')
ylabel('Altitude [km]','FontSize', 8,'FontName','Arial') % Add labels to axes and colorbar
xlabel(['Time [UT] ',num2str(y(1))],'FontSize', 8,'FontName','Arial')

%% interpolate for already existing NaNs (straight line)
idx_NaN_begin1 = find(diff(isnan(Time_new(1,:)))==1)+1;
idx_NaN_end1 = find(diff(isnan(Time_new(1,:)))==-1);
idx_NaNs1 = idx_NaN_end1 - idx_NaN_begin1 + 1;
if isempty(idx_NaNs1) == 0
    for i = 1:size(te_new,1)
        X = te_new(i,:);
        X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)), 'linear'); %linear interpolation
        te_interp1(i,:) = X;
        Y = alt_new(i,:);
        Y(isnan(Y)) = interp1(find(~isnan(Y)), Y(~isnan(Y)), find(isnan(Y)), 'linear'); %linear interpolation
        alt_interp1(i,:) = Y;
        clear X Y
    end
elseif isempty(idx_NaNs1) == 1
    te_interp1 = te_new;
    alt_interp1 = alt_new;
else disp('error')
end

%% Replace intervals by NaNs
NaN_duration = 10;%[5 10 15 20 25 30];
%tg stands for time gap
idx_start_tg = 50:100:length(alt_new);

alt_poi_tg = alt_new;
ne_poi_tg = ne_new;
te_poi_tg = te_new;
ti_poi_tg = ti_new;
Time_poi_tg = Time_new;

for i = 1:length(idx_start_tg)
    alt_poi_tg(:,idx_start_tg(i):idx_start_tg(i)+NaN_duration) = NaN; 
    ne_poi_tg(:,idx_start_tg(i):idx_start_tg(i)+NaN_duration) = NaN; 
    te_poi_tg(:,idx_start_tg(i):idx_start_tg(i)+NaN_duration) = NaN; 
    ti_poi_tg(:,idx_start_tg(i):idx_start_tg(i)+NaN_duration) = NaN; 
    Time_poi_tg(:,idx_start_tg(i):idx_start_tg(i)+NaN_duration) = NaN; 
    % Time_datetime_poi_tg(:,idx_start_tg:idx_start_tg+NaN_duration) = NaT;
end

%% Check: plot period of interest with gap
f1=figure();
colormap jet;

subplot(3,1,1);
pcolor(Time_new(1,:)',alt_poi_tg,log10(ne_poi_tg)),shading flat;
% imagesc(Time',alt,log10(ne)) %does not work, because needs a vector for Y instead of a matrix (not sure how to do this with a changing altitude) --> this function should make it easier to have empty parts for missing data.
% surface(Time(1,:)',alt,log10(ne),'EdgeColor','none');%,shading flat; %different way to plot --> works

caxis([10 11.4]); % give the limits of the colourbar
colorbar; % Add colorbar to the right of the plot
set(get(colorbar,'label'),'FontSize', 8,'string','Electron density (m^{-3})')
xlim(xLimits); 
ylim(yLimits);
datetick('x',13,'keeplimits')
ylabel('Altitude [km]','FontSize', 8,'FontName','Arial') %Add labels to axes and colorbar
xlabel(['Time [UT] ',num2str(y(1)),'-',num2str(m(1)),'-',num2str(d(1))],'FontSize', 8,'FontName','Arial')
 
subplot(3,1,2)
pcolor(Time_new(1,:)',alt_poi_tg,te_poi_tg),shading flat;
caxis([0 3000]); % Limits to colorbar
colorbar; % Add colorbar to the right of the plot
set(get(colorbar,'label'),'FontSize', 8,'string','Electron temperature (K)') % Add labels to axes and colorbar
xlim(xLimits); 
ylim(yLimits); %ylim([min(alt,[],'all') 300])
datetick('x',13,'keeplimits')
ylabel('Altitude [km]','FontSize', 8,'FontName','Arial') % Add labels to axes and colorbar
xlabel(['Time [UT] ',num2str(y(1))],'FontSize', 8,'FontName','Arial') % Add labels to axes and colorbar

subplot(3,1,3)
pcolor(Time_new(1,:)',alt_poi_tg,ti_poi_tg),shading flat;
caxis([0 3000]); % Limits to colorbar
colorbar; % Add colorbar to the right of the plot
set(get(colorbar,'label'),'FontSize', 8,'string','Ion temperature (K)')
xlim(xLimits);
ylim(yLimits) ;
datetick('x',13,'keeplimits')
ylabel('Altitude [km]','FontSize', 8,'FontName','Arial') % Add labels to axes and colorbar
xlabel(['Time [UT] ',num2str(y(1))],'FontSize', 8,'FontName','Arial')

%% interpolate for NaNs (straight line)
idx_NaN_begin = find(diff(isnan(Time_poi_tg(1,:)))==1)+1;
idx_NaN_end2 = find(diff(isnan(Time_poi_tg(1,:)))==-1);
idx_NaNs = idx_NaN_end2 - idx_NaN_begin + 1;
for i = 1:size(te_poi_tg,1)
    X = te_poi_tg(i,:);
    X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)), 'linear'); %linear interpolation
    te_interp(i,:) = X;
    Y = alt_poi_tg(i,:);
    Y(isnan(Y)) = interp1(find(~isnan(Y)), Y(~isnan(Y)), find(isnan(Y)), 'linear'); %linear interpolation
    alt_interp(i,:) = Y;
    clear X Y
end

%% Spectrogram with and without data gap (_tg and no_tg)
E_region_lim = 110; %upper limit E-region/lower limit F-region
F_region_lim = 400; %upper limit F-region

for i = 1:size(Time_datetime_poi,2)-1
    delta_T(i) = Time_datetime_poi(1,i+1) - Time_datetime_poi(1,i);
end
T = round(seconds(median(delta_T))); %Time period over which EISCAT data is averaged
Fs = 1/T;  
f = linspace(0,Fs/2,500); %number of bins in y axis (number of frequency bins)

for i = 1:length(alt_interp)%(alt_interp)
%     te_av_tg(1,i) = mean(te_interp(find(alt_interp(:,i) < E_region_lim),i));% Average over E region
%     te_av_tg(2,i) = mean(te_interp(find(alt_interp(:,i) > E_region_lim & alt_interp(:,i) < F_region_lim),i));% Average over F region
%     te_av_no_tg(1,i) = mean(te_poi(find(alt_poi(:,i) < E_region_lim),i));% Average over E region
%     te_av_no_tg(2,i) = mean(te_poi(find(alt_poi(:,i) > E_region_lim & alt_poi(:,i) < F_region_lim),i));% Average over F region
    te_av_tg(1,i) = nanmean(te_interp(find(alt_interp(:,i) < E_region_lim),i));% Average over E region
    te_av_tg(2,i) = nanmean(te_interp(find(alt_interp(:,i) > E_region_lim & alt_interp(:,i) < F_region_lim),i));% Average over F region
    te_av_no_tg(1,i) = nanmean(te_interp1(find(alt_interp1(:,i) < E_region_lim),i));% Average over E region
    te_av_no_tg(2,i) = nanmean(te_interp1(find(alt_interp1(:,i) > E_region_lim & alt_interp1(:,i) < F_region_lim),i));% Average over F region
end

if length(te_av_tg) > 3600/T%length(test_Te{i}) > 3600/T
    N_1 = 60; %to have a window length of 60 data points (\approx 1 hour)
else %time period is shorter than 1 hour
    N_1 = length(te_av_tg);
end
windowOverlap = N_1 - 10;
[S_Eregion_tg,F_Eregion_tg,T_Eregion_tg,P_Eregion_tg] = spectrogram(te_av_tg(1,:),hanning(N_1),windowOverlap,f,Fs,'yaxis'); %Calculate for E-region              
[S_Fregion_tg,F_Fregion_tg,T_Fregion_tg,P_Fregion_tg] = spectrogram(te_av_tg(2,:),hann(N_1),windowOverlap,f,Fs,'yaxis'); %Calculate for F-region   
[S_Eregion_no_tg,F_Eregion_no_tg,T_Eregion_no_tg,P_Eregion_no_tg] = spectrogram(te_av_no_tg(1,:),hanning(N_1),windowOverlap,f,Fs,'yaxis'); %Calculate for E-region              
[S_Fregion_no_tg,F_Fregion_no_tg,T_Fregion_no_tg,P_Fregion_no_tg] = spectrogram(te_av_no_tg(2,:),hann(N_1),windowOverlap,f,Fs,'yaxis'); %Calculate for F-region   

%determine which unit of time is needed for the x-axis (only works if start
%time is at the hour
lab = 'Time (secs)';
if T_Fregion_tg(end) > 60*60
   T_Fregion_tg = T_Fregion_tg/(60*60);
   T_Fregion_no_tg = T_Fregion_no_tg/(60*60);
   lab = 'Time (hours)';
elseif T_Fregion_tg(end) > 60
   T_Fregion_tg = T_Fregion_tg/60;
   T_Fregion_no_tg = T_Fregion_no_tg/60;
   lab = 'Time (mins)';
end
TT_Fregion_datetime = Time_datetime_poi(1,1) + hours(T_Fregion_tg);
formatOut = 'HH:MM';
formatOut2 = 'yyyy-mm-dd';
% xlabels = cellstr(TT_Fregion_datetime,formatOut);
%%
%time gap spectrogram
figure()
h = imagesc(datenum(TT_Fregion_datetime),F_Fregion_tg,10*log10(P_Fregion_tg));%datenum(TT_Fregion_datetime),F_Fregion,10*log10(P_Fregion));
% h = imagesc(T_Fregion_tg,F_Fregion_tg,10*log10(P_Fregion_tg));%datenum(TT_Fregion_datetime),F_Fregion,10*log10(P_Fregion));
datetick('x',formatOut,'keeplimits')%,'keepticks')
xlabel(strcat(['Time (UT) (',datestr(TT_Fregion_datetime(1),formatOut2),')']))
ylabel('Frequency (Hz)')
title(strcat(['F-region; Time gaps ',num2str(NaN_duration),' min']))
colorbar;
set(gca,'YDir','normal');
set(get(colorbar,'label'),'FontSize', 12,'string','(dB/Hz)')
caxis([0 95])

%no time gap spectrogram
figure()
h = imagesc(datenum(TT_Fregion_datetime),F_Fregion_no_tg,10*log10(P_Fregion_no_tg));%datenum(TT_Fregion_datetime),F_Fregion,10*log10(P_Fregion));
% h = imagesc(T_Fregion_no_tg,F_Fregion_no_tg,10*log10(P_Fregion_no_tg));%datenum(TT_Fregion_datetime),F_Fregion,10*log10(P_Fregion));
datetick('x',formatOut,'keeplimits')%,'keepticks')
xlabel(strcat(['Time (UT) (',datestr(TT_Fregion_datetime(1),formatOut2),')']))
ylabel('Frequency (Hz)')
title(['F-region; NO time gap'])
colorbar;
set(gca,'YDir','normal');
set(get(colorbar,'label'),'FontSize', 12,'string','(dB/Hz)')
caxis([0 95])

%plot difference time gap and no time gap
figure()
h = imagesc(datenum(TT_Fregion_datetime),F_Fregion_no_tg,10*log10(P_Fregion_no_tg)-10*log10(P_Fregion_tg));%datenum(TT_Fregion_datetime),F_Fregion,10*log10(P_Fregion));
% h = imagesc(T_Fregion_no_tg,F_Fregion_no_tg,10*log10(P_Fregion_no_tg)-10*log10(P_Fregion_tg));%datenum(TT_Fregion_datetime),F_Fregion,10*log10(P_Fregion));
datetick('x',formatOut,'keeplimits')%,'keepticks')
xlabel(strcat(['Time (UT) (',datestr(TT_Fregion_datetime(1),formatOut2),')']))
ylabel('Frequency (Hz)')
title(['F-region; Difference'])
colorbar;
set(gca,'YDir','normal');
set(get(colorbar,'label'),'FontSize', 12,'string','(dB/Hz)')
