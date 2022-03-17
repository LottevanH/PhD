clear all
close all

addpath('C:\Github\PhD\functions');

%% specify which filename to process
% filename = 'C:\Users\lotte\OneDrive - Universitetssenteret på Svalbard AS (1)\Svalbard\PhD\Data\ULF_waves_events_currently_unexamined\2017_12_18\EISCAT_2017-12-18_bella_60@vhf.hdf5';
% filename = 'C:\Users\lotte\OneDrive - Universitetssenteret på Svalbard AS (1)\Svalbard\PhD\Data\ULF_waves_events_currently_unexamined\1998_12_20\analysed_01_04\analysed_1998_12_20_01_till_04_60s_integration\hdf5_folder\EISCAT_1998-12-20_gup3c_60@32m\EISCAT_1998-12-20_gup3c_60@32m.hdf5';
filename = 'C:\Users\lotte\OneDrive - Universitetssenteret på Svalbard AS (1)\Svalbard\PhD\Data\Paper_Lisa_2017\EISCAT_2007-12-27_ipy_60@42m.hdf5';

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

%% Plot complete time period
n1 = Time(1,1); 
n2 = Time(2,end); 
xLimits=[n1 n2];
yLimits=[min(alt,[],'all'); 300];
f1=figure();
colormap jet;

subplot(3,1,1);
pcolor(Time(1,:)',alt,log10(ne)),shading flat;
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
pcolor(Time(1,:)',alt,te),shading flat;
caxis([0 3000]); % Limits to colorbar
colorbar; % Add colorbar to the right of the plot
set(get(colorbar,'label'),'FontSize', 8,'string','Electron temperature (K)') % Add labels to axes and colorbar
xlim(xLimits); 
ylim(yLimits); %ylim([min(alt,[],'all') 300])
datetick('x',13,'keeplimits')
ylabel('Altitude [km]','FontSize', 8,'FontName','Arial') % Add labels to axes and colorbar
xlabel(['Time [UT] ',num2str(y(1))],'FontSize', 8,'FontName','Arial') % Add labels to axes and colorbar

subplot(3,1,3)
pcolor(Time(1,:)',alt,ti),shading flat;
caxis([0 3000]); % Limits to colorbar
colorbar; % Add colorbar to the right of the plot
set(get(colorbar,'label'),'FontSize', 8,'string','Ion temperature (K)')
xlim(xLimits);
ylim(yLimits) ;
datetick('x',13,'keeplimits')
ylabel('Altitude [km]','FontSize', 8,'FontName','Arial') % Add labels to axes and colorbar
xlabel(['Time [UT] ',num2str(y(1))],'FontSize', 8,'FontName','Arial')
 
%% Time period of interest
n11 = datenum(datetime(2007,12,27,15,00,00)); %in case of limited time
n22 = datenum(datetime(2007,12,27,17,30,00)); %in case of limited time
% n11 = datenum(datetime(2017,12,18,02,00,00)); %in case of limited time
% n22 = datenum(datetime(2017,12,18,07,00,00)); %in case of limited time
% n11 = datenum(datetime(1998,12,20,01,00,00)); %in case of limited time
% n22 = datenum(datetime(1998,12,20,04,00,00)); %in case of limited time
% n11 = n1;
% n22 = n2;

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
        elseif (i > 1) && (i < length(idx_time_gap)) %%check this one!
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
    Time_datetime_new(:,:) = Time_datetime_poi(:,:);
end

addpath('C:\Github\PhD\functions\Inpaint_nans')
% alt_new = inpaint_nans(alt_new); %interpolate over the NaN values so that these data points are taken into account as well later on

%% Plot with time gaps
xLimits=[n11 n22];
% yLimits=[min(alt,[],'all'); 300];
yLimits=[100 350];
f1=figure();
colormap jet;

subplot(3,1,1);
pcolor(Time_new(1,:)',alt_new,log10(ne_new)),shading flat;
caxis([10 11.4]);
colorbar;
set(get(colorbar,'label'),'FontSize', 8,'string','Electron density (m^{-3})')
xlim(xLimits);
ylim(yLimits) ;
datetick('x',13,'keeplimits')
ylabel('Altitude [km]','FontSize', 8,'FontName','Arial')
xlabel(['Time [UT] ',num2str(y(1)),'-',num2str(m(1)),'-',num2str(d(1))],'FontSize', 8,'FontName','Arial')
 
subplot(3,1,2)
pcolor(Time_new(1,:)',alt_new,te_new),shading flat;
caxis([0 3000]);
colorbar;
set(get(colorbar,'label'),'FontSize', 8,'string','Electron temperature (K)')
xlim(xLimits);
ylim(yLimits);
datetick('x',13,'keeplimits')
ylabel('Altitude [km]','FontSize', 8,'FontName','Arial')
xlabel(['Time [UT] ',num2str(y(1))],'FontSize', 8,'FontName','Arial')

subplot(3,1,3)
pcolor(Time_new(1,:)',alt_new,ti_new),shading flat;
caxis([0 3000]);
colorbar;
set(get(colorbar,'label'),'FontSize', 8,'string','Ion temperature (K)')
xlim(xLimits);
ylim(yLimits);
datetick('x',13,'keeplimits')
ylabel('Altitude [km]','FontSize', 8,'FontName','Arial')
xlabel(['Time [UT] ',num2str(y(1))],'FontSize', 8,'FontName','Arial')

%% Split into time periods (problem: would rather interpolate over a data gap of say < 5 minutes --> now not done)
idx_NaN_end = find(diff(isnan(Time_new(1,:)))==1);
idx_NaN_end1 = [idx_NaN_end length(Time_new)];
idx_NaN_start = find(diff(isnan(Time_new(1,:)))==-1)+1;
idx_NaN_start1 = [1 idx_NaN_start];


for i = 1:length(idx_NaN_start1)
    alt1{i} = alt_new(:,idx_NaN_start1(i):idx_NaN_end1(i));
    ne1{i} = ne_new(:,idx_NaN_start1(i):idx_NaN_end1(i));
    te1{i} = te_new(:,idx_NaN_start1(i):idx_NaN_end1(i));
    Time1{i} = Time_new(:,idx_NaN_start1(i):idx_NaN_end1(i));
    Time_datetime1{i} = Time_datetime_new(:,idx_NaN_start1(i):idx_NaN_end1(i));
end

%% interpolate for NaNs (straight line)
idx_NaN_begin = find(diff(isnan(Time_new(1,:)))==1)+1;
idx_NaN_end2 = find(diff(isnan(Time_new(1,:)))==-1);
idx_NaNs = idx_NaN_end2 - idx_NaN_begin + 1;
for i = 1:size(te_new,1)
    X = te_new(i,:);
    X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)), 'linear'); %linear interpolation
    te_interp(i,:) = X;
    Y = alt_new(i,:);
    Y(isnan(Y)) = interp1(find(~isnan(Y)), Y(~isnan(Y)), find(isnan(Y)), 'linear'); %linear interpolation
    alt_interp(i,:) = Y;
    clear X Y
end
% for i = 1:length(idx_NaNs)
%     for j = 1:idx_NaNs(i)
%         Time_new(1,idx_NaN_begin+j-1) = Time_new(1,idx_NaN_begin(i)-1)+(Time_new(1,idx_NaN_end(i)+1)-Time_new(1,idx_NaN_begin(i)-1))/idx_NaNs(i)*j;
%         Time_new(2,idx_NaN_begin+j-1) = Time_new(2,idx_NaN_begin(i)-1)+(Time_new(2,idx_NaN_end(i)+1)-Time_new(2,idx_NaN_begin(i)-1))/idx_NaNs(i)*j;
%     end
% end

% addpath('C:\Github\PhD\functions\Inpaint_nans')
% % alt1 = inpaint_nans(alt_new); %interpolate over the NaN values so that these data points are taken into account as well later on
% X = Time_new(1,:);
% X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)), 'linear'); 


%% spectrogram of interpolated values
E_region_lim = 110;
F_region_lim = 400;
low_lim = 220;
high_lim = 280;
keep('E_region_lim','F_region_lim','te_interp','S1','F1','T1','P1','low_lim','high_lim','ne1','Time1','Time_datetime1','te1','alt1','m','delta_T','test_Te','test_Ne','N_1','N_Te','windowOverlap') %to make sure that all variables underneath are only calculated for the correct time period
% alt_lim = low_lim:5:high_lim;
% for i = 1:size(alt_lim,2)-1 %looking at different altitudes
%     alt_lim_mid(i) = (alt_lim(i)+alt_lim(i+1))/2;
% end
for i = 1:length(alt_interp)
    te_av(1,i) = mean(te_interp(find(alt_interp(:,i) < E_region_lim),i));% Average over E region
    te_av(2,i) = mean(te_interp(find(alt_interp(:,i) > E_region_lim & alt_interp(:,i) < F_region_lim),i));% Average over F region
end

T = round(seconds(median(delta_T))); %Time period over which EISCAT data is averaged
Fs = 1/T;  
f = linspace(0,Fs/2,500); %number of bins in y axis (number of frequency bins)
for i = 2:length(alt_lim)-1
    test_Te{i} = te_interp(alt_interp>alt_lim(i) & alt_interp<alt_lim(i+1));
    if size(test_Te{i},1) > size(te_interp,2)/2 %Make sure that only the altitudes with a signal are used
        N_Te(i) = length(test_Te{i});
        if length(test_Te{i}) > 3600/T
            N_1(i) = 60; %to have a window length of 60 data points (\approx 1 hour)
        else %time period is shorter than 1 hour
            N_1(i) = length(test_Te{i});
        end
        windowOverlap(i) = N_1(i) - 10;
        [S,F,TT,P] = spectrogram(test_Te{i},hanning(N_1(i)),windowOverlap(i),f,Fs,'yaxis');              
        S1{i} = S;
        F1{i} = F;
        T1{i} = TT;
        P1{i} = P;

        figure()
        lab = 'Time (secs)';
        if TT(end) > 60*60
           TT = TT/(60*60);
           lab = 'Time (hours)';
        elseif TT(end) > 60
           TT = TT/60;
           lab = 'Time (mins)';
        end
        h = imagesc(TT,F,10*log10(P));%.*F);%,[-156.5 95.8] );
        xlabel(lab)
        ylabel('Frequency (Hz)')
        title(strcat(['Alt ',num2str(alt_lim(i)),'-',num2str(alt_lim(i+1)),' km; Time ',datestr(Time_datetime_new(1),'HH:MM:SS'),'-',datestr(Time_datetime_new(end),'HH:MM:SS')]))
        colorbar;
        set(gca,'YDir','normal');
        set(get(colorbar,'label'),'FontSize', 12,'string','(dB/Hz)')
    else ;
    end
end



%% spectrogram 
% 
% low_lim = 250;
% high_lim = 270;
% for m = 1:length(alt1) %m represents the number of time periods without gaps in them
%     keep('S1','F1','T1','P1','low_lim','high_lim','ne1','Time1','Time_datetime1','te1','alt1','m','delta_T','test_Te','test_Ne','N_1','N_Te','windowOverlap') %to make sure that all variables underneath are only calculated for the correct time period
%     alt_lim = low_lim:5:high_lim;
%     for i = 1:size(alt_lim,2)-1 %looking at different altitudes
%         alt_lim_mid(i) = (alt_lim(i)+alt_lim(i+1))/2;
%     end
%     Te_lim1 = te1{m};
%     Ne_lim1 = ne1{m};
%     Time_vector = Time1{m};
% 
% %     %divide per hour?
% %     delta_T1 = round(3600/64);
% %     delta_T_array = 1:delta_T1:size(Te_lim1,2);
% 
%     % T = 1/Fs;
%     T = round(seconds(median(delta_T))); %Time period over which EISCAT data is averaged
%     Fs = 1/T; 
% %     windowOverlap = [50]; %want to have 10 or 15 minutes (in an hour). So  
%     f = linspace(0,Fs/2,500); %number of bins in y axis (number of frequency bins)
%     for i = 2:length(alt_lim)-1
%         test_Te{m,i} = Te_lim1(alt1{m}>alt_lim(i) & alt1{m}<alt_lim(i+1));
%         test_Ne{m,i} = Ne_lim1(alt1{m}>alt_lim(i) & alt1{m}<alt_lim(i+1));
%         if size(test_Te{m,i},1) > size(Te_lim1,2)/2 %Make sure that only the altitudes with a signal are used
%             
%             N_Te(m,i) = length(test_Te{m,i});
%             if length(test_Te{m,i}) > 3600/T
%                 N_1(m,i) = 60; %to have a window length of 60 data points (\approx 1 hour)
%             else %time period is shorter than 1 hour
%                 N_1(m,i) = length(test_Te{m,i});
%             end
%             windowOverlap(m,i) = N_1(m,i) - 10;
% %             N_1(m,i) = 2^(nextpow2(N_Te(m,i))-4);%round(N_Te(i)/8);%2^nextpow2(N_Te(i)); %same for Te and Ne
% %             f{i} = Fs*(0:(N_1(i)/2))/N_1(i); %frequency domain f
%             
% %                 figure()
% %                 spectrogram(test_Te{m,i},hanning(N_1(m,i)),windowOverlap(l),f,Fs,'yaxis')
% %                 title(strcat(['Alt ',num2str(alt_lim(i)),'-',num2str(alt_lim(i+1)),' km; Time ',datestr(Time_datetime1{m}(1),'HH:MM:SS'),'-',datestr(Time_datetime1{m}(end),'HH:MM:SS')]))
% %                 [S,F,TT,P] = spectrogram(test_Te{m,i},hanning(N_1(m,i)),windowOverlap(m,i),f,Fs,'yaxis');
%                 [S,F,TT,P] = spectrogram(test_Te{m,i},hanning(N_1(m,i)),windowOverlap(m,i),f,Fs,'yaxis');              
%                 S1{m,i} = S;
%                 F1{m,i} = F;
%                 T1{m,i} = TT;
%                 P1{m,i} = P;
%                 
%                 
%                 figure()
%                 lab = 'Time (secs)';
%                 if TT(end) > 60*60
%                    TT = TT/(60*60);
%                    lab = 'Time (hours)';
%                 elseif TT(end) > 60
%                    TT = TT/60;
%                    lab = 'Time (mins)';
%                 end
%                 h = imagesc(TT,F,10*log10(P));%.*F);%,[-156.5 95.8] );
%                 xlabel(lab)
%                 ylabel('Frequency (Hz)')
%                 title(strcat(['Alt ',num2str(alt_lim(i)),'-',num2str(alt_lim(i+1)),' km; Time ',datestr(Time_datetime1{m}(1),'HH:MM:SS'),'-',datestr(Time_datetime1{m}(end),'HH:MM:SS')]))
%                 colorbar;
%                 set(gca,'YDir','normal');
%                 set(get(colorbar,'label'),'FontSize', 12,'string','(dB/Hz)')
%         else ;
%         end
%     end
% end
% 
