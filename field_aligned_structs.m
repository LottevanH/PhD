%for low elevation scans, not possible to integrate in altitude (over E and
%F region), because of the spatial differences. So need to evaluate
%altitudes separately.
%examples of low elevation modes: 
% - ESR: gup3c, tau0, hilde, folke(?)
% - VHF: bella

%Examples of field-aligned modes (ESR: elevation ~81.5; UHF/VHF elevation ~77.4):
% - VHF/UHF: Beata (sometimes in scanning mode), arc1
% - ESR: 

clear all
close all

addpath('C:\Github\PhD\functions');


%% specify which filename to process
% filename = 'C:\data\1998-12-20_gup3c_20@32m\hdf5_folder\EISCAT_1998-12-20_gup3c_20@32m\EISCAT_1998-12-20_gup3c_20@32m.hdf5';
% filename = 'C:\Users\lotte\OneDrive - Universitetssenteret på Svalbard AS (1)\Svalbard\PhD\Data\ULF_waves_events_currently_unexamined\1998_12_20\analysed_01_04\analysed_1998_12_20_01_till_04_60s_integration\hdf5_folder\EISCAT_1998-12-20_gup3c_60@32m\EISCAT_1998-12-20_gup3c_60@32m.hdf5';
% filename = 'C:\Users\lotte\OneDrive - Universitetssenteret på Svalbard AS (1)\Svalbard\PhD\Data\ULF_waves_events_currently_unexamined\1998_12_20\analysed_01_04\analysed_1998_12_20_01_till_04_120s_integration\hdf5_folder\EISCAT_1998-12-20_gup3c_120@32m\EISCAT_1998-12-20_gup3c_120@32m.hdf5';
filename = 'C:\Users\lotte\OneDrive - Universitetssenteret på Svalbard AS (1)\Svalbard\PhD\Data\Test_data\EISCAT_2008-06-04_beata_60@uhfa.hdf5';
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
vi=par2D(:,:,6);
Time_datetime = datetime(Time,'ConvertFrom','datenum'); 
sys_temp = par1D(:,4);
az = par1D(:,1);
el = par1D(:,2);
peak_power = par1D(:,3);
%% Plot complete time period
n1 = Time(1,1); 
n2 = Time(2,end); 
xLimits=[n1 n2];
yLimits=[min(alt,[],'all'); 400];
% yLimits=[min(alt,[],'all'); 600];
f1=figure();
colormap jet;

subplot(3,1,1);
pcolor(Time(1,:)',alt,log10(ne)),shading flat;
caxis([10 12]); % give the limits of the colourbar
colorbar; % Add colorbar to the right of the plot
set(get(colorbar,'label'),'FontSize', 8,'string','Electron density (m^{-3})')
xlim(xLimits); 
ylim(yLimits);
datetick('x',13,'keeplimits')
ylabel('Altitude [km]','FontSize', 8,'FontName','Arial') %Add labels to axes and colorbar
xlabel(['Time [UT] ',num2str(y(1)),'-',num2str(m(1)),'-',num2str(d(1))],'FontSize', 8,'FontName','Arial')
 
subplot(3,1,2)
pcolor(Time(1,:)',alt,te),shading flat;
caxis([0 4000]); % Limits to colorbar
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
%analyse from one hour before start interesting event
% n11 = datenum(datetime(2017,12,18,01,00,00)); %in case of limited time
% n22 = datenum(datetime(2017,12,18,08,00,00)); %in case of limited time
n11 = n1;
n22 = n2;

idx_time_start = find(Time(1,:) > n11,1);
idx_time_end = find(Time(1,:) < n22,1,'last');
 
poi.alt = alt(:,idx_time_start:idx_time_end); %poi stands for: period of interest
poi.ne = ne(:,idx_time_start:idx_time_end);
poi.te = te(:,idx_time_start:idx_time_end);
poi.ti = ti(:,idx_time_start:idx_time_end);
poi.vi = vi(:,idx_time_start:idx_time_end);
poi.Time = Time(:,idx_time_start:idx_time_end);
poi.Time_datetime = Time_datetime(:,idx_time_start:idx_time_end);

%% Sort per altitude
% MORE EFFICIENT WITH STRUCT STRUCTURE
variables = {'alt','ne','te','ti','vi','Time','Time_datetime'}; %for struct setup

%find median altitudes of all the rows to determine what kind of altitude
%spacing is needed
median_alt_row = median(poi.alt,2,'omitnan');
for i = 1:length(median_alt_row)-1
alt_lim1(i+1) = (median_alt_row(i) + median_alt_row(i+1))/2;
end
alt_lim1(1) = median_alt_row(1)-(alt_lim1(2)-median_alt_row(1));
alt_lim1(length(median_alt_row)+1) = 2*median_alt_row(end) - alt_lim1(length(median_alt_row));

for i = 1:length(variables) 
    poi_sorted.(variables{i}) = NaN(size(poi.(variables{i})));
end

for k = 1:size(poi.alt,1)
    [row{k},col{k}] = find(poi.alt > alt_lim1(k) & poi.alt < alt_lim1(k+1));
    for i = 1:length(variables)-2 %do not want to reshuffle the times
        for j = 1:length(col{k})
            poi_sorted.(variables{i})(k,col{k}(j)) = poi.(variables{i})(row{k}(j),col{k}(j));
        end
    end
    for i = [6,7] %Time and Time_datetime
        poi_sorted.(variables{i}) = poi.(variables{i});
    end
end

 
%% Add NaNs in case of time gaps without detrending

for i = 1:size(poi_sorted.Time_datetime,2)-1
    delta_T(i) = poi_sorted.Time_datetime(1,i+1) - poi_sorted.Time_datetime(1,i);
end
idx_time_gap = find(delta_T > median(delta_T)+seconds(15)); %find time gaps if there are any
for j = 1:length(variables)
    if length(idx_time_gap) == 1
        i = 1;
        number_of_NaN(i) = round((delta_T(idx_time_gap(i))+delta_T(idx_time_gap(i) + 1 ))/median(delta_T));
        NaN_data.(variables{j})(:,1:idx_time_gap(i)-1) = poi_sorted.(variables{j})(:,1:idx_time_gap(i)-1);
        if isdatetime(poi_sorted.(variables{j})) == 1
            NaN_data.(variables{j})(:,idx_time_gap(i):idx_time_gap(i)+number_of_NaN(i)-1) = NaT;
        else
            NaN_data.(variables{j})(:,idx_time_gap(i):idx_time_gap(i)+number_of_NaN(i)-1) = NaN;
        end
        NaN_data.(variables{j})(:,idx_time_gap(i)+number_of_NaN(i):length(poi_sorted.Time_datetime)+number_of_NaN-2) = poi_sorted.(variables{j})(:,idx_time_gap(i)+2:end);

    elseif length(idx_time_gap) > 1 %more than one time gap
        for i = 1:length(idx_time_gap)
            if i == 1
                number_of_NaN(i) = round((delta_T(idx_time_gap(i))+delta_T(idx_time_gap(i) + 1))/median(delta_T));
                NaN_data.(variables{j})(:,1:idx_time_gap(i)-1) = poi_sorted.(variables{j})(:,1:idx_time_gap(i)-1);
                if isdatetime(poi_sorted.(variables{j})) == 1
                    NaN_data.(variables{j})(:,idx_time_gap(i):idx_time_gap(i)+number_of_NaN(i)-1) = NaT;
                else
                    NaN_data.(variables{j})(:,idx_time_gap(i):idx_time_gap(i)+number_of_NaN(i)-1) = NaN;
                end
                NaN_data.(variables{j})(:,idx_time_gap(i)+number_of_NaN(i):idx_time_gap(i+1)-1+number_of_NaN(i)-2) = poi_sorted.(variables{j})(:,idx_time_gap(i)+2:idx_time_gap(i+1)-1);
            elseif (i > 1) && (i < length(idx_time_gap)) 
                number_of_NaN(i) = round((delta_T(idx_time_gap(i))+delta_T(idx_time_gap(i) + 1 ))/median(delta_T));
                if isdatetime(poi_sorted.(variables{j})) == 1
                    NaN_data.(variables{j})(:,idx_time_gap(i)+sum(number_of_NaN(1:i-1)-i):idx_time_gap(i)+sum(number_of_NaN(1:i))-i) = NaT;
                else
                    NaN_data.(variables{j})(:,idx_time_gap(i)+sum(number_of_NaN(1:i-1)-i):idx_time_gap(i)+sum(number_of_NaN(1:i))-i) = NaN;
                end
                NaN_data.(variables{j})(:,idx_time_gap(i)+sum(number_of_NaN(1:i))-i:idx_time_gap(i+1)-1+sum(number_of_NaN(1:i))-2-i) = poi_sorted.(variables{j})(:,idx_time_gap(i)+2:idx_time_gap(i+1)-1);
            elseif i == length(idx_time_gap)
                number_of_NaN(i) = round((delta_T(idx_time_gap(i))+delta_T(idx_time_gap(i) + 1 ))/median(delta_T));
                if isdatetime(poi_sorted.(variables{j})) == 1
                    NaN_data.(variables{j})(:,idx_time_gap(i)+sum(number_of_NaN(1:i-1)-i):idx_time_gap(i)+sum(number_of_NaN(:))-1-i) = NaT;
                else
                    NaN_data.(variables{j})(:,idx_time_gap(i)+sum(number_of_NaN(1:i-1)-i):idx_time_gap(i)+sum(number_of_NaN(:))-1-i) = NaN;
                end
                NaN_data.(variables{j})(:,idx_time_gap(i)+sum(number_of_NaN)-i:length(poi_sorted.Time_datetime)+sum(number_of_NaN)-2-i) = poi_sorted.(variables{j})(:,idx_time_gap(i)+2:end);
            else
                disp('Something is wrong')
            end            
        end

    else %so no time gaps
        NaN_data.(variables{j}) = poi_sorted.(variables{j});
    end
end

%% Plot with time gaps
xLimits=[n11 n22];
% yLimits=[min(alt,[],'all'); 300];
yLimits=[100 350];
f1=figure();
colormap jet;

subplot(4,1,1);
pcolor(NaN_data.Time(1,:)',NaN_data.alt,log10(NaN_data.ne)),shading flat;
caxis([10 11.8]);
colorbar;
set(get(colorbar,'label'),'FontSize', 8,'string','Electron density (m^{-3})')
xlim(xLimits);
ylim(yLimits) ;
datetick('x',13,'keeplimits')
ylabel('Altitude [km]','FontSize', 8,'FontName','Arial')
xlabel(['Time [UT] ',num2str(y(1)),'-',num2str(m(1)),'-',num2str(d(1))],'FontSize', 8,'FontName','Arial')
 
subplot(4,1,2)
pcolor(NaN_data.Time(1,:)',NaN_data.alt,NaN_data.te),shading flat;
caxis([0 3500]);
colorbar;
set(get(colorbar,'label'),'FontSize', 8,'string','Electron temperature (K)')
xlim(xLimits);
ylim(yLimits);
datetick('x',13,'keeplimits')
ylabel('Altitude [km]','FontSize', 8,'FontName','Arial')
xlabel(['Time [UT] ',num2str(y(1))],'FontSize', 8,'FontName','Arial')

subplot(4,1,3)
pcolor(NaN_data.Time(1,:)',NaN_data.alt,NaN_data.ti),shading flat;
caxis([0 1500]);
colorbar;
set(get(colorbar,'label'),'FontSize', 8,'string','Ion temperature (K)')
xlim(xLimits);
ylim(yLimits);
datetick('x',13,'keeplimits')
ylabel('Altitude [km]','FontSize', 8,'FontName','Arial')
xlabel(['Time [UT] ',num2str(y(1))],'FontSize', 8,'FontName','Arial')

subplot(4,1,4)
% pcolor(NaN_data.Time(1,:)',NaN_data.alt,NaN_data.vi),shading flat;
% caxis([-100 100]);
% colorbar;
% set(get(colorbar,'label'),'FontSize', 8,'string','ion velocity [m/s]')
semilogy(NaN_data.Time(1,:)',az)
hold on
semilogy(NaN_data.Time(1,:)',el)
semilogy(NaN_data.Time(1,:)',sys_temp)
semilogy(NaN_data.Time(1,:)',peak_power)
xlim(xLimits);
% ylim(yLimits);
datetick('x',13,'keeplimits')
% ylabel('Altitude [km]','FontSize', 8,'FontName','Arial')
legend('az','el','Tsys (K)','peak power (kW)')
xlabel(['Time [UT] ',num2str(y(1))],'FontSize', 8,'FontName','Arial')

%% Average over E and F region
E_lower_lim = 80;
E_high_lim = 110;
F_high_lim = 350;
for i = 2:5
    poi_averaged.(variables{i}).E_region = nanmedian(NaN_data.(variables{i})(find(NaN_data.alt(:,1) > E_lower_lim & NaN_data.alt(:,1) < E_high_lim),:));
    poi_averaged.(variables{i}).F_region = nanmedian(NaN_data.(variables{i})(find(NaN_data.alt(:,1) > E_high_lim & NaN_data.alt(:,1) < F_high_lim),:));
end
for i = 6:7
    poi_averaged.(variables{i}) = NaN_data.(variables{i});
end

%% High pass filter (subtract median value of sliding window of 10 minutes lenght for Pc5)
str_region = {'E_region','F_region'};

for i = 2:5 %only apply the high pass filter for ne, te, ti and vi (NOT for alt, Time and Time_datetime)
    for j = 1:length(str_region)
        for k = 5:length(poi_averaged.(variables{i}).(str_region{j}))-5
            filtered_data.(variables{i}).(str_region{j})(k) = poi_averaged.(variables{i}).(str_region{j})(k) - median(poi_averaged.(variables{i}).(str_region{j})(k-4:k+5),'omitnan');
        end
        for k = 1:4
            filtered_data.(variables{i}).(str_region{j})(k) = poi_averaged.(variables{i}).(str_region{j})(k) - median(poi_averaged.(variables{i}).(str_region{j})(1:10),'omitnan');
        end
        for k = length(poi_averaged.(variables{i}).(str_region{j}))-4:length(poi_averaged.(variables{i}).(str_region{j}))
            filtered_data.(variables{i}).(str_region{j})(k) = poi_averaged.(variables{i}).(str_region{j})(k) - median(poi_averaged.(variables{i}).(str_region{j})(length(poi_averaged.(variables{i}).(str_region{j}))-9:end),'omitnan');
        end
    end
end

for i = setdiff(2:length(variables),2:5) %exclude 2:5
    filtered_data.(variables{i}) = poi_averaged.(variables{i});
end

%% interpolate for NaNs (straight line)
idx_NaN_begin = find(diff(isnan(NaN_data.Time(1,:)))==1)+1;
idx_NaN_end2 = find(diff(isnan(NaN_data.Time(1,:)))==-1);
idx_NaNs = idx_NaN_end2 - idx_NaN_begin + 1;
for j = 2:5
    for i = 1:length(str_region)
            X = filtered_data.(variables{j}).(str_region{i});
            X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)), 'linear'); %linear interpolation
            interp_data.(variables{j}).(str_region{i}) = X;
            clear X 
    end
end

for j = 6
    X = filtered_data.(variables{j});
    X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)), 'linear'); %linear interpolation
    interp_data.(variables{j}) = X;
    clear X 
end
    

%% Input signal to FFT
parameters_of_interest = [2,3,4];%3; %setdiff(2:length(variables),6:7);
low_lim = 200;
high_lim = 320;
row_interest = find(median_alt_row > low_lim & median_alt_row < high_lim);
median_alt_row1 = median_alt_row(row_interest);
alt_str = strcat(num2str(round(median_alt_row1(:))),'km');
alt_str1 = {'one','two','three','four'};
for i = 1:length(str_region)
    for j = parameters_of_interest
        [input_pks{i,j},input_locs{i,j}] = findpeaks(interp_data.(variables{j}).(str_region{i}),interp_data.Time(1,:));%,'MinPeakHeight',max(powerspectrum{i,j},[],'all')-10,'MinPeakProminence',7); 
        [input_pks1{i,j},input_locs1{i,j}] = findpeaks(-interp_data.(variables{j}).(str_region{i}),interp_data.Time(1,:));%,'MinPeakHeight',max(powerspectrum{i,j},[],'all')-10,'MinPeakProminence',7); 
        [input_locs_tot{i,j}, input_idx{i,j}] = sort([input_locs{i,j} input_locs1{i,j}]);
        input_pks_tot1{i,j} = [input_pks{i,j} -input_pks1{i,j}];
        input_pks_tot{i,j} = input_pks_tot1{i,j}(input_idx{i,j});
        amplitude{i,j} = diff(input_pks_tot{i,j});
        pks_time = datetime(input_locs{i,j},'ConvertFrom','datenum');
        min_time = datetime(input_locs1{i,j},'ConvertFrom','datenum');
        delta_T_pks.(alt_str1{i}).(variables{j}) = diff(pks_time);
        delta_T_min.(alt_str1{i}).(variables{j}) = diff(min_time);
        figure()
        plot(interp_data.Time(1,:),interp_data.(variables{j}).(str_region{i}))
        title(strcat([variables{j},'; ',str_region{i},'; ',datestr(NaN_data.Time_datetime(1),'yyyy-mm-dd')]))
        datetick('x',13,'keeplimits')
        xlim([n11 n22])
        hold on
%         plot(input_locs{i,j},input_pks{i,j},'*r')
%         plot(input_locs1{i,j},-input_pks1{i,j},'*b')
        plot(input_locs_tot{i,j},input_pks_tot{i,j},'*g')
%         plot(input_locs_tot{i,j}(1:end-1),abs(amplitude{i,j}),'k')
        clear pks_time min_time
    end
end


%% spectrogram of interpolated values for different altitudes
T = round(seconds(mean(delta_T))); %Time period over which EISCAT data is averaged
Fs = 1/T;  
f = linspace(0,Fs/2,500); %number of bins in y axis (number of frequency bins)
% detrnd = 2; %order of detrending (quadratical for 
window_length_seconds = 3600; 
window_overlap_seconds = 600;

for i = 1:length(str_region)%length(suitable_rows)%length(alt_lim)-1
    for j = parameters_of_interest % 1:length(variables)
        if length(interp_data.(variables{j}).(str_region{i})) > window_length_seconds/T %meaning the period of interest is larger than 1 hour.
            N_1 = round(window_length_seconds/T);%60; %to have a window length of 60 data points (\approx 1 hour)
        else %time period is shorter than 1 hour
            N_1 = length(interp_data.(variables{j}).(str_region{i}));
        end
        windowOverlap = N_1 - round(window_overlap_seconds/T);%overlap of 10 minutes; 10;
%         [S_Eregion{i},F_Eregion{i},T_Eregion{i},P_Eregion{i}] = spectrogram(te_alt{i},hanning(N_1),windowOverlap,f,Fs,'yaxis'); %Calculate for E-region              
        [S_Fregion{i,j},F_Fregion{i,j},T_Fregion{i,j},P_Fregion{i,j}] = spectrogram(interp_data.(variables{j}).(str_region{i}),hamming(N_1),windowOverlap,f,Fs,'yaxis'); %Calculate for F-region   
        idx_freq_pc5 = find(F_Fregion{i,j}> 0.6e-3 & F_Fregion{i,j}<7e-3);
        powerspectrum{i,j} = 10*log10(P_Fregion{i,j}(idx_freq_pc5,:));
        max_power(i,j) = max(powerspectrum{i,j},[],'all');
        min_power(i,j) = min(powerspectrum{i,j},[],'all');

    end
end
% idx_freq_pc5{i,j} = find(F_Fregion{i,j}> 0.6e-3 & F_Fregion{i,j}<7e-3); %between 0.6 and 7.0 mHz


for i = 1:length(str_region)%length(suitable_rows)%length(alt_lim)-1
    for j = parameters_of_interest % 1:length(variables)
        for k = 5:size(powerspectrum{i,j},2)-5
            [pks{i,j,k},locs{i,j,k}] = findpeaks(powerspectrum{i,j}(:,k), F_Fregion{i,j}(idx_freq_pc5),'MinPeakHeight',max(powerspectrum{i,j},[],'all')-10,'MinPeakProminence',7); %IDEA: ONLY USE THE MAXIMUM POWER WITHIN A 3 HOUR INTERVAL AROUND THE POINT K
        end
        %determine which unit of time is needed for the x-axis
        lab = 'Time (secs)';
        TT_Fregion = T_Fregion{i,j};
        if TT_Fregion(end) > 60*60
           T_Fregion{i,j} = T_Fregion{i,j}/(60*60);
           lab = 'Time (hours)';
        elseif TT_Fregion(end) > 60
           T_Fregion{i,j} = T_Fregion{i,j}/60;
           lab = 'Time (mins)';
        end
        clear TT_Fregion
        TT_Fregion_datetime{i,j} = NaN_data.Time_datetime(1,1) + hours(T_Fregion{i,j});
        formatOut = 'HH:mm:ss';
        Time1 = datenum(TT_Fregion_datetime{i,j});
        %%
        integrated_power = bandpower(powerspectrum{i,j},F_Fregion{i,j}(idx_freq_pc5),'psd');
        integrated_power1 = bandpower(P_Fregion{i,j}(idx_freq_pc5,:),F_Fregion{i,j}(idx_freq_pc5),'psd'); %bandpower(powerspectrum{i,j},F_Fregion{i,j}(idx_freq_pc5),[0.6e-3 7e-3],'psd');
%%
%         integrated_power{i,j} = sum(powerspectrum{i,j});

        figure()
        h = imagesc(datenum(TT_Fregion_datetime{i,j}),F_Fregion{i,j}(idx_freq_pc5,:)*1000,powerspectrum{i,j});%10*log10(P_Fregion{i,j}));%datenum(TT_Fregion_datetime),F_Fregion,10*log10(P_Fregion));
        xlabel('Time')
        ylabel('Frequency (mHz)')
        ylim([0 1000/120])
%         title(strcat(['Altitude: ',num2str(alt_lim(i)),'-',num2str(alt_lim(i+1)),'; ',datestr(Time_datetime_new(1),'yyyy-mm-dd')]))
        title(strcat([variables{j},'; ',str_region{i},'; ',datestr(NaN_data.Time_datetime(1),'yyyy-mm-dd')]))
        colorbar;
        colormap parula
        caxis([max(max_power(:,j))-70 max(max_power(:,j))])
        set(gca,'YDir','normal');
        set(get(colorbar,'label'),'FontSize', 12,'string','(dB/Hz)')
        datetick('x',13,'keeplimits')
%         xlim([datenum(datetime(n11,'ConvertFrom','datenum')+seconds(window_length_seconds))
%         datenum(datetime(n22,'ConvertFrom','datenum')-seconds(window_length_seconds))])
%         %only works in case the period of interest exceeds 2 hours
        hold on 
        for k = 5:size(powerspectrum{i,j},2)-5 
            if isempty(locs{i,j,k}) == 0
                plot(Time1(k),locs{i,j,k}*1000,'.r')
            else ;
            end
        end
%         clear Time1
            
    end
end
%% plot the peaks for all altitudes
% color_alt = jet(length(median_alt_row1));
for j = parameters_of_interest
    figure()
    hold on
    ylim([0 1000/120])
%     xlim([datenum(datetime(n11,'ConvertFrom','datenum')+seconds(window_length_seconds)) datenum(datetime(n22,'ConvertFrom','datenum')-seconds(window_length_seconds))])
    datetick('x',13,'keeplimits')
    title(strcat([variables{j},'; ',datestr(NaN_data.Time_datetime(1),'yyyy-mm-dd')]))
    ylabel('frequency [mHz]')
    xlabel('Time')
    box on
    for i = 1:size(alt_str,1)
        plot(n11,0,'*','color',color_alt(i,:))
    end
    legend(alt_str)
    for i = 1:length(median_alt_row1)
        for k = 5:size(P_Fregion{i,j},2)-5
            if isempty(locs{i,j,k}) == 0
                plot(Time1(k),locs{i,j,k}*1000,'*','color',color_alt(i,:),'HandleVisibility','off')
            else ;
            end
        end
    end
end

%% plot the peaks for all variables per altitude
color_alt1 = jet(4);
for i = 1:length(median_alt_row1)
    figure()
    box on
    hold on
    ylim([0 1000/120])
    xlim([datenum(datetime(n11,'ConvertFrom','datenum')+seconds(window_length_seconds)) datenum(datetime(n22,'ConvertFrom','datenum')-seconds(window_length_seconds))])
    datetick('x',13,'keeplimits')
    title(strcat(['Altitude: ',alt_str(i,:),'; ',datestr(NaN_data.Time_datetime(1),'yyyy-mm-dd')]))
    for j = 1:length(parameters_of_interest)
        plot(n11,0,'*','color',color_alt1(parameters_of_interest(j)-1,:))
    end
    legend(variables{parameters_of_interest})%2:5})
    for j = parameters_of_interest%2:5
        for k = 1:size(P_Fregion{i,j},2)-5
            if isempty(locs{i,j,k}) == 0
                plot(Time1(k),locs{i,j,k}*1000,'*','color',color_alt1(j-1,:),'HandleVisibility','off')
            else ;
            end
        end
    end
end

%% try to find at least 
%% get magnitude and phase FFT
% Sa = abs(S_Fregion);
% phi = angle(S_Fregion);
% PDS = 10*log10(P_Fregion);
% 
% [r, c] = find(Sa >= 70);
% Fr = F_Fregion(r);
% Tc = T_Fregion(c)';
% FT = [Tc  Fr];
% [C, ia, ic] = unique(FT(:,1));                                              % Find Unique Times
% for k1 = 1:size(C,1)                                                        % Create Cell Array By Time
%     FrqTime{k1} = FT(FT(:,1) == C(k1),:);                                   % Time & Frequency Cell
% end
% original_f = [697  770  852  941  1209  1336  1477];                        % DTMF Frequencies
% dtmf_dcd = [1 5; 1 6; 1 7; 2 5; 2 6; 2 7; 3 5; 3 6; 3 7; 4 5; 4 6; 4 7];    % Combination Codes w.r.t. ‘original_f’
% nbr_map = ['1' '2' '3' '4' '5' '6' '7' '8' '9' '*' '0' '#'];                % Number Key Map
% for k1 = 1:size(C,1)
%     freq_dist = abs(bsxfun(@minus, FrqTime{k1}(:,2), original_f));          % Distance Of ‘FrqTime’ Frequencies From ‘original_f’ Frequencies
%     [~,freq_pos(:,k1)] = min(freq_dist,[],2);                               % Frequency Positions Of ‘FrqTime’ In ‘original_f’
%     num_pad(k1) = nbr_map(ismember(dtmf_dcd, freq_pos(:,k1)', 'rows'));     % Map To Number Key Pad
% end

%%
% figure()
% h = imagesc(TT_Fregion,F_Fregion,A);%phi);%datenum(TT_Fregion_datetime),F_Fregion,10*log10(P_Fregion));
% xlabel('Time')
% ylabel('Phase')
% title(strcat(['Phase F-region; ',datestr(Time_datetime_new(1),'yyyy-mm-dd')]))
% colorbar;
% set(gca,'YDir','normal');
% set(get(colorbar,'label'),'FontSize', 12,'string','(dB/Hz)')
% xtcks = xticks;
% [~,pos_xtcks] = intersect(TT_Fregion,xtcks); %to find out which values correspond to the xtick locations
% set(gca, 'XTick',xticks, 'XTickLabel', xlabels(pos_xtcks)); %to put datetime ticks on the x-axis

%% spectrogram of interpolated values per individual altitude
% low_lim = 220;
% high_lim = 280;
% keep('E_region_lim','F_region_lim','te_interp','S1','F1','T1','P1','low_lim','high_lim','ne1','Time1','Time_datetime1','te1','alt1','m','delta_T','N_1','N_Te','windowOverlap') %to make sure that all variables underneath are only calculated for the correct time period
% 
% alt_lim = low_lim:5:high_lim;
% for i = 1:size(alt_lim,2)-1 %looking at different altitudes
%     alt_lim_mid(i) = (alt_lim(i)+alt_lim(i+1))/2;
% end
% 
% T = round(seconds(median(delta_T))); %Time period over which EISCAT data is averaged
% Fs = 1/T;  
% f = linspace(0,Fs/2,500); %number of bins in y axis (number of frequency bins)
% for i = 2:length(alt_lim)-1
%     test_Te{i} = te_interp(alt_interp>alt_lim(i) & alt_interp<alt_lim(i+1));
%     if size(test_Te{i},1) > size(te_interp,2)/2 %Make sure that only the altitudes with a signal are used
%         N_Te(i) = length(test_Te{i});
%         if length(test_Te{i}) > 3600/T
%             N_1(i) = 60; %to have a window length of 60 data points (\approx 1 hour)
%         else %time period is shorter than 1 hour
%             N_1(i) = length(test_Te{i});
%         end
%         windowOverlap(i) = N_1(i) - 10;
%         [S,F,TT,P] = spectrogram(test_Te{i},hanning(N_1(i)),windowOverlap(i),f,Fs,'yaxis');              
%         S1{i} = S;
%         F1{i} = F;
%         T1{i} = TT;
%         P1{i} = P;
% 
%         figure()
%         lab = 'Time (secs)';
%         if TT(end) > 60*60
%            TT = TT/(60*60);
%            lab = 'Time (hours)';
%         elseif TT(end) > 60
%            TT = TT/60;
%            lab = 'Time (mins)';
%         end
%         h = imagesc(TT,F,10*log10(P));%.*F);%,[-156.5 95.8] );
%         xlabel(lab)
%         ylabel('Frequency (Hz)')
%         title(strcat(['Alt ',num2str(alt_lim(i)),'-',num2str(alt_lim(i+1)),' km; Time ',datestr(Time_datetime_new(1),'HH:MM:SS'),'-',datestr(Time_datetime_new(end),'HH:MM:SS')]))
%         colorbar;
%         set(gca,'YDir','normal');
%         set(get(colorbar,'label'),'FontSize', 12,'string','(dB/Hz)')
%     else ;
%     end
% end
% 
% %% Split into time periods (problem: would rather interpolate over a data gap of say < 5 minutes --> now not done)
% idx_NaN_end = find(diff(isnan(Time_new(1,:)))==1);
% idx_NaN_end1 = [idx_NaN_end length(Time_new)];
% idx_NaN_start = find(diff(isnan(Time_new(1,:)))==-1)+1;
% idx_NaN_start1 = [1 idx_NaN_start];
% 
% 
% for i = 1:length(idx_NaN_start1)
%     alt1{i} = alt_new(:,idx_NaN_start1(i):idx_NaN_end1(i));
%     ne1{i} = ne_new(:,idx_NaN_start1(i):idx_NaN_end1(i));
%     te1{i} = te_new(:,idx_NaN_start1(i):idx_NaN_end1(i));
%     Time1{i} = Time_new(:,idx_NaN_start1(i):idx_NaN_end1(i));
%     Time_datetime1{i} = Time_datetime_new(:,idx_NaN_start1(i):idx_NaN_end1(i));
% end
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

%% Plot complete time period for geomagnetic latitude 

% [lat,lon,glon,glat] = convert_hdf_to_lat(filename);
% 
% n1 = Time(1,1); 
% n2 = Time(2,end); 
% % n1 = datenum(datetime(2017,12,18,02,00,00));
% % n2 = datenum(datetime(2017,12,18,07,00,00));
% xLimits=[n1 n2];
% % yLimits=[70 72]; %[min(glat,[],'all'), max(glat,[],'all')]; %[min(alt,[],'all'); 300];
% yLimits = [min(glat,[],'all'), max(glat,[],'all')];
% f1=figure();
% colormap jet;
% 
% subplot(3,1,1);
% pcolor(Time(1,:)',glat,log10(ne)),shading flat;
% % imagesc(Time',alt,log10(ne)) %does not work, because needs a vector for Y instead of a matrix (not sure how to do this with a changing altitude) --> this function should make it easier to have empty parts for missing data.
% % surface(Time(1,:)',alt,log10(ne),'EdgeColor','none');%,shading flat; %different way to plot --> works
% 
% caxis([10 11.4]); % give the limits of the colourbar
% colorbar; % Add colorbar to the right of the plot
% set(get(colorbar,'label'),'FontSize', 8,'string','Electron density (m^{-3})')
% xlim(xLimits); 
% ylim(yLimits);
% datetick('x',13,'keeplimits')
% ylabel('Geomagnetic Latitude','FontSize', 8,'FontName','Arial') %Add labels to axes and colorbar
% xlabel(['Time [UT] ',num2str(y(1)),'-',num2str(m(1)),'-',num2str(d(1))],'FontSize', 8,'FontName','Arial')
% 
% % subplot(3,1,1);
% % pcolor(Time(1,:)',glat,vi),shading flat;
% % % imagesc(Time',alt,log10(ne)) %does not work, because needs a vector for Y instead of a matrix (not sure how to do this with a changing altitude) --> this function should make it easier to have empty parts for missing data.
% % % surface(Time(1,:)',alt,log10(ne),'EdgeColor','none');%,shading flat; %different way to plot --> works
% % 
% % caxis([0 500]); % give the limits of the colourbar
% % colorbar; % Add colorbar to the right of the plot
% % set(get(colorbar,'label'),'FontSize', 8,'string','Electron density (m^{-3})')
% % xlim(xLimits); 
% % ylim(yLimits);
% % datetick('x',13,'keeplimits')
% % ylabel('Geomagnetic Latitude','FontSize', 8,'FontName','Arial') %Add labels to axes and colorbar
% % xlabel(['Time [UT] ',num2str(y(1)),'-',num2str(m(1)),'-',num2str(d(1))],'FontSize', 8,'FontName','Arial')
%  
% subplot(3,1,2)
% pcolor(Time(1,:)',glat,te),shading flat;
% caxis([0 3000]); % Limits to colorbar
% colorbar; % Add colorbar to the right of the plot
% set(get(colorbar,'label'),'FontSize', 8,'string','Electron temperature (K)') % Add labels to axes and colorbar
% xlim(xLimits); 
% ylim(yLimits); %ylim([min(alt,[],'all') 300])
% datetick('x',13,'keeplimits')
% ylabel('Geomagnetic Latitude','FontSize', 8,'FontName','Arial') % Add labels to axes and colorbar
% xlabel(['Time [UT] ',num2str(y(1))],'FontSize', 8,'FontName','Arial') % Add labels to axes and colorbar
% 
% subplot(3,1,3)
% pcolor(Time(1,:)',glat,ti),shading flat;
% caxis([0 1500]); % Limits to colorbar
% colorbar; % Add colorbar to the right of the plot
% set(get(colorbar,'label'),'FontSize', 8,'string','Ion temperature (K)')
% xlim(xLimits);
% ylim(yLimits) ;
% datetick('x',13,'keeplimits')
% ylabel('Geomagnetic Latitude','FontSize', 8,'FontName','Arial') % Add labels to axes and colorbar
% xlabel(['Time [UT] ',num2str(y(1))],'FontSize', 8,'FontName','Arial')

%% Plot complete time period for geographic latitude 

% n1 = Time(1,1); 
% n2 = Time(2,end); 
% % n1 = datenum(datetime(2017,12,18,02,00,00));
% % n2 = datenum(datetime(2017,12,18,07,00,00));
% xLimits=[n11 n22];
% % yLimits=[70 72]; %[min(glat,[],'all'), max(glat,[],'all')]; %[min(alt,[],'all'); 300];
% yLimits = [min(lat,[],'all'), max(lat,[],'all')];
% yLimits = [71 75];
% f1=figure();
% colormap jet;
% 
% subplot(3,1,1);
% pcolor(Time(1,:)',lat,log10(ne)),shading flat;
% % imagesc(Time',alt,log10(ne)) %does not work, because needs a vector for Y instead of a matrix (not sure how to do this with a changing altitude) --> this function should make it easier to have empty parts for missing data.
% % surface(Time(1,:)',alt,log10(ne),'EdgeColor','none');%,shading flat; %different way to plot --> works
% 
% caxis([10 11.4]); % give the limits of the colourbar
% colorbar; % Add colorbar to the right of the plot
% set(get(colorbar,'label'),'FontSize', 8,'string','Electron density (m^{-3})')
% xlim(xLimits); 
% ylim(yLimits);
% datetick('x',13,'keeplimits')
% ylabel('Geographic Latitude','FontSize', 8,'FontName','Arial') %Add labels to axes and colorbar
% xlabel(['Time [UT] ',num2str(y(1)),'-',num2str(m(1)),'-',num2str(d(1))],'FontSize', 8,'FontName','Arial')
% 
% subplot(3,1,2)
% pcolor(Time(1,:)',lat,te),shading flat;
% caxis([0 3000]); % Limits to colorbar
% colorbar; % Add colorbar to the right of the plot
% set(get(colorbar,'label'),'FontSize', 8,'string','Electron temperature (K)') % Add labels to axes and colorbar
% xlim(xLimits); 
% ylim(yLimits); %ylim([min(alt,[],'all') 300])
% datetick('x',13,'keeplimits')
% ylabel('Geographic Latitude','FontSize', 8,'FontName','Arial') % Add labels to axes and colorbar
% xlabel(['Time [UT] ',num2str(y(1))],'FontSize', 8,'FontName','Arial') % Add labels to axes and colorbar
% 
% subplot(3,1,3)
% pcolor(Time(1,:)',lat,ti),shading flat;
% caxis([0 1500]); % Limits to colorbar
% colorbar; % Add colorbar to the right of the plot
% set(get(colorbar,'label'),'FontSize', 8,'string','Ion temperature (K)')
% xlim(xLimits);
% ylim(yLimits) ;
% datetick('x',13,'keeplimits')
% ylabel('Geographic Latitude','FontSize', 8,'FontName','Arial') % Add labels to axes and colorbar
% xlabel(['Time [UT] ',num2str(y(1))],'FontSize', 8,'FontName','Arial')
