clear all
close all

filename_txt_10sec = 'C:\Users\lotte\OneDrive - Universitetssenteret på Svalbard AS (1)\Svalbard\PhD\Data\Test_data\Pilipenko2014\image_20031031.txt';
% filename_txt_10sec = 'C:\Users\lotte\OneDrive - Universitetssenteret på Svalbard AS (1)\Svalbard\PhD\Data\ULF_waves_events_currently_unexamined\2017_12_18\magnetometer\image_20171218.txt\image_20171218.txt';
% filename_Lisa2017 = 'C:\Users\lotte\OneDrive - Universitetssenteret på Svalbard AS (1)\Svalbard\PhD\Data\Paper_Lisa_2017\magnetometer\image_20071227.txt';
%% txt
% fid=fopen(filename_txt_10sec);
% cdata=textscan(fid,'%s','delimiter','\n', 'HeaderLines', 1 );
% fclose(fid);
fid = fopen(filename_txt_10sec);
% fid = fopen(filename_Lisa2017);
line1 = fgets(fid);
fclose(fid);
header = split(line1)';
data = readtable(filename_txt_10sec);
% data = readtable(filename_Lisa2017);
time = datetime([data{:,1} data{:,2} data{:,3} data{:,4} data{:,5} data{:,6}],'InputFormat','yyyy MM dd HH mm ss');
mag_data.time = time;
idx_stations = find(cellfun('size',header,2) == 3);
station = header(idx_stations);
stations = station(1:3:end);
for i = 1:length(stations)
    mag_data.(stations{i}).X = table2array(data(:,4+3*i));
    mag_data.(stations{i}).Y = table2array(data(:,5+3*i));
    mag_data.(stations{i}).Z = table2array(data(:,6+3*i));
end

% X, Y and Z are geographic north, east and downward components. The unit
% is nT. Missing values are 99999.9

%% High pass filter (subtract median value of sliding window of 10 minutes lenght for Pc5)
% for i = 1:length(stations)
%     for j = 5:length(mag_data.(stations{i}).dbn_geo)-5
%         mag_data.(stations{i}).filtered_dbn_geo(j) = mag_data.(stations{i}).dbn_geo(j) - median(mag_data.(stations{i}).dbn_geo(j-4:j+5),'omitnan');
%     end
% end

%% Bandpass filter between 1 and 10 mHz
delta_T = diff(mag_data.time);
T = round(seconds(median(delta_T)));
Fs = 1/T;
for i = 1:length(stations)
    mag_data.(stations{i}).detrend_X = detrend(mag_data.(stations{i}).X,1);
    mag_data.(stations{i}).detrend_Y = detrend(mag_data.(stations{i}).Y,1);
    mag_data.(stations{i}).detrend_Z = detrend(mag_data.(stations{i}).Z,1);
    mag_data.(stations{i}).bandpass_X = bandpass(mag_data.(stations{i}).detrend_X,[1e-3 10e-3],Fs);
    mag_data.(stations{i}).bandpass_Y = bandpass(mag_data.(stations{i}).detrend_Y,[1e-3 10e-3],Fs);
    mag_data.(stations{i}).bandpass_Z = bandpass(mag_data.(stations{i}).detrend_Z,[1e-3 10e-3],Fs);
end

%% Time period of interest
% n1 = datetime(2017,12,18,02,00,00);
% n2 = datetime(2017,12,18,07,00,00);
% n1 = datetime(2007,12,27,15,00,00);
% n2 = datetime(2007,12,27,17,30,00);
n1 = datetime(2003,10,31,11,00,00); %in case of limited time
n2 = datetime(2003,10,31,12,00,00); %in case of limited time
n11 = find(ismember(mag_data.time,n1)==1);
n22 = find(ismember(mag_data.time,n2)==1);

poi.time = mag_data.time(n11:n22);
for i = 1:length(stations)
    poi.(stations{i}).X = mag_data.(stations{i}).X(n11:n22);
    poi.(stations{i}).bandpass_X = mag_data.(stations{i}).bandpass_X(n11:n22);
    poi.(stations{i}).Y = mag_data.(stations{i}).Y(n11:n22);
    poi.(stations{i}).bandpass_Y = mag_data.(stations{i}).bandpass_Y(n11:n22);
    poi.(stations{i}).Z = mag_data.(stations{i}).Z(n11:n22);
    poi.(stations{i}).bandpass_Z = mag_data.(stations{i}).bandpass_Z(n11:n22);
end

%% Plot original and filtered signals
figure()
hold on
for i = 1:length(stations)
plot(poi.time,poi.(stations{i}).bandpass_X)
end
legend(stations{:})
title('Bandpass filtered X-component magnetometer')
% legend('original','bandpass filter')

%% Spectrogram
f = linspace(0,Fs/2,181);%500);
windowLength = 3600/T;
windowOverlap = windowLength - 600/T;
L = length(poi.time);
for i = 1:length(stations) %[5 10 11 12]%
%     [S{i},F{i},TT{i},P_x{i}] = spectrogram(poi.(stations{i}).filtered_dbn_geo,hamming(windowLength),windowOverlap,f,Fs,'yaxis');
    [S_x{i},F_x{i},TT_x{i},P_x{i}] = spectrogram(poi.(stations{i}).bandpass_X,hamming(windowLength),windowOverlap,f,Fs,'yaxis');
    [S_y{i},F_y{i},TT_y{i},P_y{i}] = spectrogram(poi.(stations{i}).bandpass_Y,hamming(windowLength),windowOverlap,f,Fs,'yaxis');
    [S_z{i},F_z{i},TT_z{i},P_z{i}] = spectrogram(poi.(stations{i}).bandpass_Z,hamming(windowLength),windowOverlap,f,Fs,'yaxis');
%     figure()
%     spectrogram(poi.(stations{i}).bandpass_X,hamming(windowLength),windowOverlap,f,Fs,'yaxis')
%     ylim([0 1000/120])
%     caxis([0 60])
%     [S_Y{i},F_Y{i},TT_Y{i},P_Y{i}] = spectrogram(poi.(stations{i}).bandpass_Y,hamming(windowLength),windowOverlap,f,Fs,'yaxis');
    powerspectrum_x{i} = 10*log10(P_x{i});
    powerspectrum_y{i} = 10*log10(P_y{i});
    powerspectrum_z{i} = 10*log10(P_z{i});
%     f_power{i} = (0:L-1)*(Fs/L);
%     power{i} = abs(fft_X{i}).^2/L;
    TT_datetime{i} = poi.time(1) + seconds(TT_x{i});
%     powerspectrum_Y{i} = 10*log10(P_Y{i});
%     TT_datetime_Y{i} = poi.time(1) + seconds(TT_Y{i});

    %find peaks
    for k = 5:size(powerspectrum_x{i},2)-5
        [pks_x{i,k},locs_x{i,k}] = findpeaks(powerspectrum_x{i}(:,k), F_x{i},'MinPeakHeight',max(powerspectrum_x{i}(:,k-4:k+4),[],'all')-10,'MinPeakProminence',7); %IDEA: ONLY USE THE MAXIMUM POWER WITHIN A 3 HOUR INTERVAL AROUND THE POINT K
        [pks_y{i,k},locs_y{i,k}] = findpeaks(powerspectrum_y{i}(:,k), F_y{i},'MinPeakHeight',max(powerspectrum_y{i}(:,k-4:k+4),[],'all')-10,'MinPeakProminence',7); %IDEA: ONLY USE THE MAXIMUM POWER WITHIN A 3 HOUR INTERVAL AROUND THE POINT K
        [pks_z{i,k},locs_z{i,k}] = findpeaks(powerspectrum_z{i}(:,k), F_z{i},'MinPeakHeight',max(powerspectrum_z{i}(:,k-4:k+4),[],'all')-10,'MinPeakProminence',7); %IDEA: ONLY USE THE MAXIMUM POWER WITHIN A 3 HOUR INTERVAL AROUND THE POINT K  
    end
    Time1 = datenum(TT_datetime{i});
    
    figure()
    h = imagesc(datenum(TT_datetime{i}),F_x{i}*1000,powerspectrum_x{i}); 
    xlabel('Time')
    ylabel('Frequency (mHz)')
    ylim([0 1000/120])
%     caxis([0 60])
    datetick('x',13,'keeplimits')
    set(gca,'YDir','normal');
    set(get(colorbar,'label'),'FontSize', 12,'string','(dB/Hz)')
    title(strcat(stations{i},'; X-component'))
    

    hold on 
    for k = 5:size(powerspectrum_x{i},2)-5 
        if isempty(locs_x{i,k}) == 0
            plot(Time1(k),locs_x{i,k}*1000,'.r')
        else ;
        end
    end
    
    figure()
    h = imagesc(datenum(TT_datetime{i}),F_y{i}*1000,powerspectrum_y{i}); 
    xlabel('Time')
    ylabel('Frequency (mHz)')
    ylim([0 1000/120])
%     caxis([0 60])
    datetick('x',13,'keeplimits')
    set(gca,'YDir','normal');
    set(get(colorbar,'label'),'FontSize', 12,'string','(dB/Hz)')
    title(strcat(stations{i},'; Y-component'))
    

    hold on 
    for k = 5:size(powerspectrum_y{i},2)-5 
        if isempty(locs_y{i,k}) == 0
            plot(Time1(k),locs_y{i,k}*1000,'.r')
        else ;
        end
    end

    figure()
    h = imagesc(datenum(TT_datetime{i}),F_z{i}*1000,powerspectrum_z{i}); 
    xlabel('Time')
    ylabel('Frequency (mHz)')
    ylim([0 1000/120])
%     caxis([0 60])
    datetick('x',13,'keeplimits')
    set(gca,'YDir','normal');
    set(get(colorbar,'label'),'FontSize', 12,'string','(dB/Hz)')
    title(strcat(stations{i},'; Z-component'))
    

    hold on 
    for k = 5:size(powerspectrum_z{i},2)-5 
        if isempty(locs_z{i,k}) == 0
            plot(Time1(k),locs_z{i,k}*1000,'.r')
        else ;
        end
    end
    clear Time1
    
%     figure()
%     imagesc(datenum(TT_datetime{i}),F{i}*1000,20*log10(abs(S{i}(:,:,1))))
%     set(gca,'YDir','normal');
%     datetick('x',13,'keeplimits')
%     ylim([0 1000/120])
%     set(get(colorbar,'label'),'FontSize', 12,'string','Magnitude [dB]')
%     caxis([0 70])
%     title(stations{i})

end

%% find peaks

%% FFT

for i = 1:length(stations)
    fft_X{i} = fft(poi.(stations{i}).bandpass_X);
    fft_P2{i} = abs(fft_X{i}/L);
    fft_P1{i} = fft_P2{i}(1:L/2+1);
    fft_P1{i}(2:end-1) = 2*fft_P1{i}(2:end-1);
    fft_f{i} = Fs*(0:(L/2))/L;
    fft_Y{i} = fft(poi.(stations{i}).bandpass_Y);
    fft_P2_Y{i} = abs(fft_Y{i}/L);
    fft_P1_Y{i} = fft_P2_Y{i}(1:L/2+1);
    fft_P1_Y{i}(2:end-1) = 2*fft_P1_Y{i}(2:end-1);
    fft_f_Y{i} = Fs*(0:(L/2))/L;
    
    figure()
    plot(fft_f{i},fft_P1{i})
    title(strcat(stations{i},'; Single-Sided Amplitude Spectrum'))
    hold on
    plot(fft_f_Y{i},fft_P1_Y{i})
%     hold on
%     plot(f_fft_X1{i},mx{i})
% %     plot(f_power{i},power{i})
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
    legend('X','Y')
    xlim([0 5e-3])
end

%% STFT

for i = 1:length(stations)
    figure()
    stft(poi.(stations{i}).bandpass_X,Fs,'Window',hamming(windowLength),'OverlapLength',windowOverlap,'FrequencyRange','onesided');
    ylim([0 1000/120])
    caxis([0 70])
    %for all channels (1 = X, 2 = Y, 3 = Z)
    [S_STFT{i},F_STFT{i},T_STFT{i}] = stft([poi.(stations{i}).bandpass_X'; poi.(stations{i}).bandpass_Y'; poi.(stations{i}).bandpass_Z']',Fs,'Window',hamming(windowLength),'OverlapLength',windowOverlap,'FrequencyRange','onesided');
    figure()
    imagesc(T_STFT{i}/3600,F_STFT{i}*1000,20*log10(abs(S_STFT{i}(:,:,1))))
    set(gca,'YDir','normal');
    ylim([0 1000/120])
    set(get(colorbar,'label'),'FontSize', 12,'string','Magnitude [dB]')
    caxis([0 70])
end

%% plot multiple frequency spectra on top of each other
% for j = 15%1:size(powerspectrum{1},2)
%     figure()
%     hold on
%     for i = 1:length(stations)%strmatch('BJN',stations)%[5 10 11 12]%1:length(stations) %
%         plot(F{i}*1000,powerspectrum{1,i}(:,j))
%     end
%     legend(stations(2))%stations{5}, stations{10}, stations{11}, stations{12})
%     xlabel('Frequency [mHz]')
%     ylabel('Power [dB]')
%     xlim([0 8.3])
%     title(strcat('Time: ',datestr(TT_datetime{5}(j))))
% end



