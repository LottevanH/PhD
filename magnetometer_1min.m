clear all
close all

filename_csv = 'C:\Users\lotte\OneDrive - Universitetssenteret på Svalbard AS (1)\Svalbard\PhD\Data\ULF_waves_events_currently_unexamined\2017_12_18\magnetometer\20220512-13-47-supermag.csv';

%% csv 1 minute data
mag_data = readtable(filename_csv);

%% Sort the data by station
stations = unique(mag_data.IAGA);
for i = 1:length(stations)
    myStruct.(stations{i}) = mag_data(ismember(mag_data.IAGA,stations{i}),:); %create a struct in which the magnetometer data is saved per station.
    myStruct.(stations{i}).Date_UTC = datetime(myStruct.(stations{i}).Date_UTC,'InputFormat','yyyy-MM-dd''T''HH:mm:ss'); %make sure the time is in the datetime format
end

%% High pass filter (subtract median value of sliding window of 10 minutes lenght for Pc5)
for i = 1:length(stations)
    for j = 5:length(myStruct.(stations{i}).dbn_geo)-5
        myStruct.(stations{i}).filtered_dbn_geo(j) = myStruct.(stations{i}).dbn_geo(j) - median(myStruct.(stations{i}).dbn_geo(j-4:j+5),'omitnan');
    end
end

%% Bandpass filter between 1 and 8 mHz
delta_T = diff(myStruct.TRO.Date_UTC);
T = round(seconds(median(delta_T)));
Fs = 1/T;
for i = 1:length(stations)
    myStruct.(stations{i}).bandpass_dbn_geo = bandpass(myStruct.(stations{i}).dbn_geo,[1e-3 8e-3],Fs);
end

%% Time period of interest
n1 = datetime(2017,12,18,02,00,00);
n2 = datetime(2017,12,18,07,00,00);
n11 = find(ismember(myStruct.(stations{i}).Date_UTC,n1)==1);
n22 = find(ismember(myStruct.(stations{i}).Date_UTC,n2)==1);

for i = 1:length(stations)
    poi.(stations{i}).Date_UTC = myStruct.(stations{i}).Date_UTC(n11:n22);
    poi.(stations{i}).dbn_geo = myStruct.(stations{i}).dbn_geo(n11:n22);
    poi.(stations{i}).filtered_dbn_geo = myStruct.(stations{i}).filtered_dbn_geo(n11:n22);
    poi.(stations{i}).bandpass_dbn_geo = myStruct.(stations{i}).filtered_dbn_geo(n11:n22);
end

%% Plot original and filtered signals
figure()
plot(poi.TRO.Date_UTC,poi.TRO.dbn_geo)
hold on
plot(poi.TRO.Date_UTC,poi.TRO.filtered_dbn_geo)
plot(poi.TRO.Date_UTC,poi.TRO.bandpass_dbn_geo)
legend('original','high pass filter','bandpass filter')

%% FFT
f = linspace(0,Fs/2,500);
windowLength = 60;
windowOverlap = 50;
for i = 1:length(stations) %[5 10 11 12]%
%     [S{i},F{i},TT{i},P{i}] = spectrogram(poi.(stations{i}).filtered_dbn_geo,hamming(windowLength),windowOverlap,f,Fs,'yaxis');
    [S{i},F{i},TT{i},P{i}] = spectrogram(poi.(stations{i}).bandpass_dbn_geo,parzenwin(windowLength),windowOverlap,f,Fs,'yaxis');
    powerspectrum{i} = 10*log10(P{i});
    TT_datetime{i} = poi.(stations{i}).Date_UTC(1) + seconds(TT{i});
    figure()
    h = imagesc(datenum(TT_datetime{i}),F{i}*1000,powerspectrum{i}); xlabel('Time')
    ylabel('Frequency (mHz)')
    ylim([0 1000/120])
    caxis([0 60])
    datetick('x',13,'keeplimits')
    set(gca,'YDir','normal');
    set(get(colorbar,'label'),'FontSize', 12,'string','(dB/Hz)')
    title(stations{i})
end
%% plot multiple frequency spectra on top of each other
for j = 15%1:size(powerspectrum{1},2)
    figure()
    hold on
    for i = 1:length(stations)%strmatch('BJN',stations)%[5 10 11 12]%1:length(stations) %
        plot(F{i}*1000,powerspectrum{1,i}(:,j))
    end
    legend(stations(2))%stations{5}, stations{10}, stations{11}, stations{12})
    xlabel('Frequency [mHz]')
    ylabel('Power [dB]')
    xlim([0 max(F{12}*1000)])
    title(strcat('Time: ',datestr(TT_datetime{5}(j))))
end




