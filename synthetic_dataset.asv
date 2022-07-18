clear all
close all

%% generate synthetic dataset
Fs = 1/60; %sampling frequency
dt = 1/Fs; %seconds per sample
StopTime = 3*60*60; %seconds
t = (0:dt:StopTime);
%chosen frequencies
f1=2e-3; %4 different frequencies f1, f2, f3, f4
f2=3.3e-3;
f3=5e-3;
f4=6.2e-3;

%random frequencies in the Pc5 range
f11 = 1e-3 + 6e-3*rand(4,1); %random number between 1e-3 and 7e-3

noise_amp = randi(5)*50; %random number between 1 and 5 for noise amplitude, multiplied by 50 to get a realistic amplitude for EISCAT data
noise=noise_amp* rand(1, length(t)); %creates a noise array equal of the length of t where the noise amplitude (1, 2, 3, 4 or 5) is multiplied by a random number between 0 and 1
y1 = sin(f1*t) + sin(f2*t) + sin(f3*t) + sin(f4*t) + noise;
y2 = sum(sin(f11.*t)) + noise;

figure(1)
plot(t/3600,y1,'r')
hold on
plot(t/3600,y2,'k')
xlabel('Time [hr]')
ylabel('Amplitude')
title('Synthetic input')

%% optional: add NaNs
dur_time_gap = randi([1,20]);
loc_time_gap = randi([1,length(t)-dur_time_gap]);
y1_time_gap = y1;
y1_time_gap(loc_time_gap:loc_time_gap+dur_time_gap) = NaN;
y2_time_gap = y2;
y2_time_gap(loc_time_gap:loc_time_gap+dur_time_gap) = NaN;

%% interpolate over NaNs
X = y1_time_gap;
X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)), 'linear'); %linear interpolation
y1_interp = X;
Y = y2_time_gap;
Y(isnan(Y)) = interp1(find(~isnan(Y)), Y(~isnan(Y)), find(isnan(Y)), 'linear'); %linear interpolation
y2_interp = Y;            
clear X Y

figure(2)
plot(t/3600,y1_interp,'r')
hold on
plot(t/3600,y2_interp,'k')
xlabel('Time [hr]')
ylabel('Amplitude')
title('Synthetic input')

%% to do
%spectrogram of complete synthetic dataset
%spectrogram of interpolated synthetic dataset
%lomb-scargle of complete synthetic dataset
%lomb-scargle of synthetic dataset with time gaps

%% spectrogram complete synthetic dataset
f = linspace(0,Fs/2,500); %number of bins in y axis (number of frequency bins)
window_length_seconds = 3600; 
window_overlap_seconds = 600;
N_1 = round(window_length_seconds/dt);
[S,F,T,P] = spectrogram(y1),hamming(N_1),windowOverlap,f,Fs,'yaxis'); %Calculate for F-region   
idx_freq_pc5 = find(F> 0.6e-3 & F<7e-3);
powerspectrum = 10*log10(P(idx_freq_pc5,:));
max_power(i,j) = max(powerspectrum{i,j},[],'all');
min_power(i,j) = min(powerspectrum{i,j},[],'all');
