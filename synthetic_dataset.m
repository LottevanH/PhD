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

noise=5*rand(1, length(t)); %creates a noise array equal of the length of t where the noise amplitude (1, 2, 3, 4 or 5) is multiplied by a random number between 0 and 1
phase_shift = 2*pi*rand(4,1);
y0 = sin(2*pi*f1*t+phase_shift(1)) + sin(2*pi*f2*t+phase_shift(2)) + sin(2*pi*f3*t+phase_shift(3)) + sin(2*pi*f4*t+phase_shift(4)); %no noise
y1 = sin(2*pi*f1*t+phase_shift(1)) + sin(2*pi*f2*t+phase_shift(2)) + sin(2*pi*f3*t+phase_shift(3)) + sin(2*pi*f4*t+phase_shift(4)) + noise;

% 
% figure(1)
% plot(t/3600,y0,'k')
% % hold on
% % plot(t/3600,y1,'r')
% % hold on
% % plot(t/3600,y2,'k')
% xlabel('Time [hr]')
% ylabel('Amplitude')
% title('Synthetic input')

%% optional: add NaNs
dur_time_gap = 20;%randi([1,20]);
loc_time_gap = randi([1,length(t)-dur_time_gap]);
y1_time_gap = y1;
y1_time_gap(loc_time_gap:loc_time_gap+dur_time_gap) = NaN;
y0_time_gap = y0;
y0_time_gap(loc_time_gap:loc_time_gap+dur_time_gap) = NaN;

%% interpolate over NaNs
X = y1_time_gap;
X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)), 'linear'); %linear interpolation
y1_interp = X;
Y = y0_time_gap;
Y(isnan(Y)) = interp1(find(~isnan(Y)), Y(~isnan(Y)), find(isnan(Y)), 'linear'); %linear interpolation
y0_interp = Y;            
clear X Y

% figure(2)
% plot(t/3600,y1_interp,'r')
% hold on
% plot(t/3600,y2_interp,'k')
% xlabel('Time [hr]')
% ylabel('Amplitude')
% title('Synthetic input')

%% to do
%spectrogram of complete synthetic dataset
%spectrogram of interpolated synthetic dataset
%lomb-scargle of complete synthetic dataset
%lomb-scargle of synthetic dataset with time gaps

%% spectrogram complete synthetic dataset
f = linspace(1e-3,7e-3,50);%(0,Fs/2,500); %number of bins in y axis (number of frequency bins)
window_length_seconds = 3600; 
window_overlap_seconds = 600;
N_1 = round(window_length_seconds/dt);
windowOverlap = N_1 - round(window_overlap_seconds/dt);
[S,F,T,P] = spectrogram(y1,hamming(N_1),windowOverlap,f,Fs,'yaxis');   
% [S,F,T,P] = spectrogram(y1,rectwin(N_1),windowOverlap,f,Fs,'yaxis');   %promising windows: rectwin, bartlett, gausswin, kaiser, tukeywin(N_1,0.25) (or 0.5)
idx_freq_pc5 = find(F> 0.6e-3 & F<7e-3);
powerspectrum = 10*log10(P(idx_freq_pc5,:));
max_power = max(powerspectrum,[],'all');
min_power = min(powerspectrum,[],'all');

fig = figure();
h = imagesc(T/3600,F(idx_freq_pc5,:)*1000,powerspectrum);%10*log10(P_Fregion{i,j}));%datenum(TT_Fregion_datetime),F_Fregion,10*log10(P_Fregion));
xlabel('Time [hr]')
ylabel('Frequency (mHz)')
ylim([0 1000/120])
title('Synthetic dataset no time gap')
colorbar;
colormap parula
caxis([max_power-20 max_power])
set(gca,'YDir','normal');
set(get(colorbar,'label'),'FontSize', 12,'string','(dB/Hz)')

%% spectrogram of interpolated synthetic dataset
[S_interp,F_interp,T_interp,P_interp] = spectrogram(y1_interp,hamming(N_1),windowOverlap,f,Fs,'yaxis');   
% [S,F,T,P] = spectrogram(y1,rectwin(N_1),windowOverlap,f,Fs,'yaxis');   %promising windows: rectwin, bartlett, gausswin, kaiser, tukeywin(N_1,0.25) (or 0.5)
idx_freq_pc5 = find(F_interp> 0.6e-3 & F_interp<7e-3);
powerspectrum_interp = 10*log10(P_interp(idx_freq_pc5,:));
max_power = max(powerspectrum_interp,[],'all');
min_power = min(powerspectrum_interp,[],'all');

fig = figure();
h = imagesc(T_interp/3600,F_interp(idx_freq_pc5,:)*1000,powerspectrum_interp);%10*log10(P_Fregion{i,j}));%datenum(TT_Fregion_datetime),F_Fregion,10*log10(P_Fregion));
xlabel('Time [hr]')
ylabel('Frequency (mHz)')
ylim([0 1000/120])
title(strcat('Synthetic dataset interpolated; Duration time gap: ',num2str(dur_time_gap)))
colorbar;
colormap parula
caxis([max_power-50 max_power])
set(gca,'YDir','normal');
set(get(colorbar,'label'),'FontSize', 12,'string','(dB/Hz)')

%% Lomb-Scargle dataset without datagaps
%need to calculate plomb of a sliding window 
% note: with plomb, you can calculate pth (power-level threshold). A preak
% with a value larger than pth has a probability pdvec of being a true
% signal peak and not the result of random fluctuations (pth is not
% available if fvec is used)
df = 1e-4;
fvec = 1e-3:df:7e-3;
[pxx,f] = plomb(y1,t,fvec); %plomb of the entire synthetic dataset

% figure()
% line(f,pxx)

%now: evolving plomb with a certain window overlap and window length.
w = N_1;
pxx_evol = zeros(length(fvec),length(t));
f_evol   = zeros(length(fvec),length(t));

for i = w/2+1 : length(t)-w/2 
    [pxx_evol(:,i),f_evol(:,i)] = plomb(y1(i-w/2:i+w/2),t(i-w/2:i+w/2),fvec);
end
F_EVOL = f_evol(:,i);
powerspectrum_ls = 10*log10(pxx_evol);
max_power = max(powerspectrum_ls,[],'all');
min_power = min(powerspectrum_ls,[],'all');

figure()
imagesc(t/3600,F_EVOL*1000,powerspectrum_ls)
xlabel('Time [hr]')
ylabel('Frequency (mHz)')
% ylim([0 1000/120])
title(strcat('Synthetic dataset Lomb Scargle; No time gap'))
colorbar;
colormap parula
caxis([max_power-20 max_power])
hold on
xline(t(loc_time_gap)/3600,'r','LineWidth',1)
xline(t(loc_time_gap+dur_time_gap)/3600,'r','LineWidth',1)
set(gca,'YDir','normal');
set(get(colorbar,'label'),'FontSize', 12,'string','(dB/Hz)')

%% Lomb-Scargle with datagap
w = N_1;
pxx_evol = zeros(length(fvec),length(t));
f_evol   = zeros(length(fvec),length(t));

for i = w/2+1 : length(t)-w/2 
    [pxx_evol(:,i),f_evol(:,i)] = plomb(y1_time_gap(i-w/2:i+w/2),t(i-w/2:i+w/2),fvec);
end
F_EVOL = f_evol(:,i);
powerspectrum_ls = 10*log10(pxx_evol);
max_power = max(powerspectrum_ls,[],'all');
min_power = min(powerspectrum_ls,[],'all');

figure()
imagesc(t/3600,F_EVOL*1000,powerspectrum_ls)
xlabel('Time [hr]')
ylabel('Frequency (mHz)')
% ylim([0 1000/120])
title(strcat('Synthetic dataset Lomb Scargle; Duration time gap: ',num2str(dur_time_gap)))
colorbar;
colormap parula
caxis([max_power-20 max_power])
hold on
xline(t(loc_time_gap)/3600,'r','LineWidth',1)
xline(t(loc_time_gap+dur_time_gap)/3600,'r','LineWidth',1)
if loc_time_gap-w/2 > 1
    xline(t(loc_time_gap-w/2)/3600,'m','LineWidth',1)
else 
end
if loc_time_gap+dur_time_gap+w/2 < length(t)
    xline(t(loc_time_gap+dur_time_gap+w/2)/3600,'m','LineWidth',1)
else 
end
set(gca,'YDir','normal');
set(get(colorbar,'label'),'FontSize', 12,'string','(dB/Hz)')
