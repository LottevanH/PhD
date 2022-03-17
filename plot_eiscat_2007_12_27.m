clear all
close all

addpath('C:\Github\PhD\functions');

%%
[Time,par2D,par1D,rpar2D,err2D,errs2d]=load_param_hdf5('C:\Users\lotte\OneDrive - Universitetssenteret på Svalbard AS (1)\Svalbard\PhD\Data\Paper_Lisa_2017\EISCAT_2007-12-27_ipy_60@42m.hdf5');
% [Time,par2D,par1D,rpar2D,err2D,errs2d]=load_param_hdf5('C:\Svalbard\PhD\Data\Test_data\EISCAT_2019-11-20_folke_64@42mb.hdf5');
% [Time,par2D,par1D,rpar2D,err2D,errs2d]=load_param_hdf5('C:\Users\lotte\OneDrive - Universitetssenteret på Svalbard AS (1)\Svalbard\PhD\Data\ULF_waves_events_currently_unexamined\1998_12_20\EISCAT_1998-12-20_gup3_aclp_60@32m.hdf5');

[y,m,d,h,mn,s]=datevec(Time(1,:));
ut_time=h+mn/60.;
alt=par2D(:,:,2);
ne=par2D(:,:,3);
ne(find(ne<0)) = NaN; %in case of negative electron densities (physically impossible..)
te=par2D(:,:,4).*par2D(:,:,5);
ti=par2D(:,:,5);
alt_min = min(alt,[],'all');
Time_datetime = datetime(Time,'ConvertFrom','datenum'); %does something weird for 

%% Plot complete time period

n1 = Time(1,1); 
n2 = Time(2,end); 
xLimits=[n1 n2];
yLimits=[alt_min 300];
f1=figure();
colormap jet;

subplot(3,1,1);
% pcolor(datenum([y(:),m(:),d(:),h(:),mn(:),s(:)]),alt,log10(ne)),shading flat;
pcolor(Time(1,:)',alt,log10(ne)),shading flat;
% imagesc(Time',alt,log10(ne))
% surface(Time(1,:)',alt,log10(ne),'EdgeColor','none');%,shading flat;

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
 
%% Time period of interest
n11 = datenum(datetime(2007,12,27,03,00,00));
n22 = datenum(datetime(2007,12,27,06,00,00));
% n11 = datenum(datetime(2007,12,27,15,00,00));
% n22 = datenum(datetime(2007,12,27,17,30,00));
% n11 = n1;
% n22 = n2;

idx_time_start = find(Time(1,:) > n11,1);
idx_time_end = find(Time(1,:) < n22,1,'last');
 
alt_new = alt(:,idx_time_start:idx_time_end);
ne_new = ne(:,idx_time_start:idx_time_end);
te_new = te(:,idx_time_start:idx_time_end);
ti_new = ti(:,idx_time_start:idx_time_end);
Time_new = Time(:,idx_time_start:idx_time_end);
Time_datetime_new = Time_datetime(:,idx_time_start:idx_time_end);

%% Plot time period of interest
xLimits=[n11 n22];
yLimits=[alt_min 300];
f1=figure();
colormap jet;

subplot(3,1,1);
pcolor(Time_new(1,:)',alt_new,log10(ne_new)),shading flat;
% surface(Time(1,:)',alt,log10(ne),'EdgeColor','none');%,shading flat;

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

%% Split in time periods if there are time gaps
max_duration = 1.5; %number of hours of which the FFT window can maximally exist
%FFT requires evenly-spaced data without data gaps
for i = 1:size(Time_datetime_new,2)-1
delta_T(i) = Time_datetime_new(1,i+1) - Time_datetime_new(1,i);
end
idx_time_gap = find(delta_T > median(delta_T)+seconds(1));
if isempty(idx_time_gap) == 0
    time_duration(1) = minutes(Time_new(1,idx_time_gap(1)) - Time_new(1,1))*60*24;
    % time_duration1(1) = Time_datetime_new(1,idx_time_gap(1)) - Time_datetime_new(1,1);
    for i = 2:size(idx_time_gap,2)
        time_duration(i) = minutes(Time_new(1,idx_time_gap(i)) - Time_new(1,idx_time_gap(i-1)+1))*60*24;
    %     time_duration1(i) = Time_datetime_new(1,idx_time_gap(i)) - Time_datetime_new(1,idx_time_gap(i-1)+1);
    end
    time_duration(size(time_duration,2)+1) = minutes(Time_new(1,end) - Time_new(1,idx_time_gap(end)+1))*60*24;
    time_duration.Format = 'hh:mm:ss';

% introduce extra time gaps if the time span between gaps exceeds a certain
% period (for example one hour)
for i = 1:size(time_duration,2)
    if time_duration(i) > hours(max_duration)%hours(1)
        number_extra_gaps(i) = floor(hours(time_duration(i)));
    else ;
    end
end
        
idx_time_gap1 = zeros(1,sum(number_extra_gaps(:))+size(idx_time_gap,2));
i = 1;
    for j = 1:number_extra_gaps(i)
        idx_time_gap1(j) = find(Time_new(1,:) > Time_new(1,1) + 1/24*j,1);
    end
idx_time_gap1(number_extra_gaps(i)+i) = idx_time_gap(i);
for i = 2:size(time_duration,2)-1
    if number_extra_gaps(i) == 0
        idx_time_gap1(sum(number_extra_gaps(1:i))+i) = idx_time_gap(i);
    else
        for j = 1:number_extra_gaps(i)
            idx_time_gap1(j + i - 1 + sum(number_extra_gaps(1:i-1))) = find(Time_new(1,:) > Time_new(1,idx_time_gap(i-1)) + 1/24*j,1);
            idx_time_gap1(sum(number_extra_gaps(1:i))+i) = idx_time_gap(i);
        end
    end
end
i = size(time_duration,2);
    if number_extra_gaps(i) == 0
        ;
    else
        for j = 1:number_extra_gaps(i)
            idx_time_gap1(j + i - 1 + sum(number_extra_gaps(1:i-1))) = find(Time_new(1,:) > Time_new(1,idx_time_gap(i-1)) + 1/24*j,1);
        end
    end
    
else 
    % introduce extra time gaps if the time span between gaps exceeds a certain
    % period (for example one hour)
    time_duration(1) = minutes(Time_new(1,end) - Time_new(1,1))*60*24;
    for i = 1:size(time_duration,2)
        if time_duration(i) > hours(max_duration)%hours(1)
            number_extra_gaps(i) = floor(hours(time_duration(i)));
        else ;
        end
    end
        for j = 1:number_extra_gaps(i)
            idx_time_gap1(j) = find(Time_new(1,:) > Time_new(1,1) + 1/24*j,1);
        end
end


% time_duration1(size(time_duration1,2)+1) = Time_datetime_new(1,end) - Time_datetime_new(1,idx_time_gap(end)+1);

if isempty(idx_time_gap1) == 0
    alt_new1{1} = alt_new(:,1:idx_time_gap1(1)-1); %-1 to get rid of the last data point before the break, since that one is usually shorter
    te_new1{1} = te_new(:,1:idx_time_gap1(1)-1);
    Time_new1{1} = Time_new(:,1:idx_time_gap1(1)-1);
    ne_new1{1} = ne_new(:,1:idx_time_gap1(1)-1);
    if size(idx_time_gap1,2) > 1
    for i = 2:size(idx_time_gap1,2)
        alt_new1{i} = alt_new(:,idx_time_gap1(i-1)+1:idx_time_gap1(i)-1);
        te_new1{i} = te_new(:,idx_time_gap1(i-1)+1:idx_time_gap1(i)-1);
        Time_new1{i} = Time_new(:,idx_time_gap1(i-1)+1:idx_time_gap1(i)-1);
        ne_new1{i} = ne_new(:,idx_time_gap1(i-1)+1:idx_time_gap1(i)-1);
    end
    alt_new1{i+1} = alt_new(:,idx_time_gap1(end)+1:end);
    te_new1{i+1} = te_new(:,idx_time_gap1(end)+1:end);
    Time_new1{i+1} = Time_new(:,idx_time_gap1(end)+1:end);
    ne_new1{i+1} = ne_new(:,idx_time_gap1(end)+1:end);
    elseif size(idx_time_gap1,2) == 1
        alt_new1{2} = alt_new(:,idx_time_gap1(1)+1:end);
        te_new1{2} = te_new(:,idx_time_gap1(1)+1:end);
        Time_new1{2} = Time_new(:,idx_time_gap1(1)+1:end);
        ne_new1{2} = ne_new(:,idx_time_gap1(1)+1:end);
    end
else
    alt_new1{1} = alt_new;
end

%make sure that time intervals are max 1 (or 0.5?) hour


%% FFT to find any frequency peaks

% low_lim = 200;
% high_lim = 220;
% for m = 1:size(alt_new1,2) %m represents the number of time periods without gaps in them
%     keep('low_lim','high_lim','ne_new1','Time_new1','te_new1','alt_new1','m','delta_T') %to make sure that all variables underneath are only calculated for the correct time period
%     alt_lim = low_lim:5:high_lim;
%     for i = 1:size(alt_lim,2)-1
%         alt_lim_mid(i) = (alt_lim(i)+alt_lim(i+1))/2;
%     end
%     Te_lim1 = te_new1{m};
%     Ne_lim1 = ne_new1{m};
%     Time_vector = Time_new1{m};
% 
%     %divide per hour?
%     delta_T1 = round(3600/64);
%     delta_T_array = 1:delta_T1:size(te_new1{m},2);
% 
%     % T = 1/Fs;
%     T = round(seconds(median(delta_T))); %Time period over which EISCAT data is averaged
%     Fs = 1/T; 
%   
%     for i = 1:size(alt_lim,2)-1
%         test_Te{i} = Te_lim1(alt_new1{m}>alt_lim(i) & alt_new1{m}<alt_lim(i+1));
%         test_Ne{i} = Ne_lim1(alt_new1{m}>alt_lim(i) & alt_new1{m}<alt_lim(i+1));
%         
%         if size(test_Te{i},1) > size(te_new1{m},2)/2 %Make sure that only the altitudes with a signal are used
%             test_detrend_Te{i} = detrend(test_Te{i},1); %to remove the linear trend that shows up as a huge peak around 0 Hz 
%             test_detrend_Ne{i} = detrend(test_Ne{i},1); %to remove the linear trend that shows up as a huge peak around 0 Hz 
%             
%             %Fast Fourier Transform for electron temperature
%             NFFT_Te(i) = 2^nextpow2(size(test_Te{i},1)); %next power of 2 from length of Te (in order to improve the performance of fft)
%             X_Hann_Te{i} = test_detrend_Te{i}.*hanning(size(test_Te{i},1));%hamming(size(test_Te{i},1));
%             FFT1_Te{i} = fft(X_Hann_Te{i},NFFT_Te(i)); %%FFT transform
% %             FFT1_Te{i} = fft(test_detrend_Te{i},NFFT_Te(i)); %%FFT transform
%             P2_Te{i} = abs(FFT1_Te{i}/NFFT_Te(i)); %two sided spectrum (amplitude)
%             P1_Te{i} = P2_Te{i}(1:NFFT_Te(i)/2+1); %single sided spectrum
%             P1_Te_a = P1_Te{i};
%             P1_Te_a(2:end-1) = 2*P1_Te_a(2:end-1); %multiply since spectrum is symmetric and combining both sides into one side to get the original amplitude
%             P1_Te_1{i} = P1_Te_a;
%             f_Te{i} = Fs*(0:(NFFT_Te(i)/2))/NFFT_Te(i); %frequency domain f
%             
%             % Fast Fourier Transform for electron density
%             NFFT_Ne(i) = 2^nextpow2(size(test_Ne{i},1)); %next power of 2 from length of Te (in order to improve the performance of fft)
%             X_Hann_Ne{i} = test_detrend_Ne{i}.*hann(size(test_Ne{i},1));%hamming(size(test_Ne{i},1));
%             FFT1_Ne{i} = fft(X_Hann_Ne{i},NFFT_Ne(i)); %%FFT transform
%             P2_Ne{i} = abs(FFT1_Ne{i}/NFFT_Ne(i)); %two sided spectrum (amplitude)
%             P1_Ne{i} = P2_Ne{i}(1:NFFT_Ne(i)/2+1); %single sided spectrum
%             P1_Ne_a = P1_Ne{i};
%             P1_Ne_a(2:end-1) = 2*P1_Ne_a(2:end-1); %multiply since spectrum is symmetric and combining both sides into one side to get the original amplitude
%             P1_Ne_1{i} = P1_Ne_a;
%             f_Ne{i} = Fs*(0:(NFFT_Ne(i)/2))/NFFT_Ne(i); %frequency domain f
%             
%             % Plot the original and Fourier transformed signal
% %             figure()
% %             subplot(2,2,1)
% %             plot(datetime(Time_vector(1,:),'ConvertFrom','datenum'),test_Te{i},datetime(Time_vector(1,:),'ConvertFrom','datenum'),test_detrend_Te{i},datetime(Time_vector(1,:),'ConvertFrom','datenum'),test_Te{i}-test_detrend_Te{i})
% %             title('Original signal of Te, detrend version and trend')
% %             legend('original','detrended version','trend')
% %             xlabel('Time')
% %             ylabel('Te (K)')
% %             
% %             subplot(2,2,2)
% %             plot(f_Te{i},P1_Te_a)
% %             title('Single-Sided Amplitude Spectrum of Te(t)')
% %             xlim([-inf inf])
% %             xlabel('f (Hz)')
% %             ylabel('|P1(f)|')
% %             
% %             subplot(2,2,3)
% %             plot(datetime(Time_vector(1,:),'ConvertFrom','datenum'),test_Ne{i},datetime(Time_vector(1,:),'ConvertFrom','datenum'),test_detrend_Ne{i},datetime(Time_vector(1,:),'ConvertFrom','datenum'),test_Ne{i}-test_detrend_Ne{i})
% %             title('Original signal of Ne, detrend version and trend')
% %             legend('original','detrended version','trend')
% %             xlabel('Time')
% %             ylabel('Te (K)')
% % 
% %             subplot(2,2,4)
% %             plot(f_Ne{i},P1_Ne_a)
% %             title('Single-Sided Amplitude Spectrum of Ne(t)')
% %             xlim([-inf inf])
% %             xlabel('f (Hz)')
% %             ylabel('|P1(f)|')
%             
%             clear P1_Te_a P1_Ne_a
%         else ;
%         end
%     end
%       
%     % Plot altitude FFT signals stacked 
%     figure()
%     suptitle(['FFT from ' datestr(datetime(Time_vector(1,1),'ConvertFrom','datenum')) ' to ' datestr(datetime(Time_vector(1,end),'ConvertFrom','datenum'))])
%     subplot(2,1,1)
%     hold on
%     for i = 1:size(alt_lim,2)-2
%         if isempty(test_Te{i}) == 0
%             plot(f_Te{i},P1_Te_1{i})
%             idx_alt(i) = find(isempty(test_Te{i}) == 0);
%         else ;
%         end 
%     end
% %     title('Single-Sided Amplitude Spectrum of Te(t)')
%     xL=xlim;
%     yL=ylim;
%     text(0.99*xL(2),0.99*yL(2),"Te",'HorizontalAlignment','right','VerticalAlignment','top')
%     xlim([-inf inf])
%     xlabel('f (Hz)')
%     ylabel('|P1(f)|')
%     legend(strsplit(num2str(alt_lim_mid(find(idx_alt(:)==1)))))
%     
%     subplot(2,1,2)
%     hold on
%     for i = 1:size(alt_lim,2)-2
%         if isempty(test_Ne{i}) == 0
%             plot(f_Ne{i},P1_Ne_1{i})
%         else 
%         end;
%     end
% %     title('Single-Sided Amplitude Spectrum of Ne(t)')
%     xL=xlim;
%     yL=ylim;
%     text(0.99*xL(2),0.99*yL(2),"Ne",'HorizontalAlignment','right','VerticalAlignment','top')
%     xlim([-inf inf])
%     xlabel('f (Hz)')
%     ylabel('|P1(f)|')
%     legend(strsplit(num2str(alt_lim_mid(find(idx_alt(:)==1)))))
%     
% %     clear Time_vector
% end

%% Try with spectrogram (to use window overlap)

low_lim = 200;
high_lim = 220;
for m = 1:size(alt_new1,2) %m represents the number of time periods without gaps in them
    keep('low_lim','high_lim','ne_new1','Time_new1','te_new1','alt_new1','m','delta_T') %to make sure that all variables underneath are only calculated for the correct time period
    alt_lim = low_lim:5:high_lim;
    for i = 1:size(alt_lim,2)-1
        alt_lim_mid(i) = (alt_lim(i)+alt_lim(i+1))/2;
    end
    Te_lim1 = te_new1{m};
    Ne_lim1 = ne_new1{m};
    Time_vector = Time_new1{m};

    %divide per hour?
    delta_T1 = round(3600/64);
    delta_T_array = 1:delta_T1:size(te_new1{m},2);

    % T = 1/Fs;
    T = round(seconds(median(delta_T))); %Time period over which EISCAT data is averaged
    Fs = 1/T; 
    windowOverlap = [10, 20, 30]; 
  
    for i = 1:size(alt_lim,2)-1
        test_Te{i} = Te_lim1(alt_new1{m}>alt_lim(i) & alt_new1{m}<alt_lim(i+1));
        test_Ne{i} = Ne_lim1(alt_new1{m}>alt_lim(i) & alt_new1{m}<alt_lim(i+1));
        if size(test_Te{i},1) > size(te_new1{m},2)/2 %Make sure that only the altitudes with a signal are used
            
            N_Te(i) = length(test_Te{i});
            N_1(i) = N_Te(i)/8;%2^nextpow2(N_Te(i)); %same for Te and Ne
%             f{i} = Fs*(0:(N_1(i)/2))/N_1(i); %frequency domain f
            f{i} = linspace(0,Fs/2,100);

            for l = 1:length(windowOverlap)
                figure()
%                 spectrogram(test_Te{i},hanning(N_1(i)),[],f{i},Fs,'yaxis')
                spectrogram(test_Te{i},hanning(N_1(i)),[],[],[],'yaxis')
                title(['window overlap = ' num2str(l)])
            end
        else ;
        end
    end
end
