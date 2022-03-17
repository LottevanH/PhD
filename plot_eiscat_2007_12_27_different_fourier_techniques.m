%% pwelch function or plomb (in case of <5 min time interval)
% 
% low_lim = 200;
% high_lim = 220;
% 
% num_segments = 1;
% for m = 1:size(alt_new1,2) %m represents the number of time periods without gaps in them
% %     keep('low_lim','high_lim','ne_new1','Time_new1','te_new1','alt_new1','m','delta_T','num_segments') %to make sure that all variables underneath are only calculated for the correct time period
%     alt_lim = low_lim:5:high_lim;
%     for i = 1:size(alt_lim,2)-1
%         alt_lim_mid(i) = (alt_lim(i)+alt_lim(i+1))/2;
%     end
%     Te_lim1 = te_new1{m};
%     Ne_lim1 = ne_new1{m};
%     Time_vector = Time_new1{m};
% 
%     T = round(seconds(median(delta_T))); %Time period over which EISCAT data is averaged
%     Fs = 1/T; 
%     
%     for i = 1:size(alt_lim,2)-1
%         
%         test_Te{i} = Te_lim1(alt_new1{m}>alt_lim(i) & alt_new1{m}<alt_lim(i+1));
%         test_Ne{i} = Ne_lim1(alt_new1{m}>alt_lim(i) & alt_new1{m}<alt_lim(i+1));
%         if length(test_Te{i}) > length(te_new1{m})/2
%         [row,col] = find(alt_new1{m}>alt_lim(i) & alt_new1{m}<alt_lim(i+1));
%         test_Te_plomb{i} = Te_lim1(median(row),:);
%         test_Ne_plomb{i} = Ne_lim1(median(row),:);
%             h1 = hann(round(length(test_Ne{i})/num_segments)); %window used
%             if isempty(find(isnan(Time_vector(1,:)) == 1)) == 0 %meaning there is a (small) gap in the data
%                 figure()
%                 suptitle(['FFT from ' datestr(datetime(Time_vector(1,1),'ConvertFrom','datenum')) ' to ' datestr(datetime(Time_vector(1,end),'ConvertFrom','datenum'))])
%                 subplot(2,1,1)
%                 plomb(test_Te_plomb{i},datetime(Time_vector(1,:),'ConvertFrom','datenum'))%,h1,[],round(length(test_Ne{i})/(2*num_segments)),Fs)
%                 subplot(2,1,2)
%                 plomb(test_Ne_plomb{i},datetime(Time_vector(1,:),'ConvertFrom','datenum'))%,h1,[],round(length(test_Ne{i})/(2*num_segments)),Fs)
%             else 
%                 test_detrend_Te{i} = detrend(test_Te{i},1); %to remove the linear trend that shows up as a huge peak around 0 Hz  detrending only works if no NANs in array
%                 test_detrend_Ne{i} = detrend(test_Ne{i},1); %to remove the linear trend that shows up as a huge peak around 0 Hz 
%                 figure()
%                 suptitle(['FFT from ' datestr(datetime(Time_vector(1,1),'ConvertFrom','datenum')) ' to ' datestr(datetime(Time_vector(1,end),'ConvertFrom','datenum'))])
%                 subplot(2,1,1)
%                 pwelch(test_detrend_Te{i},h1,[],round(length(test_Ne{i})/(2*num_segments)),Fs)
%                 subplot(2,1,2)
%                 pwelch(test_detrend_Ne{i},h1,[],round(length(test_Ne{i})/(2*num_segments)),Fs)
%             end
%         else ;
%         end
%     end
% end


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
%         if length(test_Te{i}) > size(te_new1{m},2)/2 %Make sure that only the altitudes with a signal are used
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
%         if length(test_Te{i}) > size(te_new1{m},2)/2 %isempty(test_Te{i}) == 0
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
%         if length(test_Ne{i}) > size(ne_new1{m},2)/2
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