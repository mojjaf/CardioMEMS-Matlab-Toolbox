function [ filtered_signal ] = fft_filter( signal, sampling_freq, high_pass, low_pass )
%FFT_FILTER Summary of this function goes here
%   Detailed explanation goes here
%   [ filtered_signal ] = fft_filter( signal, sampling_freq, high_pass,  low_pass )
    plot_flag = 0;
    
    Fs = sampling_freq;
    len = length(signal);
    
    
    
    filter_start = round((high_pass*len)/Fs) +1;
    filter_stop  = round((low_pass*len)/Fs) +1;
    
    signal_fft = fft(signal);
    
    freq_array = 0:Fs/len:Fs-(Fs/len);
            
    signal_fft_filt=zeros(len,1);
    signal_fft_filt(filter_start:filter_stop) = signal_fft(filter_start:filter_stop);
    filtered_signal = real(ifft(signal_fft_filt))*2;
    
    if plot_flag
        figure, plot(signal), title('orginal signal')
        figure, plot(freq_array,abs(signal_fft)), title('fft')
        figure, plot(filtered_signal), title('filtered signal')
    end
end

