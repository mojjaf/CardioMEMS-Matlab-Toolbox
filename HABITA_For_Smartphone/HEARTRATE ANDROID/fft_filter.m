function [ filtered_signal ] = fft_filter( signal, sampling_freq, high_pass, low_pass, options )
%FFT_FILTER: Does filtering with fft
%  
%   Inputs:
%       signal (will convert it to double)
%       sampling_freq -in Hz
%       high_pass -cut off frequency (Hz) for high pass
%       low_pass -cut off frequency (Hz) for low pass
%
    plot_flag = 0;
    
    Fs = sampling_freq;
    len = length(signal);
    remove_dc = 1;
    
    if exist ('options')
        if strfind(options,'keep dc')
            remove_dc = 0;
        end
    end
        
    
    filter_start = round((high_pass*len)/Fs) +1;
    filter_stop  = round((low_pass*len)/Fs) +1;
    
    signal_fft = fft(double(signal));
    
    signal_fft_filt=zeros(len,1);
    signal_fft_filt(filter_start:filter_stop) = signal_fft(filter_start:filter_stop);
    if remove_dc
        signal_fft_filt(1) = 0;
    end
    filtered_signal = real(ifft(signal_fft_filt))*2;
    
    if plot_flag
        freq_array = 0:Fs/len:Fs-(Fs/len);
        figure, plot(signal), title('orginal signal')
        figure, plot(freq_array,abs(signal_fft)), title('fft')
        figure, plot(filtered_signal), title('filtered signal')
    end
end

