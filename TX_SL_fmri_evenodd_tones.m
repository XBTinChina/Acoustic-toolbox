function [ SOUNDWAVE ] = TX_SL_fmri_evenodd_tones(fs, frequencies,amplitude_gain,even_period,odd_period,seg_ramp_time,intervals,total_time,total_ramp_time,band_width,random_phase,amplitude_random_range)
% Xiangbin Teng, 14-April-2018
% version controled - github

% [ soundwave ] = TX_SL_evenodd_tones(fs, frequency_component,frequencies,amplitude_gain,even_period,odd_period,seg_ramp_time,intervals,total_time,total_ramp_time,random_phase)

%   Detailed explanation goes here
% fs = 44100
% frequencies = [225, 602, 1611, 4311, 11535]
% amplitude_gain = ones(1,length(frequencies)) .* 1
% even_period = 80 ; % in ms
% odd_period = 80 ; % in ms
% seg_ramp_time = 20; % in ms
% intervals = 20 ; %in ms
% total_time = 8000; %in ms
% total_ramp_time = 1000; %in ms
% band_width = 1/3  ; % 0 = sinusoid;  if exist and > 0, create narrowband noise 1/3 = 0.3 octave
% random_phase = 'yes';
% amplitude_random_range  =  2 ; % in dB



num_of_even_odd_cycles = floor(total_time/ (even_period + odd_period + 2 * intervals)) ;



% assign frequenices to even and odd

even_frequency =  frequencies(1:2:end);
odd_frequency =  frequencies(2:2:end);

even_gain = amplitude_gain(1:2:end);
odd_gain = amplitude_gain(2:2:end);

% generage ramping cosine functions

seg_ramp_freq = 1000 / (seg_ramp_time * 4);

period_ramp_up = sin(2*pi*seg_ramp_freq*(1/fs:(1/fs):seg_ramp_time/1000));
period_ramp_down = cos(2*pi*seg_ramp_freq*(1/fs:(1/fs):seg_ramp_time/1000));



total_ramp_freq = 1000 / (total_ramp_time * 4);

total_ramp_up = sin(2*pi*total_ramp_freq*(1/fs:(1/fs):total_ramp_time/1000));
total_ramp_down = cos(2*pi*total_ramp_freq*(1/fs:(1/fs):total_ramp_time/1000));



SOUNDWAVE = zeros(1,total_time * fs / 1000);

interval_points = zeros(1,round(fs * intervals / 1000));

for num_cycle = 1:num_of_even_odd_cycles
    
    
    
    
    if exist('random_phase')
        
        even_phase = rand(1,length(even_frequency)) * 2 * pi;
        odd_phase = rand(1,length(odd_frequency)) * 2 * pi;
    else
        even_phase = zeros(1,length(even_frequency)) ;
        odd_phase = zeros(1,length(odd_frequency)) ;
    end
    
    
    %% generate even frequencies
    temp_even_wave = zeros(length(even_frequency),round(fs * even_period / 1000));
    
    for even_f  = 1:length(even_frequency)
        
        if exist('band_width') && band_width > eps
            even_band = [ even_frequency(even_f) * 2^(-band_width/2)  even_frequency(even_f) * 2^(band_width/2) ] ;
            broadband = randn(1,fs*even_period/1000);
            w1 = even_band(1);
            w2 = even_band(2) ;
            [n,Wn] = buttord([w1 w2] * 2/fs,[w1-w1/10 w2+w2/10] * 2/fs,3,10);
            [b,a]=butter(n,Wn);
            temp_even_wave(even_f,:) = filter(b,a,broadband) * even_gain(even_f) ;
            
            
        else
            
            temp_even_wave(even_f,:) = sin(2*pi*even_frequency(even_f)*(1/fs:(1/fs):even_period/1000) + even_phase(even_f)) * even_gain(even_f) ;
            
        end
        
        
    end
    
    even_wave = mean(temp_even_wave,1);
    
    even_wave(1:length(period_ramp_up)) = even_wave(1:length(period_ramp_up)).*period_ramp_up;
    even_wave(end-length(period_ramp_down )+1:end) = even_wave(end-length(period_ramp_down )+1:end).*period_ramp_down ;
    
    even_wave = even_wave ./ max(abs(even_wave)) * 0.5;
    
    if exist('amplitude_random_range')
        even_wave = even_wave * 10^( amplitude_random_range / 20 * (rand - 0.5)/0.5);
    end
    
    %% generate odd frequencies
    
    temp_odd_wave = zeros(length(odd_frequency),round(fs * odd_period / 1000));
    
    for odd_f  = 1:length(odd_frequency)
        
        if exist('band_width') && band_width > eps
            odd_band = [ odd_frequency(odd_f) * 2^(-band_width/2)  odd_frequency(odd_f) * 2^(band_width/2) ] ;
            broadband = randn(1,fs*odd_period/1000);
            w1 = odd_band(1);
            w2 = odd_band(2) ;
            [n,Wn] = buttord([w1 w2] * 2/fs,[w1-w1/10 w2+w2/10] * 2/fs,3,10);
            [b,a]=butter(n,Wn);
            temp_odd_wave(odd_f,:) = filter(b,a,broadband) * even_gain(even_f) ;
            
        else
            
            temp_odd_wave(odd_f,:) = sin(2*pi*odd_frequency(odd_f)*(1/fs:(1/fs):odd_period/1000) + odd_phase(odd_f)) * even_gain(even_f) ;
            
        end
        
        
    end
    
    
    odd_wave = mean(temp_odd_wave,1);
    
    odd_wave(1:length(period_ramp_up)) = odd_wave(1:length(period_ramp_up)).*period_ramp_up;
    odd_wave(end-length(period_ramp_down )+1:end) = odd_wave(end-length(period_ramp_down )+1:end).*period_ramp_down ;
      odd_wave = odd_wave ./ max(abs(odd_wave)) * 0.5;
    
    if exist('amplitude_random_range')
        odd_wave = odd_wave * 10^( amplitude_random_range / 20 * (rand - 0.5)/0.5);
    end
    %%
    
    temp_cycle = [even_wave interval_points odd_wave interval_points];
    
    SOUNDWAVE((num_cycle-1)*length(temp_cycle)+1 : num_cycle*length(temp_cycle)) = temp_cycle;
    
    
    
end

SOUNDWAVE(1:length(total_ramp_up)) = SOUNDWAVE(1:length(total_ramp_up)).*total_ramp_up;
SOUNDWAVE(end-length(total_ramp_down )+1:end) = SOUNDWAVE(end-length(total_ramp_down )+1:end).*total_ramp_down ;

%audiowrite('sound_temp.wav',SOUNDWAVE,fs)