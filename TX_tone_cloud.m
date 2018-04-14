function [ tone_cloud ] = TX_tone_cloud( band_number,frequency_range,tone_duration,stimulus_length,ramp_dur,fs )
% Xiangbin Teng, 1.1.2018
% Generate tone cloud without temporal and spatial overlap


center_frequency = ERBSpace(frequency_range(1),frequency_range(2),band_number);


number_of_tones = round(stimulus_length / tone_duration);



temp = randperm(number_of_tones);
tone_sequence = mod(temp,band_number) + 1;


tone_cloud = [];
for t = 1:length(tone_sequence)
    
    [target] = TX_generate_tone( tone_duration,center_frequency(tone_sequence(t)),fs,ramp_dur,1);
    
    tone_cloud = [tone_cloud target];
    
end






end

