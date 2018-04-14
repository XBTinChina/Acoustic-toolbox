function [sound_2] = TX_generate_2tone( dur,freq,ramp_dur,interval )
%[target] = TX_generate_tone( dur,freq,ramp_dur )

fs = 44100;

sound_2 = zeros(1,dur/1000 * fs *2 + interval /1000 * fs);




target = sin(2*pi*freq*(1/fs:(1/fs):dur/1000));
am1 = sin(2*pi*250/ramp_dur*(1/fs:(1/fs):ramp_dur/1000));
am2 = cos(2*pi*250/ramp_dur*(1/fs:(1/fs):ramp_dur/1000));
target = target(1:length(am1)).*am1;
target = target(end-length(am2)+1:end).*am2;
target = target/sqrt(sum(target.^2)/length(target));

sound_2(1:length(target)) = target;
sound_2(end - length(target) + 1 : end) = target;

end

