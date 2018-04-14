function [target] = TX_generate_tone( dur,freq,fs,ramp_dur,random_phase,savefile)
%[target] = TX_generate_tone( dur,freq,ramp_dur )


target = zeros(1,dur/1000 * fs );

if exist('random_phase')
    rphase = random_phase * pi;
else
    rphase = 0;
end


target = sin(2*pi*freq*(1/fs:(1/fs):dur/1000) + rphase);


if ~isempty(ramp_dur)
    
    am1 = sin(2*pi*250/ramp_dur*(1/fs:(1/fs):ramp_dur/1000));
    am2 = cos(2*pi*250/ramp_dur*(1/fs:(1/fs):ramp_dur/1000));
    
    
    target(1:length(am1)) = target(1:length(am1)).*am1;
    target(end-length(am2)+1:end) = target(end-length(am2)+1:end).*am2;
    
end

target = target/sqrt(sum(target.^2)/length(target)) * 0.1;

if exist('savefile')
    audiowrite(['tone_' num2str(freq) 'Hz_' num2str(dur) 'ms.wav'], target,fs );
end


end

