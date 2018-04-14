function [ wideband ] = TX_wideband_noise(dur,ramp_dur,savefile )
%  TX_narrowband_noise(duration,cf,bandwidth,ramp_dur )
%

fs = 16000;
wideband = randn(1,fs*dur/1000);



am1 = sin(2*pi*250/ramp_dur*(1/fs:(1/fs):ramp_dur/1000));
am2 = cos(2*pi*250/ramp_dur*(1/fs:(1/fs):ramp_dur/1000));

wideband(1:length(am1)) = wideband(1:length(am1)).*am1;
wideband(end-length(am2)+1:end) = wideband(end-length(am2)+1:end).*am2;


wideband = wideband/sqrt(sum(wideband.^2)/length(wideband)) * 0.1;

if exist('savefile')
    audiowrite(['wideband_'  num2str(dur) 'ms.wav'], wideband,fs );
end



end

