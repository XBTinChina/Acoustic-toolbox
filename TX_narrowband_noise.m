function [ narrowband ] = TX_narrowband_noise(dur,fs,ramp_dur,bandwidth,savefile )
%  TX_narrowband_noise(duration,cf,bandwidth,ramp_dur )
%


broadband = randn(1,fs*dur/1000);

% if strcmp(bandwidth,'auditory')
%     bandwidth = 24.7*(4.37*cf/1000 + 1)
% end


w1 = bandwidth(1);
w2 = bandwidth(2) ;


stopband = 200;

[n,Wn] = buttord([w1 w2] * 2/fs,[w1-w1/10 w2+w2/10] * 2/fs,3,10)

[b,a]=butter(n,Wn);

narrowband = filter(b,a,broadband);


am1 = sin(2*pi*250/ramp_dur*(1/fs:(1/fs):ramp_dur/1000));
am2 = cos(2*pi*250/ramp_dur*(1/fs:(1/fs):ramp_dur/1000));

narrowband(1:length(am1)) = narrowband(1:length(am1)).*am1;
narrowband(end-length(am2)+1:end) = narrowband(end-length(am2)+1:end).*am2;


narrowband = narrowband/sqrt(sum(narrowband.^2)/length(narrowband)) * 0.1;





end

