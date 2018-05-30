function [ am ] = TX_AM_nbNoise( am_frq,duration,cf,bandwidth,percentage,env,savefile )
%[ am ] = TX_AM_signal( am_frq,duration,cf,bandwidth )
% example
% Whitenoise [ am ] = TX_AM_nbNoise( 2,0.5,0 )

%   Detailed explanation goes here
fs = 44100;

if  ~exist('env')
    env = sin(2 * pi * am_frq/2 * (1/fs:1/fs:duration/1000));
    
    env = abs(env) + 1 - percentage;
end

if cf == 0
    carrier = randn(1,fs * duration/1000);
elseif bandwidth == 0
    carrier =  sin(2 * pi * cf * (1/fs:1/fs:duration/1000));
else
    carrier = TX_narrowband_noise(duration,fs,20,bandwidth );
end

am = carrier .* env;

am = am * 0.1;

if exist('savefile')
    audiowrite(['AM_' num2str(am_frq) '_NbNoise_' num2str(cf) 'Hz' num2str(bandwidth) 'Hz_' num2str(duration) 'ms.wav'], am,fs );
end


end

