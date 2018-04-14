function [ am ] = TX_AM_WbNoise( am_frq,duration,percentage, fs,envelope,ref_value,savefile, num )
%[ am ] = TX_AM_signal( am_frq,duration,cf,bandwidth )
% example
% Whitenoise [ am ] = TX_AM_nbNoise( 2,0.5,0 )

%   Detailed explanation goes here

%  'warning: the sampling rate is'
%  fs
if ~exist('envelope','var')
    
     envelope = sin(2 * pi * am_frq * (1/fs:1/fs:duration/1000) -pi/2);
    
     envelope = envelope * percentage + 1 ;
    
end

carrier = randn(1,fs * duration/1000);
carrier = carrier - mean(carrier);


if ~exist('ref_value','var')
    carrier = carrier / std(carrier) * ref_value;
end
    
    
am = carrier .*  envelope;

am = am * 0.1;

if exist('savefile')
    audiowrite(['AM_' num2str(am_frq) '_WbNoise_'  'Hz_' num2str(duration) 'ms.wav'], am,fs );
end


end

