function [ signal_new ] = TX_power_refer( signal,refer,snr )
%   Correct the power of signal compared to the referrence;
%   [ signal_new ] = TX_power_refer( signal,refer );

gain = 10^(snr/20);

signal_new = signal/std(signal)*max(std(refer))*gain;

signal_new = signal_new';

end

