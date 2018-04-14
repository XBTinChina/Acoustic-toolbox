function [ signal ] = TX_AM(carrior_freq,modulater_freq,phase,depth,dur,fs )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

seq = 0: 1/fs : dur / 1000;

carrier = sin(2 * pi * carrior_freq * seq );
modulator = 0.5 + depth * 0.5 * sin(2 * pi * modulater_freq * (seq + phase /(2 * pi) * 1/modulater_freq));

signal = carrier .* modulator;

end

