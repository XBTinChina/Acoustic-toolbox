function [ waveform ] = TX_dip_train(click_dur,interval,num_click,AMP, fs)
%  Xiangbin Teng, 1.1.2018
%  Generate a click train with the amplitude of pulses equal to AMP
%  Sampling rate: fs
%  The click train ends with a click

display(fs)

interval = ceil( interval * fs / 1000 );
click_dur = ceil( click_dur * fs / 1000 );

seg_dur = interval;

waveform = ones(seg_dur,num_click);

waveform(1:click_dur,:) = AMP;

waveform =  waveform(:);


end
