

function snr=TX_SNR(Si,No)

% by TXB
% 3/23/2012
% Si :original signal
% No:noisy signal(ie. original signal + noise signal)
% snr=10*log10(sigma2(Si)/sigma2(No-Si))
% the length of signal must be the same as noise

snr=0;

Ps=sum((Si-mean(Si)).^2)/length(Si);%signal power
Pn=sum((No-Si).^2)/length(No);%noise power
snr=10*log10(Ps/Pn);





end






