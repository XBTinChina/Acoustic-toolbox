function M = generate_fmam(SampleRate,F1,F2,D1,D2,Bandwidth,NumComps,...
					noise_dur,jitter, SNR)

% generate_fmam(SampleRate,F1,F2,D1,D2,Bandwidth,NumComps,...
%					noise_dur,jitter, SNR)
%    
%  11.26.2012   Xiangbin Teng     email: xiangbin.teng@nyu.edu
%  This program is adapted from 
%        {stimtrain_synth(SampleRate,F1,F2,D1,D2,Bandwidth,...
%              NumComps,ShowComps,ShowMatchFM,ShowMatchTone,...
%              ShowEndPoints,ShowFreqs,AssignFM,WavWrite)
%              ver 01-05-03}
%
% SampleRate:	Sample rate (Hz)
% F1:		    Low-freq (Hz)
% F2:		    High-freq (Hz)
% D1:		    Mean duration of indiv sweep (ms)
% D2:		    Duration of entire stimulus (sec)
% Bandwidth:	"Formant" bandwidth (Hz)
% NumComps:	    Num of freq components (typ. > 20)
% jitter:   onset of the stimuli, input as an array, e.g. [400 600 800] ms
% SNR:  signal to noise ratio. The power of signal is constant; the power
%       of noise is changed. e.g. [-30:2:0] dB
% 
%
%e.g., generate_fmam(44100,1000,1500,132,800,200,50,2000,[400 600 800],[-30:2:0])
%
% Sound level calibration need to be done on stimuli. The power of Noise is then scaled by SNR.  
% 
%
% ------------ User-Stuff------------------------------------


FadeTime=100;  % ms

AmpEnv=1; % apply a ramp-ON and ramp-OFF to each segment
RampPeriod=0.003; % seconds; for ON/OFF AmpEnv

FiltStim=1; % LP-filter stim waveform
FiltOrder=3; % Filter order  HP-filter

NumCycles=2; % the num of cycles looked at to find the smallest zero-crossing point
			 % if this is made too big, it will generate an error for short D1's
			 % since it will produce search lengths larger than the segment length

ChirpCycles=1; % the num of chirp interpolation cycles


LowerBoundFact=0.8; 	% this factor multiplies "D1" and establishes a hard
						% lower-bound to the duration distribution

FudgeFact1=1.5; % used to produce an "Amp"-array large enough to scale "L"
				% due to variability of (random) sweep duration

FudgeFact2=2; 	% used to produce a dur-distr with num of elements greater than
				% NumSweeps given the hard lower-bound threshold

PolyDistrExp=1.5; % exponent for polynomial freq distr
FreqScale=2;      % draws the mean of the freq distr down by this factor--compensates
				  % for the componnets in the heavy, high-freq tail of the polynomial distr

PlotAnalytic=0;


WavPath='~/Desktop';  % where you put your stimuli

% -----------------------------------------------------------

NumArgs=nargin('generate_fmam');
if nargin ~= NumArgs
    help generate_fmam;
	return
end


StartTime=now;
Nyquist=SampleRate/2;
SamplePeriod=1/SampleRate;
NumSweeps=round(D2/D1);

% don't change order of the lines below!
Sigma_Dur= 30 ; % this scales the width of the dur-distr so that it's comparable across D1


SweepDuration=randn(1,FudgeFact2*NumSweeps);
Direction=(SweepDuration>0); 
Direction=Direction(randperm(length(Direction)));
SweepDuration=SweepDuration*Sigma_Dur+D1; % converts std-norm distr to duration-distr.
SweepDurationThresh=SweepDuration > (LowerBoundFact*D1); % establishes lower-bound to dur-distr
SweepDuration=SweepDuration(SweepDurationThresh);

ToneFreqsSeeds=rand(1,NumSweeps); % Controls the freq-component distribution 
ToneFreqs=ToneFreqsSeeds*(F2-F1)+F1; % converts std-norm distr to duration-distr.



x2=0;
y2=0;
amp2=0;

h=waitbar(0, 'Constructing stimuli...');
Counter=0;

if exist(WavPath, 'dir')~=7
	[FILENAME, PATHNAME] = uigetfile('*.*', 'Select Stim Directory');
	cd(PATHNAME);
else
	cd(WavPath);
end



%fprintf('\n\n');
FMStimSegTemp=[];


for n=1:NumSweeps

Counter=Counter+1;
waitbar(Counter/NumSweeps, h);


t=0:SamplePeriod:SweepDuration(n)/1000;


r=rand(size(t));
PolyDistr=(t.^PolyDistrExp);
PolyDistr=(Bandwidth*(PolyDistr./max(PolyDistr)).*r)+1;

if Direction(n)==0 % inverts Freq1 & Freq2 for downward-FMs
	Freq1=F1+(Bandwidth/2)*(randn(1,NumComps)-0.5);
	Freq2=F2+(Bandwidth/2)*(randn(1,NumComps)-0.5);
else
	Freq2=F1+(Bandwidth/2)*(randn(1,NumComps)-0.5);
	Freq1=F2+(Bandwidth/2)*(randn(1,NumComps)-0.5);
end

TF=ToneFreqs(n)+(Bandwidth/2)*(randn(1,NumComps)-0.5);

Amp=rand(1,NumComps);
InitPhase=2*pi*(rand(1,NumComps)-0.5);

for i=1:NumComps
	DeltaFreq=(Freq2(i)-Freq1(i))/(length(t)-1);
	FreqRamp=Freq1(i):DeltaFreq:(Freq2(i));
	Modulation=(2*pi*FreqRamp/2.*t);
   X=Amp(i)*sin(2*pi*Freq1(i)/2*t+InitPhase(i)+Modulation); 	% FM
	Y=Amp(i)*sin(2*pi*TF(i)*t+InitPhase(i));							% Tone (narrow-band)
   
	%Y=Amp*sin(2*pi*ToneFreqs*t+InitPhase);					% Tone (single-freq)
	%fprintf('Freq1=%5.0f  Freq2=%5.0f\n',Freq1(i),Freq2(i));
	
	if i==1
		FMStimSeg=zeros(size(X));
		ToneStimSeg=zeros(size(Y));
	end

	% -------- AM-segment envelope---------

if AmpEnv==1
	tR=linspace(0,pi/2,round(SampleRate*RampPeriod));
	XEnvON=sin(tR).^2;
	XEnvOFF=cos(tR).^2;
	XEnv=[XEnvON,ones(1,length(X)-length(XEnvON)-length(XEnvOFF)),XEnvOFF];
	X=X.*XEnv;

	YEnvON=sin(tR).^2;
	YEnvOFF=cos(tR).^2;
	YEnv=[YEnvON,ones(1,length(Y)-length(YEnvON)-length(YEnvOFF)),YEnvOFF];
	Y=Y.*YEnv;
end

	% -------------------------------------

	FMStimSeg=FMStimSeg+X;
	ToneStimSeg=ToneStimSeg+Y;
end % for i=1:NumComps

% ---- FM Phase Matching ------------------------------------------

if n > 1
	e2=round(SampleRate*NumCycles/(F1+(F2-F1)/2));
	FMStimSegBegin=FMStimSeg(1:e2);

	FM_ZeroCrossBegin=abs(diff(FMStimSegBegin/max(FMStimSegBegin)>0));
	FM_ZeroCrossIndicesBeginTemp=1:length(FM_ZeroCrossBegin);
	FM_ZerosCrossDeltasBegin=FM_ZeroCrossIndicesBeginTemp(FM_ZeroCrossBegin==1);
	FM_ZerosCrossDeltasBeginMin=min(FM_ZerosCrossDeltasBegin);
	FMStimSegBeginShort=FMStimSeg(FM_ZerosCrossDeltasBeginMin:end);
	FMStimSegBeginShort(1)=0;
	%FMStimSegBeginShort(length(FMStimSegBeginShort))=0;
	FMStimSegBeginDirection=FMStimSegBeginShort(2)-FMStimSegBeginShort(1);
	FM_ZeroCrossIFsBegin=1./(2*diff(FM_ZerosCrossDeltasBegin)*SamplePeriod);
	FM_FirstIF=FM_ZeroCrossIFsBegin(length(FM_ZeroCrossIFsBegin));

	if 0==1 % shows location of zero-crossings against FMSeg
		plot(FMStimSegBegin,'ro-');
		hold
		plot(FMStimSegBeginShort,'b-');
		plot((FM_ZeroCrossBegin==1)-1,'k+');
		grid;
		pause;
	end

	tChirp=0:SamplePeriod:(1/((FM_LastIF+FM_FirstIF)/2))*ChirpCycles;
	tChirpInitPhase=90;
	FM_InterSweepChirp=chirp(tChirp,LastIF,tChirp(length(tChirp)),FM_FirstIF,...
		'linear',tChirpInitPhase); % don't use 'quadratic' as it will change the zero-crossing of the chirp!
						
	InterSweepChirp(1)=0;
	InterSweepChirp(length(InterSweepChirp))=0;
													 
	if FMStimSegBeginDirection > 0 % flips FMStimSegBeginShort to match chirp polarity 
		FMStimSegBeginShort=-FMStimSegBeginShort;
	end
	FMStimSegTemp=[FM_InterSweepChirp,FMStimSegBeginShort];

	if FMStimSegEndDirection > 0 % flips FMStimSegTemp to match FMStimSegEndShort polarity 
		FMStimSegTemp=-FMStimSegTemp;
		FMStimSegBeginShortPlot=-FMStimSegBeginShort;
		InterSweepChirpPlot=-FM_InterSweepChirp;
	else
		FMStimSegBeginShortPlot=FMStimSegBeginShort;
		InterSweepChirpPlot=FM_InterSweepChirp;
	end


end % if n > 1

	if n==1
		FMStimSegTemp=FMStimSeg;
	end

	e1=round(length(FMStimSegTemp)-SampleRate*NumCycles/(F1+(F2-F1)/2));
	FMStimSegEnd=FMStimSegTemp(e1:end);

	FM_ZeroCrossEnd=abs(diff(FMStimSegEnd/max(FMStimSegEnd)>0));
	FM_ZeroCrossIndicesEndTemp=1:length(FM_ZeroCrossEnd);
	FM_ZerosCrossDeltasEnd=FM_ZeroCrossIndicesEndTemp(FM_ZeroCrossEnd==1);
	FM_ZerosCrossDeltasEndMax=max(FM_ZerosCrossDeltasEnd);
	FMStimSegEndShort=FMStimSegTemp(1:(FM_ZerosCrossDeltasEndMax+e1));
	FMStimSegEndShort(length(FMStimSegEndShort))=0;
	FMStimSegEndDirection=FMStimSegEndShort(length(FMStimSegEndShort))-FMStimSegEndShort(length(FMStimSegEndShort)-1);
	FM_ZeroCrossIFsEnd=1./(2*diff(FM_ZerosCrossDeltasEnd)*SamplePeriod);
	FM_LastIF=FM_ZeroCrossIFsEnd(length(FM_ZeroCrossIFsEnd));

% --------------------------------------------------------------


% ---- Tone Phase Matching ------------------------------------------

if n > 1
	e2=round(SampleRate*NumCycles/mean(TF));
	ToneStimSegBegin=ToneStimSeg(1:e2);

	ZeroCrossBegin=abs(diff(ToneStimSegBegin/max(ToneStimSegBegin)>0));
	ZeroCrossIndicesBeginTemp=1:length(ZeroCrossBegin);
	ZerosCrossDeltasBegin=ZeroCrossIndicesBeginTemp(ZeroCrossBegin==1);
	ZerosCrossDeltasBeginMin=min(ZerosCrossDeltasBegin);
	ToneStimSegBeginShort=ToneStimSeg(ZerosCrossDeltasBeginMin:end);
	ToneStimSegBeginShort(1)=0;
	%ToneStimSegBeginShort(length(ToneStimSegBeginShort))=0;
	ToneStimSegBeginDirection=ToneStimSegBeginShort(2)-ToneStimSegBeginShort(1);
	ZeroCrossIFsBegin=1./(2*diff(ZerosCrossDeltasBegin)*SamplePeriod);
	FirstIF=ZeroCrossIFsBegin(length(ZeroCrossIFsBegin));

	if 0==1 % shows location of zero-crossings against ToneSeg
		plot(ToneStimSegBegin,'ro-');
		hold
		plot(ToneStimSegBeginShort,'b-');
		plot((ZeroCrossBegin==1)-1,'k+');
		grid;
		pause;
	end

	tChirp=0:SamplePeriod:(1/((LastIF+FirstIF)/2))*ChirpCycles;
	tChirpInitPhase=90;
	InterSweepChirp=chirp(tChirp,LastIF,tChirp(length(tChirp)),FirstIF,...
		'linear',tChirpInitPhase); % don't use 'quadratic' as it will change the zero-crossing of the chirp!
						
	InterSweepChirp(1)=0;
	InterSweepChirp(length(InterSweepChirp))=0;
													 
	if ToneStimSegBeginDirection > 0 % flips ToneStimSegBeginShort to match chirp polarity 
		ToneStimSegBeginShort=-ToneStimSegBeginShort;
	end
	ToneStimSegTemp=[InterSweepChirp,ToneStimSegBeginShort];

	if ToneStimSegEndDirection > 0 % flips ToneStimSegTemp to match ToneStimSegEndShort polarity 
		ToneStimSegTemp=-ToneStimSegTemp;
		ToneStimSegBeginShortPlot=-ToneStimSegBeginShort;
		InterSweepChirpPlot=-InterSweepChirp;
	else
		ToneStimSegBeginShortPlot=ToneStimSegBeginShort;
		InterSweepChirpPlot=InterSweepChirp;
	end

	

end % if n > 1

	if n==1
		ToneStimSegTemp=ToneStimSeg;
	end

	e1=round(length(ToneStimSegTemp)-SampleRate*NumCycles/mean(TF));
	ToneStimSegEnd=ToneStimSegTemp(e1:end);

	ZeroCrossEnd=abs(diff(ToneStimSegEnd/max(ToneStimSegEnd)>0));
	ZeroCrossIndicesEndTemp=1:length(ZeroCrossEnd);
	ZerosCrossDeltasEnd=ZeroCrossIndicesEndTemp(ZeroCrossEnd==1);
	ZerosCrossDeltasEndMax=max(ZerosCrossDeltasEnd);
	ToneStimSegEndShort=ToneStimSegTemp(1:(ZerosCrossDeltasEndMax+e1));
	ToneStimSegEndShort(length(ToneStimSegEndShort))=0;
	ToneStimSegEndDirection=ToneStimSegEndShort(length(ToneStimSegEndShort))-ToneStimSegEndShort(length(ToneStimSegEndShort)-1);
	ZeroCrossIFsEnd=1./(2*diff(ZerosCrossDeltasEnd)*SamplePeriod);
	LastIF=ZeroCrossIFsEnd(length(ZeroCrossIFsEnd));

% --------------------------------------------------------------


x1=x2+1;
x2=x1+length(FMStimSegEndShort)-1;
y1=y2+1;
y2=y1+length(ToneStimSeg)-1;

FM_Train(x1:x2)=FMStimSegEndShort;
Tone_Train(y1:y2)=ToneStimSeg;

end % for n=1:NumSweeps

close(h);









SweepDuration=randn(1,FudgeFact2*NumSweeps);
Direction=(SweepDuration>0); 
Direction=Direction(randperm(length(Direction)));
SweepDuration=SweepDuration*Sigma_Dur+D1; % converts std-norm distr to duration-distr.
SweepDurationThresh=SweepDuration > (LowerBoundFact*D1); % establishes lower-bound to dur-distr
SweepDuration=SweepDuration(SweepDurationThresh);




% --- Creates AM Stim ---------------
t=0:SamplePeriod:sum(SweepDuration)/1000;
Amp=rand(1,NumComps);
InitPhase=2*pi*(rand(1,NumComps)-0.5);

ToneFreqsSeeds=rand(1,NumComps); % Controls the freq-component distribution 
ToneFreqs=ToneFreqsSeeds*(F2-F1)+F1; % converts std-norm distr to duration-distr.
NoiseFreqs=ToneFreqs+(Bandwidth/2)*(randn(1,NumComps)-0.5);

for i=1:NumComps
	Z=Amp(i)*sin(2*pi*NoiseFreqs(i)*t+InitPhase(i));	% Noise
	if i==1
		NoiseStimSeg=zeros(size(Z));
	end
	NoiseStimSeg=NoiseStimSeg+Z;
end % for i=1:NumComps
Noise_Train=NoiseStimSeg;



Am = zeros(size(Noise_Train));
point = 0;
for am = 1:round(length(SweepDuration)*2/3)
    t=0:SamplePeriod:SweepDuration(am)/1000;
    Am(point+1:point+length(t)) = sin(2*pi*(500/SweepDuration(am))*t);
    point = point + length(t);
end
AM_Train = Noise_Train .* (1+Am*1);  
%%%% Increase or decrease the degree of AM by changing '0.3'. It is helpful to maintain a moderate threshold. 


    






% ------------- TrimStim ------------------------


	TrimSamples=floor(SampleRate*D2/1000);

	if length(FM_Train)>TrimSamples
		FM_Train=FM_Train(1:TrimSamples);
      
	end
	if length(Tone_Train)>TrimSamples
		Tone_Train=Tone_Train(1:TrimSamples);
	end
	if length(AM_Train)>TrimSamples
		AM_Train=AM_Train(1:TrimSamples);
	end


% ------------- FadeOut -------------------------


	FadeSamples=floor(SampleRate*FadeTime/1000);

	FadeAmp=(FadeSamples:-1:0)/FadeSamples;  % linear fade to zero
	FM_TrainFade=FM_Train(length(FM_Train)-FadeSamples:end).*FadeAmp;	% apply fade
	FM_Train=[FM_Train(1:length(FM_Train)-FadeSamples-1),FM_TrainFade];	% replace orig seg with faded seg

	FadeAmp=(FadeSamples:-1:0)/FadeSamples;  % linear fade to zero
	Tone_TrainFade=Tone_Train(length(Tone_Train)-FadeSamples:end).*FadeAmp;	% apply fade
	Tone_Train=[Tone_Train(1:length(Tone_Train)-FadeSamples-1),Tone_TrainFade];	% replace orig seg with faded seg

	

% -----------------------------------------------


% -------------- HP-Filter Smoother --------------

	HighFilt=F2;
	[b,a]=butter(FiltOrder,HighFilt/Nyquist);

	FM_Train_Filt=filter(b,a, FM_Train);
	FM_Train=FM_Train_Filt;

	Tone_Train_Filt=filter(b,a, Tone_Train);
	Tone_Train=Tone_Train_Filt;

	

% ------------------------------------------------


% --------------------------------------





RMS_FM_Train=sqrt(mean(FM_Train.^2)); % Calcuilate RMS power values
RMS_Tone_Train=sqrt(mean(Tone_Train.^2));
RMS_AM_Train=sqrt(mean(AM_Train.^2));

FM_Train = FM_Train/RMS_FM_Train;
Tone_Train = Tone_Train/RMS_Tone_Train;
AM_Train = AM_Train/RMS_AM_Train;

% Insert sound level calibration here
% Let me know if you need it: 
% email: xiangbin.teng@nyu.edu
% 







% --------------- Add noise ----------------------
%SNR = [-30:2:0];
%noise_dur = 2000; % ms
%jitter = [400 600 800]


noise = randn(1,noise_dur/1000 * 44100);
RMS_Noise = sqrt(mean(noise.^2));
noise = noise / RMS_Noise;

% scale the signal %
for s = 1: length(SNR)
    
    RMS_signal = 10^(SNR(s)/20);
    FM_noise = noise / RMS_signal;
    Tone_noise = noise / RMS_signal;
    AM_noise = noise / RMS_signal;
    Con_noise = noise / RMS_signal;
    for j = 1:3
        stim_fm = zeros(size(noise));
        stim_fm(jitter(j)/1000 * SampleRate + 1 :jitter(j)/1000 * SampleRate + length(FM_Train)) = FM_Train;  
        stim_FM{s,j} = FM_noise + stim_fm;

        stim_tone = zeros(size(noise));
        stim_tone(jitter(j)/1000 * SampleRate + 1 :jitter(j)/1000 * SampleRate + length(Tone_Train)) = Tone_Train;  
        stim_Tone{s,j} = Tone_noise + stim_tone;
        
        stim_am = zeros(size(noise));
        stim_am(jitter(j)/1000 * SampleRate + 1 :jitter(j)/1000 * SampleRate + length(AM_Train)) = AM_Train;  
        stim_AM{s,j} = AM_noise + stim_am;
        
        stim_Con{s,j} = Con_noise; 
    end
end



% ------------ save wave file and doc ------------ %
cd(WavPath)
mkdir FM_stimuli
mkdir AM_stimuli
mkdir Tone_stimuli
mkdir Control_stimuli

cd ([WavPath '/FM_stimuli'])


fid=fopen(['params_FM.txt'], 'wt');
fprintf(fid, 'WAV-File:  %s\n',['FM_SNR_jitter.wav']);
fprintf(fid, 'Date:      %s\n\n',date);
fprintf(fid, 'Taget StimDuration (D): %5.2f sec\n', D2);
fprintf(fid, 'Mean SweepDuration (S): %5.2f ms\n', D1);
fprintf(fid, 'Num of Sweeps (NS): %i\n', NumSweeps);
fprintf(fid, 'Freq1 (Fa): %5.1f Hz\n', F1);
fprintf(fid, 'Freq1 (Fb): %5.1f Hz\n', F2);
fprintf(fid, 'Bandwidth (BW): %5.1f Hz\n', Bandwidth);
fprintf(fid, 'Sigma_Dur: %5.1f \n', Sigma_Dur);
fprintf(fid, 'NumComps: %i \n', NumComps);
fprintf(fid, 'Totoal duration %i \n', noise_dur);
fprintf(fid, ['Jitter:' num2str(jitter)]);
fprintf(fid, '\n\n----------------------------------\n');
fprintf(fid, ' Freq-Components (Hz)');
fprintf(fid, '  \n----------------------------------\n');

for j = 1:length(jitter)
    for s = 1:length(SNR)
        WAVFileName=['FM_SNR_', num2str(SNR(s)), '_jitter_', num2str(jitter(j))];
        wavwrite(stim_FM{s,j} *0.005,SampleRate,16,WAVFileName);
        
        
        
       
    end
end

cd ([WavPath '/AM_stimuli'])

fid=fopen(['params_AM.txt'], 'wt');
fprintf(fid, 'WAV-File:  %s\n',['AM_SNR_jitter.wav']);
fprintf(fid, 'Date:      %s\n\n',date);
fprintf(fid, 'Taget StimDuration (D): %5.2f sec\n', D2);
fprintf(fid, 'Mean SweepDuration (S): %5.2f ms\n', D1);
fprintf(fid, 'Num of Segments (NS): %i\n', NumSweeps);
fprintf(fid, 'Freq1 (Fa): %5.1f Hz\n', F1);
fprintf(fid, 'Freq1 (Fb): %5.1f Hz\n', F2);
fprintf(fid, 'Bandwidth (BW): %5.1f Hz\n', Bandwidth);
fprintf(fid, 'Sigma_Dur: %5.1f \n', Sigma_Dur);
fprintf(fid, 'NumComps: %i \n', NumComps);
fprintf(fid, 'Totoal duration %i \n', noise_dur);
fprintf(fid, ['Jitter:' num2str(jitter)]);
fprintf(fid, '\n\n----------------------------------\n');
fprintf(fid, ' Freq-Components (Hz)');
fprintf(fid, '  \n----------------------------------\n');


for j = 1:length(jitter)
    for s = 1:length(SNR)
        WAVFileName=['AM_SNR_', num2str(SNR(s)), '_jitter_', num2str(jitter(j))];
        wavwrite(stim_AM{s,j}*0.005,SampleRate,16,WAVFileName);
    end
end




cd ([WavPath '/Tone_stimuli'])

fid=fopen(['params_Tone.txt'], 'wt');
fprintf(fid, 'WAV-File:  %s\n',['Tone_SNR_jitter.wav']);
fprintf(fid, 'Date:      %s\n\n',date);
fprintf(fid, 'Taget StimDuration (D): %5.2f ms\n', D2);
fprintf(fid, 'Mean SweepDuration (S): %5.2f ms\n', D1);
fprintf(fid, 'Num of Segments (NS): %i\n', NumSweeps);
fprintf(fid, 'Freq1 (Fa): %5.1f Hz\n', F1);
fprintf(fid, 'Freq1 (Fb): %5.1f Hz\n', F2);
fprintf(fid, 'Bandwidth (BW): %5.1f Hz\n', Bandwidth);
fprintf(fid, 'Sigma_Dur: %5.1f ms\n', Sigma_Dur);
fprintf(fid, 'NumComps: %i \n', NumComps);
fprintf(fid, 'Totoal duration: %i ms \n', noise_dur);
fprintf(fid, ['Jitter:' num2str(jitter)]);
fprintf(fid, '\n\n----------------------------------\n');
fprintf(fid, ' Freq-Components (Hz)');
fprintf(fid, '  \n----------------------------------\n');


for j = 1:length(jitter)
    for s = 1:length(SNR)
        WAVFileName=['Tone_SNR_', num2str(SNR(s)), '_jitter_', num2str(jitter(j))];
        wavwrite(stim_Tone{s,j}*0.005,SampleRate,16,WAVFileName);
    end
end

cd ([WavPath '/Control_stimuli'])

fid=fopen(['params_Control.txt'], 'wt');
fprintf(fid, 'WAV-File:  %s\n',['Control_SNR.wav']);
fprintf(fid, 'Date:      %s\n\n',date);
fprintf(fid, 'Freq1 (Fa): %5.1f Hz\n', F1);
fprintf(fid, 'Freq1 (Fb): %5.1f Hz\n', F2);
fprintf(fid, 'Bandwidth (BW): %5.1f Hz\n', Bandwidth);

fprintf(fid, 'NumComps: %i \n', NumComps);

fprintf(fid, 'Totoal duration: %i ms \n', noise_dur);

fprintf(fid, '\n\n----------------------------------\n');
fprintf(fid, ' Freq-Components (Hz)');
fprintf(fid, '  \n----------------------------------\n');



    for s = 1:length(SNR)
        WAVFileName=['Control_SNR_', num2str(SNR(s))];
        wavwrite(stim_Con{s,j}*0.005,SampleRate,16,WAVFileName);
    end









   cd(WavPath);
	


   




