function FM_Train = object_synth_up(Direction,SampleRate,F1,F2,D1,D2,Bandwidth,NumComps,...
    AssignFM)

% stimtrain_synth(SampleRate,F1,F2,D1,D2,Bandwidth,...
%   NumComps,ShowComps,ShowMatchFM,ShowMatchTone,...
%    ShowEndPoints,ShowFreqs,AssignFM,WavWrite)
%    ver 01-05-03
%
% SampleRate:	Sample rate (Hz)
% F1:		    Low-freq (Hz)
% F2:		    High-freq (Hz)
% D1:		    Mean duration of indiv sweep (ms)
% D2:		    Duration of entire stimulus (sec)
% Bandwidth:	"Formant" bandwidth (Hz)
% NumComps:	    Num of freq components (typ. > 20)

% AssignFM:	    Assign FM waveform to workspace [0=No, 1=Yes]
% WavWrite:	    Write wav-file [0=No, 1=Yes]
%
%  e.g., stimtrain_synth(44100,1000,1500,132,5,200,50,0,0,0,0,0,0,1)
%
% *** Make sure to run stim_norm on all stim after running this routine! ***
% This routine only normalizes the three stim types within an SOA, but not
%  across SOA.

%D1 = 1/D1 * 1000;  % convert freq to time, by TXB
% ------------ User-Stuff------------------------------------

%ShortNames=1;  % uses short WAV-filename for Presentation

TrimStim=1;	% trims stim to less than or equal to Duration D2
FadeOut=0;	% linearly ramps stim (from "FadeTime" to end) to zero
FadeTime=100;  % ms

AmpEnv=1; % apply a ramp-ON and ramp-OFF to each segment
RampPeriod=0.003; % seconds; for ON/OFF AmpEnv

FiltStim=1; % LP-filter stim waveform
FiltOrder=3; % Filter order  HP-filter

NumCycles=2; % the num of cycles looked at to find the smallest zero-crossing point
% if this is made too big, it will generate an error for short D1's
% since it will produce search lengths larger than the segment length

ChirpCycles=1; % the num of chirp interpolation cycles

StimAmpEnv=0; % amp-scaling of raw noise: "0"=flat; "1"=linear-decay


DegradeStim=0; % adds noise to stim
DegradeFact=0.25; % additive factor for degrading stim


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




NumArgs=nargin('object_synth');
if nargin ~= NumArgs
    help object_synth;
    return
end


StartTime=now;
Nyquist=SampleRate/2;
SamplePeriod=1/SampleRate;
NumSweeps=round(D2/(D1/1000));

% don't change order of the lines below!
Sigma_Dur=0; % this scales the width of the dur-distr so that it's comparable across D1
SweepDuration=randn(1,FudgeFact2*NumSweeps);

SweepDuration=SweepDuration*Sigma_Dur+D1; % converts std-norm distr to duration-distr.
SweepDurationThresh=SweepDuration > (LowerBoundFact*D1); % establishes lower-bound to dur-distr
SweepDuration=SweepDuration(SweepDurationThresh);

ToneFreqsSeeds=rand(1,NumSweeps); % Controls the freq-component distribution
ToneFreqs=ToneFreqsSeeds*(F2-F1)+F1; % converts std-norm distr to duration-distr.

Linear=linspace(0,1,ceil(FudgeFact1*D2*SampleRate));
Constant=ones(1,ceil(FudgeFact1*D2*SampleRate));
if StimAmpEnv==1
    L=Linear;
else
    L=Constant;
end
Amp=fliplr(L/max(L));

x2=0;
y2=0;
amp2=0;

h=waitbar(0, 'Constructing stimuli...');
Counter=0;
% 
% if exist(WavPath, 'dir')~=7
%     [FILENAME, PATHNAME] = uigetfile('*.*', 'Select Stim Directory');
%     cd(PATHNAME);
% else
%     cd(WavPath);
% end

WAVFileName=['_cnnb', '_D', num2str(D2), '_S', num2str(D1), '_NS',...
    num2str(NumSweeps), '_Fa', num2str(F1), '_Fb', num2str(F2), '_BW',...
    num2str(Bandwidth),'_SYNTH'];
FM_WAVFileName=['FM', WAVFileName];

Short_FM_WAVFileName=['FM', num2str(D1)];

fid=fopen(['params', WAVFileName, '.txt'], 'wt');
fprintf(fid, 'WAV-File:  %s\n',[FM_WAVFileName ,'.wav']);
fprintf(fid, 'Date:      %s\n\n',date);
fprintf(fid, 'Total StimDuration (D): %5.2f sec\n', D2);
fprintf(fid, 'Mean SweepDuration (S): %5.2f ms\n', D1);
fprintf(fid, 'Num of Sweeps (NS): %i\n', NumSweeps);
fprintf(fid, 'Freq1 (Fa): %5.1f Hz\n', F1);
fprintf(fid, 'Freq1 (Fb): %5.1f Hz\n', F2);
fprintf(fid, 'Bandwidth (BW): %5.1f Hz\n', Bandwidth);
fprintf(fid, 'Sigma_Dur: %5.1f \n', Sigma_Dur);
fprintf(fid, 'NumComps: %i \n', NumComps);
fprintf(fid, '\n\n----------------------------------\n');
fprintf(fid, ' Freq-Components (Hz)');
fprintf(fid, '  \n----------------------------------\n');

%fprintf('\n\n');
FMStimSegTemp=[];


for n=1:NumSweeps
    
    Counter=Counter+1;
    waitbar(Counter/NumSweeps, h);
    
 
    t=0:SamplePeriod:SweepDuration(n)/1000;
    
    LogNSigma=0.5;
    LogNMu=0.5;
    
    r=rand(size(t));
    PolyDistr=(t.^PolyDistrExp);
    PolyDistr=(Bandwidth*(PolyDistr./max(PolyDistr)).*r)+1;
    
    if Direction(n)==0 % inverts Freq1 & Freq2 for downward-FMs
        Freq1=F1+(Bandwidth/2)*(randn(1,NumComps)-0.5);
        Freq1 = (Freq1-mean(Freq1)) + F1 - 200 + 400/NumSweeps  * ( n - 1 );
        Freq2=F2+(Bandwidth/2)*(randn(1,NumComps)-0.5);
        Freq2 = (Freq2-mean(Freq2)) + F2 - 200 + 400/NumSweeps * ( n - 1 );
        if 0==1
            [LogNMean,LogNVar]=lognstat(LogNSigma,LogNMu);
            Freq1=F1+(Bandwidth/2)*((rand(1,NumComps)+0.5).*lognrnd(LogNMu,.5,1,NumComps)-LogNMean);
            Freq2=F2+(Bandwidth/2)*((rand(1,NumComps)+0.5).*lognrnd(LogNMu,.5,1,NumComps)-LogNMean);
            
            Freq1=(F1/FreqScale)*PolyDistr;
            Freq2=(F2/FreqScale)*PolyDistr;
        end
    else
        Freq2=F1+(Bandwidth/2)*(randn(1,NumComps)-0.5);
        Freq2 = (Freq2-mean(Freq2)) + F1 - 200 + 400/NumSweeps  * ( n - 1 );
        Freq1=F2+(Bandwidth/2)*(randn(1,NumComps)-0.5);
        Freq1 = (Freq1-mean(Freq1)) + F2 - 200 + 400/NumSweeps  * ( n - 1 );
        if 0==1
            [LogNMean,LogNVar]=lognstat(LogNSigma,LogNMu);
            Freq2=F1+(Bandwidth/2)*((rand(1,NumComps)+0.5).*lognrnd(LogNMu,.5,1,NumComps)-LogNMean);
            Freq1=F2+(Bandwidth/2)*((rand(1,NumComps)+0.5).*lognrnd(LogNMu,.5,1,NumComps)-LogNMean);
            
            Freq1=(F2/FreqScale)*PolyDistr;
            Freq2=(F1/FreqScale)*PolyDistr;
        end
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
        
        
        
        if 0==1 % shows chirp
            plot(FM_InterSweepChirp)
            title(['LastIF=',num2str(LastIF),'  FirstIF=', num2str(FM_FirstIF)]);
            grid;
            pause
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
        
        
        
        if 0==1 % shows chirp
            plot(InterSweepChirp)
            title(['LastIF=',num2str(LastIF),'  FirstIF=', num2str(FirstIF)]);
            grid;
            pause
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


for fp=1:length(Freq1)
    fprintf(fid, ' %8.2f  %8.2f\n', Freq1(fp), Freq2(fp));
end

fprintf(fid, '\n------------------------------------\n');
fprintf(fid, ' Sweep-Durations (ms)');
fprintf(fid, '\n------------------------------------\n');
fprintf(fid, ' %8.2f\n', SweepDuration(1:NumSweeps));
fprintf(fid, '\n------------------------------------\n');
fprintf(fid, ' Tone Frequencies (Hz)');
fprintf(fid, '\n------------------------------------\n');
fprintf(fid, ' %8.2f\n', ToneFreqs);
fclose(fid);


%if 0==1
% --- Creates Noise Stim ---------------
t=0:SamplePeriod:D2;
Amp=rand(1,NumComps);
InitPhase=2*pi*(rand(1,NumComps)-0.5);

ToneFreqsSeeds=rand(1,NumComps); % Controls the freq-component distribution
ToneFreqs=ToneFreqsSeeds*(F2-F1)+F1; % converts std-norm distr to duration-distr.
NoiseFreqs=ToneFreqs(i)+(Bandwidth/2)*(randn(1,NumComps)-0.5);

for i=1:NumComps
    Z=Amp(i)*sin(2*pi*NoiseFreqs(i)*t+InitPhase(i));	% Noise
    if i==1
        NoiseStimSeg=zeros(size(Z));
    end
    NoiseStimSeg=NoiseStimSeg+Z;
end % for i=1:NumComps
Noise_Train=NoiseStimSeg;
% --------------------------------------
%end


FM_Train=FM_Train./(max(FM_Train)); % Normalize data
Tone_Train=Tone_Train./(max(Tone_Train));
Noise_Train=Noise_Train./(max(Noise_Train));

RMS_FM_Train=sqrt(mean(FM_Train.^2)); % Calcuilate RMS power values
RMS_Tone_Train=sqrt(mean(Tone_Train.^2));
RMS_Noise_Train=sqrt(mean(Noise_Train.^2));

RMS_Scaling_FM=1/RMS_FM_Train; % Makes all stim equal RMS power
RMS_Scaling_Tone=1/RMS_Tone_Train;
RMS_Scaling_Noise=1/RMS_Noise_Train;

FM_Train=FM_Train.*RMS_Scaling_FM.^2;
Tone_Train=Tone_Train.*RMS_Scaling_Tone.^2;
Noise_Train=Noise_Train.*RMS_Scaling_Noise.^2;

MaxStimAmp1=max(FM_Train); % Normalizes all stim to max-amp stim
MaxStimAmp2=max(Tone_Train); % Normalizes all stim to max-amp stim
MaxStimAmp3=max(Noise_Train); % Normalizes all stim to max-amp stim

MaxStimAmp=max([MaxStimAmp1,MaxStimAmp2,MaxStimAmp3]); % Normalizes all stim to max-amp stim

FM_Train=FM_Train./MaxStimAmp;
Tone_Train=Tone_Train./MaxStimAmp;
Noise_Train=Noise_Train./MaxStimAmp;


% ------------- TrimStim ------------------------

if TrimStim==1
    TrimSamples=floor(SampleRate*D2);
    
    if length(FM_Train)>TrimSamples
        FM_Train=FM_Train(1:TrimSamples);
    end
    if length(Tone_Train)>TrimSamples
        Tone_Train=Tone_Train(1:TrimSamples);
    end
    if length(Noise_Train)>TrimSamples
        Noise_Train=Noise_Train(1:TrimSamples);
    end
end % if TrimStim==1

% ------------- FadeOut -------------------------

if FadeOut==1
    FadeSamples=floor(SampleRate*FadeTime/1000);
    
    FadeAmp=(FadeSamples:-1:0)/FadeSamples;  % linear fade to zero
    FM_TrainFade=FM_Train(length(FM_Train)-FadeSamples:end).*FadeAmp;	% apply fade
    FM_Train=[FM_Train(1:length(FM_Train)-FadeSamples-1),FM_TrainFade];	% replace orig seg with faded seg
    
    FadeAmp=(FadeSamples:-1:0)/FadeSamples;  % linear fade to zero
    Tone_TrainFade=Tone_Train(length(Tone_Train)-FadeSamples:end).*FadeAmp;	% apply fade
    Tone_Train=[Tone_Train(1:length(Tone_Train)-FadeSamples-1),Tone_TrainFade];	% replace orig seg with faded seg
    
    FadeAmp=(FadeSamples:-1:0)/FadeSamples;  % linear fade to zero
    Noise_TrainFade=Noise_Train(length(Noise_Train)-FadeSamples:end).*FadeAmp;	% apply fade
    Noise_Train=[Noise_Train(1:length(Noise_Train)-FadeSamples-1),Noise_TrainFade];	% replace orig seg with faded seg
end % if FadeOut==1

% -----------------------------------------------


% -------------- HP-Filter Smoother --------------
if FiltStim==1
    HighFilt=F2;
    [b,a]=butter(FiltOrder,HighFilt/Nyquist);
    
    FM_Train_Filt=filter(b,a, FM_Train);
    FM_Train=FM_Train_Filt;
    
    Tone_Train_Filt=filter(b,a, Tone_Train);
    Tone_Train=Tone_Train_Filt;
    
    Noise_Train_Filt=filter(b,a, Noise_Train);
    Noise_Train=Noise_Train_Filt;
end
% ------------------------------------------------

if DegradeStim==1
    HighFilt=F2;
    [b,a]=butter(FiltOrder,HighFilt/Nyquist); % filter from above
    FM_AdditiveNoise=filter(b,a, rand(size(FM_Train)));
    Tone_AdditiveNoise=filter(b,a, rand(size(Tone_Train)));
    Noise_AdditiveNoise=filter(b,a, rand(size(Noise_Train)));
    
    FM_Train=FM_Train+DegradeFact*FM_AdditiveNoise;
    Tone_Train=Tone_Train+DegradeFact*Tone_AdditiveNoise;
    Noise_Train=Noise_Train+DegradeFact*Noise_AdditiveNoise;
end



% -------- mov-avg filter (DON'T USE!)---------------------
if 0==1
    FM_Train_Filt=0;
    h=waitbar(0, 'Smoothing Stim...');
    Counter=0;
    
    for k=(NumFiltSamples+1):length(FM_Train)
        FM_Train_Filt(k-NumFiltSamples)=...
            sum(FM_Train((k-NumFiltSamples):k))/(NumFiltSamples+1);
        Counter=Counter+1;
        waitbar(Counter/length(FM_Train), h);
    end
    close(h);
end %if 0==1
% -----------------------------------------------------------

if 0==1  	% use these to determine (indiv) segment bandwidth
    % only use this on Tones or Noise, not FM
    figure;
    psd(Tone_Train,[],SampleRate);
    
    [P,F]=psd(Tone_Train,[],SampleRate);
    figure;
    semilogy(F,P,'r-');
    hold
    semilogy(F,sqrt(P),'b-');
end % if 0==1


if 0==1
    figure;
    subplot(2,1,1);
    specgram(FM_Train,[],SampleRate,[],[]);
    colormap('autumn');
    title(['FM-Specgram: ', FM_WAVFileName], 'Interpreter','none');
    subplot(2,1,2);
    plot(FM_Train);
    text(0.7*length(FM_Train),0.9,['FM RMS: ',num2str(RMS_FM_Train)],'FontWeight','bold');
    set(gca,'YLim',[-1.1,1.1]);
    grid;
    title(['FM-Waveform: ', FM_WAVFileName], 'Interpreter','none');
    
    figure;
    subplot(2,1,1);
    specgram(Tone_Train,[],SampleRate,[],[]);
    colormap('autumn');
    title(['Tone-Specgram: ', Tone_WAVFileName], 'Interpreter','none');
    subplot(2,1,2);
    plot(Tone_Train);
    text(0.7*length(Tone_Train),0.9,['Tone RMS: ',num2str(RMS_Tone_Train)],'FontWeight','bold');
    set(gca,'YLim',[-1.1,1.1]);
    grid;
    title(['Tone-Waveform: ', Tone_WAVFileName], 'Interpreter','none');
    
    figure;
    subplot(2,1,1);
    specgram(Noise_Train,[],SampleRate,[],[]);
    colormap('autumn');
    title(['Noise-Specgram: ', Noise_WAVFileName], 'Interpreter','none');
    subplot(2,1,2);
    plot(Noise_Train);
    text(0.7*length(Noise_Train),0.9,['Noise RMS: ',num2str(RMS_Noise_Train)],'FontWeight','bold');
    set(gca,'YLim',[-1.1,1.1]);
    grid;
    title(['Noise-Waveform: ', Noise_WAVFileName], 'Interpreter','none');
end % if 0==1


% figure;
% subplot(3,1,1);
% specgram(FM_Train,[],SampleRate,[],[]);
% set(gca,'XTick',[]);
% xlabel('');
% colormap('autumn');
% title(['FM-Specgram: ', FM_WAVFileName], 'Interpreter','none');
% subplot(3,1,2);
% specgram(Tone_Train,[],SampleRate,[],[]);
% set(gca,'XTick',[]);
% xlabel('');
% colormap('autumn');
% %title(['Tone-Specgram: ', Tone_WAVFileName], 'Interpreter','none');
% subplot(3,1,3);
% specgram(Noise_Train,[],SampleRate,[],[]);
% colormap('autumn');
% %title(['Noise-Specgram: ', Noise_WAVFileName], 'Interpreter','none');


[PF,FF]=psd(FM_Train,[],SampleRate,[],[]);
[PT,FT]=psd(Tone_Train,[],SampleRate,[],[]);
[PN,FN]=psd(Noise_Train,[],SampleRate,[],[]);

if 0==1
    figure;
    subplot(3,1,1);
    loglog(FF,PF,'r-');
    set(gca,'XLim',[1E2,1E4]);
    set(gca,'YLim',[1E-10,1E3]);
    set(gca,'XTick',[]);
    xlabel('');
    grid;
    title(['FM-PSD: ', FM_WAVFileName], 'Interpreter','none');
    subplot(3,1,2);
    loglog(FT,PT,'r-');
    set(gca,'XTick',[]);
    set(gca,'XLim',[1E2,1E4]);
    set(gca,'YLim',[1E-10,1E3]);
    xlabel('');
    grid;
    %title(['Tone-Specgram: ', Tone_WAVFileName], 'Interpreter','none');
    subplot(3,1,3);
    loglog(FN,PN,'r-');
    set(gca,'XLim',[1E2,1E4]);
    set(gca,'YLim',[1E-10,1E3]);
    grid;
    %title(['Noise-Specgram: ', Noise_WAVFileName], 'Interpreter','none');
end


[MinPlotSize, MinPlotIndex]=min([length(FF),length(FT),length(FN)]);
if MinPlotIndex==1
    F=FF;
elseif MinPlotIndex==2
    F=FT;
    elseIF MinPlotIndex==3
    F=FN;
end;
% figure;
% loglog(F,PF,'r-', F,PT,'b-', F,PN,'k-');
% set(gca,'XLim',[1E2,1E4]);
% set(gca,'YLim',[1E-10,1E3]);
% legend('FM','Tone','Noise', 1);
% grid;
% title(['PSDs: ', WAVFileName], 'Interpreter','none');


% if PlotAnalytic==1
%     ht=hilbert(FM_Train);
%     hta=abs(ht);
%     htp=atan(real(ht)./imag(ht));
%     htpd=diff(htp);
%     
%     figure;
%     subplot(4,1,1);
%     plot(FM_Train,'b-');
%     title('FM Waveform');
%     subplot(4,1,2);
%     plot(hta,'k-');
%     title('Instantaneous Amplitude');
%     subplot(4,1,3);
%     plot(htp,'g-');
%     title('Instantaneous Phase');
%     subplot(4,1,4);
%     plot(1./htpd,'r-');
%     title('Instantaneous Frequency');
% end % if PlotAnalytic==1
% 

% figure;
% hist(SweepDuration(1:NumSweeps));
% title(['Sweep-Dur Distr.', ' (Total=',	num2str(NumSweeps), ')'],...
%     'Interpreter','none');
% ylabel([' Mean: ', num2str(round(mean(SweepDuration(1:NumSweeps)))),...
%     ' ms', '  StdDev: ', num2str(round(std(SweepDuration(1:NumSweeps)))), ' ms']);
% 
% figure;
% subplot(2,2,1);
% hist(Freq1);
% title(['F1 Freq-Comps Distr.']);
% ylabel([' Mean: ', num2str(round(mean(Freq1))),...
%     ' Hz', '  StdDev: ', num2str(round(std(Freq1))), ' Hz']);
% subplot(2,2,2);
% hist(Freq2);
% title(['F2 Freq-Comps Distr.']);
% ylabel([' Mean: ', num2str(round(mean(Freq2))),...
%     ' Hz', '  StdDev: ', num2str(round(std(Freq2))), ' Hz']);
% subplot(2,2,3);
% hist(ToneFreqs);
% title(['Tone Freq-Comps Distr.']);
% ylabel([' Mean: ', num2str(round(mean(ToneFreqs))),...
%     ' Hz', '  StdDev: ', num2str(round(std(ToneFreqs))), ' Hz']);
% subplot(2,2,4);
% hist(NoiseFreqs);
% title(['Noise Freq-Comps Distr.']);
% ylabel([' Mean: ', num2str(round(mean(NoiseFreqs))),...
%     ' Hz', '  StdDev: ', num2str(round(std(NoiseFreqs))), ' Hz']);


%sound(FM_Train,SampleRate,16);
%pause(max(size(FM_Train))*SamplePeriod);

%sound(Tone_Train,SampleRate,16);
%pause(max(size(Tone_Train))*SamplePeriod);

%sound(Noise_Train,SampleRate,16);


fprintf('\n');
fprintf('Total StimDuration (D): %5.2f sec\n', D2);
fprintf('Mean SweepDuration (S): %5.2f ms\n', D1);
fprintf('Num of Sweeps (NS): %i\n', NumSweeps);
fprintf('Freq1 (Fa): %5.1f Hz\n', F1);
fprintf('Freq1 (Fb): %5.1f Hz\n', F2);
fprintf('Bandwidth (BW): %5.1f Hz\n', Bandwidth);
fprintf('Sigma_Dur: %5.1f \n', Sigma_Dur);
fprintf('NumComps: %i \n\n', NumComps);





if AssignFM==1
    assignin('base','FM',FM_Train);
end

EndTime=now;

fprintf('\nTotal Run-Time: %5.2f min\n\n', (60*24)*(EndTime-StartTime));

