ranseed = 2;
fs = 44100;
res = [ 0  ]
fundamental = 200;
gamma = [ 1 2];
numHarmonics = 20;

for g = gamma
    for resolution = res
        
        
        
        pointsInSpectrum=2^17;
        amp=1./((1:pointsInSpectrum-1).^g);
        rand('seed',ranseed);
        
        % Get frequency component
        phase=rand(1,pointsInSpectrum-1)*2*pi;
        x=amp.*(cos(phase)+i*sin(phase));
        
        % Converse into time series
        s=[0 x zeros(1,resolution * pointsInSpectrum-1) fliplr(conj(x)) ];
        
        
        f=real(ifft(s));
        % normalize
        f=f-mean(f);
        f=f/sqrt(var(f));
        
        
        global fvect
        sMod=zeros(size(f)); % initialise the output signal to zeros
        f=f-min(f)+1;
        phaseFact=2*pi*fundamental/fs;
        while max(f) > 2,
            f=f./2;
        end;
        fvect=[];
        while min(f) < numHarmonics,
            phaseAngle=0;
            % fadeFact is used to fade a tonal component in or out if it leaves or
            % re-enters the passband of our stimulus in the course of its frequency
            % modulation
            fadeFact=fadeFactor(f,numHarmonics);
            % ampCorrections is used to apply an amplitude correction to the pure
            % tone component if necessary to equalise sound transducer output.
            % If you want to incorporate amplitude corrections, you must implement
            % a findAmpCorrections function which returns a vector of appropriate
            % correction in dB for the given input vector of frequencies.
            % If no such function is found, then we assume that no amplitude
            % corrections are necessary, i.e. the transducer is flat.
            try
                ampCorrections=findAmpCorrections(f*fundamental);
            catch
                ampCorrections=ones(1,length(f));
            end;
            for ii=1:length(f),
                e=f(ii);
                if and( (e > 1) , (e < numHarmonics))
                    phaseStep=e*phaseFact;
                    phaseAngle=phaseAngle+phaseStep;
                    sMod(ii)=sMod(ii)+fadeFact(ii)*ampCorrections(ii)*sin(phaseAngle);
                    %sMod(ii)=sMod(ii)+sin(phaseAngle);
                else
                    phaseAngle=0;
                end;
            end;
            f=f.*2^(1/2);  % take a 3rd octave step
            fvect=[fvect min(f)];
        end;
        length(f)/44100
        % finally scale the sound to have an RMS value of unity
        rmsVal=sqrt(mean(sMod.^2));
        sMod=sMod./rmsVal;
        
        
        name = ['stim_' num2str(g) '_' num2str(resolution) '.wav']
        
        wavwrite(sMod(1:5*44100)*0.01,44100,name);
        % subplot(2,1,2)
        % [b,f,t]=specgram(sMod,1024,44100);
        % pcolor(t,f/1000,(abs(b))); shading flat; set(gca,'yscale','log');
        % set(gca,'ytick',[0.5 1 2 4 8 16],'yticklabel',[0.5 1 2 4 8 16]);
        % xlabel('seconds'); ylabel('kHz');
        
        
    end
end
