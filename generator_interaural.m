%% generate temporal sequence

FudgeFact1 = 1.5;
SampleRate = 44100;
LowerBoundFact = 0.7;


%%%%% generate theta sequence
dur = 5;
seg_dur = 190;
sigma_dur = 30;
NumSweeps=round(dur/(seg_dur/1000));
sweep_dur = randn(1,2*NumSweeps)


Direction=(sweep_dur>0)
Direction=Direction(randperm(length(Direction)))
SweepDuration=sweep_dur*sigma_dur+seg_dur
SweepDurationThresh=SweepDuration > (LowerBoundFact*seg_dur);
SweepDuration=SweepDuration(SweepDurationThresh);

theta_sequence = SweepDuration(1:NumSweeps)

%%%%% generate alpha sequence
dur = 5;
seg_dur = 100;
sigma_dur = 15;
NumSweeps=round(dur/(seg_dur/1000));
sweep_dur = randn(1,2*NumSweeps)


Direction=(sweep_dur>0)
Direction=Direction(randperm(length(Direction)))
SweepDuration=sweep_dur*sigma_dur+seg_dur
SweepDurationThresh=SweepDuration > (LowerBoundFact*seg_dur);
SweepDuration=SweepDuration(SweepDurationThresh);

alpha_sequence = SweepDuration(1:NumSweeps)



%%%%% generate alpha sequence
dur = 5;
seg_dur = 25;
sigma_dur = 3;
NumSweeps=round(dur/(seg_dur/1000));
sweep_dur = randn(1,2*NumSweeps)


Direction=(sweep_dur>0)
Direction=Direction(randperm(length(Direction)))
SweepDuration=sweep_dur*sigma_dur+seg_dur
SweepDurationThresh=SweepDuration > (LowerBoundFact*seg_dur);
SweepDuration=SweepDuration(SweepDurationThresh);

gamma_sequence = SweepDuration(1:NumSweeps)

%% generate stimulus

%%%%% theta
theta_noise = [];
for segment = 1:length(theta_sequence)
    noise_seg = zeros(floor(theta_sequence(segment)*SampleRate/1000),1);
    noise_raw = randn(floor(theta_sequence(segment)*SampleRate/1000),1);
    noise_seg(:,1) = noise_raw;noise_seg(:,2) = noise_raw;
    if rem(segment,2) > eps
    noise_seg = randn(floor(theta_sequence(segment)*SampleRate/1000),2);
    noise_seg = TX_ormatrix(noise_seg);
    end
    scalor = repmat(std(noise_seg),length(noise_seg),1);
    theta_noise = [theta_noise;noise_seg./scalor];
    
end
theta_noise = theta_noise * 0.3;
wavwrite(theta_noise,44100,'theta_noise')



%%%%% alpha
alpha_noise = [];
for segment = 1:length(alpha_sequence)
    noise_seg = zeros(floor(alpha_sequence(segment)*SampleRate/1000),1);
    noise_raw = randn(floor(alpha_sequence(segment)*SampleRate/1000),1);
    noise_seg(:,1) = noise_raw;noise_seg(:,2) = noise_raw;
    if rem(segment,2) > eps
    noise_seg = randn(floor(alpha_sequence(segment)*SampleRate/1000),2);
    noise_seg = TX_ormatrix(noise_seg);
    end
    scalor = repmat(std(noise_seg),length(noise_seg),1);
    alpha_noise = [alpha_noise;noise_seg./scalor];
    
end
alpha_noise = alpha_noise * 0.3;
wavwrite(alpha_noise,44100,'alpha_noise')






%%%%% gamma
gamma_noise = [];
for segment = 1:length(gamma_sequence)
    noise_seg = zeros(floor(gamma_sequence(segment)*SampleRate/1000),1);
    noise_raw = randn(floor(gamma_sequence(segment)*SampleRate/1000),1);
    noise_seg(:,1) = noise_raw;noise_seg(:,2) = noise_raw;
    if rem(segment,2) > eps
    noise_seg = randn(floor(gamma_sequence(segment)*SampleRate/1000),2);
    noise_seg = TX_ormatrix(noise_seg);
    end
    scalor = repmat(std(noise_seg),length(noise_seg),1);
    gamma_noise = [gamma_noise;noise_seg./scalor];
    
end
gamma_noise = gamma_noise * 0.3;
wavwrite(gamma_noise,44100,'gamma_noise')



%%%%% control
control_noise = [];
for segment = 1:length(theta_sequence)
    noise_seg = zeros(floor(theta_sequence(segment)*SampleRate/1000),1);
    noise_raw = randn(floor(theta_sequence(segment)*SampleRate/1000),1);
    noise_seg(:,1) = noise_raw;noise_seg(:,2) = noise_raw;
    scalor = repmat(std(noise_seg),length(noise_seg),1);
    control_noise = [control_noise;noise_seg./scalor];
    
end
control_noise = control_noise * 0.3;
wavwrite(control_noise,44100,'control_noise')









