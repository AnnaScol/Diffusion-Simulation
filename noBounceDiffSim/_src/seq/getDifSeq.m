function result_struct = getDifSeq(pulsedurE, pulsedurR, TE, sDelta, lDelta, bVal, dt)
% resulting struct contains gradient amplitudes and rf pulse for sequence
    
    
    % allocate memory 
    nTimeSteps = round(TE/dt) + round(pulsedurE/dt/2);
    
    result_struct.G  = zeros(3,nTimeSteps);
    result_struct.RF = zeros(1,nTimeSteps);

    %%% RF Pulses %%%

    % get RF Excitation Waveform
    rfStepsE = round(pulsedurE/dt);
    rfPulseE = apodize_sinc_rf(rfStepsE,3,pi/2,dt); %B1+ in Tesla
    
    % get RF Refocusing Waveform
    rfStepsR = round(pulsedurR/dt);
    rfPulseR = apodize_sinc_rf(rfStepsR,3,pi,dt); %B1+ in Tesla
    
    % get Refocusing pulse start and end time 
    refPulseStaTime = rfPulseE/2 + TE/dt/2 - rfStepsR/2;
    refPulseEndTime = refPulseStaTime + rfStepsR -1;
    
    % prepare final RF time series
    result_struct.RF(1,              1:       rfStepsE)        = rfPulseE;
    result_struct.RF(1,round(refPulseStaTime:refPulseEndTime)) = rfPulseR;
    
    %%% Diffusion Pulses %%%
    
    % get gradient amplitude
    DifGradAmp = findGValues((bVal),lDelta,sDelta);
%     fprintf("Gamp %d\n",DifGradAmp);
    nGradSteps = round(sDelta/dt);
    
    % get second gradient start and end time
    secGradStartTime = refPulseEndTime;
    
    % get First pulse start and end time
    firGradStartTime = secGradStartTime - round(lDelta/dt);
  
    
    % prepare final grad time series
    result_struct.G(3,round(firGradStartTime:(firGradStartTime+nGradSteps-1))) = DifGradAmp;
    result_struct.G(3,round(secGradStartTime:(secGradStartTime+nGradSteps-1))) = DifGradAmp;
 
end