function result_struct = MRI_Sequence(pulsedurE, pulsedurR, sdelta, bdelta, G, nTimesteps, dt)
% resulting struct contains gradient amplitudes and rf pulse for sequence
        
    % Generate and Display MRI Sequence
    % RF Excitation Waveform
    rfStepsE = round(1:(pulsedurE/(dt)));
    rfPulse(rfStepsE) = apodize_sinc_rf(length(rfStepsE),3,pi/2,dt); %B1+ in Tesla
    
    %RF Refocusing pules
    rfStepsR = round(1:(pulsedurR/(dt)));
    rfPulseR = apodize_sinc_rf(length(rfStepsR),3,pi,dt); %B1+ in Tesla
    
    %%% Diffusion Pulses %%%
    % First diffusion gradient
    diffusionGradient1_loc = round((1:(sdelta/dt))+ pulsedurE/dt);
    gradAmp(3,diffusionGradient1_loc) =  G; %Z gradients in Tesla per meter
    % Diffusion pulse 2
    diffusionGradient2_loc = round((1:(sdelta/dt)) + pulsedurE/dt + pulsedurR/dt + bdelta/dt);
    gradAmp(3,diffusionGradient2_loc) =  G; %Z gradients in Tesla per meter
    %%% RF excitation pulse %%%
    TE = diffusionGradient2_loc(end)*dt;
    rfPulse(round(TE/2/dt) + length(rfStepsE)/2 + rfStepsR) = rfPulseR;
    
    %%% OUTPUT %%%
    result_struct.gradAmp = gradAmp;
    result_struct.rfPulse = [rfPulse zeros(1,(diffusionGradient2_loc(end)-length(rfPulse))+1)];
    result_struct.END_diffusionGradient2_loc = diffusionGradient2_loc(end);
    result_struct.diffusionGradient1_loc = diffusionGradient1_loc;
    result_struct.diffusionGradient2_loc = diffusionGradient2_loc;
end