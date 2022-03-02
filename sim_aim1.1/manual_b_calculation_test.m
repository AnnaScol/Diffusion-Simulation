%% Spin Echo Line Scan
% TODO:
% Vectorize the mirror bounce 

clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

dt    = 1*10^-6; 
gamma = 2*pi*42.577*10^6;

TE = 20*(10^-3);

%Allocate the memory needed
nTimeSteps  = round(TE/dt)  +20000;%70ms
rfPulse     = zeros(1,nTimeSteps); %variable to hold a RF waveform
gradAmp     = zeros(3,nTimeSteps); %variable to hold a gradient waveform
time        = zeros(1,nTimeSteps); %variable to hold the time points

% %% Parameters %% %

% Diffusion Gradients
ldelta = 10e-3;        
sdelta = 5.0e-3; 
D = 3e-9; %m^2/ms

nSpins = 500;
% G = ([0,5,10,15,25,30,40,50,70,80,100]*1e-3); %mT
G = (258.98*1e-3); %mT


nG = length(G);
mFinalVect = zeros(nG,2); %variable to hold the final magnetization calculated for each position

b = gamma^2*G.^2*sdelta^2*(ldelta-sdelta/3);%s/m

T1 = 1500*(10^-3);
T2 = 1000*(10^-3); 

% constraint_radii = ([1 2.5 5 7.5 10 12.5 15 20 1.0e6]*1e-6); 
constraint_radii = ([20]*1e-6); 

for i=1:nTimeSteps %i starts at 1 go's to 15000
    time(i)    = i*dt;                       %Time in seconds
end

% RF Excitation Waveform
pulsedurE = 0.001; % duration of the RF in s

%RF Refocusing pules
pulsedurR = 0.001; % duration of the RF in s

%%% Diffusion Pulses %%%
% First diffusion gradient
diffusionGradient1_loc = round((1:(sdelta/dt))+ pulsedurE/dt);
gradAmp(3,diffusionGradient1_loc) =  G(1); %Z gradients in Tesla per meter

% Diffusion pulse 2
diffusionGradient2_loc = round((1:(sdelta/dt)) + pulsedurE/dt + pulsedurR/dt + ldelta/dt);
gradAmp(3,diffusionGradient2_loc) =  G(1); %Z gradients in Tesla per meter

%%
% %% Checking B-Values %% %
b_test = zeros(1,length(b));
for g_values = 1:nG
    
    moment = 0;
    int1 = zeros(1,round(TE/dt));
    int2 = zeros(1,round(TE/dt));
    int3 = zeros(1,round(TE/dt));

    gradAmp(3,diffusionGradient1_loc) =  G(g_values); 
    gradAmp(3,diffusionGradient2_loc) =  -G(g_values); % negative moment due to being after 180deg pulse
    
    for t = 1:round(TE/dt)
        int1(t) = ((dt) * gradAmp(3,t));
        
    end
    
    for t = 1:round(TE/dt)
        if (t > 1)
            
            int2(t) = (dt)*int1(t) + (sum(int1(1:(t-1)))*(length(1:(t-1))*dt));
        else
            int2(t) = (dt)*int1(t) + int1(1)*(dt);
        end

    end
    
    for t = 1:round(TE/dt)
        if (t > 1)
            int3(t) = (dt)*int2(t)+(sum(int2(1:(t-1)))*(length(1:(t-1))*dt));
        else
            int3(t) = (dt)*int2(t)+int2(1)*(dt);
        end
    end
    
    
    moment_sum = sum(int3);
    
    b_test(g_values) = (gamma^2)*(moment_sum^2 * dt^2); %*1e-6 uses mm measurments
end
fprintf("actual b\n")
disp(b_test*1.0e-6);
fprintf("formula b\n")
disp(b*1.0e-6);
