%% Spin Echo Line Scan
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

dt    = 10^-5; 
gamma = 2*pi*42.577*10^6;


%Allocate the memory needed
nTimeSteps  = 7000;%500ms
rfPulse     = zeros(1,nTimeSteps); %variable to hold a RF waveform
gradAmp     = zeros(3,nTimeSteps); %variable to hold a gradient waveform
adc         = zeros(1,nTimeSteps); %variable to hold a gradient waveform
time        = zeros(1,nTimeSteps); %variable to hold the time points

% %% Parameters %% %

% Diffusion Gradients
ldelta = 0.030; %ms
sdelta = 0.015; %ms
D = 2.0e-12;    %m^2/s

nSpins = 5000;
start_position = zeros(1,nSpins);
% dz = (randn(1,nSpins)-0.5)*sqrt(2*D*(ldelta-sdelta/3)); 

G = ([0,5,10,15,20,25,30,35,40,50,60,70,80,100,150]*1e-3); %T

nG = length(G);
mFinalVect = zeros(nG,2); %variable to hold the final magnetization calculated for each position

b = gamma^2*G.^2*sdelta^2*(ldelta-sdelta/3);

T1 = 850*(10^-3);
T2 = 60*(10^-3); 
TE = 45*(10^-3);

constraint_radius = 100;

% %%   %% %
for i=1:nTimeSteps %i starts at 1 go's to 15000
    time(i)    = i*dt;                       %Time in seconds
end

% Generate and Display MRI Sequence
% %% RF %% %
% RF Excitation Waveform
pulsedurE = 0.001; % duration of the RF in s
rfStepsE = round(1:(pulsedurE/(dt)));
rfPulse(rfStepsE) = apodize_sinc_rf(length(rfStepsE),6,pi/2,dt); %B1+ in Tesla

% First diffusion gradient
diffusionGradient1_loc = round((1:(sdelta/dt))+ pulsedurE/dt);
% diffusionGradient1_loc = rfStepsE(end)+round(1:(sdelta/dt));
gradAmp(3,diffusionGradient1_loc) =  G(2); %Z gradients in Tesla per meter

%RF Refocusing pules
pulsedurR = 0.001; % duration of the RF in s
rfStepsR = round(1:(pulsedurR/(dt)));
rfPulseR = apodize_sinc_rf(length(rfStepsR),6,pi,dt); %B1+ in Tesla
rfPulse(round(TE/2/dt) +length(rfStepsE)/2 + rfStepsR) = rfPulseR;

% diffusion pulse 2
diffusionGradient2_loc = round((1:(sdelta/dt)) + pulsedurE/dt + pulsedurR/dt + ldelta/dt);
gradAmp(3,diffusionGradient2_loc) =  G(2); %Z gradients in Tesla per meter


% %% PLOTTING %% %
figure
subplot(3,1,1); plot(time,rfPulse,'k-','LineWidth',2);title('RF Pulse'); 
xlabel('time (s)'), ylabel('|B_{1}^{+}| (T)');grid on;

subplot(3,1,2); plot(time,phase(rfPulse),'k-','LineWidth',2);title('RF Pulse Phase'); 
xlabel('time (s)'), ylabel('|B_{1}^{+}| (T)');grid on;

subplot(3,1,3); plot(time,gradAmp(3,:),'b-','LineWidth',2);title('Slice Select Gradient');
xlabel('time (s)'), ylabel('G_{z}(T/m)');grid on;

%% Perform seuquence

for i = 1:nG
    gradAmp(3,diffusionGradient1_loc) =  G(i); %Z gradients in Tesla per meter
    gradAmp(3,diffusionGradient2_loc) =  G(i); %Z gradients in Tesla per meter
    
    mT = zeros(1,nSpins);
    mZ = ones(1,nSpins);
    
    %starting spin locations
    deltaZ = zeros(3,nSpins);
    [mT,mZ] =  bloch(dt,deltaZ,0,T1,T2,mT,mZ); 
    
    for j = 2:nTimeSteps
        dz = (-1 + (2)*rand(3,nSpins)).*sqrt(2*D*(ldelta-sdelta/3));
        deltaZ = dz+deltaZ; %adding a random step to each spin in each timestep between diffusion grads
        dB0 = gradAmp(:,j).*deltaZ;
        
        [mT,mZ] =  bloch(dt,dB0,rfPulse(j),T1,T2,mT,mZ); 
        
        % condition is true if it has reached the read point
        if (j > diffusionGradient2_loc(end))
            mFinalVect(i,:) = [mean(mT,'all'), mean(mZ,'all')];  
            break;
        end

    end
   
    disp(['b-value = ' num2str(round(b(i)*1e-6))]);
    
end


figure;
subplot(3,1,1);plot(b*1e-6,abs(mFinalVect(:,1)),'o-','LineWidth',2);
xlabel('b (s/mm^2)'), ylabel('|M_{xy}|')
title('Absolute of Diffusion Attenuation');

subplot(3,1,2);plot(b*1e-6,angle(mFinalVect(:,1)),'o-','LineWidth',2);
xlabel('b (s/mm^2)'), ylabel('\angleM_{xy} '); ylim([-pi pi])
title('Angle of Diffusion Attenuation')


subplot(3,1,3);
for spin = 1:nSpins
    plot3(deltaZ(1,spin),deltaZ(2,spin),deltaZ(3,spin),'x');
    hold on; 
end
title('Final location of particle')

%% FUNCTIONS %%

function result = InBounds(r,x,y,z)
    origin = [0 0 0];
    mag = sqrt((abs(x)-origin(1))^2 + (abs(y)-origin(2))^2 + (abs(z)-origin(3))^2);
    if ( mag < r)
        result = 1; % (x,y) is inside the circle
    elseif ( mag == r)
        result = 2; % (x,y) is on the circle
    elseif ( mag > r)
        result = 3; % (x,y) is outside the circle
    end
end
