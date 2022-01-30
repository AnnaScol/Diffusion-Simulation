%% Spin Echo Line Scan
%do a 1d , get log plot
%map to a paper
%try different strating points
%e^{b*distanct traveled}
%finding the phase change
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

dt    = 10^-5; 
gamma = 2*pi*42.577*10^6;


%Allocate the memory needed
nTimeSteps  = 7000;%70ms
rfPulse     = zeros(1,nTimeSteps); %variable to hold a RF waveform
gradAmp     = zeros(3,nTimeSteps); %variable to hold a gradient waveform
adc         = zeros(1,nTimeSteps); %variable to hold a gradient waveform
time        = zeros(1,nTimeSteps); %variable to hold the time points

% %% Parameters %% %

% Diffusion Gradients
ldelta = 0.050; %ms
sdelta = 0.015; %ms
D = 2.0e-7;    %m^2/s

nSpins = 1000;
start_position = zeros(1,nSpins);

% G = ([50]*1e-3); %mT

G = ([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,90,100,110]*1e-3); %mT

nG = length(G);
mFinalVect = zeros(nG,2); %variable to hold the final magnetization calculated for each position

b = gamma^2*G.^2*sdelta^2*(ldelta-sdelta/3);

T1 = 8500*(10^-3);
T2 = 8000*(10^-3); 
TE = 60*(10^-3);

%make size soma cell and brain cell sizes

% constraint_radius = 2.0e-5;m
constraint_radius = 5;

% %%   %% %
for i=1:nTimeSteps %i starts at 1 go's to 15000
    time(i)    = i*dt;                       %Time in seconds
end

% Generate and Display MRI Sequence
% %% RF %% %
% RF Excitation Waveform
pulsedurE = 0.001; % duration of the RF in s
rfStepsE = round(1:(pulsedurE/(dt)));
rfPulse(rfStepsE) = apodize_sinc_rf(length(rfStepsE),3,pi/2,dt); %B1+ in Tesla

% First diffusion gradient
diffusionGradient1_loc = round((1:(sdelta/dt))+ pulsedurE/dt);
gradAmp(3,diffusionGradient1_loc) =  G(1); %Z gradients in Tesla per meter

%RF Refocusing pules
pulsedurR = 0.001; % duration of the RF in s
rfStepsR = round(1:(pulsedurR/(dt)));
rfPulseR = apodize_sinc_rf(length(rfStepsR),3,pi,dt); %B1+ in Tesla
rfPulse(round(TE/2/dt) +length(rfStepsE)/2 + rfStepsR) = rfPulseR;

% diffusion pulse 2
diffusionGradient2_loc = round((1:(sdelta/dt)) + pulsedurE/dt + pulsedurR/dt + ldelta/dt);
gradAmp(3,diffusionGradient2_loc) =  G(1); %Z gradients in Tesla per meter

location = zeros(3,nTimeSteps);

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
    r = zeros(3,nSpins);
    [mT,mZ] =  bloch(dt,r,0,T1,T2,mT,mZ); 
    
    for j = 2:nTimeSteps
        
        %%% Bounds checking %%%
        rnd = (-1 + 2*rand(3,nSpins));
        rnd = rnd(:,:)./sqrt(sum(rnd(:,:).*rnd(:,:)));
        
%         dz =  rnd.*sqrt(2*D*(ldelta-sdelta/3));
        dz = rnd*D;
        r = dz+r; %adding a random step to each spin in each timestep between diffusion grads
%         deltaZ = InBounds(constraint_radius,temp(1,:),temp(2,:),temp(3,:), dz(1,:),dz(2,:),dz(3,:));
        %%%%%%%%%%%%%%%%%%%%%%%
        
        dB0 = gradAmp(:,j).*r; 
        [mT,mZ] =  bloch(dt,dB0,rfPulse(j),T1,T2,mT,mZ); 
        
        % condition is true if it has reached the read point
        if (j > diffusionGradient2_loc(end))
            mFinalVect(i,:) = [mean(mT,'all'), mean(mZ,'all')];  
            break;
        end

    end
    
    disp(['b-value = ' num2str(round(b(i)*1e-6))]);
end
%0.02
%% PLOTTING %%
% Make it log values --> -(1/b)*ln(S_{DWI}/S_{b}) (dont need th 1/b maybe)
%make sure to scale both axis
%
figure;
% subplot(3,1,1);plot(b*1e-6,-log(abs(mFinalVect(:,1))./b'),'o-','LineWidth',2);
% xlabel('b (s/mm^2)'), ylabel('ln(|M_{xy}|/b)');grid on;
% title('Absolute of Diffusion Attenuation');

subplot(2,1,1);plot(b*1e-6,abs(mFinalVect(:,1)),'o-','LineWidth',2);
xlabel('b (s/mm^2)'), ylabel('|M_{xy}|');grid on;
title('Absolute of Diffusion Attenuation');

% subplot(3,1,2);plot(b*1e-6,angle(mFinalVect(:,1)),'o-','LineWidth',2);
% xlabel('b (s/mm^2)'), ylabel('\angleM_{xy} '); ylim([-pi pi]);grid on;
% title('Angle of Diffusion Attenuation')


subplot(2,1,2);
for spin = 1:nSpins
    plot3(r(1,spin),r(2,spin),r(3,spin),'x');
    hold on; 
end
xlabel('x'), ylabel('y'),zlabel('z');
title('Final location of particle')

%% FUNCTIONS %%

function result = InBounds(r,x1,y1,z1,dx,dy,dz)

    origin = [0 0 0];
    mag = sqrt((abs(x1)-origin(1)).^2 + (abs(y1)-origin(2)).^2 + (abs(z1)-origin(3)).^2);
    
    for i = 1:length(mag)
        
        if ( mag(i) > r) %outside circle
            x1(i) = x1(i)-dx(i); 
            y1(i) = y1(i)-dy(i); 
            z1(i) = z1(i)-dz(i);  % (x,y) is outside the circle and reverse step
        end
        
    end
    
    result = [x1;y1;z1];
end
