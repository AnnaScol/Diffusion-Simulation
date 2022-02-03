%% Spin Echo Line Scan


clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

dt    = 10*10^-6; 
gamma = 2*pi*42.577*10^6;


%Allocate the memory needed
nTimeSteps  = 7000*10^-5/dt;%70ms
rfPulse     = zeros(1,nTimeSteps); %variable to hold a RF waveform
gradAmp     = zeros(3,nTimeSteps); %variable to hold a gradient waveform
time        = zeros(1,nTimeSteps); %variable to hold the time points

% %% Parameters %% %

% Diffusion Gradients
ldelta = 0.040; %ms
sdelta = 0.020; %ms
D = 1e-10;    % (1.0e-6)cm^2/s --> (1.0e-10)m^2/s

nSpins = 500;
G = ([0,5,7.5,10,12.5,15,17.5,20,22,25,30,35,40,45,50,55,60,65,70,75,80,100,125,150]*1e-3); %mT


nG = length(G);
mFinalVect = zeros(nG,2); %variable to hold the final magnetization calculated for each position

b = gamma^2*G.^2*sdelta^2*(ldelta-sdelta/3);

T1 = 1000*(10^-3);
T2 = 1000*(10^-3); 
TE = 65*(10^-3);

% Make size soma cell and brain cell sizes - vary 0.005-0.1mm
constraint_radii = ([1 2.5 5 7.5 10 12.5 15 20 1.0e6]*1e-6); 
% %%   %% %
for i=1:nTimeSteps %i starts at 1 go's to 15000
    time(i)    = i*dt;                       %Time in seconds
end

% Generate and Display MRI Sequence
% %% RF %% %
% RF Excitation Waveform
pulsedurE = 0.002; % duration of the RF in s
rfStepsE = round(1:(pulsedurE/(dt)));
rfPulse(rfStepsE) = apodize_sinc_rf(length(rfStepsE),3,pi/2,dt); %B1+ in Tesla

% First diffusion gradient
diffusionGradient1_loc = round((1:(sdelta/dt))+ pulsedurE/dt);
gradAmp(3,diffusionGradient1_loc) =  G(3); %Z gradients in Tesla per meter

%RF Refocusing pules
pulsedurR = 0.002; % duration of the RF in s
rfStepsR = round(1:(pulsedurR/(dt)));
rfPulseR = apodize_sinc_rf(length(rfStepsR),3,pi,dt); %B1+ in Tesla
rfPulse(round(TE/2/dt) +length(rfStepsE)/2 + rfStepsR) = rfPulseR;

% diffusion pulse 2
diffusionGradient2_loc = round((1:(sdelta/dt)) + pulsedurE/dt + pulsedurR/dt + ldelta/dt);
gradAmp(3,diffusionGradient2_loc) =  G(3); %Z gradients in Tesla per meter

location = zeros(3,nTimeSteps);

%%% Plotting the data points %%%
Coords = zeros(length(constraint_radii),3,nSpins,nTimeSteps); %particle start loc is assume 0,0,0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% PLOTTING %% %
figure
subplot(2,1,1); plot(time,rfPulse,'k-','LineWidth',2);title('RF Pulse'); 
xlabel('time (s)'), ylabel('|B_{1}^{+}| (T)');grid on;

subplot(2,1,2); plot(time,gradAmp(3,:),'b-','LineWidth',2);title('Slice Select Gradient');
xlabel('time (s)'), ylabel('G_{z}(T/m)');grid on;

%%
nSet = 1;
z0 = zeros(3,nSpins);
dz = (randn(3,nSpins)-0.5)*sqrt(2*D*(ldelta-sdelta/3)); % randn has mean 0.5 and std(1)


%% Simulate just the last points S(b) Value

for j = 1:nG
    

    % simulation
    z = z0+dz;
    time_to_record = TE+pulsedurE/2;
    
    res = [0;0;G(j)]*gamma.*z*time_to_record;
    
    mFinalVect(j,1) = mean(res,'all');
           
    
    disp(['b = ' num2str(round(b(j)*1e-6)) 's/mm^2']);
end
figure;
plot(b,1-abs(mFinalVect(:,1)),'.-');
