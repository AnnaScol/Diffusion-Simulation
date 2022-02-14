%% Spin Echo Line Scan
% TODO:
% Vectorize the mirror bounce 

clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

dt    = 10*10^-6; 
gamma = 2*pi*42.577*10^6;


%Allocate the memory needed
% nTimeSteps  = 7000*10^-5/dt;%70ms
nTimeSteps  = 3000*10^-5/dt;%70ms
rfPulse     = zeros(1,nTimeSteps); %variable to hold a RF waveform
gradAmp     = zeros(3,nTimeSteps); %variable to hold a gradient waveform
adc         = zeros(1,nTimeSteps); %variable to hold a gradient waveform
time        = zeros(1,nTimeSteps); %variable to hold the time points

% %% Parameters %% %

% Diffusion Gradients
ldelta = 0.011; %ms was 0.04
sdelta = 0.003; %ms was 0.02
D = 3e-9; %m^2/ms

nSpins = 2000;
G = ([0,5,7.5,10,12.5,15,17.5,20,22,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110]*1e-3)*15; %mT
% G = ([0,5,10,15,25,30,40,50,70,80,100,125,150]*1e-3); %mT
% G = ([0,5,10,15,25,30,40,50,70,80,100]*1e-3)*15; %mT


nG = length(G);
mFinalVect = zeros(nG,2); %variable to hold the final magnetization calculated for each position

b = gamma^2*G.^2*sdelta^2*(ldelta-sdelta/3);

T1 = 1500*(10^-3);
T2 = 1000*(10^-3); 
TE = 20*(10^-3);

% constraint_radii = ([1 2.5 5 7.5 10 12.5 15 20 1.0e6]*1e-6); 
constraint_radii = ([2 5 10 20 40 80]*1e-6); 

for i=1:nTimeSteps %i starts at 1 go's to 15000
    time(i)    = i*dt;                       %Time in seconds
end

% Generate and Display MRI Sequence
% RF Excitation Waveform
pulsedurE = 0.001; % duration of the RF in s
rfStepsE = round(1:(pulsedurE/(dt)));
rfPulse(rfStepsE) = apodize_sinc_rf(length(rfStepsE),3,pi/2,dt); %B1+ in Tesla

%RF Refocusing pules
pulsedurR = 0.001; % duration of the RF in s
rfStepsR = round(1:(pulsedurR/(dt)));
rfPulseR = apodize_sinc_rf(length(rfStepsR),3,pi,dt); %B1+ in Tesla
rfPulse(round(TE/2/dt) +length(rfStepsE)/2 + rfStepsR) = rfPulseR;

%%% Diffusion Pulses %%%
% First diffusion gradient
diffusionGradient1_loc = round((1:(sdelta/dt))+ pulsedurE/dt);
gradAmp(3,diffusionGradient1_loc) =  G(3); %Z gradients in Tesla per meter

% Diffusion pulse 2
diffusionGradient2_loc = round((1:(sdelta/dt)) + pulsedurE/dt + pulsedurR/dt + ldelta/dt);
gradAmp(3,diffusionGradient2_loc) =  G(3); %Z gradients in Tesla per meter

%%% Plotting the data points %%%
Coords = zeros(length(constraint_radii),3,nSpins,nTimeSteps); %particle start loc is assume 0,0,0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% PLOTTING %% %
figure
subplot(2,1,1); plot(time,rfPulse,'k-','LineWidth',2);title('RF Pulse'); 
xlabel('time (s)'), ylabel('|B_{1}^{+}| (T)');grid on;

subplot(2,1,2); plot(time,gradAmp(3,:),'b-','LineWidth',2);title('Slice Select Gradient');
xlabel('time (s)'), ylabel('G_{z}(T/m)');grid on;

% %% Get location matrix --> 3 x nSpins dimensions
nSet = 40;
SetsOfSpins = nSet;
store_final_vect = zeros(SetsOfSpins,length(b),length(constraint_radii));
n    = 3;
%% %%%%%%%%%%%%% COORDS SETUP %%%%%%%%%%%%%
ReTry_RandomWalkPathGenerator(nSet,nSpins,constraint_radii,nTimeSteps,diffusionGradient2_loc,D,n,dt)

%% %%%%%%%%%%%%% Perform sequence %%%%%%%%% 

for SpinSet = 1:nSet
    
    fprintf("\n\n Set %d/%d \n\n",SpinSet,SetsOfSpins);
    file_path = "C:\Users\s4427550\3D_Coords";%3D_Coords/Coords
    load(sprintf("%s%d.mat",file_path,SpinSet));

    
    for radius_bounds = 1:length(constraint_radii)

        disp("Starting sequence");

        for i = 1:nG

            j = 1;
            gradAmp(3,diffusionGradient1_loc) =  G(i); %Z gradients in Tesla per meter
            gradAmp(3,diffusionGradient2_loc) =  G(i); %Z gradients in Tesla per meter

            mT = zeros(3,nSpins);
            mZ = ones(3,nSpins);

            %starting spin locations
            [mT,mZ] =  bloch(dt,([squeeze(Coords(radius_bounds,1,:,1)),...
                                  squeeze(Coords(radius_bounds,2,:,1)), ...
                                  squeeze(Coords(radius_bounds,3,:,1))]'),0,T1,T2,mT,mZ); 

            for j = 2:nTimeSteps

                dB0 = gradAmp(:,j)'*([squeeze(Coords(radius_bounds,1,:,j)),...
                                      squeeze(Coords(radius_bounds,2,:,j)), ...
                                      squeeze(Coords(radius_bounds,3,:,j))]'); 

                [mT,mZ] =  bloch(dt,dB0,rfPulse(j),T1,T2,mT,mZ); 

                % condition is true if it has reached the read point
                if (j > diffusionGradient2_loc(end))
                    mFinalVect(i,:) = [mean(mT,'all'), mean(mZ,'all')];  
                    break;
                end

            end

            disp(['b-value = ' num2str(round(b(i)*1e-6))]);
        end


        store_final_vect(SpinSet,:,radius_bounds) = (mFinalVect(:,1));
    end
        
end


%%
save('mat_store/6_Large_final_vect.mat','store_final_vect')
%% Sum all results 
final_res = squeeze(sum((store_final_vect),1)./SetsOfSpins);
save('mat_store/compare_result2_40x2000spins.mat','final_res')
%% Plot the results
nFig = 4;
plot_signal_vs_b_results(length(constraint_radii), final_res, b, nFig)
