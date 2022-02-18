%% Spin Echo Line Scan
% TODO:
% Vectorize the mirror bounce 

clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

dt    = 1*10^-6; 
gamma = 2*pi*42.577*10^6;

% %% Parameters %% %

% Diffusion Gradients
% ldelta = 0.003; %ms was 0.04
sdelta = 0.003; %ms was 0.02
D = 3e-9; %m^2/ms
nSpins = 50;
G = [120]*1e-3; %mT
nG = length(G);
%the result of bDelta is only larger than sDelta at b-values > 200
b_values = [0 50 100 150 200 250 300 400 500 750 1000 1250 1500 1750 2000 2500 3000 4000 ...
                5000 6000 7000 8000 9000 10000]/1e-6;
bDelta_values = solveBforDelta('b', sdelta, b_values ,G);
nB = length(b_values);
mFinalVect = zeros(nG,2); %variable to hold the final magnetization calculated for each position

T1 = 1500*(10^-3);
T2 = 1000*(10^-3); 

%Allocate the memory needed
nTimeSteps  = round(bDelta_values(end)/dt + 2*sdelta/dt);%70ms
rfPulse     = zeros(1,nTimeSteps); %variable to hold a RF waveform
gradAmp     = zeros(3,nTimeSteps); %variable to hold a gradient waveform
adc         = zeros(1,nTimeSteps); %variable to hold a gradient waveform
time        = zeros(1,nTimeSteps); %variable to hold the time points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

constraint_radii = ([3 4 5 6 8 10]*1e-6); 

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
%%% Diffusion Pulses %%%
% First diffusion gradient
diffusionGradient1_loc = round((1:(sdelta/dt))+ pulsedurE/dt);
gradAmp(3,diffusionGradient1_loc) =  G; %Z gradients in Tesla per meter
% Diffusion pulse 2
diffusionGradient2_loc = round((1:(sdelta/dt)) + pulsedurE/dt + pulsedurR/dt + bDelta_values(5)/dt);
gradAmp(3,diffusionGradient2_loc) =  G; %Z gradients in Tesla per meter
%%% RF excitation pulse %%%
TE = diffusionGradient2_loc(end)*dt;
rfPulse(round(TE/2/dt) + length(rfStepsE)/2 + rfStepsR) = rfPulseR;

%%% Plotting the data points %%%
Coords = zeros(length(constraint_radii),3,nSpins,nTimeSteps); %particle start loc is assume 0,0,0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% PLOTTING first b_value/bDelta sequence %% %
figure
subplot(2,1,1); plot(time,rfPulse,'k-','LineWidth',2);title('RF Pulse'); 
xlabel('time (s)'), ylabel('|B_{1}^{+}| (T)');grid on;
subplot(2,1,2); plot(time,gradAmp(3,:),'b-','LineWidth',2);title('Slice Select Gradient');
xlabel('time (s)'), ylabel('G_{z}(T/m)');grid on;
%%
% %% Get location matrix --> 3 x nSpins dimensions
nSet = 20;
SetsOfSpins = nSet;
store_final_vect = zeros(SetsOfSpins,length(bDelta_values),length(constraint_radii));
n    = 3;
%% %%%%%%%%%%%%% COORDS SETUP %%%%%%%%%%%%%
%TODO: need to change the diffusion gradient 2 location since it will
%change with bdelta
ReTry_RandomWalkPathGenerator(nSet,nSpins,constraint_radii,nTimeSteps,diffusionGradient2_loc,D,n,dt)
%% %%%%%%%%%%%%% Perform sequence %%%%%%%%% 

for SpinSet = 1:nSet
    
    fprintf("\n\n Set %d/%d \n\n",SpinSet,SetsOfSpins);
    file_path = "3D_Coords/Coords";%3D_Coords/Coords
    load(sprintf("%s%d.mat",file_path,SpinSet));

    
    for radius_bounds = 1:length(constraint_radii)

        disp("Starting sequence");

        for i = 1:nB

            j = 1;
            %% need to change this to cycle through Bdelta and rf variations
            result_struct = MRI_Sequence(pulsedurE, pulsedurR, sdelta,  bDelta_values(i), G,nTimeSteps, dt);
            
           
            gradAmp(3,1:length(result_struct.gradAmp(3,:))) =  result_struct.gradAmp(3,:);
            rfPulse = result_struct.rfPulse;
            %%
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
                if (j > result_struct.END_diffusionGradient2_loc)
                    mFinalVect(i,:) = [mean(mT,'all'), mean(mZ,'all')];  
                    break;
                end

            end

            disp(['b-value = ' num2str(round(b_values(i)*1e-6))]);
        end


        store_final_vect(SpinSet,:,radius_bounds) = (mFinalVect(:,1));
    end
        
end


%%
save('mat_store/testing.mat','store_final_vect')
%% Sum all results 
final_res = squeeze(sum((store_final_vect),1)./SetsOfSpins);
%% Plot the results
nFig = 4;
plot_signal_vs_b_results(length(constraint_radii), final_res, b_values, nFig)