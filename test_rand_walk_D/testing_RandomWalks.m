clear all; close all; % clean up
addpath(genpath('./_src'));

% Simulation parameter
nSpins            = 100;
duration          = 10.0e-3; %s
movements_per_sec = 1.0e6;   %1um steps
nTimeSteps     = duration*movements_per_sec;
dt = 1/movements_per_sec;    %s

% Water Diffusion Parameter
D  = 3.0e-9; %free water
nD = 3;

%%%%%%%%%%% % TEST WALKS %%%%%%%%%%%%%
% rndWalks1 is twice as fast as rndWalks2
% for 100 spins, averaged over 100 iterations rndWalks1 = 0.06s and rndWalks 2 = 0.12s

coords = rndWalks1(D,nSpins,nTimeSteps,dt); %constant time, distance from distribution
% coords = rndWalks2(D,nSpins,nTimeSteps,dt); %constant length normalized from normal distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Plot Distance Travelled %%%%%%%%%%%%
magnitudeFromOrigin = sqrt(squeeze(dot(coords(:,:,end),coords(:,:,end),1)));
figure(1); 
histObject = histogram(magnitudeFromOrigin, 30);
caption    = sprintf('Distribution of %d Final Distances', nSpins);
title(caption); xlabel('Distance'); ylabel('Count');
grid on; hold off;

%%%%%%%%%%%% CALCULATE D %%%%%%%%%%%%
dr_squared = squeeze(dot(coords(:,:,end),coords(:,:,end),1));              
MSD        = mean(dr_squared);

final_D = MSD/(2*nD*dt); %m^2/1us
final_D = (final_D* (1/nTimeSteps)); %convert to be in seconds
fprintf("\nFinal Calculated Diffusivity %d\n",final_D);

% Predicted Mean distance
predMean = sqrt(2*D*3*duration);
fprintf("predicted mean: %d, actual mean %d\n",predMean,mean(magnitudeFromOrigin));