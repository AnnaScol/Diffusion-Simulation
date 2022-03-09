clear all; close all; % clean up
addpath(genpath('./_src'));

% Simulation parameter
nSpins            = 1000;
duration          = 0.05; %sec
movements_per_sec = 1.0e6; %1um steps
nTimeSteps     = duration*movements_per_sec;
dt = 1/movements_per_sec;
t  = movements_per_sec*dt;

% Water Diffusion Parameter
D  = 3.0e-9; %free water
nD = 3;

%%%%%%%%%%%% TEST WALKS %%%%%%%%%%%%%
coords = rndWalks(D,nSpins,nTimeSteps,dt);
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

final_D = MSD/(2*nD*dt*nTimeSteps);
%convert to be in seconds
% final_D = (final_D*(4));
fprintf("Final D %d\n",final_D);



