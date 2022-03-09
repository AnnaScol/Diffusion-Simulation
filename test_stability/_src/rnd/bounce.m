D = 3e-9; %m^2/s
cellRadius = 5.0e-6; %m
nTimeSteps = 40000;
dt = 1.0e-6; %time steps
nSpins = 2; 
coords = getRndWalks(D,cellRadius,nSpins,nTimeSteps,dt);
