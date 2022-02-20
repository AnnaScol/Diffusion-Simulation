clear all; clc; close all; % clean up

% Diffusion Gradients
ldelta = 0.004; %ms was 0.011
sdelta = 0.002; %ms was 0.002

% Parameters 
dt     = 1*10^-6; 
gamma  = 2*pi*42.577*10^6;
D      = 3e-9; %m^2/ms
T1     = 1500*(10^-3);
T2     = 1000*(10^-3); 
TE     = 20*(10^-3);

nSpins      = 500;
nSet        = 1;
dim         = 3; %dimensions in random walk

%Allocate the memory needed
nTimeSteps  = round(ldelta/dt + 2*sdelta/dt)+10;%70ms
rfPulse     = zeros(1,nTimeSteps); %variable to hold a RF waveform
gradAmp     = zeros(3,nTimeSteps); %variable to hold a gradient waveform

% Generate gradient amplitude from a set of b (s/mm^2) values 
b = [0 50 100 150 200 250 300 400 500 750 1000 1250 1500 1750 2000 2500 3000 4000 5000 6000 7000 8000 9000 10000]; % s/mm^2
G = findGValues((b/1e-6),gamma,ldelta,sdelta);
nG = length(G);

constraint_radii = ([2 3 4 4.5 5 6 8 10]*1e-6); 
%% %%%%%%%%%%%%% CHECK SEQUENCE %%%%%%%%%%%%%
pulsedurE = 0.001; %s
pulsedurR = 0.001; %s
result_struct = MRI_Sequence(pulsedurE, pulsedurR, sdelta, ldelta, G(2), nTimeSteps, dt);

%% %%%%%%%%%%%%% COORDS SETUP %%%%%%%%%%%%%
rndWalkPathGenerator(nSet,nSpins,constraint_radii,nTimeSteps,result_struct.END_diffusionGradient2_loc,D,dim,dt)

%% %%%%%%%%%%%%% Perform sequence %%%%%%%%% 
store_final_vect = zeros(nSet,length(b),length(constraint_radii));
mFinalVect = zeros(nG,2); %variable to hold the final magnetization calculated for each position
diffusionGradient1_loc = result_struct.diffusionGradient1_loc;
diffusionGradient2_loc = result_struct.diffusionGradient2_loc;
rfPulse = result_struct.rfPulse;
gradAmp = result_struct.gradAmp;

tic
for SpinSet = 1:nSet
    
    fprintf("\n\n Set %d/%d \n\n",SpinSet,nSet);
    file_path = "3D_Coords/Coords";%3D_Coords/Coords
    load(sprintf("%s%d.mat",file_path,SpinSet));
 
    for radius_bounds = 1:length(constraint_radii)
        disp("Starting sequence");
        spinLocs = squeeze(Coords(radius_bounds,:,:,:));
       
        mFinalVect = GVec_simulateMRISequence(G,gradAmp,rfPulse,T1,T2,diffusionGradient1_loc,diffusionGradient2_loc,spinLocs,nSpins,dt);
        disp(['b-value = ' num2str(round(b(1:nG)))]);

        store_final_vect(SpinSet,:,radius_bounds) = (mFinalVect(:,1));
    end
        
end
toc
str = sprintf("raw_%dx%dSpins_vec.mat",nSet,nSpins);
save(sprintf('mat_store/%s.mat', str),'store_final_vect');
%% Sum all results 
final_res = squeeze(sum((store_final_vect),1)./nSet);
str = sprintf("averaged_%dx%dSpins.mat",nSet,nSpins);
save(sprintf('mat_store/%s',str),'final_res');
%% Plot the results
nFig = 4;
plot_signal_vs_b_results(length(constraint_radii), final_res, b, nFig)