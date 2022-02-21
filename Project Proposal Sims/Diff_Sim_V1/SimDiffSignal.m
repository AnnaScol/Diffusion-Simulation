clear all; clc; close all; % clean up
addpath(genpath('./_src'));

% Sim Parameters
dt     = 1.0e-6;        % s
% gamma  = 2*pi*42.577e6; % rad/s/T
nSpins = 1000;
nSet   = 1;


% Seq Parameters 
pulsedurE = 0.001; %s
pulsedurR = 0.001; %s
TE     = 55.0e-3;       % s
ldelta = 40.0e-3;        % s was 0.011
sdelta = 10.0e-3;        % s was 0.002
bVal   = [0 50 100 200 500 1000 1250 1500 1750 2000 2500 3000 4000 ...
            5000 6000 7000 8000 9000 10000 12500 15000 20000]*1.0e-6;


% Tissue Paramaters
T1     = 1.5;           % s  (not used)
T2     = 1.0;           % s  (not used)
D      = 3.0e-9;        % m^2/s
Radius = 10.0e-6;        % m



res = zeros(length(bVal),2,nSet);

for b_idx = 1:length(bVal) 
    
        %% %%%%%%%%%%%%% GET SEQUENCE %%%%%%%%%%%%%%%

        Sequence_Struct = getDifSeq(pulsedurE, pulsedurR, TE, sdelta, ldelta, bVal(b_idx), dt);

        %% %%%%%%%%%%%%% CHECK SEQUENCE %%%%%%%%%%%%%

%         dispSequence(Sequence_Struct);
%         pause

     %% %%%%%%%%%%%%% GET RND WALKS %%%%%%%%%%%%% 
    for iSet = 1:nSet
        spinCor = getRandomWalks(D,Radius,nSpins,length(Sequence_Struct.RF),dt);%nSet,nSpins,constraint_radii,nTimeSteps,Sequence_Struct.END_diffusionGradient2_loc,D,dim,dt)

        % fubar = sqrt(dot(coord,coord,1));
        % max(squeeze(fubar(:,:,end)))

        %% %%%%%%%%%%%%% Perform sequence %%%%%%%%% 
        res(b_idx,:,iSet) = simulateMRISequence(Sequence_Struct,T1,T2,spinCor,dt);
    end
end

%% %%%%%%%%%%%%% Final Plotting %%%%%%%%% 
final_res = squeeze(sum((res),3)./nSet);

for b_idx = 1:length(bVal) 
    figure(1); 
    plot(bVal(b_idx)/1e-6,abs(final_res(b_idx,1)),'.-','LineWidth',2);
    xlabel('b (s/mm^2)'), ylabel('|M_{xy}|');grid on;
    title('Absolute of Diffusion Attenuation');
    hold on
end