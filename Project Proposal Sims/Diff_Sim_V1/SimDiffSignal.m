clear all; clc; close all; % clean up
addpath(genpath('./_src'));

% Sim Parameters
dt     = 1.0e-6;        % s
% gamma  = 2*pi*42.577e6; % rad/s/T
nSpins = 500;
nSet   = 1;


% Seq Parameters 
pulsedurE = 0.002; %s
pulsedurR = 0.002; %s
TE     = 40.0e-3;       % s
sdelta = 0.01e-3;        % s was 0.020
ldelta = 20.0e-3;        % s was 0.010
bVal   = [0 500 1000 1500 2000 3000 4000 5000 6000 8000 10000 20000]*1.0e-6;

% Tissue Paramaters
T1     = 1.5;           % s  (not used)
T2     = 1.0;           % s  (not used)
D      = 3.0e-9;        % m^2/s
Radius = 1000.0e-6;        % m

res = zeros(length(bVal),2,nSet);
%%
for b_idx = 1:length(bVal) 
    fprintf("B Index %d\n", b_idx);
        %% %%%%%%%%%%%%% GET SEQUENCE %%%%%%%%%%%%%%%
        Sequence_Struct = getDifSeq(pulsedurE, pulsedurR, TE, sdelta, ldelta, bVal(b_idx)/1.0e-6, dt);

        %% %%%%%%%%%%%%% CHECK SEQUENCE %%%%%%%%%%%%%
%         dispSequence(Sequence_Struct);
        
     %% %%%%%%%%%%%%% GET RND WALKS %%%%%%%%%%%%% 
    for iSet = 1:nSet
        spinCor = getRandomWalks(D,Radius,nSpins,length(Sequence_Struct.RF),dt);

        % fubar = sqrt(dot(coord,coord,1));
        % max(squeeze(fubar(:,:,end)))

        %% %%%%%%%%%%%%% Perform sequence %%%%%%%%% 
        res(b_idx,:,iSet) = simulateMRISequence(Sequence_Struct,T1,T2,spinCor,dt);
    end
end

%% %%%%%%%%%%%%% Final Plotting %%%%%%%%% 
final_res = squeeze(sum((res),3)./nSet);

figure; 
plot(bVal/1e-6,abs(final_res(:,1)),'.-','LineWidth',2);
xlabel('b (s/mm^2)'), ylabel('|M_{xy}|');grid on;
title('Absolute of Diffusion Attenuation');