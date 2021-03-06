clear all; clc; close all; % clean up
addpath(genpath('./_src'));

% Sim Parameters
dt     = 20.0e-6;        % s
nSpins = 10000;
nSet   = 5;


% Seq Parameters 
pulsedurE = 0.001; %s
pulsedurR = 0.001; %s
TE     = 55.0e-3;       % s
ldelta = 25e-3;        % s was 0.011
sdelta = 5.0e-3;        % s was 0.002
bVal   = [0,100,200,300,500,600,750,1000,1500,2000,3000,4000,5000,6000,7000,8000,10000,12500,15000,20000]; %s/mm^2

% Tissue Paramaters
T1     = 1.5;           % s  (not used)
T2     = 1.0;           % s  (not used)
D      = 3.0e-9;        % m^2/s
Radii = [3.5e-6,4.0e-6,5.0e-6,8.0e-6];        % s/m

%remove radius loop and do only one at a time and check

%%
res = zeros(nSet,length(bVal),length(Radii));
for iSet = 1:nSet 
    fprintf('\n------------Set %d/%d ------------\n', iSet,nSet);
    for rIdx = 1:length(Radii)
        fprintf('Radius %d/%d\n', rIdx,length(Radii));
        %% %%%%%%%%%%%%% GET RND WALKS %%%%%%%%%%%%% 
        spinCor = getRndWalks(D,Radii(rIdx),nSpins,(round(TE/dt) + round(pulsedurE/dt/2)),dt);
        
        for bIdx = 1:length(bVal)
            %% %%%%%%%%%%%%% GET SEQUENCE %%%%%%%%%%%%%%%
            Sequence_Struct = getDifSeq(pulsedurE, pulsedurR, TE, sdelta, ldelta, bVal(bIdx)/1e-6, dt);
            %check b-result
            
            %% %%%%%%%%%%%%% CHECK SEQUENCE %%%%%%%%%%%%%
            %dispSequence(Sequence_Struct);
            
            %% %%%%%%%%%%%%% Perform sequence %%%%%%%%%%%
            res(iSet,bIdx,rIdx) = simulateMRISequence(Sequence_Struct,T1,T2,spinCor,dt);
            
        end
    end
end

%% %%%%%%%%%%%%% Final Plotting %%%%%%%%% 
final_res = squeeze(sum(res,1)./nSet);
plot_signal_vs_b_results(final_res, bVal);
