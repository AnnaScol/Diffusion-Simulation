clear all; clc; close all; % clean up
addpath(genpath('./_src'));
addpath(genpath('./_matStoreBounce'));
addpath(genpath('./_matStoreNoBounce'));

% vdt = [1.0e-5 1.0e-6 5.0e-6];
vdt = [1.0e-5]; %only done this up to 15 so far


% Seq Parameters 
TE     = 20.0e-3;       % s
ldelta = 10e-3;        % s was 0.011
bVal   = [0,100,200,300,500,600,750,1000,1500,2000,3000,4000,5000,6000,7000,8000,10000,12500,15000,20000]; %s/mm^2
nbVal  = length(bVal);

numTrials  = 16;
testCases = [500 1000 2000 3000 4000 5000 10000 25000 50000 100000];
nCase      = length(testCases);

% Tissue Paramaters
Radii  = [3.5e-6,4.0e-6,5.0e-6,8.0e-6];        % s/m
nRadii = length(Radii);

% Pre-allocated Memory
allData      = zeros(length(vdt),nCase,numTrials,nbVal,nRadii);
nSpinData    = zeros(numTrials,nbVal,nRadii);
stdSpinData  = zeros(length(vdt),nCase,nbVal,nRadii);
varnSpinData = zeros(length(vdt),nCase,nbVal,nRadii);

for test_dt_ind = 1:length(vdt)
   
    dt = vdt(test_dt_ind); 

    for stability_idx  = 1:length(testCases)
        for test_idx = 1:numTrials 
            
            %load data for all numTrials, for all nSpin test cases of 20x4
            allData(test_dt_ind,stability_idx,test_idx,:,:) = load(sprintf("./_matStoreBounce/averageRes_%d_%d_dt%d",stability_idx,test_idx,test_dt_ind)).final_res;
            
        end
        
        %find standard deviation and variance within each nSpin case
        nSpinData(:,:,:) = squeeze(allData(test_dt_ind,stability_idx,:,:,:));
        
        stdSpinData(test_dt_ind,stability_idx,:,:)  = squeeze(std(abs(nSpinData),0,1));
        varnSpinData(test_dt_ind,stability_idx,:,:) = squeeze(var(abs(nSpinData),0,1));
        %% Plot single std and variance case 
%         plot_std_vs_bVals(squeeze(stdSpinData(test_dt_ind,stability_idx,:,:)), bVal, squeeze(varnSpinData(test_dt_ind,stability_idx,:))); 

    end %end of test case loop
end

%% Plot Std Dev. vs nSpins for Single dt over the 10 nSpin Cases
BVALS = 2;

stdSpinData1 = stdSpinData(1,:,:,:);
meanSTD = squeeze(mean(squeeze(stdSpinData1),BVALS));

%convert to normalized log data
logTestCases = log(testCases);
finalRes     = (meanSTD./norm(meanSTD));
vVar         = (squeeze(mean(squeeze(varnSpinData(1,:,:,:)),BVALS))./...
                    norm(squeeze(mean(squeeze(varnSpinData(1,:,:,:)),BVALS))));

figure(1); hold on;

for ii = 1:size(meanSTD,2)
    errorbar(logTestCases,finalRes(:,ii),vVar(:,ii),'LineWidth',1);
end


set(gca,'xticklabels',{'500', '1000', '2000', '3000', '4000', '5000', '10000', '25000', '50000', '100000'});

xlabel('nSpins'), ylabel('norm(Std Dev. |M_{xy}|)');

grid on; title(sprintf('Stability across %d trials for nSpins',numTrials));
% legend(sprintf("3.5*10^{-6}  var = %d",vVar(1)),sprintf("4*10^{-6}  var = %d",vVar(2)), ...
%             sprintf("5*10^{-6} var = %d",vVar(3)),sprintf("8*10^{-6}  var = %d",vVar(4)))
legend("3.5*10^{-6}","4*10^{-6}","5*10^{-6}","8*10^{-6}")


hold off;

