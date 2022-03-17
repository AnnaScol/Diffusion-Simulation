clear all; clc; %close all; % clean up
addpath(genpath('./_src'));
addpath(genpath('./_matStoreBounce'));
addpath(genpath('./_matStoreNoBounce'));
addpath(genpath('./dataAnalysis'));

% Seq Parameters 
TE     = 20.0e-3;       % s
ldelta = 10e-3;        % s was 0.011
bVal   = [0,100,200,300,500,600,750,1000,1500,2000,3000,4000,5000,6000,7000,8000,10000,12500,15000,20000]; %s/mm^2
nbVal  = length(bVal);

numTrials  = 18;
testCases  = [500 1000 2000 3000 4000 5000 10000 25000 50000 100000];
nCase      = length(testCases);

% Tissue Paramaters
Radii  = [3.5e-6,4.0e-6,5.0e-6,8.0e-6];        % s/m
nRadii = length(Radii);

% Pre-allocated Memory
vdt = [5e-6];
allData      = zeros(length(vdt),nCase,numTrials,nbVal,nRadii);
nSpinData    = zeros(numTrials,nbVal,nRadii);
stdSpinData  = zeros(length(vdt),nCase,nbVal,nRadii);
varnSpinData = zeros(length(vdt),nCase,nbVal,nRadii);


test_dt_ind = 1;
for stability_idx  = 1:length(testCases)
    for test_idx = 1:numTrials 

        %load data for all numTrials, for all nSpin test cases of 20x4
        allData(test_dt_ind,stability_idx,test_idx,:,:) = load(sprintf("./_matStoreNoBounce/averageRes_%d_%d_dt2",stability_idx,test_idx)).final_res;

    end

    %find standard deviation and variance within each nSpin case
    nSpinData(:,:,:) = squeeze(allData(test_dt_ind,stability_idx,:,:,:));

    stdSpinData(test_dt_ind,stability_idx,:,:)  = squeeze(std(abs(nSpinData),0,1));
    varnSpinData(test_dt_ind,stability_idx,:,:) = squeeze(var(abs(nSpinData),0,1));
    %% Plot single std and variance case 
%         plot_std_vs_bVals(squeeze(stdSpinData(test_dt_ind,stability_idx,:,:)), bVal, squeeze(varnSpinData(test_dt_ind,stability_idx,:))); 

end %end of test case loop



save("./dataAnalysis/stdData-NoBounce-dt2",'stdSpinData');
save("./dataAnalysis/varData-NoBounce-dt2",'varnSpinData');
% save("./dataAnalysis/varData-NoBounce-dt1",'varnSpinData');

%% Calculate and plot difference of dt = 1e-5s for bounce, no bounce
trialsDT1 = 30;
rows = 2;

stdSpinData_NB  = load("./dataAnalysis/stdData-NoBounce-dt1").stdSpinData;
varnSpinData_NB = load("./dataAnalysis/varData-NoBounce-dt1").varnSpinData;

stdSpinData_B   = load("./dataAnalysis/stdData-Bounce-dt1").stdSpinData;
varnSpinData_B  = load("./dataAnalysis/varData-Bounce-dt1").varnSpinData;

diffStdSpinData = stdSpinData_B-stdSpinData_NB;
diffvarnSpinData = varnSpinData_B-varnSpinData_NB;

figure(1); subplot(rows,3,1); hold on;
pltAnalysis(stdSpinData_NB,testCases,trialsDT1,"10us - No Bounce")
hold off;

subplot(rows,3,2); hold on;
pltAnalysis(stdSpinData_B,testCases,trialsDT1, "10us - Bounce")
hold off;

subplot(rows,3,3); hold on;
pltAnalysis(diffStdSpinData,testCases,trialsDT1,"Difference")
hold off;

trialsDT2 = 18;
stdSpinData_NB  = load("./dataAnalysis/stdData-NoBounce-dt2").stdSpinData;
varnSpinData_NB = load("./dataAnalysis/varData-NoBounce-dt2").varnSpinData;

stdSpinData_B   = load("./dataAnalysis/stdData-Bounce-dt2").stdSpinData;
varnSpinData_B  = load("./dataAnalysis/varData-Bounce-dt2").varnSpinData;

diffStdSpinData = stdSpinData_B-stdSpinData_NB;
diffvarnSpinData = varnSpinData_B-varnSpinData_NB;

subplot(2,3,4); hold on;
pltAnalysis(stdSpinData_NB,testCases,trialsDT2,"5us - No Bounce")
hold off;

subplot(2,3,5); hold on;
pltAnalysis(stdSpinData_B,testCases,trialsDT2,"5us - Bounce")
hold off;

subplot(2,3,6); hold on;
pltAnalysis(diffStdSpinData,testCases,trialsDT2,"Difference")
hold off;


