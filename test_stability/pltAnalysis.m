function pltAnalysis(stdSpinData,testCases,criteria)
% plots the coefficient of variation

    BVALS = 2;

    stdSpinData1 = stdSpinData(1,:,:,:);
%     squeeze(mean(squeeze(stdSpinData1),BVALS));
    STD = squeeze(std(squeeze(stdSpinData1),0,BVALS));

    %convert to normalized log data
    logTestCases = testCases;
    finalRes     = (STD./squeeze(mean(squeeze(stdSpinData1),BVALS)));
    errBar       = (STD);%(squeeze(mean(squeeze(varSpinData(1,:,:,:)),BVALS)));


    for ii = 1:size(STD,2)
        errorbar(logTestCases,finalRes(:,ii),errBar(:,ii),'LineWidth',1);
    end
%     set(gca,'xticklabels',{'500', '1000', '2000', '3000', '4000', '5000', '10000', '25000', '50000', '100000'});
    xlabel('nSpins'), ylabel('(Std Dev. |M_{xy}|)/mean(|M_{xy}|)');
    grid on; title(sprintf('%s - Stability across 30 trials for nSpins',criteria));
    legend("3.5*10^{-6}","4*10^{-6}","5*10^{-6}","8*10^{-6}")
end