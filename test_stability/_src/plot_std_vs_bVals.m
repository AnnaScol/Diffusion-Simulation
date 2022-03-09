function plot_std_vs_bVals(final_results, b,vVar)
     figure; hold on;
%      vVar = zeros(1,size(final_results,2));
     
    for ii = 1:size(final_results,2)
        plot(b,final_results(:,ii),'.-','LineWidth',2);
        xlabel('b (s/mm^2)'), ylabel('Std Dev. of |M_{xy}|');grid on;
        title('Standard Deviation across 14 trials');
    end
    legend(sprintf("3.5*10^{-6}  var = %d",vVar(1)),sprintf("4*10^{-6}  var = %d",vVar(2)), ...
                sprintf("5*10^{-6} var = %d",vVar(3)),sprintf("8*10^{-6}  var = %d",vVar(4)))
    hold off;
end
