function plot_analysis(final_results, b)
    hold on;
    for ii = 1:size(final_results,2)
        plot(b,abs(final_results(:,ii)),'.-','LineWidth',2);
        xlabel('b (s/mm^2)'), ylabel('|M_{xy}|');grid on;
    end
    legend("2.5*1e-6","3.5*1e-6","4*1e-6","5*1e-6")
    hold off;
end
