function plot_signal_vs_b_results(nRadii, final_results, b, nFig)

    for i = 1:nRadii
        figure(nFig);
        subplot(2,1,1);
        log_vals = -log(abs(final_results(:,i))./(b*1e-6)');
        log_vals(1) = 1;
        plot(b,log_vals,'.-','LineWidth',2);
        xlabel('b (s/mm^2)'), ylabel('ln(|M_{xy}|/b)');grid on;
        title('Log Diffusion Attenuation (Spins 3/11ms)');
        hold on
        drawnow

        figure(nFig);
        subplot(2,1,2);
        plot(b,abs(final_results(:,i)),'.-','LineWidth',2);
        xlabel('b (s/mm^2)'), ylabel('|M_{xy}|');grid on;
        title('Absolute of Diffusion Attenuation (Spins 3/11ms)');
        hold on
        drawnow
    end
    % figure(2);
    subplot(2,1,1);
    legend("1*1e-6","5*1e-6","10*1e-6","20*1e-6","40*1e-6","80*1e-6")
    subplot(2,1,2);
    legend("1*1e-6","5*1e-6","10*1e-6","20*1e-6","40*1e-6","80*1e-6")

end