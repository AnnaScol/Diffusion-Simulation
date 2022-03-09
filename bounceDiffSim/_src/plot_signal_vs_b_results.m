function plot_signal_vs_b_results(final_results, b)
     figure; hold on;
    for ii = 1:size(final_results,2)
        plot(b,abs(final_results(:,ii)),'.-','LineWidth',2);
        xlabel('b (s/mm^2)'), ylabel('|M_{xy}|');grid on;
        title('Absolute of Diffusion Attenuation');
    end
    legend("3.5*1e-6","4*1e-6","5*1e-6","8*1e-6")
    hold off;
end
