function plot_signal_vs_b(final_results, b,nfig)
    figure(nfig); hold on;
    for ii = 1:size(final_results,2)
        plot(b,abs(final_results(:,ii)),'.-','LineWidth',2);
        xlabel('b (s/mm^2)'), ylabel('|M_{xy}|');grid on;
        title('|Diffusion Attenuation|');
    end
    legend("3.5*1e-6","4*1e-6","5*1e-6","8*1e-6")
    hold off;
end

