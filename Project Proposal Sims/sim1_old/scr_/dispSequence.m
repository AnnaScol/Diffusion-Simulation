function dispSequence(gradAmp,rfPulse)
    figure;
    subplot(4,1,1);plot(rfPulse,'k-','LineWidth',2);title('RF Pulse');
    xlabel('time (s)'), ylabel('|B_{1}^{+}| (T)');grid on;
    subplot(4,1,2);plot(gradAmp(1,:),'r-','LineWidth',2); title('Read Gradient');
    xlabel('time (s)'), ylabel('G_{r}(T/m)');grid on;
    subplot(4,1,3);plot(gradAmp(2,:),'g-','LineWidth',2); title('Phase Gradient')
    xlabel('time (s)'), ylabel('G_{p}(T/m)');grid on;
    subplot(4,1,4);plot(gradAmp(3,:),'b-','LineWidth',2); title('Slice Select Gradient')
    xlabel('time (s)'), ylabel('G_{z}(T/m)');grid on;
end