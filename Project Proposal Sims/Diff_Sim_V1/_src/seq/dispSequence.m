function dispSequence(Sequence_Struct)
    figure;
    subplot(4,1,1);plot(Sequence_Struct.RF,'k-','LineWidth',2);title('RF Pulse');
    xlabel('time (s)'), ylabel('|B_{1}^{+}| (T)');grid on;
    subplot(4,1,2);plot(Sequence_Struct.G(1,:),'r-','LineWidth',2); title('Read Gradient');
    xlabel('time (s)'), ylabel('G_{r}(T/m)');grid on;
    subplot(4,1,3);plot(Sequence_Struct.G(2,:),'g-','LineWidth',2); title('Phase Gradient')
    xlabel('time (s)'), ylabel('G_{p}(T/m)');grid on;
    subplot(4,1,4);plot(Sequence_Struct.G(3,:),'b-','LineWidth',2); title('Slice Select Gradient'); %ylim([0,5.0e-3]);
    xlabel('time (s)'), ylabel('G_{z}(T/m)');grid on;
end