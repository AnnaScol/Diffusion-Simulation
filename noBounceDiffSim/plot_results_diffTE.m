clear all; clc;

bVal   = [0,100,200,300,500,600,750,1000,1500,2000,3000,4000,5000,6000,7000,8000,10000,12500,15000,20000]; %s/mm^2


folder_path = "./_data/diffTE_dt10";
d25 = load(sprintf('%s/final_res_Delta_25_50000spin.mat',folder_path)).final_res;
d20 = load(sprintf('%s/final_res_Delta_20_50000spin.mat',folder_path)).final_res;
d15 = load(sprintf('%s/final_res_Delta_15_50000spin.mat',folder_path)).final_res;
d10 = load(sprintf('%s/final_res_Delta_10_50000spin.mat',folder_path)).final_res;


figure(4);
subplot(2,4,1);plot_analysis(d25, bVal);
title('Delta = 25ms, dt = 10ms');

subplot(2,4,2);plot_analysis(d20, bVal);
title('Delta = 20ms, dt = 10ms');

subplot(2,4,3);plot_analysis(d15, bVal);
title('Delta = 15ms, dt = 10ms');

subplot(2,4,4);plot_analysis(d10, bVal);
title('Delta = 10ms, dt = 10ms');

%% Find coefficient of variation between data
mean_sig = mean(cat(3,d25,d20,d15,d10),3);
std_sig  = std(cat(3,d25,d20,d15,d10),0, 3);

coefVar = std_sig./mean_sig;
coefVar = coefVar./max(coefVar);

subplot(2,4,5); 
errorbar(bVal,abs(coefVar(:,1)),std_sig(:,1),'k.-'); 
xlabel('b (s/mm^2)'), ylabel('2.5 Norm Coefficient of Var.');grid on;
title("Var between \Delta 's, R = 3.5um");

subplot(2,4,6); 
errorbar(bVal,abs(coefVar(:,2)),std_sig(:,2),'k.-'); hold on;
xlabel('b (s/mm^2)'), ylabel('Norm Coefficient of Var.');grid on;
title("Var between \Delta 's, R = 4um");

subplot(2,4,7); 
errorbar(bVal,abs(coefVar(:,3)),std_sig(:,3),'k.-'); hold on;
xlabel('b (s/mm^2)'), ylabel('Norm Coefficient of Var.');grid on;
title("Var between \Delta 's, R = 5um");

subplot(2,4,8); 
errorbar(bVal,abs(coefVar(:,4)),std_sig(:,4),'k.-'); hold on;
xlabel('b (s/mm^2)'), ylabel('Norm Coefficient of Var.');grid on;
title("Var between \Delta 's, R = 8um");

%%
folder_path = "./_data/diffTE_dt20";
d25 = load(sprintf('%s/final_res_Delta_25_50000spin.mat',folder_path)).final_res;
d20 = load(sprintf('%s/final_res_Delta_20_50000spin.mat',folder_path)).final_res;
d15 = load(sprintf('%s/final_res_Delta_15_50000spin.mat',folder_path)).final_res;
d10 = load(sprintf('%s/final_res_Delta_10_50000spin.mat',folder_path)).final_res;

figure(5);
subplot(2,4,1);plot_analysis(d25, bVal);
title('Delta = 25ms, dt = 20ms');

subplot(2,4,2);plot_analysis(d20, bVal);
title('Delta = 20ms, dt = 20ms');

subplot(2,4,3);plot_analysis(d15, bVal);
title('Delta = 15ms, dt = 20ms');

subplot(2,4,4);plot_analysis(d10, bVal);
title('Delta = 10ms, dt = 20ms');

%% Find coefficient of variation between data
mean_sig2 = mean(cat(3,d25,d20,d15,d10),3);
std_sig2  = std(cat(3,d25,d20,d15,d10),0, 3);

coefVar2 = std_sig2./mean_sig2;
coefVar2 = coefVar2./max(coefVar2);

subplot(2,4,5); 
errorbar(bVal,abs(coefVar(:,1)),std_sig(:,1),'k.-'); 
xlabel('b (s/mm^2)'), ylabel('Norm Coefficient of Var.');grid on;
title("Var between \Delta 's, R = 3.5um");

subplot(2,4,6); 
errorbar(bVal,abs(coefVar(:,2)),std_sig(:,2),'k.-'); hold on;
xlabel('b (s/mm^2)'), ylabel('Norm Coefficient of Var.');grid on;
title("Var between \Delta 's, R = 4.5um");

subplot(2,4,7); 
errorbar(bVal,abs(coefVar(:,3)),std_sig(:,3),'k.-'); hold on;
xlabel('b (s/mm^2)'), ylabel('Norm Coefficient of Var.');grid on;
title("Var between \Delta 's, R = 5um");

subplot(2,4,8); 
errorbar(bVal,abs(coefVar(:,4)),std_sig(:,4),'k.-'); hold on;
xlabel('b (s/mm^2)'), ylabel('Norm Coefficient of Var.');grid on;
title("Var between \Delta 's, R = 8um");

%% Plot diff and CV
figure(6)

subplot(3,4,1); 
errorbar(bVal,abs(coefVar(:,1)),std_sig(:,1),'k.-'); 
xlabel('b (s/mm^2)'), ylabel('Normalised CV');grid on;
title('R = 3.5um - dt = 10ms');

subplot(3,4,2); 
errorbar(bVal,abs(coefVar(:,2)),std_sig(:,2),'k.-'); hold on;
xlabel('b (s/mm^2)'), ylabel('Normalised CV');grid on;
title('R = 4um - dt = 10ms');

subplot(3,4,3); 
errorbar(bVal,abs(coefVar(:,3)),std_sig(:,3),'k.-'); hold on;
xlabel('b (s/mm^2)'), ylabel('Normalised CV');grid on;
title('R = 5um - dt = 10ms');

subplot(3,4,4); 
errorbar(bVal,abs(coefVar(:,4)),std_sig(:,4),'k.-'); hold on;
xlabel('b (s/mm^2)'), ylabel('Normalised CV');grid on;
title('R = 8um - dt = 10ms');

subplot(3,4,5); 
errorbar(bVal,abs(coefVar2(:,1)),std_sig2(:,1),'k.-'); 
xlabel('b (s/mm^2)'), ylabel('Normalised CV');grid on;
title('R = 3.5um - dt = 20ms');

subplot(3,4,6); 
errorbar(bVal,abs(coefVar2(:,2)),std_sig2(:,2),'k.-'); hold on;
xlabel('b (s/mm^2)'), ylabel('Normalised CV');grid on;
title('R = 4um - dt = 20ms');

subplot(3,4,7); 
errorbar(bVal,abs(coefVar2(:,3)),std_sig2(:,3),'k.-'); hold on;
xlabel('b (s/mm^2)'), ylabel('Normalised CV');grid on;
title('R = 5um - dt = 20ms');

subplot(3,4,8); 
errorbar(bVal,abs(coefVar2(:,4)),std_sig2(:,4),'k.-'); hold on;
xlabel('b (s/mm^2)'), ylabel('Normalised CV');grid on;
title('R = 8um - dt = 20ms');


subplot(3,4,9); 
errorbar(bVal,abs(coefVar2(:,1)-coefVar(:,1)),std_sig2(:,1)-std_sig(:,1),'k.-'); 
xlabel('b (s/mm^2)'), ylabel('Normalised CV');grid on;
title('R = 3.5um Diff');ylim([0 0.008])

subplot(3,4,10); 
errorbar(bVal,abs(coefVar2(:,2)-coefVar(:,2)),std_sig2(:,2)-std_sig(:,2),'k.-'); hold on;
xlabel('b (s/mm^2)'), ylabel('Normalised CV');grid on;
title('R = 4um Diff');ylim([0 0.008])

subplot(3,4,11); 
errorbar(bVal,abs(coefVar2(:,3)-coefVar(:,3)),std_sig2(:,3)-std_sig(:,3),'k.-'); hold on;
xlabel('b (s/mm^2)'), ylabel('Normalised CV');grid on;
title('R = 5um Diff'); ylim([0 0.1])

subplot(3,4,12); 
errorbar(bVal,abs(coefVar2(:,4)-coefVar(:,4)),std_sig2(:,4)-std_sig(:,4),'k.-'); hold on;
xlabel('b (s/mm^2)'), ylabel('Normalised CV');grid on;
title('R = 8um Diff');ylim([0 1])
