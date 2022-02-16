clear all; 
clc;

string1 = "compare_result1_40x2000spins.mat";%"compare_result1.mat"
string2 = "compare_result2_40x2000spins.mat";%"compare_result2.mat"

res1 = load(string1);
res2 = load(string2);

ldelta = 0.011; %ms was 0.04
sdelta = 0.003; %ms was 0.02

G = ([0,5,7.5,10,12.5,15,17.5,20,22,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,150]*1e-3)*15; %mT
nG = length(G);
mFinalVect = zeros(nG,2); %variable to hold the final magnetization calculated for each position
gamma = 2*pi*42.577*10^6;
b = gamma^2*G.^2*sdelta^2*(ldelta-sdelta/3);


for i = 1:size(res1.final_res,2) %gets the number of radii values
    sig_difference = (res1.final_res(:,i))  - (res2.final_res(:,i));
    logSig_difference = (-log(abs(res1.final_res(:,i))./(b*1e-6)'))  - (-log(abs(res2.final_res(:,i))./(b*1e-6)'));
    
    figure(4);
    subplot(2,1,1);
    log_vals = (-log(abs(res1.final_res(:,i))./(b*1e-6)'))  - (-log(abs(res2.final_res(:,i))./(b*1e-6)')) ;
    log_vals(1) = 1;
    plot(b*1e-6,log_vals,'.-','LineWidth',2);
    xlabel('b (s/mm^2)'), ylabel('ln(|M_{xy}|/b)');grid on;
    title('Log Diffusion Attenuation Difference between Two Sets of 4000 Spins');
    hold on
    drawnow

    figure(4);
    subplot(2,1,2);
    plot(b*1e-6,abs(sig_difference),'.-','LineWidth',2);
    xlabel('b (s/mm^2)'), ylabel('|M_{xy}|');grid on;
    title('Absolute of Diffusion Attenuation Difference between Two Sets of 4000 Spins');
    hold on
    drawnow
end
% figure(2);
subplot(2,1,1);
legend("2*1e-6","5*1e-6","10*1e-6","20*1e-6","40*1e-6","80*1e-6")
subplot(2,1,2);
legend("2*1e-6","5*1e-6","10*1e-6","20*1e-6","40*1e-6","80*1e-6")


%% Cumulative difference between signals
cumulative_mean = zeros(1,size(res1.final_res,2));

for i = 1:size(res1.final_res,2) %gets the number of radii values
    sig_difference = (res1.final_res(:,i))  - (res2.final_res(:,i));
    cumulative_sum = cumsum(abs(sig_difference));
    cumulative_mean(i) = mean(cumulative_sum);
  
    figure(5);
    plot(b*1e-6,cumulative_sum,'.-','LineWidth',2);
    xlabel('b (s/mm^2)'), ylabel('|M_{xy}|');grid on;
    title('Absolute of Diffusion Attenuation Difference between Two Sets of 4000 Spins');
    hold on
    drawnow
end
% figure(2);
legend("2*1e-6","5*1e-6","10*1e-6","20*1e-6","40*1e-6","80*1e-6")
disp(cumulative_mean);


