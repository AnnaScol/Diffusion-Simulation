clear all;

num_sets = 7;
nB = 17;
all_data = zeros(num_sets,nB,4);

bVal   = [0,100,200,300,500,600,750,1000,1500,2000,3000,4000,5000,6000,7000,8000,10000,12500,15000,20000]; %s/mm^2

bounce_data = load("70000spins.mat");
old_data = load("fina_res70000old.mat");
%%
figure(1);
subplot(1,3,1);plot_signal_vs_b(bounce_data.final_res, bVal, 1);
subplot(1,3,2);plot_signal_vs_b(old_data.final_res_old, bVal, 1);
subplot(1,3,3);
figure(1); hold on;
diff_res = bounce_data.final_res - old_data.final_res_old;
for ii = 1:size(diff_res,2)
    plot(bVal,abs(diff_res(:,ii)),'.-','LineWidth',2);
    xlabel('b (s/mm^2)'), ylabel('|M_{xy}|');grid on;
    title('Difference Diffusion Attenuation');
end
legend("3.5*1e-6","4*1e-6","5*1e-6","8*1e-6")
hold off;



