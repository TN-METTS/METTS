% generate figures for specific heat and susceptibility from saved result. 
%
%
% Written by M.Kim (Dec.07,2022)


steps = [50,50,30,30,20]; % sampling number 
SNs = [30,30,30,30,20]; % sampling number 
Ts = [0.5, 0.4, 0.3, 0.2, 0.1];
Nkeeps = [40, 40, 40, 50, 60];
dts=[1/25,1/25,1/20,1/20,1/15];

Cvs = zeros(1, 5);
Cv_errors = zeros(1, 5);
chis = zeros(2, 5);
chi_errors = zeros(2, 5);

N=100;

for itN=(1:5) 
    T = Ts(itN);  Nkeep = Nkeeps(itN); dt = dts(itN);
    SN = SNs(itN);    step = steps(itN);
    load([PATH, sprintf('specific_heat_T=%0.2f_SN=%d_dt=%0.4f_Nkeep=%d.mat', T, SN, dt, Nkeep)]);
    load([PATH, sprintf('susceptibility_T=%0.2f_SN=%d_dt=%0.4f_Nkeep=%d.mat', T, step, dt, Nkeep)], 'mz', 'mz2', 'mz_MPO', 'mz2_MPO')
    [Cvs(itN), Cv_errors(itN)] =  cal_specific_heat(E, E2, N, T);
    [chis(1,itN), chi_errors(1, itN)] = cal_susceptibility(mz(:, 2:2:end), mz2(:, 2:2:end), N, T);
    [chis(2,itN), chi_errors(2, itN)] = cal_susceptibility(mz_MPO(:, 5:end), mz2_MPO(:, 5:end), N, T);
end 

figure; 
hold on ;
errorbar(Ts, Cvs, Cv_errors,  '-o', 'MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue', 'Color', 'blue');
hold off;
set(gca, 'FontSize', 13, 'LineWidth', 1);
legend({'METTS'})
%grid on;
xlabel('Temperature T');
ylabel('Specific Heat C_{v}');
savefig([PATH, sprintf('specific_heat_last_SN=%d_dt=%0.4f_Nkeep=%d.fig', SN, dt, Nkeep)])


figure; 
hold on ;
errorbar(Ts, chis(2, :), chi_errors(2, :),  '-o', 'MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red', 'Color', 'red');
%errorbar(Ts, chis(1, :), chi_errors(1, :),  '-o', 'MarkerSize',5,'MarkerEdgeColor','green','MarkerFaceColor','green', 'Color', 'green');

hold off;
set(gca, 'FontSize', 13, 'LineWidth', 1);
legend({ 'MPO'}) %'direct measurement',
%grid on;
xlabel('Temperature T');
ylabel('suscetbility \chi');
savefig([PATH, sprintf('susceptibility_last_SN=%d_dt=%0.4f_Nkeep=%d.fig', SN, dt, Nkeep)])


