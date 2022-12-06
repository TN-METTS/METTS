SN = 10 ; % sampling number 
Ts = [0.1, 0.2, 0.3, 0.4, 0.5];
Cvs = zeros(1, 5);
Cv_errors = zeros(1, 5);

dt = 1/25;
Nkeeps = [40, 40, 40, 50, 60, 60];
N=100;
for itN=(1:5) 
    T = Ts(itN);  Nkeep = Nkeeps(itN);
    load(sprintf('../result/specific_heat_T=%0.2f_SN=%d_dt=%0.4f_Nkeep=%d.mat', T, SN, dt, Nkeep));

    SN = size(E, 1);
    step_num = size(E, 2);
    Cvs(itN) = (mean(E2(:, 3:step_num),  'all')-mean(E(:, 3:step_num), 'all').^2)./N;
    Cv_errors(itN) = Cvs(itN).^2.*sqrt(2/(SN*step_num-1))./T^2; 
    Cvs(itN) = Cvs(itN)./(T^2);
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