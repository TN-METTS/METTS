
SN=500;
% load data 
load(sprintf('../result/E_measurement_SN=%d.mat', SN));

% calculate std 
z_only_std= std(z_only, 0, 1);
random_std= std(random, 0, 1);
max_mixed_std= std(max_mixed, 0, 1);

Step_num = 10; 
% plot the result 
figure; 
hold on ;
errorbar((1:Step_num), mean(z_only, 1), z_only_std./sqrt(SN),'-s', 'MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red', 'Color', 'red');
errorbar((1:Step_num), mean(random, 1), random_std./sqrt(SN), '-o', 'MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue', 'Color', 'blue');
errorbar((1:Step_num), mean(max_mixed, 1), max_mixed_std./sqrt(SN), '-^', 'MarkerSize',5,'MarkerEdgeColor','green','MarkerFaceColor','green', 'Color', 'green');
hold off;
set(gca, 'FontSize', 13, 'LineWidth', 1);
legend({'Z Only', 'Random', 'Maximally Mixed'})
grid on;
xlabel('Step Number');
ylabel('Energy per Site');