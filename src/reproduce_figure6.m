% This script reproduces Figure 6 in the paper. 
% compare convergence of Energy according to the basis chosed for the CPS collapse
tobj = tic2;
MaxCoreNum= 32; % this should be revised for the computer executed 

N = 50; % number of sites in a chain
beta = 1; 
SN =500; % sampling number

%log
print_log = false;

% Imaginary time evolution parameters
Nkeep = 40; % bond dimension
dt = 1/20; % discrete time step size

% Local operators
[S,I] = getLocalSpace('Spin',1);
basis = zeros(size(I, 1), size(I, 2), 2);
basis( :, :, 1)= [exp(2*pi/3*1i), -1, -exp(-2*pi/3*1i); 1, -1, -1; exp(-2*pi/3*1i), -1, -exp(2*pi/3*1i)]/sqrt(3);
[basis(:, :, 2), ~] = eig(S(:,:,2)); %basis of Sz 

%Set gates
% nearest-neighbor interaction terms
H = cell(1, N-1);
H(:) = {J*contract(S,3,3,permute(conj(S),[2 1 3]),3,3)};
%       2      4      [legs 1 and 2 are for site n;
%       |      |       legs 3 and 4 are for site n+1]
%      [ Hs{n}  ]
%       |      |
%       1      3


% generate operators to get expectation value
% Set H as MPO
Hloc = cell(5,5);
Hloc(:) = {zeros(size(I))};
Hloc{1,1} = I;
Hloc{2,1} = S(:,:,1);
Hloc{3,1} = S(:,:,2);
Hloc{4,1} = S(:,:,3);
Hloc{5,2} = J*S(:,:,1)';
Hloc{5,3} = J*S(:,:,2)';
Hloc{5,4} = J*S(:,:,3)';
Hloc{end,end} = I;
Hloc = cell2mat(reshape(Hloc,[1 1 size(Hloc,1) size(Hloc,2)]));

H_MPO = cell(1,N);
H_MPO(:) = {Hloc};
H_MPO{1} = H_MPO{1}(:,:,end,:); % choose the last components of the left leg
H_MPO{end} = H_MPO{end}(:,:,:,1); % choose the first components of the right leg

[H_MPO, lognorm_H_MPO] = MPO_canonForm(H_MPO, 0, [], 1e-8);
fprintf("For N=%d, spin=%d, beta=%1.3f, dt=%0.4f, Heisenberg chain, calculation started\n", N, 1, beta, dt);

Step_num = 10; 
Ms=cell(1,SN);

% case 1 
z_only = zeros(SN, Step_num);
for itS=(1:SN)
    % save initial M for case 2 and case 3 
    M = cell(1,N);
    for itN = (1:N)
        M{itN} = permute(rand(1, 3),[1 3 2]);
    end
    %normalize
    M = canonForm(M, 0, [], []);        
    M = tmp_CPS_collapse(M, basis( :, :,2),  print_log);
    Ms{itS}=M;

    for step= (1:Step_num)
        % Use CPS_collapse to generate orthonormal basis 
        [M, isright] = TS_1D(M, H, Nkeep, dt, beta/2, print_log);
        z_only(itS, step) = exp_val(M, H_MPO);
        M = tmp_CPS_collapse(M, basis( :, :,2),  print_log);

    end 
end 
z_only = z_only./N.*exp(lognorm_H_MPO);
z_only_std= std(z_only, 0, 1);
disptime('z only ended')

% case 2 
random = zeros(SN, Step_num);
for itS=(1:SN)
    M=Ms{itS};
    for step= (1:Step_num)
        % Use CPS_collapse to generate orthonormal basis 
        [M, isright] = TS_1D(M, H, Nkeep, dt, beta/2, print_log);
        random(itS, step) = exp_val(M, H_MPO);
        rand_basis=basis_random_axis(S);
        M = tmp_CPS_collapse(M, rand_basis,  print_log);
    end 
end 
random = real(random)./N.*exp(lognorm_H_MPO);
random_std= std(random, 0, 1);
disptime('random ended')
% case 3 
max_mixed = zeros(SN, Step_num);
for itS=(1:SN)
    M=Ms{itS};
    for step= (1:Step_num)
        % Use CPS_collapse to generate orthonormal basis 
        [M, isright] = TS_1D(M, H, Nkeep, dt, beta/2, print_log);
        max_mixed(itS, step) = exp_val(M, H_MPO);
        M = tmp_CPS_collapse(M, basis( :, :, mod(step+1, 2)+1),  print_log);

    end 
end 
max_mixed = real(max_mixed)./N.*exp(lognorm_H_MPO);
max_mixed_std= std(max_mixed, 0, 1);
disptime('maximaly mixed ended')

toc2(tobj,'-v');
chkmem;
% save result 
save(sprintf('../result/E_measurement_SN=%d.mat', SN), 'z_only', 'max_mixed', 'random')
disptime('Save succeed');


% % plot the result 
% figure; 
% hold on ;
% errorbar((1:Step_num), mean(z_only, 1), z_only_std./sqrt(SN),'-s', 'MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red', 'Color', 'red');
% errorbar((1:Step_num), mean(random, 1), random_std./sqrt(SN), '-o', 'MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue', 'Color', 'blue');
% errorbar((1:Step_num), mean(max_mixed, 1), max_mixed_std./sqrt(SN), '-^', 'MarkerSize',5,'MarkerEdgeColor','green','MarkerFaceColor','green', 'Color', 'green');
% hold off;
% set(gca, 'FontSize', 13, 'LineWidth', 1);
% legend({'Z Only', 'Random', 'Maximally Mixed'})
% grid on;
% xlabel('Step Number');
% ylabel('Energy per Site');
% savefig(sprintf('../result/Figure6_SN=%d.fig', SN))



function rand_basis=basis_random_axis(S)
    dim =size(S, 3);
    axis = rand( dim, 1);
    axis = axis./norm(axis(:));
    new_S=zeros(size(S));
    new_S(:,:,3)=S(:,:,2);
    new_S(:,:,1)=(S(:,:,1)+S(:,:,3))/sqrt(2);
    new_S(:,:,2)=(S(:,:,1)-S(:,:,3))/sqrt(2)/1i;
    mat = contract(new_S, 3, 3, axis, 2, 1);
    [rand_basis, ~] = eig((mat+mat')/2);
end 
