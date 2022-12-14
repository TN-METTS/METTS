% system parameter

%MaxCoreNum= 32; % this should be revised for the computer executed 

J = 1; % coupling strength
N = 100; % number of sites in a chain
%T = 0.3; % temperature
%Ts = [0.5, 0.4,  0.2, 0.1]; % erase 0.3
Ts = [0.5];
SN = 20; % sampling number

%log
print_log = false;

% Imaginary time evolution parameters
Nkeep = 30; % bond dimension
Skeep = 1e-16;
dt = 1/20; % discrete time step size

% Local operators
[S,I] = getLocalSpace('Spin',1);
basis = zeros(size(I, 1), size(I, 2), 2);
basis( :, :, 1) =[exp(2*pi/3*1i), -1, -exp(-2*pi/3*1i); 1, -1, -1; exp(-2*pi/3*1i), -1, -exp(2*pi/3*1i)]'/sqrt(3);
[basis(:, :, 2), D] = eig(S(:,:,2)); %basis of Sz 
D = diag(D);

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
Nsweep = 3; 
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

%[H_MPO, lognorm_H_MPO] = MPO_canonForm(H_MPO, 0, [], Skeep);
%H_MPO{1} = H_MPO{1}.*exp(lognorm_H_MPO);

for T=Ts
    beta = 1/T; 
    
    fprintf("For N=%d, spin=%d, beta=%1.3f, dt=%0.4f, Heisenberg chain, calculation started\n", N, 1, beta, dt);
    tobj = tic2;

    % Initialize state with random number
    M = cell(1,N);
    for itN = (1:N)
        M{itN} = permute(rand(1, 3),[1 3 2]);
    end
    
    %normalize
    M = canonForm(M, 0, [], []);
    
    % malloc space to save result
    mag_z= zeros(SN,1);
    mag2_z = zeros(SN,1);
    mag_m = zeros(SN,1);
    mag2_m = zeros(SN,1);
    % time evolution with TS_1D
    for itS=(1:SN)
        %[M, ~] = TS_1D(M, H, Nkeep, dt, beta/2, print_log);
        size(M{end})
        M = time_evolution(H_MPO, beta/2, M, 100, 0, Nkeep, Skeep);
        [M, idxes] = tmp_CPS_collapse(M, basis( :, :, 1), print_log);
        m = D(idxes)';
        mag_m(itS) = sum(m, 'all');
        mag2_m(itS) = sum(m' * m, 'all');
    
        [M, ~] = TS_1D(M, H, Nkeep, dt, beta/2, print_log);
        [M, idxes] = tmp_CPS_collapse(M, basis( :, :, 2), print_log);
        m = D(idxes)';
        mag_z(itS) = sum(m, 'all');
        mag2_z(itS) = sum(m' * m, 'all');

        if mod(itS, 10)==0
            disptime(sprintf('Iteration %d ended', itS));
        end 
    end 
    
    cutoff = 4; 
    chi_z = (mean(mag2_z(cutoff:end))-mean(mag_z(cutoff:end))^2)./(N*T);
    chi_m = (mean(mag2_m(cutoff:end))-mean(mag_m(cutoff:end))^2)./(N*T);
    % chi_error = (std(mag2_z(cutoff:end)-mag_z(cutoff:end).^2))/(N*T)/sqrt(SN-cutoff);
    
    fprintf('T = %0.3f, susceptibility=%0.4f\n', T, chi_z);
    % fprintf('susceptibility in [%0.4f, %0.4f] \n', chi-chi_error, chi+chi_error);
    
    
    toc2(tobj,'-v');
    chkmem;
    
    save(sprintf('../result/susceptibility_taylor_T=%0.2f_SN=%d_dt=%0.4f.mat', T, SN, dt), 'mag2_z', 'mag_z', "mag_m", "mag2_m")
end
