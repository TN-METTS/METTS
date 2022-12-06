% system parameter
tobj = tic2;
%MaxCoreNum= 32; % this should be revised for the computer executed 

J = 1; % coupling strength
N = 100; % number of sites in a chain
T = 0.5; % temperature
% T = [0.5, 0.4, 0.3, 0.2, 0.1]
SN = 10; % sampling number

%log
print_log = false;

% Imaginary time evolution parameters
Nkeep = 20; % bond dimension
Skeep = 1e-10;
beta = 1/T; 
dt = 1/25; % discrete time step size

% Local operators
[S,I] = getLocalSpace('Spin',1);
basis = zeros(size(I, 1), size(I, 2),3);
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

H_MPO = MPO_canonForm(H_MPO, 0, [], Skeep);
H2_MPO = mtimes_MPO (H_MPO, H_MPO, Nkeep,Nsweep);
%[H2_MPO, norm]= zipup_algo(H_MPO, H_MPO, Nkeep, Nsweep, Skeep);

fprintf("For N=%d, spin=%d, beta=%1.3f, dt=%0.4f, Heisenberg chain, calculation started\n", N, 1, beta, dt);

%       |2                   |2          
%    3  |    4            3  |    4       
% ----H_MPO{n}-----   ----H2_MPO{n}-----
%       |                    |      
%       |1                   |1   

Sz = S(:,:,2);
Szloc=cell(2, 2);
Szloc{1, 1}=I;
Szloc{1, 2}=zeros(size(I));
Szloc{2, 1}=Sz;
Szloc{2, 2} = I; 
Szloc = cell2mat(reshape(Szloc, [1, 1, size(Szloc, 1), size(Szloc, 2)]));
Sz_MPO = cell(1, N);
Sz_MPO(:) = {Szloc};
Sz_MPO{1} = Sz_MPO{1}(:,:,end,:); % choose the last components of the left leg
Sz_MPO{end} = Sz_MPO{end}(:,:,:,1); % choose the first components of the right leg
Sz_MPO = MPO_canonForm(Sz_MPO, 0, [], Skeep);
Sz2_MPO = mtimes_MPO(Sz_MPO, Sz_MPO, Nkeep,Nsweep);

MSzloc=cell(2, 2);
MSzloc{1, 1}=I;
MSzloc{1, 2}=zeros(size(I));
MSzloc{2, 1}=basis( :, :, 1)'*Sz*basis( :, :, 1);
MSzloc{2, 2} = I; 
MSzloc = cell2mat(reshape(MSzloc, [1, 1, size(MSzloc, 1), size(MSzloc, 2)]));
MSz_MPO = cell(1, N);
MSz_MPO(:) = {MSzloc};
MSz_MPO{1} = MSz_MPO{1}(:,:,end,:); % choose the last components of the left leg
MSz_MPO{end} = MSz_MPO{end}(:,:,:,1); % choose the first components of the right leg
MSz_MPO = MPO_canonForm(MSz_MPO, 0, [], Skeep);
MSz2_MPO = mtimes_MPO(MSz_MPO, MSz_MPO, Nkeep, Nsweep);


% Initialize state with random number
M = cell(1,N);
for itN = (1:N)
    M{itN} = permute(rand(1, 3),[1 3 2]);
end

%normalize
M = canonForm(M, 0, [], []);

% Use CPS_collapse to generate orthonormal basis 
M = tmp_CPS_collapse(M, basis( :, :,2),  print_log);

% malloc space to save result
mag = zeros(SN,1);
mag2 = zeros(SN,1);
E = zeros(SN, 1);
E2 = zeros(SN, 1);

% time evolution with TS_1D
for itS=(1:SN)
    [M, isright] = TS_1D(M, H, Nkeep, dt, beta/2, print_log);
%     if mod(itS+1, 2)
%         mag(itS) = real(exp_val(M, MSz_MPO, ~isright));
%         mag2(itS) = real(exp_val(M, MSz2_MPO, ~isright));
%     else
%         mag(itS) = real(exp_val(M, Sz_MPO, ~isright));
%         mag2(itS) = real(exp_val(M, Sz2_MPO, ~isright));
%     end     
    E(itS) = real(exp_val(M, H_MPO, ~isright));
    E2(itS) = real(exp_val(M, H2_MPO, ~isright));
    %M = CPS_collapse(M, size(M{1}, 3), print_log);
    [M, idxes] = tmp_CPS_collapse(M, basis( :, :, mod(itS+1, 2)+1), print_log);
    mag(itS) = sum(idxes-2, 'all');
    mag2(itS) = sum( (idxes-2)' * (idxes-2), 'all');
    if mod(itS, 10)==0
        disptime(sprintf('Iteration %d ended', itS));
    end 
end 

cutoff = 4; 
chi = (mean(mag2(cutoff:end))-mean(mag(cutoff:end))^2)./(N*T);
chi_error = (std(mag2(cutoff:end)-mag(cutoff:end).^2))/(N*T)/sqrt(SN-cutoff);

specific_heat = (mean(E2(cutoff:end))-mean(E(cutoff:end))^2)./(N*T^2);
specific_heat_error = (std(E2(cutoff:end)-E(cutoff:end).^2))/(N*T^2)/sqrt(SN-cutoff);

fprintf('T = %0.3f, susceptibility=%0.4f, specific heat = %0.4f\n', T, chi, specific_heat);
fprintf('susceptibility in [%0.4f, %0.4f] specific heat in [%0.4f, %0.4f]\n', chi-chi_error, chi+chi_error,specific_heat-specific_heat_error, specific_heat+specific_heat_error);


toc2(tobj,'-v');
chkmem;

save(sprintf('../result/mag2_T=%0.2f_SN=%d.mat', T, SN), 'mag2')
save(sprintf('../result/mag_T=%0.2f_SN=%d.mat', T, SN), 'mag')
save(sprintf('../result/E2_T=%0.2f_SN=%d.mat', T, SN), 'E2')
save(sprintf('../result/E_T=%0.2f_SN=%d.mat', T, SN), 'E')

