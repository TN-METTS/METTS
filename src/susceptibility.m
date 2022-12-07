%MaxCoreNum= 32; % this should be revised for the computer executed 

J = 1; % coupling strength
N = 100; % number of sites in a chain
Ts = [0.5, 0.4, 0.3, 0.2, 0.1];
SNs = [2, 2, 2, 2, 2]; % sampling number
steps=[50, 50, 40, 30, 20];
%step_num = 50; 

%log
print_log = false;

% Imaginary time evolution parameters
Nkeeps = [40, 40, 40, 50, 60];% bond dimension
Skeep = 1e-10;
dts=[1/25,1/25,1/20,1/20,1/15];

% Local operators
[S,I] = getLocalSpace('Spin',1);
basis = zeros(size(I, 1), size(I, 2),3);
% Maximally mixed states 
basis( :, :, 1)= [exp(2*pi/3*1i), 1,  exp(-2*pi/3*1i); -1, -1,-1; -exp(-2*pi/3*1i), -1, -exp(2*pi/3*1i)]/sqrt(3);
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

[H_MPO, lognorm_H_MPO] = MPO_canonForm(H_MPO, 0, [], Skeep);
H2_MPO = mtimes_MPO (H_MPO, H_MPO, 40, Nsweep);

%       |2                   |2          
%    3  |    4            3  |    4       
% ----H_MPO{n}-----   ----H2_MPO{n}-----
%       |                    |      
%       |1                   |1   

% generate S^z_tot operator in MPO form 
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
[Sz_MPO, lognorm_Sz_MPO] = MPO_canonForm(Sz_MPO, 0, [], Skeep);

Sz2_MPO = mtimes_MPO(Sz_MPO, Sz_MPO, 40, Nsweep);



for test = (5:5)%(1:5)
    dt = dts(test);
    T = Ts(test);
    Nkeep = Nkeeps(test);
    SN = SNs(test); step_num = steps(test);
    beta = 1/T; 
    disptime(sprintf('For N=%d, spin=%d, beta=%1.3f, dt=%0.4f, Heisenberg chain, calculation started', N, 1, beta, dt));
    tobj = tic2;
    % malloc space to save result
    mz = zeros(SN, step_num);
    mz2 = zeros(SN, step_num);
    mz_MPO = zeros(SN, step_num);
    mz2_MPO = zeros(SN, step_num);
    for itS=(1:SN)
        % Initialize state with random number
        M = cell(1,N); 
        for itN = (1:N)
            M{itN} = permute(rand(1, 3),[1 3 2]);
        end
        %normalize
        M = canonForm(M, 0, [], []);
        
        % Use CPS_collapse to generate orthonormal basis 
        M = tmp_CPS_collapse(M, basis( :, :,2),  print_log);
        
        for step = (1:step_num)
            % time evolution with TS_1D
            [M, isright] = TS_1D(M, H, Nkeep, dt, beta/2, print_log);
            mz_MPO(itS, step) = real(exp_val(M, Sz_MPO, ~isright)).*exp(lognorm_Sz_MPO);
            mz2_MPO(itS, step) = real(exp_val(M, Sz2_MPO, ~isright)).*exp(2*lognorm_Sz_MPO);    
        
            [M, idxes] = tmp_CPS_collapse(M, basis( :, :, mod(step+1, 2)+1), print_log);
            m = D(idxes); 
            mz(itS, step) = sum(m, 'all');
            mz2(itS, step) = mz(itS, step)^2; %sum(m*m', 'all');
        
        if mod(step, 5)==0
            [chi, chi_error] = cal_susceptibility(mz(1:itS, 2:2:end), mz2(1:itS, 2:2:end), N, T);
            [chi_MPO, chi_MPO_error] = cal_susceptibility(mz_MPO(1:itS, 3:end), mz2_MPO(1:itS, 3:end), N, T);
            
            disptime(sprintf('T = %0.3f, susceptibility = %0.4f in [%0.4f, %0.4f]', T,  chi(end), chi(end)-chi_error(end), chi(end)+chi_error(end)));            
            disptime(sprintf('T = %0.3f, susceptibility = %0.4f in [%0.4f, %0.4f] by MPO', T,  chi_MPO(end), chi_MPO(end)-chi_MPO_error(end), chi_MPO(end)+chi_MPO_error(end)));            
            disptime(sprintf('Iteration %d ended', step));
            save([PATH, sprintf('tmp_susceptibility_T=%0.2f_SN=%d_dt=%0.4f_Nkeep=%d.mat', T, step, dt, Nkeep)], 'mz', 'mz2', 'mz_MPO', 'mz2_MPO')
        
        end 
        end
    end 

    [chi, chi_error] = cal_susceptibility(mz(1:itS, 2:2:end), mz2(1:itS, 2:2:end), N, T);
    [chi_MPO, chi_MPO_error] = cal_susceptibility(mz_MPO(1:itS, 1:end), mz2_MPO(1:itS, 1:end), N, T);
            
    disptime(sprintf('T = %0.3f, susceptibility = %0.4f in [%0.4f, %0.4f]', T,  chi(end), chi(end)-chi_error(end), chi(end)+chi_error(end)));            
    disptime(sprintf('T = %0.3f, susceptibility = %0.4f in [%0.4f, %0.4f] by MPO', T,  chi_MPO(end), chi_MPO(end)-chi_MPO_error(end), chi_MPO(end)+chi_MPO_error(end)));            
            
    
    toc2(tobj,'-v');
    chkmem;
    save([PATH, sprintf('susceptibility_T=%0.2f_SN=%d_dt=%0.4f_Nkeep=%d.mat', T, step_num, dt, Nkeep)], 'mz', 'mz2', 'mz_MPO', 'mz2_MPO')

    disptime('Save succeed')
end 


