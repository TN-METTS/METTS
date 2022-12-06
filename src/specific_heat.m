% system parameter

%MaxCoreNum= 32; % this should be revised for the computer executed 

J = 1; % coupling strength
N = 100; % number of sites in a chain
%T = 0.5; % temperature
Ts = [0.5, 0.4, 0.3, 0.2, 0.1];
SN = 100; % sampling number
step_num = 6; 
%log
print_log = false;

% Imaginary time evolution parameters
Nkeeps = [40, 40, 50, 60, 60]; % bond dimension
Skeep = 1e-10;
dt = 1/25; % discrete time step size

% Local operators
[S,I] = getLocalSpace('Spin',1);
basis = zeros(size(I, 1), size(I, 2),3);
basis( :, :, 1)= [exp(2*pi/3*1i), 1,  exp(-2*pi/3*1i); -1, -1,-1; -exp(-2*pi/3*1i), -1, -exp(2*pi/3*1i)]/sqrt(3);
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

[H_MPO, lognorm_H_MPO] = MPO_canonForm(H_MPO, 0, [], Skeep);
%H_MPO{1} = H_MPO{1}.*exp(lognorm_H_MPO);
H2_MPO = mtimes_MPO (H_MPO, H_MPO, 25, Nsweep);
%[H2_MPO, norm]= zipup_algo(H_MPO, H_MPO, Nkeep, Nsweep, Skeep);


%       |2                   |2          
%    3  |    4            3  |    4       
% ----H_MPO{n}-----   ----H2_MPO{n}-----
%       |                    |      
%       |1                   |1   



for test = (1:5)
    
    T = Ts(test);
    Nkeep = Nkeeps(test);
    beta = 1/T; 
    disptime(sprintf('For N=%d, spin=%d, beta=%1.3f, dt=%0.4f, Heisenberg chain, calculation started\n', N, 1, beta, dt));
    tobj = tic2;
    % malloc space to save result
    E = zeros(SN, step_num);
    E2 = zeros(SN, step_num);
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
            E(itS, step) = real(exp_val(M, H_MPO, ~isright)).*exp(lognorm_H_MPO);
            E2(itS, step) = real(exp_val(M, H2_MPO, ~isright)).*exp(2*lognorm_H_MPO);    
        
            [M, idxes] = tmp_CPS_collapse(M, basis( :, :, mod(step+1, 2)+1), print_log);
            
        end 
%         [M, isright] = TS_1D(M, H, Nkeep, dt, beta/2, print_log);
%         E(itS) = real(exp_val(M, H_MPO, ~isright)).*exp(lognorm_H_MPO);
%         E2(itS) = real(exp_val(M, H2_MPO, ~isright)).*exp(2*lognorm_H_MPO);    
        
        if mod(itS, 5)==0
            [Cv, Cv_error] =  cal_specific_heat(E(1:itS, :), E2(1:itS, :), N, T);
            
            fprintf('T = %0.3f, specific heat = %0.4f\n', T,  Cv(end));
            fprintf('specific heat in [%0.4f, %0.4f]\n',Cv(end)-Cv_error(end), Cv(end)+Cv_error(end));
            disptime(sprintf('Iteration %d ended', itS));
        end 
    end 

    [Cv, Cv_error] =  cal_specific_heat(E(1:end, :), E2(1:end, :), N, T);
    
    fprintf('T = %0.3f, specific heat = %0.4f\n', T,  Cv(end));
    fprintf('specific heat in [%0.4f, %0.4f]\n',Cv(end)-Cv_error(end), Cv(end)+Cv_error(end));
    
    
    toc2(tobj,'-v');
    chkmem;
    save(sprintf('../result/specific_heat_T=%0.2f_SN=%d_dt=%0.4f_Nkeep=%d.mat', T, SN, dt, Nkeep), 'E', 'E2')
    disptime('Save succeed')
    break; 
end 

function [Cv, Cv_error] = cal_specific_heat(E, E2, N, T)
    SN = size(E, 1);
    Cv = (mean(E2, 1)-mean(E, 1).^2);
    Cv_error = Cv.^2.*sqrt(2/(SN-1))./(N*T^2);%(std(E2, 1)-std(E,1).^2)./(N*T^2)./sqrt(SN-1);
    Cv = Cv./(N*T^2);
    
    
end 

