function equality = Test_TS_1D()

% This function test sweeping function "MPO_multiplication_sweep" 
% B = MPO_multiplication_sweep(W, A, Nkeep, Nsweep);
% generate random MPS, MPO and compare whether the function works well. 
% 
% Returns 
%   Equality(boolean) : whether the output of the function equals expected output 
% 
% Written by M.Kim (Nov.29,2022)

% system parameter
J = 1; % coupling strength
L = 4; % number of sites in a chain
T = 0.1; 

% Imaginary time evolution parameters
Nkeep = 30; % bond dimension
beta = 1/T; 
dt = 1/1000; % discrete time step size

% Local operators
[S,I] = getLocalSpace('Spin',1);

%Set Hamiltonian in MPO 
% nearest-neighbor interaction terms
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

% full chain
Hs = cell(1,L);
Hs(:) = {Hloc};
Hs{1} = Hs{1}(:,:,end,:); % choose the last components of the left leg
Hs{end} = Hs{end}(:,:,:,1); % choose the first components of the right leg
%       |2            [legs 1 and 2 are for site n;
%    3  |    4         legs 3 and 4 are for site n+1]
% ----- Hs{n}-----
%       |      
%       |1      
H = cell(1,L-1);
H(:) = {J*contract(S,3,3,permute(conj(S),[2 1 3]),3,3)};
%       2      4      [legs 1 and 2 are for site n;
%       |      |       legs 3 and 4 are for site n+1]
%      [ H{n}  ]
%       |      |
%       1      3


% Initialize state 

M = cell(1,L);
for itN = (1:L)
    if itN <= (L/2)
        M{itN} = permute([1,0,0],[1 3 2]);
    else
        M{itN} = permute([0,0,1],[1 3 2]);
    end
end

% Set exact value 
% Test function
[ts, M, Ovals, EE,dw] = TS_1D(M, H, Hs, Nkeep, dt, beta/2, true);
[M, S] = canonForm(M, 0, [], []);
S
figure;
plot(ts, Ovals)
tol = 2^(-13);
if abs(Ovals(end)-Ovals(end-1))<tol
    equality = true; 
else 
    equality = false;

end

