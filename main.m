
clear

%Set Hamiltonian in MPO 

% system parameter
J = 1; % coupling strength
L = 50; % number of sites in a chain

% Local operators
[S,I] = getLocalSpace('Spin',1);

% nearest-neighbor interaction terms
Hs = cell(1,L-1);
Hs(:) = {J*contract(S,3,3,permute(conj(S),[2 1 3]),3,3)};
%       2      4      [legs 1 and 2 are for site n;
%       |      |       legs 3 and 4 are for site n+1]
%      [ Hs{n}  ]
%       |      |
%       1      3

M = cell(1,L);
for itN = (1:L)
    if itN <= (L/2)
        M{itN} = permute([1,0,0],[1 3 2]);
    else
        M{itN} = permute([0,0,1],[1 3 2]);
    end
end

% tDMRG parameters
Nkeep = 20; % bond dimension
dt = 1/20; % discrete time step size
tmax = 20; % maximum time

% operator to measure magnetization
Sz = S(:,:,2);

% tDMRG
[ts,M,Ovals,EE,dw] = TS_1D(M,Hs,Sz,Nkeep,dt,tmax);