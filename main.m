clear
% system parameter
J = 1; % coupling strength
L = 4; % number of sites in a chain
T = 0.05; % temperature
SN =50; % sampling number

%log
print_log = false;

% Imaginary time evolution parameters
Nkeep = 20; % bond dimension
beta = 1/T; 
dt = 1/20; % discrete time step size

% Local operators
[S,I] = getLocalSpace('Spin',1);

%Set gates
% nearest-neighbor interaction terms
H = cell(1,L-1);
H(:) = {J*contract(S,3,3,permute(conj(S),[2 1 3]),3,3)};
%       2      4      [legs 1 and 2 are for site n;
%       |      |       legs 3 and 4 are for site n+1]
%      [ Hs{n}  ]
%       |      |
%       1      3

% full chain
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

Hs = cell(1,L);
Hs(:) = {Hloc};
Hs{1} = Hs{1}(:,:,end,:); % choose the last components of the left leg
Hs{end} = Hs{end}(:,:,:,1); % choose the first components of the right leg
%       |2            [legs 1 and 2 are for site n;
%    3  |    4         legs 3 and 4 are for site n+1]
% ----- Hs{n}-----
%       |      
%       |1      

% Initialize state 

M = cell(1,L);
for itN = (1:L)
    if itN <= (L/2)
        M{itN} = permute(rand(1, 3),[1 3 2]);
    else
        M{itN} = permute(rand(1, 3),[1 3 2]);
    end
end

M = CPS_collapse(M, size(M{1}, 3));

% operator to measure magnetization
Sz = S(:,:,2);
operator = Hs;
E = zeros(1, SN);
% TS_1D
for itS=(1:SN)
    [ts,M,Ovals,EE,dw] = TS_1D(M, H, [], Nkeep, dt, beta/2, print_log);
    E(itS) = exp_val(M, operator);
    M = CPS_collapse(M, size(M{1}, 3));%,print_log);
end 
mean(E)



