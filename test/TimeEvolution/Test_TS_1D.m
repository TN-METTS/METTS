function equality = Test_TS_1D()

% This function test time evolution for nearest-neighbor interaction which is implemented ad TS_1D. 
% Also, it tests exp_val function 
% B = MPO_multiplication_sweep(W, A, Nkeep, Nsweep);
% 
% Returns 
%   Equality(boolean) : whether the output of the function equals expected output 
% 
% Written by M.Kim (Nov.30,2022)

% system parameter

J = 1; % coupling strength
L = 5; % number of sites in a chain
T = 0.01; 

% Imaginary time evolution parameters
Nkeep = 30; % bond dimension
beta = 1/T; 
dt = 1/30; % discrete time step size

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

% get ground state energy from exact diagonalization 
expected_output = get_EG(Hs, I);

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
    % Set random input
    M{itN} = permute(rand(1, 3),[1 3 2]);
end
M = canonForm(M, 0, [], []); % normalize

% Set exact value 
% Test function
M = TS_1D(M, H, Nkeep, dt, beta/2, false);
output = exp_val(M, Hs);


tol = 2^(-13);
if abs(output-expected_output)<tol
    equality = true; 
else 
    equality = false;
    disptime(['Max difference : ', sprintf('%f\n', abs(output-expected_output))]);

end
end

function EG=get_EG(Hs, I)
   Hs_tot = 1; % initialize
   L = numel(Hs);
    for itN = (1:L)
        Hs_tot = contract(Hs_tot,2*itN,2*itN,Hs{itN},4,3);
    end
    % permute the left- and rightmost legs to the end
    Hs_tot = permute(Hs_tot,[(2:2*L+2) 1]);
    
    % merge the incoming legs into a thick incoming leg;
    % merge the outgoing legs into a thick outgoing leg
    Hs_tot = permute(Hs_tot,[(1:2:2*L) (2:2:2*L)]);
    Hs_tot = reshape(Hs_tot,(size(I,1)^L)*[1 1]);

    D = eig(Hs_tot);
    EG = min(D);
end 
