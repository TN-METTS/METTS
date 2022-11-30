function [B, C] = zipup_algo(W, A, Nkeep)
% < Description >
%
% [B, C] = zipup_algo(W, A)
% 
% It implemented zipup algorithm introduced in the paper 
% [ E. M. Stoudenmire and S. R. White, New J. Phys. 12, 055026 (2010).]
% The detail diagram is in Figure 5. in the paper. 
% 
% < Input >
% W: [1 X L cell array] MPO to be multiplied
%    leg order: bottom-top-left-right
% A: [1 X L cell array] MPS to be multiplied
%    leg order: left-right-bottom
% 
% Nkeep : [numeric] The maximum bond dimension for the result MPS B.
%
% < Output >
% B : [1 x L cell array] Multiplication of A and W using zipup algorithm. 
%     Each tensor B{n} follows the same leg order
%     leg order: left-right-bottom(physical)
% C: [numeric] Norm of the state. 
%
% Written by M.Kim (Nov.30,2022)


N = numel(A);

% sanity check
if N ~= numel(W)
    error('ERR: Length of two input W(MPO), A(MPS) do not match');
end

% Bring W into right-canonical form
% To make W as canonical form, merge physical legs and use canonForm
% function. 

% Isometries for merging physical legs 
Aloc = cell(1,N); % isometries for merging the bottom and top legs of MPO tensors
for itN = (1:N)
    Aloc{itN} = getIdentity(W{itN},1,W{itN},2);
    W{itN} = contract(W{itN},4,[1 2],Aloc{itN},3,[1 2]);
end

% Use canonForm for MPS
[W, norm_] = canonForm(W,0,[],[]);
W{1} = W{1}.*norm_;

% Bring back to rank-4
for itN = (1:N)
    W{itN} = contract(W{itN},3,3,conj(Aloc{itN}),3,3,[3 4 1 2]);
end

% Bring A into right-canonical form 
[A, norm_] = canonForm(A, 0, [], []);
A{1}  = A{1}.*norm_;

% malloc mem for return 
B = cell(1, N);

C = 1; 
for itN=(1:N)
    CWA = contract(C, 3, 3, A{itN}, 3, 1); % CD CWR AR AP
    CWA = contract(CWA, 4, [2, 4], W{itN}, 4, [3 2] , [1 3 4 2]); % CD, WD, WR, AR
    % Use twice of Nkeep according to the paper
    [U, s, Vd] = svdTr(CWA, 4, [1 2], 2*Nkeep, []);
    B{itN} = permute(U, [1 3 2]); % C1, D, WD
    C = contract(diag(s), 2, 2, Vd, 3, 1); % D, WR, AR
end 
% For the last iteration, C is just number due to the dummy index 


end

