function [B, C] = zipup_algo(W, A, Nkeep, Skeep)
% < Description >
%
% [B, C] = zipup_algo(W, A, Nkeep, Nsweep, Skeep)
% 
% It implemented zipup algorithm introduced in the paper 
% [ E. M. Stoudenmire and S. R. White, New J. Phys. 12, 055026 (2010).]
% The detail diagram is in Figure 5. in the paper. 
% 
% < Input >
% W: [1 X L cell array] MPO to be multiplied
%    leg order: bottom-top-left-right
% A: [1 X L cell array] MPS/MPO to be multiplied
%    [MPS]leg order: left-right-bottom(physical)
%    [MPO]leg order: bottom-top-left-right
% 
% Nkeep : [numeric] The maximum bond dimension for the result MPS B.
% Skeep: [numeric] Minimum magnitude of the singular value to keep at SVD
% < Output >
% B : [1 x L cell array] Multiplication of A and W using zipup algorithm. 
%     Each tensor B{n} follows the same leg order
%     [MPS]leg order: left-right-bottom(physical)
%     [MPO]leg order: bottom-top-left-right
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
[W, norm_] = canonForm(W,0,[],[]);% No truncation for initialization step 
W{1} = W{1}.*norm_;

% Bring back to rank-4
for itN = (1:N)
    W{itN} = contract(W{itN},3,3,conj(Aloc{itN}),3,3,[3 4 1 2]);
end

if numel(size(A{2}))==3 % A is MPS form  
    % Bring A into right-canonical form 
    [A, norm_] = canonForm(A, 0, [], []); % No truncation for initialization step 
    A{1}  = A{1}.*norm_;
    % Use twice of Nkeep according to the paper
    [B, C] =zipup_MPS(W, A, 2*Nkeep, Skeep);

elseif numel(size(A{2}))==4 % A is MPO form 
    % Convert A into MPS form 
    % Isometries for merging physical legs 
    loc = cell(1,N); % isometries for merging the bottom and top legs of MPO tensors
    for itN = (1:N)
        loc{itN} = getIdentity(A{itN},1,A{itN},2);
        A{itN} = contract(A{itN},4,[1 2],loc{itN},3,[1 2]);
    end
        [A, norm_] = canonForm(A, 0, [], []); % No truncation for initialization step 
        % Bring back to rank-4
    for itN = (1:N)
        A{itN} = contract(A{itN},3,3,conj(loc{itN}),3,3,[3 4 1 2]);
    end
    A{1}  = A{1}.*norm_;
    [B, C] =zipup_MPO(W, A, 2*Nkeep, Skeep);
    
end
end

function [B, C] =zipup_MPS(W, A, Nkeep, Skeep)
    % malloc mem for return 
    B = cell(1, N);
    C = 1; 
    for itN=(1:N)
        CWA = contract(C, 3, 3, A{itN}, 3, 1); % CD CWR AR AP
        CWA = contract(CWA, 4, [2, 4], W{itN}, 4, [3 2] , [1 3 4 2]); % CD, WD, WR, AR
        [U, s, Vd] = svdTr(CWA, 4, [1 2], Nkeep, Skeep);
        B{itN} = permute(U, [1 3 2]); % C1, D, WD
        C = contract(diag(s), 2, 2, Vd, 3, 1); % D, WR, AR
    end 
    % For the last iteration, C is just number due to the dummy index 
end 

function [B, C] =zipup_MPO(W, A, Nkeep, Skeep)
    % malloc mem for return 
    B = cell(1, N);
    C = 1; 
    for itN=(1:N)
        CWA = contract(C, 3, 3, A{itN}, 4, 3); % CD CWR AD AU AR
        CWA = contract(CWA, 5, [2, 3], W{itN}, 4, [3 2] , [1 2 4 5 3]); % CD, AU, WD, WR, AR
        
        [U, s, Vd] = svdTr(CWA, 5, [1 2,3], Nkeep, Skeep);
        B{itN} = permute(U, [3 2 1 4]); % CD, AU, WD, D -> WD AU CD D 
        C = contract(diag(s), 2, 2, Vd, 3, 1); % D, WR, AR
    end 
    % For the last iteration, C is just number due to the dummy index 
end 