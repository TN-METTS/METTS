function B = MPO_multiplication_sweep(W, A, Init_B, Nkeep, Nsweep)
% < Description >
%
% B = MPO_multiplication_sweep(W, A, Init_B, Nkeep, Nsweep)
%
% This function is revised from function "C = mtimes_MPO
% (B,A,Nkeep,Nsweep)" written by S.Lee 
% It implemented fitting algorithm for multiplying an MPO(W) times an MPS(A).
% Calculate B = WA by fitting algorithm with/without initial guess 
%                                      
%     1    2 1    2            2         1    2 1    2            2 
%    --B{1}-- B{2}--.... -B{N}--        --A{1}----A{2}---.... - A{N}--
%       |       |           |     =         |3     |3            |3        
%       |3      |3          |3              *      *             *
%                                           |2     |2            |2                                      
%                                       --W{1}----W{2}---.... - W{N}--
%                                         3 | 4   3 | 4         3 | 4 
%                                           |1      |1            |1 
% Note that it start sweep from right to left. 
%
% The detail diagram is in Figure 4 in 
% [ E. M. Stoudenmire and S. R. White, New J. Phys. 12, 055026 (2010).]
% 
% < Input >
% W: [1 X L cell array] MPO to be multiplied
%    leg order: bottom-top-left-right
% A: [1 X L cell array] MPS to be multiplied
%    leg order: left-right-bottom
% Init_B: [1 X L cell array] Initial guess of B=WA
%    leg order: left-right-bottom
% Nkeep : [numeric] The maximum bond dimension for the result MPS B.
% Nsweep : [numeric] Number of round trips of sweeping. There will be
%       Nsweep pairs of left-to-right sweep and right-to-left sweep. That
%       is, the first sweep is left-to-right and the last is right-to-left.
%
% < Output >
% B : [1 x L cell array] Multiplication of A and W. 
%     Each tensor B{n} follows the same leg order
%     leg order: left-right-bottom(physical)
%
% Written by M.Kim (Nov.30,2022)
N = numel(A);
BWA = cell(1, N+2);
BWA{1} = 1; BWA{end} = 1; 

if isempty(Init_B)
    B = A; 
    B = canonForm(B, N, [], []);
end 

% sanity check
if N ~= numel(W)
    error('ERR: Length of two input MPO and MPS do not match.');
end

% Initialize BWA. Contract from left
% First sweep: right to left 
for itN = (1:N)
    % exchange left and right since the contract_BWA function is
    % implemented for leftwise BWA
    BWA{itN+1} = contract_BWA(BWA{itN}, B{itN}, W{itN}, A{itN});
end 

for itS = (1:Nsweep)
    for itN=((N-1):-1:1) % right-to-left sweep
        L = BWA{itN}; Al = A{itN}; Wl = W{itN};
        R = BWA{itN+3}; Ar = A{itN+1}; Wr = W{itN+1};
        T = contract_all(L, Al, Wl, R, Ar, Wr); % This is B{itN} * B{itN+1}
        % calculate B{itN}, B{itN+1} using svd, move orthogonality center
        % to B{itN}
        [U, S, Vd] = svdTr(T, 4, [1 2], Nkeep, []);
        B{itN+1} = permute(Vd,[1 3 2]); 
        B{itN}=contract(U, 3, 3, diag(S), 2, 1, [1 3 2]); 
        % update right BWA(R)
        % exchange left, right legs since the function is implemented
        % for leftwise BWA
        BWA{itN+2} = contract_BWA(BWA{itN+3},  permute(B{itN+1}, [2 1 3]), ...
        permute(W{itN+1}, [1 2 4 3]), permute(A{itN+1}, [2 1 3]) );
    end 

    for itN=(1:(N-2)) % left-to-right sweep
        L = BWA{itN}; Al = A{itN}; Wl = W{itN};
        R = BWA{itN+3}; Ar = A{itN+1}; Wr = W{itN+1};
        T = contract_all(L, Al, Wl, R, Ar, Wr); % This is B{itN} * B{itN+1}
        % calculate B{itN}, B{itN+1} using svd, move orthogonality center
        % to B{itN+1}
        % leg order left-right-bottom(physical)
        [U, S, Vd] = svdTr(T, 4, [1 2], Nkeep, []);
        B{itN} = permute(U, [1 3 2]);
        B{itN+1} = contract(diag(S), 2, 2, Vd, 3, 1, [1 3 2]);
        % update left BWA(L)
        BWA{itN+1} =  contract_BWA(BWA{itN}, B{itN}, W{itN}, A{itN});
    end 
    
end 
end

% These subfunctions are implemented for contraction 
% Contractions are decribed with diagram 

function T = contract_WA(T, W, A)
%        
%        
%       1     2
%   /----- A ----
%   |      | 
%   | 3    | 3
%   |      | 2
%   |  2 3 |   4
%   T ---- W ----
%   |      |
%   |      | 1
%   | 1
%   |
%   \----
    T = contract(T,3, 3, A, 3, 1); % T1 T2 A2 A3
    
    T = contract(T, 4, [2,4], W, 4, [3,2], [1 2 4 3]);
    
    % return with leg order T1, A2, W4, W1
end 

function T = contract_BWA(T, B, W, A)
%        
%        
%             2
%   /----- A ----             3    
%   |      |                  
%   |      |                  
%   |      |                   
%   |      |  3 
%   T ---- W ----    ->       2
%   |      |
%   |      | 4
%   | 1    | 3
%   |  1   |   2
%   \----- B ----             1   
    T = contract_WA(T, W, A);
    T = contract(conj(B), 3, [1, 3], T, 4, [1 4], [1 3 2]);
    
end 

function T = contract_all(L, Al, Wl, R, Ar, Wr)

%
% 
%   /--1-- Al -2*1--Ar -2----\
%   |      |        |        |
%   |3     |3       | 3      | 3
%   |      |2       | 2      |
%   |  2*3 |   4*3  |  4*2   |
%   L ---- Wl ----- Wr ----- R
%   |      |        |        |
%   |      |1       | 1      |
%   | 1                      | 1
%   |                        |
%   \----                ----/
    L = contract_WA(L, Wl, Al);
    R = contract_WA(R, permute(Wr, [1 2 4 3]), permute(Ar, [2 1 3]));
    T = contract(L, 4,[2 3], R, 4, [2 3], [1 2 4 3]);
end 
