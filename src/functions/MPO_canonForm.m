function [W, lognorm] = MPO_canonForm(W,id, Nkeep, Skeep)
% Make canonical form of MPO using MPS canonForm function 
% 
% 
% < Input >
% W: [1 X L cell array] MPO
%    leg order: bottom-top-left-right
% id: [numeric] Index for the bond connecting the tensor W{id} and W{id+1}
%               if id ==0(L) right(left) canonical form 
% Nkeep, Skeep: for cononForm function 
% < Output >
% W: [1 X L cell array] MPO
%    leg order: bottom-top-left-right
%
% Written by M.Kim (Dec.01,2022)
% Revised by M.Kim (Dec.01,2022); add truncation
    N = numel(W);
    % Isometries for merging physical legs 
    Aloc = cell(1,N); % isometries for merging the bottom and top legs of MPO tensors
    for itN = (1:N)
        Aloc{itN} = getIdentity(W{itN},1,W{itN},2);
        W{itN} = contract(W{itN},4,[1 2],Aloc{itN},3,[1 2]);
    end
    
    % Use canonForm for MPS
    [W, norm_] = canonForm(W, id, Nkeep,Skeep);% No truncation 
    lognorm = log(norm_);
%     if id==0
%         W{1} = W{1}.*norm_;
%     else
%         W{id} = W{id}.*norm_;
%     end 
    
    % Bring back to rank-4
    for itN = (1:N)
        W{itN} = contract(W{itN},3,3,conj(Aloc{itN}),3,3,[3 4 1 2]);
    end

end

