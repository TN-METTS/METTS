
function M = tmp_CPS_collapse(M, V, print_log)
% CPS = CPS_collapse(M, isodd, print_log) 
%
% Obtain CPS by collapsing METTS; collapsed state is determined at each site randomly
% 
%
%
% Input :
% M : [1 x N cell array] Matrix product state to collapse 
%
%    1      2   1      2         1        2
%   ---M{1}---*---M{2}---* ... *---M{end}---
%       |          |                 |
%       ^3         ^3                ^3
% V: [local dim X local dim] : new basis to collapse 
% 
% print_log [boolean] : whether print time, used memory or not
%
% M : [1 x N cell array] collapsed matrix product space (to classical
% product space)
%
if print_log
    tobj = tic2;
end

    
N = numel(M); % the number of sites

% % convert basis 
% for it=(1:N)
%     M{it} = contract(M{it}, 3, 3, V, 2, 1);
% end 

dim = size(M{1}, 3); % dimension of local Hilbert space 
C = cell(1,dim); % set for pure state projector

% generate pure state projector
for it = (1:dim)
    Projector = zeros(1,dim);
    Projector(1,it) = 1;
    C{it} = V'*diag(Projector)*V;
end

prob = zeros(1,dim); % probability for each pure state

[M,s,~] = canonForm(M,1,[],0); % site-canonical form with orthogonality center site 1
M{1} = contract(M{1}, 3, 2, diag(s), 2, 1, [1 3 2]);

for it = (1:N)
    % calculate probability 
    for it2 = (1:dim)
        Tmp = contract(M{it}, 3, 3, C{it2}, 2, 2);
        prob(it2) = contract(Tmp, 3, (1:3), conj(M{it}), 3,(1:3));
    end

%     fprintf("At site %d| prob %0.4f, %0.4f, %0.4f, sum = %0.4f, norm = %0.4f\n", it, prob(1),prob(2),prob(3), sum(prob), contract(M{it}, 3, (1:3), conj(M{it}), 3,(1:3)));
    
    prob = prob./sum(prob, 'all'); %sum(prob, 'all') is 1 but divide with it to remove noise. 
    
    % choose 1 state with projected probability
    rn = rand(1);
    probC = cumsum(prob);
    
    %Set collapse axis 
    for it2=(1:dim)
        if rn<=probC(it2)
            idx = it2; 
            break; 
        end 
    end 

    % Collapse the state with C 
    % divide sqrt(prob(idx)) to ensure norm = 1 
    M{it} = contract(M{it}, 3, 3, C{idx}, 2, 2)./sqrt(prob(idx));
        
    % move on to the next site (it+1) via SVD
    if it<N
        MM = contract(M{it}, 3, 2, M{it+1}, 3, 1);
        [U,S,Vd] = svdTr(MM, 4, [1 2],[], dim*eps);
        M{it} = permute(U,[1 3 2]);
        M{it+1} = contract(diag(S),2,2,Vd,3,1);
    end    
    % CPS noise ? 
end

% % restore basis 
% for it=(1:N)
%     M{it} = contract(M{it}, 3, 3, V, 2, 1);
% end

if print_log
    toc2(tobj,'-v');
    chkmem;
end

end


