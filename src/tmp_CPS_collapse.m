function CPS = CPS_collapse(M,dim,print_log)
% CPS = CPS_collapse(M,dim) 
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
%
% dim : [integer] dimension of local Hilbert space (i.e. the dimension of physical
% leg)
% print_log [boolean] : whether print time, used memory or not
%
% CPS : [1 x N cell array] collapsed matrix product space (to classical
% product space)
%
if print_log
    tobj = tic2;
end

N = numel(M); % the number of sites

C = cell(1,dim); % set for pure state projector

Proj = cell(1,dim); % set for pure state projector
CPS = cell(1, N); % collapsed state
idx = zeros(1, N); % axis to collapse 
prob = zeros(1,dim); % probability for each pure state
T = cell(1, dim);
% generate pure state projector
for it = (1:dim)
    null = zeros(1,dim);
    null(1,it) = 1;
    C{it} = diag(null);
    Proj{it} = null;
end

[M,s,~] = canonForm(M,1,[],0); % site-canonical form with orthogonality center site 1

M{1} = contract(M{1}, 3, 2, diag(s), 2, 1, [1 3 2]);
for it = (1:N)
    for it2 = (1:dim)
        T{it2} = contract(M{it}, 3, 3, Proj{it2}, 2, 2);
        prob(it2) = contract(T{it2}, 3, (1:3), conj(T{it2}), 3,(1:3));
    end
    % choose 1 state with projected probability
    prob = prob./sum(prob, 'all');
    rn = rand(1);
    probC = cumsum(prob);
    %Set collapse axis 
    for it2=(1:dim)
        if rn<=probC(it2)
            idx(it) = it2; 
            break; 
        end 
    end 
    CPS{it} = contract(M{it}, 3, 3, C{idx(it)}, 2, 2);

    contract(T{idx(it)}, 3, 2, M{it+1}, 3, 1);
    % move on to the next site (it+1) via SVD
    if it<N
        [U,S,Vd] = svdTr(M{it},3,[1 3],[], dim*eps);
        M{it} = permute(U,[1 3 2]);
        S = contract(diag(S),2,2,Vd,2,1);
        M{it+1} = contract(S,2,2,M{it+1},3,1)./sqrt(prob(idx(it)));
    end    
    
    % CPS noise ? 
end
% collapse 
% for it=(1:N)
%     CPS{it} = contract(M{it}, 3, 3, C{idx(it)}, 2, 2);
% end 
CPS = M;


if print_log
    toc2(tobj,'-v');
    chkmem;
end

end


