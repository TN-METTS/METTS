function M = CPS_collapse(M,dim,print_log)
% M = CPS_collapse(M,dim) 
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
% M : [1 x N cell array] collapsed matrix product space (to classical
% product space)
%
if print_log
    tobj = tic2;
end

N = numel(M); % the number of sites

C = cell(1,dim); % set for pure state projector

prob = zeros(1,dim); % probability for each pure state

for it = (1:dim)
    null = zeros(1,dim);
    null(1,it) = 1;
    C{it} = diag(null);
end


[M,S,~] = canonForm(M,1,[],0); % site-canonical form with orthogonality center site 1
M{it1} = contract(M,3,2,diag(S),2,1,[1 3 2]);


for it = (1:N-1)
    %  obtain probability
    nrm = sqrt(contract(M{it},3,[1 2 3],conj(M{it}),3,[1 2 3]));
    M{it} = M{it}/nrm; % normalize MPS
    
    for it2 = (1:dim)
       prob(it2) = trace(updateLeft([],[],M{it},C{it2},2,M{it}));
    end
    % choose 1 state with projected probability
    prob = round(prob*10000);
    prob = cumsum(prob);
    idx = randi(prob(end));
    for it2 = (1:dim)
       if idx <= prob(it2)
            c = it2;
           break;
       end
    end
    
    M{it} = contract(M{it},3,3,C{c},2,2);
    
    % move on to the next site (it+1) via SVD
    
    [U,S,Vd] = svdTr(M{it},3,[1 3],[],dim*eps);
    
    M{it} = permute(U,[1 3 2]);
    S = contract(diag(S),2,2,Vd,2,1);
    M{it+1} = contract(S,2,2,M{it+1},3,1);
    
    % CPS noise ? 
end

if print_log
    toc2(tobj,'-v');
    chkmem;
end

end


