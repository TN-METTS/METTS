function M = CPS_collapse(M,dim)
%[] = CPS_collapse(M,dim) 
%
% Obtain CPS by collapsing METTS
% 
% Input :
% M : MPS
% dim : dimension of local Hilbert space
%
%
%
%

N = numel(M); % the number of sites

C = cell(1,dim); % pure state projector

prob = zeros(1,dim); % probability for each pure state

for it = (1:dim)
    null = zeros(1,dim);
    null(1,it) = 1;
    C(it) = diag(null);
end


[M,~,~] = canonform(M,N,[],0);
[M,~,~] = canonform(M,1,[],0); % site-canonical form with orthogonality center site 1



for it = (1:N-1)
    %  obtain probability
    nrm = sqrt(contract(M{it},3,[1 2 3],conj(M{it}),3,[1 2 3]));
    M{it} = M{it}/nrm; % normalize MPS
    
    for it2 = (1:dim)
       prob(it2) = trace(updateLeft([],[],M{it},C(it2),2,M{it}));
    end
    % choose 1 state with projected probability
    prob = round(prob*10000);
    idx = randi(sum(prob));
    
    for it2 = (1:dim)
       if idx <= prob(it2)
            c = it2;
           break;
       end
    end
    
    new = M{it};
    new(:,:,c) = 1;
    new(:,:,[(1:1:c-1),(c+1:1:dim)]) = 0;
    
    % move on to the next site (it+1) via SVD
    
    [U,S,Vd] = svdTr(M{it},3,[1 3],[],0);
    
    M{it} = permute(U,[1 3 2]);
    S = contract(diag(S),2,2,Vd,2,1);
    M{it+1} = contract(S,2,2,M{it+1},3,1);
end


end

