function M = CPS_collapse(M,dim)
%[] = CPS_collapse(M,dim) 
%
% Obtain CPS by collapsing METTS
% 
% Input :
% M : (1,N) cell MPS
% dim : dimension of local Hilbert space (i.e. the dimension of physical
% leg)
%
%
%
%
tobj = tic2;


N = numel(M); % the number of sites

C = cell(1,dim); % set for pure state projector

prob = zeros(1,dim); % probability for each pure state

for it = (1:dim)
    null = zeros(1,dim);
    null(1,it) = 1;
    C{it} = diag(null);
end


[M,~,~] = canonForm(M,N,[],0);
[M,~,~] = canonForm(M,1,[],0); % site-canonical form with orthogonality center site 1



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
    prob(end)
    for it2 = (1:dim)
       if idx <= prob(it2)
            c = it2;
           break;
       end
    end
    
    M{it} = contract(M{it},3,3,C{c},2,2);
    
    % move on to the next site (it+1) via SVD
    
    [U,S,Vd] = svdTr(M{it},3,[1 3],[],0);
    
    M{it} = permute(U,[1 3 2]);
    S = contract(diag(S),2,2,Vd,2,1);
    M{it+1} = contract(S,2,2,M{it+1},3,1);
end

toc2(tobj,'-v');

end

