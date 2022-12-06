function K_epsilon = propagator(Hs, p, epsilon, Nkeep, Nsweep)
minus_eH = cell(size(Hs));

for it=(1:numel(Hs))
    
    minus_eH{it} = - epsilon.* Hs{it};
end 

W = minus_eH;
for itN=(1:p)
%     W = zipup_algo(Hs, add_MPO_scalar( 1/(p - itN+1), W)); 
    W = add_MPO_scalar(1/(p - itN+1), W);
    W = mtimes_MPO (W, minus_eH, Nkeep,Nsweep);% W = - eH * (1+ 1/(p-itN+1) * W)
    
end 

K_epsilon = add_MPO_scalar(1, W); % 1+ W_p
end


function MPO = add_MPO_scalar(a, M)
    % calculate MPO = 1 + a * M 
    N = numel(M);
    MPO = cell(1,N);
    MPO{1} = cat(4, getIdentity(M{1},2) , a.*M{1});
    for itN = (2:N-1)
        sz = [size(M{itN}),ones(1,4-ndims(M{itN}))];
        MPO{itN} = cat(3, ...
            cat(4,getIdentity(M{itN},2),zeros([sz(1:2) 1 sz(4)])), ...
            cat(4,zeros([sz(1:2) sz(3)]),M{itN}));
    end
    MPO{N} = cat(3,getIdentity(M{1},2),M{N});
end 
