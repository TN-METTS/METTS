function K_epsilon = propagator(Hs, p, epsilon)

W = - epsilon.* Hs;
for itN=(1:p)
    W = zipup_algo(Hs, add_MPO_scalar(1, 1/(p - itN+1), W));
    W = W .*(-epsilon);
end 

K_epsilon = add_MPO_Scalar(1, W);
end


function MPO = add_MPO_scalar(a, M)
    % calculate MPO = 1 + a * M 
    MPO = cell(1,N);
    MPO{1} = cat(4,getIdentity(M{1},2) , a*M{1});
    for itN = (2:N-1)
        sz = [size(M{itN}),ones(1,4-ndims(M{itN}))];
        MPO{itN} = cat(3, ...
            cat(4,getIdentity(M{itN},2),zeros([sz(1:2) 1 sz(4)])), ...
            cat(4,zeros([sz(1:2) sz(3)]),M{itN}));
    end
    MPO{N} = cat(3,getIdentity(M{1},2),M{N});
end 
