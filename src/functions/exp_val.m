function vals = exp_val(M, operator, isleft)
% Calculating expectation value of local operator (operator) for given MPS
% M. 
% This function is revised from the subfunction of t_DMRG function by S.Lee
% < Input >
% M: [cell] Input MPS.
% operator : [cell] 
%           if number of element of operator  == N : operator in MPO form.
%           else : use local operator
% isleft : [logical] If true, it means that the MPS M is in left-canonical
%       form. Otherwise, right-canonical form.
%
% < Output >
% vals : [numeric/vector] Expectation value of operator in MPO form. 
%
%
% Written by M.Kim (Nov.30,2022)

N = numel(M);
if numel(operator)==N % operator in MPO form 
    MM = 1; 
    vals=1;

    for itN = (1:N)
        T2 = contract(vals,3,3,M{itN},3,1); %L1 L2 MR MD
        T1 = contract(MM,3,3,M{itN},3,1); %L1 L2 MR MD
        T2 = contract(T2,4,[2, 4],operator{itN},4,[3 2]); % L1 MR OD OR

        vals = contract(conj(M{itN}),3,[1, 3],T2,4,[1 3], [1 3 2]); % MdR OR MR
        MM = contract(conj(M{itN}),3,[1,3],T1,4,[1,4]);
    end 
    vals = vals/MM;
elseif numel(operator) == 1  
    vals = zeros(1,N);
    O = operator{1};
    MM = 1; % contraction of bra/ket tensors
    if isleft % left-normalized
        for itN = (N:-1:1)
            T = permute(M{itN},[2 1 3]); % permute left<->right to make use of updateLeft
            T2 = contract(MM,2,2,T,3,1);
            T2 = contract(T2,3,3,O,2,2);
            vals(itN) = contract(conj(T),3,(1:3),T2,3,(1:3));
            MM = updateLeft(MM,2,T,[],[],T);
        end
    else % right-normalized
        for itN = (1:N)
            T2 = contract(MM,2,2,M{itN},3,1);
            T2 = contract(T2,3,3,O,2,2);
            vals(itN) = contract(conj(M{itN}),3,(1:3),T2,3,(1:3));
            MM = updateLeft(MM,2,M{itN},[],[],M{itN});
        end 
    end
end 

end