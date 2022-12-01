function state = time_evolution(Hs, tau, state,  n, p, Nkeep, Skeep)
%TIME_EVOLUTION 이 함수의 요약 설명 위치
%   자세한 설명 위치

% get logarithmic discrete times. 
% Sum of epsilons is tau up to error tau*(1/2)^n
epsilons = 2^(-(1:n)).*tau;

for itn=(1:n)
    K_epsilon = propagator(Hs, p, epsilons{itn});
    state = zipup_algo(K_epsilon, state, Nkeep, Skeep); 
end 

end

