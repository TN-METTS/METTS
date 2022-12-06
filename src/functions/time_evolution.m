function [state, epsilons] = time_evolution(Hs, tau, state,  num_step, p, Nkeep, Skeep)
%TIME_EVOLUTION 이 함수의 요약 설명 위치
%   자세한 설명 위치

% get logarithmic discrete times. 
% Sum of epsilons is tau up to error tau*(1/2)^n
%epsilons = 2.^(-(1:num_step)).*tau;


Nsweep = 3; 
K_epsilon = propagator(Hs, p, tau/num_step, Nkeep, Nsweep);
for itn=(1:num_step)
    itn
    state = zipup_algo(K_epsilon, state, Nkeep, Skeep);
    
    size(state{50})
end 

end

