function [chi, chi_error] = cal_susceptibility(mz, mz2, N, T)
    SN = size(mz, 1);
    step_num = size(mz, 2);
    chi = (mean(mz2(:, 2:step_num), 'all')-mean(mz(:, 2:step_num), 'all').^2)./N;
    chi_error = chi.^2.*sqrt(2/(SN*(step_num-1)-1))./T;
    chi = chi./T;
end 