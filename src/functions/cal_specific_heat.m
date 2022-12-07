function [Cv, Cv_error] = cal_specific_heat(E, E2, N, T)
    SN = size(E, 1);
    step_num = size(E, 2);
    Cv = (mean(E2(:, 3:step_num),  'all')-mean(E(:, 3:step_num), 'all').^2)./N;
    Cv_error = Cv.^2.*sqrt(2/(SN*(step_num-2)-1))./T^2;
    Cv = Cv./(T^2);
end 