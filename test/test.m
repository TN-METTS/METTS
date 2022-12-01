% This is a script of executing test functions. 
% Written by M.Kim (Nov.29,2022)
disptime('==============================================================\n')
disp('Test time evolution algorithm')

% Test initial guess 
disp('Test initial guess');
equality = Test_initial_guess(); 
if equality 
    disp('Test succeed');
else 
    disp('Test failed for initial guess!')
end 

disp('Test sweeps with fitting function')
equality= Test_sweeps();
if equality 
    disp('Test succeed')
else 
    disp('Test failed for the sweeps with the fitting function!!')
end 

disp('Test time evolution with nearest neighbor Hamiltonian')
equality = Test_TS_1D(); 

if equality 
    disp('Test succeed')
else 
    disp('Test failed for TS_1D function!!')
end 

disptime('==============================================================\n')
disptime('Test CPS algorithm')
