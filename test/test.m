% This is a script of executing test functions. 
% Written by M.Kim (Nov.29,2022)
disptime("==============================================================\n")
disptime("Test time evolution algorithm\n")

% Test initial guess 
disptime("Test initial guess");
equality = Test_initial_guess(); 
if equality 
    disp("Test succeed");
else 
    disp("Test failed for initial guess!")
end 

disp("Test sweeps with fitting function")
equality= Test_sweeps();
if equality 
    disp("Test succeed")
else 
    disp("Test failed for the sweeps with the fitting function!!")
end 

disp("Test time evolution with Hamiltonian")
equality = Test_time_evolution(); 

if equality 
    disp("Test succeed")
else 
    disp("Test failed for time evolution with Hamiltonian!!")
end 

disptime("==============================================================\n")
disptime("Test CPS algorithm\n")
