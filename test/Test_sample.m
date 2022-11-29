function Equality= Test_sample() 
% This is a template of test functions. 
% The test function tests a sample function with input. 
% It compares the output of the function with expected output and return the equality. 
% 
% Returns 
%   Equality(boolean) : whether the output of the function equals expected output 
% 
% Written by M.Kim (Nov.29,2022)

input = []; % Generate sample input
output = func(input);
expected_output = []; % write expected output of the sample input

Equality = (output == expected_output); 
end

