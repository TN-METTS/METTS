clear

% Set path 
PATH = 'C:/Users/IQSL/Documents/MATLAB/TN/FinalProject/result/';

% To calculate, set isSaved = false
isSaved = true; 

% solve (c) : reproduce figure 6 
if ~isSaved 
    reproduce_figure6;
end 
generate_figure6;

%Calculate susceptibility and specific heat 
if ~isSaved 
    specific_heat; %should be revised 
    susceptibility;

end 

% generate figure 
generate_specific_heat_susceptibility_figure;


