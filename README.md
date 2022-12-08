# METTS
This is the final exam of Tensor Network course in SNU. 2022    
Minimally Entangled Typical Quantum States(METTS) at Finite Temperature algorithms 

## Overview 
This small project calculates expectation value of observables using METTS.    
First, it calculates the energy and compare the convergence using the ensemble.    
Second, it calculates the magnetic susceptibility and specific heat using METTS.    
For the detail, refer the [paper](https://iopscience.iop.org/article/10.1088/1367-2630/12/5/055026).    



## Execution 

### Setting local path 
In the main.m, change PATH to your PATH of result/ directory in your local computer.    
For example,    
```
PATH = 'C:/Users/User/Documents/MATLAB/FinalProject/result/';
```
Also, execute startup.m file to add functions in searchspace. 
```
startup.m
```

### Calculation and Visualization
To calculate specific heat and susceptibility using METTS, set boolean variable named "isSaved" in main.m as "false".   
Else, set as "true"    
```
% To calculate, set isSaved = false
isSaved = true; 
```
You can execute calculation and visualization using following command in the terminal.    

```
main.m 
```
