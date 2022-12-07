# METTS

Minimally Entangled Typical Quantum States(METTS) at Finite Temperature    

## Overview 
This small project calculates expectation value of observables using METTS.    
First, it calculates the energy and compare the convergence using the ensemble.    
Second, it calculates the magnetic susceptibility and specific heat using METTS.    
For the detail, refer the paper(https://iopscience.iop.org/article/10.1088/1367-2630/12/5/055026).    



## Execution 

### Setting local path 
In main.m, startip.m files you shouled revise the PATH.    
For the main.m, change PATH to your PATH of result/ directory in your local computer.    
Also, in startup.m, revise the path to the current project.    

### Calculation and Visualization
To calculate specific heat and susceptibility using METTS, set boolean variable named "isSaved" in main as "false".   
Else, set as "true"    
You can execute calculation and visualization using following command in the terminal.    

```
main.m 
```
