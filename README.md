# GeneticCrn
Simulation source code for the paper ["Genetic Algorithm Aided Transmit Power Control in Cognitive Radio Networks"](https://ieeexplore.ieee.org/abstract/document/6849663).

## Abstract
We address the power control problem in cognitive radio networks where secondary users exploit spatial spectrum opportunities without causing unacceptable interference to primary users. An optimization problem is formulated aiming at maximizing the utility of secondary users and to ensure the QoS for both primary and secondary users. To solve the power allocation problem a genetic algorithm is developed, and two fitness functions are proposed. The first is oriented towards minimizing the total transmit power consumption of the secondary network. The second is a multi-objective function and is oriented to the joint optimization of total capacity and transmit power consumption of the secondary network. Results show a near-optimum performance of the genetic algorithm aided power control scheme based on the multi-objective fitness function.

## Getting Started
In order to run the simulations you need Matlab 2015a or higher and a C compiler compatible with the installed Matlab version. From the command line type:
```bash
git clone https://github.com/raikel/GeneticCrn
```
Open Matlab and add the source directory `src` (and all its subfolders) to the Matlab search path. In the Matlab workspace, open the directory `src\lib\mex` and type in the command window:
```bash
compile
```
This will compile all the source `mex` files. To run the simulation with default parameters values, type in the Matlab command window:
```Matlab
stats = netsim()
```

## Cite
Please cite this work as
> Bordón, R. B., Sánchez, S. M., Fernandez, E. M., Souza, R. D., & Alves, H. (2014, June). Genetic algorithm aided transmit power control in cognitive radio networks. In 2014 9th International Conference on Cognitive Radio Oriented Wireless Networks and Communications (CROWNCOM) (pp. 61-66). IEEE.
