[![DOI](https://zenodo.org/badge/556244409.svg)](https://zenodo.org/doi/10.5281/zenodo.10230141)
# ctpToMstoch
This C++ code was used to produce the small mutation limit results in the paper T. Lenaerts, M. Saponara, J.M. Pacheco and F.C. Santos (2023) Evolution of a theory of mind. 
The functions for producing each result are part of the main file, which needs to be edited and recompiled for each result.

## Contents

- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Compilation Guide](#compilation-guide)
- [Demo Guide](#demo)
- [Results](#demo)
- [Citation](#citation)
# Repo Contents

All *.CPP and *.HPP files:
- `rangen.h/rangen.cpp`
- `ctpgame.hpp/ctpgame.cpp`
- `ctpdata.hpp/ctpdata.cpp`
- `strategy.hpp/strategy.cpp`
- `kernel.hpp/kernel.cpp`
- `main.cpp`

# System Requirements

## Hardware Requirements

This code was built on a standard Macbook pro (Monterey 12.6) with 
- RAM: 8 GB
- CPU: 2.4 GHz Quad-Core Intel Core i5

## Software Requirements

OS Requirements

The package compilation was tested on OS X (laptop) and Unix (cluster) operating systems. 

To compile you will need
- GSL (tested for version 2.6 Unix platform and version 2.7 on mac os X) 
- GCC (version gcc-9.3) on Unix platform and clang on Mac os X.

# Compilation Guide
A CMAKE file is provide. Version 3.16.4 was used on the Unix platform. 
Compile by running first
```
cmake .
```
 to create the makefiles.  Once that works run 
```
make
```
which should produce an executable `ctpstoch`. 


# Demo
When executing `ctpstoch` without modifications in the `main.cpp`, the program will produce for an ICG with L=4 the results for β=0.31 and ε=0.19. Changes to the parameters will produce also the L=6 results.

In `main.cpp` there are 13 parameters that can be set.  These are
- `length` The length of the ICG (default value is 4).
- `first` The payoff player 1 gets in the first step (default value is 0.4).
- `second` The payoff player 2 gets in the first step (default value is 0.1).
- `factor` The growth of the resource at each step (default value is 2.0).
- `levels` Reasoning levels of the individuals (default value is 4).
- `maxlevel` Maximum reasoning level in the population (default value is 4).
- `epsilon` Reasoning error for each individual (default value is 0.18).
- `betas` Selection strength in the stochastic evolutionary dynamic (default value is 0.3).
- `repeats` Number of repeats to calculate the fitness of each individual in the Monte Carlo sampling process (default value is 50000).
- `psize` Population size (default value is 500).
- `mut` Mutation probability (default value is 0.0).
- `scaling` Tuning of difference in transitions between reasoning levels (default value is 1.0).
- `cost` Cost associated with each additional reasoning level (default value is 0.0).
- `approach` type of reasoning process; 0)inertia, 1)payoff conditional and 2) unconditional


In `main.cpp` there are 7 functions that provide the results needed to reproduce the data in the figures:
- `runAnalytical` : this function determines the stationary distribution and other data for a given β and ε (see Figure 2 for instance). 
- `runBeta` : this function determines the stationary distribution and other data for a β within a range and a fixed ε. 
- `runEpislon` : this function determines the stationary distribution and other data for a given β and a range of ε values.
- `runBoth` : this function determines for a range of β- and ε-values the RMSD with the experimental data and average belief, reasoning, action and misbelief results (see Figure 1A for instance).
- `runNetwork` : this function generates the networks as visualized in Extended Data Figure 5.
- `runGradient` : this function produces the gradient of selection for specific values of  β and ε.
- `RunGradientEpsilon` : this function produces multiple gradients of selections for different ε-values.

# Output 
The raw output files as well as the R files needed to reproduce the figures can be found in the folder `results`.

The R-scripts can simply be executed in the same folder as the raw data files.

# Citation

When using or extending this code, please cite :

T. Lenaerts, M. Saponara, J.M. Pacheco and F.C. Santos (2023) Evolution of a Theory of Mind ...
