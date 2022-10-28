# ctpToMstoch
This C++ code was used to produce the small mutation limit results in the paper T. Lenaerts, J.M. Pacheco and F.C. Santos (2022) Evolution of a theory of mind. 
The functions for producing each result are part of the main file, which needs to be edited and recompiled for each result.

## Contents

- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Compilation Guide](#compilation-guide)
- [Running Guide](#running-guide)
- [Citation](#citation)
# Repo Contents

All *.CPP and *.HPP files 

# System Requirements

## Hardware Requirements

This code was built on a standard Macbook pro with 

- RAM: 8+ GB
- CPU: 4 cores, 2.4 GHz

## Software Requirements

OS Requirements

The package compilation was tested on OS X and Linux operating systems. 
- Linux: Ubuntu 16.04
- Mac OSX: 12.6 
- Windows:

To compile you will need
- GSL (tested for version 2.6 linux platform and version 2.7 on mac os X) 
- GCC (tested with gcc-9.3 om linux platform and clang on mac os X)

# Compilation Guide
A CMAKE file is provide. Version 3.16.4 was used on the linux platform. 
Compile by running first
```
cmake .
```
 to create the makefiles.  Once that works run 
```
make
```
which should produce an executable `ctpstoch`. 


# Running guide
When executing `ctpstoch` wihtou modifications in the `main.cpp`, the program will produce for an ICG with L=4 the results for β=0.3 and ε=0.18.


