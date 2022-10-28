# ctpToMstoch
This C++ code was used to produce the small mutation limit results in the paper T. Lenaerts, J.M. Pacheco and F.C. Santos (2022) Evolution of a theory of mind. 
The functions for producing each result are part of the main file, which needs to be edited and recompiled for each result.

## Contents

- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Compilation Guide](#compilation-guide)
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
GSL (tested for version 2.6 with gcc 9.3 on linux platform, and version 2.7 on mac os X) 

# Compilation Guide
A CMAKe file is provide.
Compile by running first, to create the makefiles 
```
cmake .
```
Once that works run 
```
make
```


