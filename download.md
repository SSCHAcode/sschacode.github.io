---
layout: page
title: Download and Installation
---

## Download

In order to use the SSCHA code you will need to download both the  *CellConstructor* package and the *python-sscha* package. These packages can be downloaded directly from the following github pages:

> [**github.com/SSCHAcode/CellConstructor**](https://github.com/SSCHAcode/CellConstructor)
>
> [**github.com/SSCHAcode/python-sscha**](https://github.com/SSCHAcode/python-sscha)

In case you want to use the force field calculator used in the [SnTe tutorial](http://sscha.eu/Tutorials/SnSe/), you will need to download as well the *F3ToyModel* package from this link:

> [**github.com/SSCHAcode/F3ToyModel**](https://github.com/SSCHAcode/F3ToyModel) 

## Requirements for installation

In order to install both *CellConstructor* and *python-sscha* packages, you need to install previously python and several python packages. As *python-sscha* depends on *CellConstructor*, the former will not work unless the latter is installed. As the SSCHA code also is partly written in Fortran, you will aslo need a Fortran compiler as well as Lapack and Blas librearies. 

The full list of dependencies to install *CellConstructor* and *python-sscha* packages is:
1. python
2. numpy
3. matplotlib
3. Lapack
4. Blas
5. gfortran (or any fortran compiler)
6. Atomic Simulation Environment (ASE)
7. SPGLIB

All the needed python dependencies can be easily installed for *CellConstructor* or *python-sscha* by simply running
```
pip install -r requirements.txt
```
in the folder of each package.


## Installation

To install the *CellConstructor* and *python-sscha* packages it is recommended to use the last version of anaconda-python2, which cames with all the updated numpy and matplotlib packages already compiled to work in parallel. Moreover the last version of matplotlib will allow the user to modify the plots after they are produced.

Once all the dependencies have been installed, the *CellConstructor* and *python-sscha* codes can be easily installed from the command line as:
```
python setup.py install
```
This command must be executed in the directory where the setup.py script is, insice the *CellConstructor* and *python-sscha* folders. If you installed python in a system directory, administration rights may be requested (add a sudo before the command). 

Installing *CellConstructor* and *python-sscha* in clusters may me more tricky and one needs to adapt the setup.py to the cluster characteristics. For instance, if you use the intel compiler, you need to delete the lapack linking from the setup.py and include -mkl. Note that you must force to use the same linker as the one used for the compilation. For example, specific setup.py scripts are provided with the distribution to install *CellConstructor* easily in FOSS or INTEL clusters. 

## Installation through pip

Alternatively, both *CellConstructor* and *python-sscha* can be installed through pip simply as:
```
pip install CellConstructor
``` 
and
```
pip install python-sscha 
```
