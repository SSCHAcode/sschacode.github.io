---
layout: page
title: Installation
---

Here we provide the installation guide of the SSCHA code as well as other related codes distributed within the [SSCHA GitHub](https://github.com/SSCHAcode) repository. 

# 1. Requirements

Most of the codes require a fortran or C compiler and MPI configured. Here we install all the requirements to properly setup the SSCHA code. To properly compile and install the SSCHA code, you need a fortran compiler and LAPACK/BLAS available.

On Debian-based Linux distribution, all the software required is installed with (Tested on Quantum Mobile and ubuntu 20.04):
```
sudo apt update
sudo apt install libblas-dev liblapack-dev liblapacke-dev gfortran openmpi-bin
```
Note that some of the names of the libraries may change slightly in different linux versions or on MacOS.

## 1.1 Python installation

SSCHA is a Python library and program. Most linux distribution come with python already installed, however, for a performance boost, it is usually better to use the python distribution provided by the anaconda environment. Anaconda python can be downloaded and installed from [www.anaconda.com/download](https://www.anaconda.com/download).

Once you installed the software, at the beginning of your terminal you should see a
```
(base) $
```
The (base) identifies the current anaconda environment. It may be necessary to restart the terminal after the installation to have Anaconda start properly.

If you are **not** using anaconda but the default python from the linux distribution, it may be necessary to install the python header files to correctly compile the SSCHA extension
```
sudo apt install python-dev
```

## 1.2 Install python packages

Most of the code can be easily installed with the following pip command:
```
pip install ase spglib
```
The Atomic Simulation Environment (ASE) is employed to read and write structure files. SPGLIB is used to recognize the space-group and perform symmetry anaylisys. 

# 2. SSCHA

Once the prerequisites have been installed, python-sscha can be downloaded and installed with
```
pip install cellconstructor python-sscha
```

Alternatively, it is possible to use the most recent version from the [SSCHA GitHub](https://github.com/SSCHAcode) repository, under CellConstructor and python-sscha repositories. The installation is performed in this case with
```
python setup.py install
```

## 2.1 Personalize the compiler

If you have multiple compilers installed, and want to force pip to employ a specific fortran compiler, you can specify its path in the FC environment variable. Remember that the compiler employed to compile the code should match with the linker, indicated in the LDSHARED variable.

For example
```
FC=gfortran LDSHARED=gfortran pip install cellconstructor python-sscha
```
For the development version of the code, subtitute the pip call with the python setup.py install.

## 2.2 Running the testsuite

To be sure everything is working, you can run the testsuite. Make sure to install the pytest package with
```
pip install pytest
```
Then run the testsuite with
```
cellconstrutor_test.py
```
If it works without errors, then the code has been correctly installed.

# 3. TDSCHA

As for the SSCHA code, also TDSCHA (time-dependent SCHA) is distributed on PyPi
```
pip install tdscha
```
Alternatively, it can be downloaded from [GitHub](https://github.com/SSCHAcode/tdscha).

To install the GitHub code that enables the MPI parallelization also without the JULIA speedup, you can use:
```
git clone https://github.com/SSCHAcode/tdscha.git
cd tdscha
MPICC=mpicc python setup.py install
```
where mpicc is a valid mpi c compiler (the specification of MPICC can be dropped, but parallelization will not be available aside for the julia mode discussed below).

## 3.1 JULIA speedup enhancement

The TDSCHA code exploits JULIA to speedup the calculation by a factor of 10x-15x with the same number of processors.

To have it working, download and install julia from [https://julialang.org/downloads/](https://julialang.org/downloads/). Alternatively, to install julia on linux we can employ juliaup:
```
curl -fsSL https://install.julialang.org | sh
```
Hit enter when asked to install julia.

To use julia, either open a new terminal, or hit:
```
source ~/.bashrc
```
Then, open a terminal and type `julia`. Inside the julia prompt, type `]`. The prompt should change color and display the julia version ending with `pkg>`

Install the required julia libraries
```
pkg> add SparseArrays, LinearAlgebra, InteractiveUtils, PyCall
```
This should install the required libraries. Press backspace to return to the standard julia prompt and exit with
```
julia> exit()
```
Then, install the python bindings for julia with
```
pip install julia
```
Now, you should be able to exploit the julia speedup in the TDSCHA calculations. It is not required to install julia before TDSCHA, it can also be done in a later moment.

## 3.2 MPI Parallelization

MPI parallelization is not necessary, however you may like to configure it in practical calculations to further speedup the code. For production runs, it is suggested to combine the mpi parallelization with the julia speedup.

The TDSCHA code exploits the mpi parallelization using mpi4py. This assumes that you have a MPI C compiler installed. This is done by installing the library `openmpi-bin`, installed in the requirements.

You can now install mpi4py
```
pip install mpi4py
```
The parallelization is automatically enabled in the julia version and if mpi4py is available. However, to run the parallel code without the julia speedup, you need to recompile the code from the github repository as (not the version installed with pip)
```
MPICC=mpicc python setup.py install
```
Make sure that at the end of the installation no error is displayed, and the write PARALLEL ENVIRONMENT DECTECTED SUCCESFULLY is displayed. Note that, if using the julia enhanced version, the last command is not required, and you can install only mpi4py.

# 4. Install qe-5.1.0_elph

In order to install this old version of Quantum Espresso, which is tuned to allow the combination of electron-phonon matrix elements with SSCHA dynamical matrices, follow these instructions:
```
git clone https://github.com/SSCHAcode/qe-5.1.0_elph.git
cd qe-5.1.0_elph
./configure
make all
```

It may happen that the compilation fails with a message like
```
Error: Rank mismatch between actual argument at...
```
In this case you need to edit the make.sys file with the following command
```
sed -i "s/FFLAGS         = -O3 -g/FFLAGS         = -O3 -g -fallow-argument-mismatch/g" make.sys
```
and rerun
```
make all
```
again.

# 5. F3ToyModel installation

F3ToyModel is a force-field that can mimic the physics of ferroelectric transitions in FCC lattices. All preresquisites are met with the SSCHA installation.

The code for this force-field can be downloaded from the SSCHA GitHub page [here](https://github.com/SSCHAcode/F3ToyModel) with the command:
```
git clone https://github.com/SSCHAcode/F3ToyModel.git
```
Now enter the F3ToyModel directory and install with:
```
python setup.py install
```
