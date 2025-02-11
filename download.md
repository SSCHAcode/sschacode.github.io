---
layout: page
title: Installation
---

Here we provide the installation guide of the SSCHA code as well as other related codes distributed within the [SSCHA GitHub](https://github.com/SSCHAcode) repository. 

# 1. Easy installation through Anaconda/Mamba

The SSCHA code comes as a python library, with computationally intense part speedup with C, Fortran and Julia. The easiest way to install is through Anaconda ([how to install anaconda](https://www.anaconda.com/download)).

If anaconda is too big, you can alternatively install [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html), a much faster and lightweight alternative to ``conda``. Then, replace ``conda`` with ``micromamba`` in the following commands.

```
conda create -n sscha -c conda-forge python=3.11 gfortran libblas lapack openmpi julia openmpi-mpicc pip numpy scipy spglib=2.2 setuptools=64
conda activate sscha
pip install ase julia mpi4py
pip install cellconstructor python-sscha tdscha
```


This is the safest and best way to install the SSCHA. The first line creates a new pristine Python environment with all the required libraries to compile the source code. The second line activates the newly installed environment. Then, the third command installs the additional dependencies, and the last line compiles and installs the SSCHA code.

To use the SSCHA, you must activate the Python environment with the following:

```
conda activate sscha
```

This installation method should also work on clusters and computers with custom configurations. Remember to activate the ``sscha`` environment even in your submission scripts on clusters.

You must ensure Julia's dependencies are correctly set up to activate Julia's speedup on the SSCHA minimization. To do this, run the following line:

```
python -c 'import julia; julia.install()'
```

Note: this command may fail if you are using micromamba. To solve the issue, you need to manually specify the binary location of micromamba to julia:

```
export CONDA_JL_CONDA_EXE=$HOME/.local/bin/micromamba
```
Replacing ``$HOME/.local/bin/micromamba`` with the path to the micromamba binary if you changed the default.
To make it work after the next login, add the environment variable to the init script
```
echo "export CONDA_JL_CONDA_EXE=$HOME/.local/bin/micromamba" >> $HOME/.bashrc
```

To configure Julia PyCall to work with anaconda (or micromamba), open a Julia shell, typing ``Julia ``. Enter in the package manager by typing ``]``. You should see your prompt turning into a ``pkg>``. Then build the conda extension and compile PyCall.
```
build Conda
add PyCall
```


## Troubleshooting

New Python versions and numpy dropped the support for the automatic Fortran compilation required by cell constructor and python-sscha.
If your installation errors with something similar to
```
ModuleNotFoundError: No module named 'distutils.msvccompiler'
```
Run 
```
pip install --force-reinstall setuptools==64
```
then retry to install cellconstructor and python-sscha.

# 2. Installing without Anaconda 

If you do not have anaconda to handle your dependencies, you need to compile the code manually.

You need a FORTRAN and C compiler with MPI configured. Here, we install all the requirements to properly set up the SSCHA code. To properly compile and install the SSCHA code, you need a FORTRAN compiler and LAPACK/BLAS available.

On Debian-based Linux distribution, all the software required is installed with (Tested on Quantum Mobile and ubuntu 20.04):
```
sudo apt update
sudo apt install libblas-dev liblapack-dev liblapacke-dev gfortran openmpi-bin
```
Note that some of the names of the libraries may change slightly in different Linux versions or on MacOS.

## Python installation

Up to version 1.4 of SSCHA, it supports only python <= 3.11. If you are using the default Python in the system, install the development header files. On Ubuntu, they can be installed with:

```
sudo apt install python-dev
```

If you use Anaconda or Micromamba, they are automatically installed.

## Prerequisites

### Python libraries
We need to install the prerequisites with pip:
```
pip install ase spglib==2.2
```

### Julia speedup
The SSCHA benefits from having Julia installed in the system. If present,
it will be automatically used to speed up the calculation.

To install julia, refer to the official website [julialang.org/downloads/](https://julialang.org/downloads/)
Alternatively, to install Julia on linux we can employ Juliaup:
```
curl -fsSL https://install.julialang.org | sh
```
Hit enter when asked to install Julia.

Then, install the Python bindings for Julia with
```
pip install julia
```


The tdscha extension to compute Raman and IR requires additional julia packages that can be installed within a julia terminal. Update your configuration to have access to the newly installed julia
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
Now, you should be able to exploit the julia speedup in the TDSCHA calculations. It is not required to install julia before TDSCHA, it can also be done in a later moment.



## Compiling SSCHA

Once the prerequisites have been installed, python-sscha can be downloaded and installed with
```
pip install cellconstructor python-sscha
```

Alternatively, it is possible to use the most recent version from the [SSCHA GitHub](https://github.com/SSCHAcode) repository, under CellConstructor and python-sscha repositories. The installation is performed in this case with

```
pip install .
```

### Personalize the compiler

If you have multiple compilers installed, and want to force pip to employ a specific fortran compiler, you can specify its path in the FC environment variable. Remember that the compiler employed to compile the code should match with the linker, indicated in the LDSHARED variable.

For example
```
FC=gfortran LDSHARED=gfortran pip install cellconstructor python-sscha
```
For the development version of the code, substitute the pip call with the python setup.py install.

### Running the test suite

To be sure everything is working, you can run the test suite. Make sure to install the pytest package with
```
pip install pytest
```
Then run the test with
```
cellconstrutor_test.py
```
If it works without errors, then the code has been correctly installed.

## TDSCHA installation

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

###  MPI Parallelization for TDSCHA

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
Make sure that no error is displayed at the end of the installation and that the write PARALLEL ENVIRONMENT DECTECTED SUCCESFULLY is displayed. Note that if using the Julia enhanced version, the last command is not required, and you can install only mpi4py.

# 3. Install qe-5.1.0_elph

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

# 4. F3ToyModel installation

F3ToyModel is a force-field that can mimic the physics of ferroelectric transitions in FCC lattices. All preresquisites are met with the SSCHA installation.

The code for this force-field can be downloaded from the SSCHA GitHub page [here](https://github.com/SSCHAcode/F3ToyModel) with the command:
```
git clone https://github.com/SSCHAcode/F3ToyModel.git
```
Now enter the F3ToyModel directory and install with:
```
pip install .
```
