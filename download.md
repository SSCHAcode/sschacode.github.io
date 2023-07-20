---
layout: page
title: Installation
---

Here we provide the installation guide of the SSCHA code as well as other related codes distributed within the [SSCHA github](https://github.com/SSCHAcode) repository. 

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

Alternatively, it is possible to use the most recent version from the github repository [https://github.com/SSCHAcode](https://github.com/SSCHAcode), under CellConstructor and python-sscha repositories. The installation is performed in this case with
```
python setup.py install
```

### Personalize the compiler[¶](#personalize-the-compiler "Permalink to this headline")

If you have multiple compilers installed, and want to force pip to employ a specific fortran compile, you can specify its path in the FC environment variable. Remember that the compiler employed to compile the code should match with the linker, indicated in the LDSHARED variable.

For example

FC=gfortran LDSHARED=gfortran pip install cellconstructor python-sscha

For the development version of the code, subtitute the pip call with the python setup.py install.

### Running the testsuite[¶](#running-the-testsuite "Permalink to this headline")

To be sure everything is working, you can run the testsuite. Make sure to install the pytest package with

pip install pytest

Then run the testsuite with

cellconstrutor\_test.py

If it works without errors, then the code has been correctly installed.

TDSCHA[¶](#tdscha "Permalink to this headline")
-----------------------------------------------

As for the SSCHA code, also TDSCHA is distributed on PyPi

pip install tdscha

Alternatively, the code to compute Raman and IR spectrum can be downloaded from GitHub at \[[https://github.com/SSCHAcode/tdscha](https://github.com/SSCHAcode/tdscha)\]

To install the github code, that enables the MPI parallelization also without the JULIA speedup, you can use:

git clone https://github.com/SSCHAcode/tdscha.git
cd tdscha
MPICC=mpicc python setup.py install

where mpicc is a valid mpi c compiler (the specification of MPICC can be dropped, but parallelization will not be available aside for the julia mode discussed below).

### JULIA speedup enhancement[¶](#julia-speedup-enhancement "Permalink to this headline")

The TDSCHA code exploits JULIA to speedup the calculation by a factor of 10x-15x with the same number of processors.

To have it working, download and install julia from \[[https://julialang.org/downloads/](https://julialang.org/downloads/)\]. Alternatively, to install julia on linux we can employ juliaup:

curl -fsSL https://install.julialang.org | sh

Hit enter when asked to install julia.

To use julia, either open a new terminal, or hit:

source ~/.bashrc

Then, open a terminal and type `julia`. Inside the julia prompt, type `]` The prompt should change color and display the julia version ending with `pkg>`

Install the required julia libraries

pkg> add SparseArrays, LinearAlgebra, InteractiveUtils, PyCall

This should install the required libraries. press backspace to return to the standard julia prompt and exit with

julia> exit()

Then, install the python bindings for julia with

pip install julia

Now, you should be able to exploit the julia speedup in the TDSCHA calculations. It is not required to install julia before TDSCHA, it can also be done in a later moment.

### MPI Parallelization[¶](#mpi-parallelization "Permalink to this headline")

MPI parallelization is not necessary for the tutorial, however you may like to configure it in practical calculation to further speedup the code. For production runs, it is suggested to combine the mpi parallelization with the julia speedup.

The TDSCHA code exploits the mpi parallelization using mpi4py, This assumes that you have a MPI C compiler installed. This is done by installing the library `openmpi-bin` which we installed in the requirements.

You can now install mpi4py

pip install mpi4py

The parallelization is automatically enabled in the julia version and if mpi4py is available. However, to run the parallel code without the julia speedup, you need to recompile the code from the github repository as (not the version installed with pip)

MPICC=mpicc python setup.py install

e sure that at the end of the installation no error are displayed, and the write PARALLEL ENVIRONMENT DECTECTED SUCCESFULLY is displayed. Note that, if using the julia enhanced version, the last command is not required, and you can install only mpi4py.

Install qe-5.1.0\_elph[¶](#install-qe-5-1-0-elph "Permalink to this headline")
------------------------------------------------------------------------------

In order to install this old version of Quantum Espresso, which is tuned to allow the combination of electron-phonon matrix elements with SSCHA dynamical matrices, follow these instructions:

git clone https://github.com/SSCHAcode/qe-5.1.0\_elph.git
cd qe-5.1.0\_elph
./configure
make all

It may happen that the compilation fails with a message like

Error: Rank mismatch between actual argument at...

In this case you need to edit the make.sys file with the following command

sed -i "s/FFLAGS         = -O3 -g/FFLAGS         = -O3 -g -fallow-argument-mismatch/g" make.sys

and rerun

make all

again.

EPIq[¶](#epiq "Permalink to this headline")
-------------------------------------------

The EPIq code is hosted in a git repository. The last stable version can be downloaded [here.](https://gitlab.com/the-epiq-team/epiq/-/releases/EPIq-1.0)

> Once the source code has been downloaded, unzip the archive and enter the epiq main folder (`cd epiq`). EPIq has very few prerequisites:
> 
> *   `BLAS` and `LAPACK` libraries.
>     
> *   Any MPI fortran compiler ( e.g. `mpif90` for `openmpi` ).
>     
> 
> > Then compile _EPIq_. Enter in the source directory and run `make` as:

cd epiq
make all

In some cases (like in quantum mobile), the compilation may fail. If it fails with error:

gfortran: error: unrecognized command line option ‘-fallow-argument-mismatch’; did you mean ‘-Wno-argument-mismatch’?

This can be fixed replacing `-fallow-argument-mismatch` with `-Wno-argument-mismatch` in the make.sys file. This can be done automatically with the following command:

sed -i 's/-fallow-argument-mismatch/-Wno-argument-mismatch/g' make.sys

Then run again `make all`.

If everything went smoothly, an executable file named epiq.x will be created in the `bin` folder. If the compilation was not successful, this probabily means that the `configure` could not find the necessary libraries/compiler. You should manually modify the make.sys file in order to correctly locate them.

fermisurfer installation[¶](#fermisurfer-installation "Permalink to this headline")
-----------------------------------------------------------------------------------

Fermisurfer is a program for the visualization of Fermi surface resolved physical quantities. First, install preresquisites with:

sudo apt-get install -y libwxgtk3.0-gtk3-dev

Download fermisurfer [here](https://osdn.net/projects/fermisurfer/releases/71529) and extract the tar archive:

tar -xsf fermisurfer-2.1.0.tar.gz

Finally, enter the fermisurfer directory and install with:

./configure
 make
 sudo make install

F3ToyModel installation[¶](#f3toymodel-installation "Permalink to this headline")
---------------------------------------------------------------------------------

F3ToyModel is a force-field that can mimic the physics of ferroelectric transitions in FCC lattices. All preresquisites are met with the SSCHA installation.

The code for this force-field can be downloaded from the SSCHA github [here](https://github.com/SSCHAcode/F3ToyModel) with the command:

> git clone https://github.com/SSCHAcode/F3ToyModel.git

Now enter the F3ToyModel directory and install with:

> python setup.py install

### Related Topics

*   [Documentation overview](index.html)
    *   Previous: [Welcome to Tutorial SSCHA School’s documentation!](index.html "previous chapter")
    *   Next: [Hands-on-session 1 - First SSCHA simulations: free energy and structural relaxations](tutorial_01_first_simulations.html "next chapter")

### Quick search

 

$('#searchbox').show(0);

©2023, Lorenzo Monacelli. | Powered by [Sphinx 4.2.0](http://sphinx-doc.org/) & [Alabaster 0.7.12](https://github.com/bitprophet/alabaster) | [Page source](_sources/installation.rst.txt)


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

To install the code, a fortran compiler is required. We recommend gfortran. 
For example, on Ubuntu 20.04, the fortran prerequisites may be installed with:
```
sudo apt-get install libblas-dev liblapack-dev liblapacke-dev gfortran
```



## Installation

To install the *CellConstructor* and *python-sscha* packages it is recommended to use the last version of anaconda-python2, which cames with all the updated numpy and matplotlib packages already compiled to work in parallel. Moreover the last version of matplotlib will allow the user to modify the plots after they are produced.

Once all the dependencies have been installed, the *CellConstructor* and *python-sscha* codes can be easily installed from the command line as:
```
python setup.py install
```
This command must be executed in the directory where the setup.py script is, insice the *CellConstructor* and *python-sscha* folders. If you installed python in a system directory, administration rights may be requested (add a sudo before the command). 

Installing *CellConstructor* and *python-sscha* in clusters may me more tricky and one needs to adapt the setup.py to the cluster characteristics. For instance, if you use the intel compiler, you need to delete the lapack linking from the setup.py and include -mkl. Note that you must force to use the same linker as the one used for the compilation. For example, specific setup.py scripts are provided with the distribution to install *CellConstructor* easily in FOSS or INTEL clusters. 

### Installation through pip

Alternatively, both *CellConstructor* and *python-sscha* can be installed through pip simply as:
```
pip install CellConstructor
``` 
and
```
pip install python-sscha 
```


### Installation through docker
If you are not able to compile the code, or you want to use it on a cluster, where compilation could be cumbersome, we provide a docker container with the SSCHA code already compiled and installed.

You can download it from the docker hub. To run the docker command you need Docker already installed.
```
docker pull mesonepigreco/python-sscha
```
If you get an error of permissions, you need to add your user to the docker group. 

This can be done with the commands:
```
sudo usermod -aG docker $USER
newgrp docker
```

Most clusters provide docker through a module load command, please, ask the cluster maintainer how to run a docker container on your favorite HPC system.

Once the container is installed, you can access it with

```
docker run -it mesonepigreco/python-sscha
```

the previous command opens a new shell with the python-sscha code installed.
To share the content of the current directory with the container, execute the command with

```
docker run -it -v $PWD:/root mesonepigreco/python-sscha
```

This loads the content of the local directory inside the home (/root) of the container with python-sscha.



