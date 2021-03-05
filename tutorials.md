---
layout: page
title: Tutorials
---

In order to get familiar with the code we strongly suggest that before attempting a real calculation you run the tutorials provided here.

We suggest that you start to run the tutorials in the following order:

1. [**Lead Telluride**](http://sscha.eu/Tutorials/Tutorial_PbTe/)

    In this tutorial you will learn how you can set up a calculation from scratch starting from a CIF file. 

    It is illustrated how the ASE package can be used to prepare the input files to run the initial Born-Oppenheimer structural relaxation and harmonic phonon calculation. It is also illustrated how the output of these calculations can be used to start the SSCHA free energy minimization by creating first the ensemble and, second, calculating the energies, forces, and stress tensors for them. It is described how the latter calculations can be performed in three different ways: by using ASE to perform the calculations locally, by running the DFT calculations manually locally or in a cluster, and setting up an automatic submission to a cluster.

2. [**Lead Telluride structural instability**](http://sscha.eu/Tutorials/StructuralInstability/)

    In this tutorial you will understand how to calculate the free energy Hessian used to determine the stability of the system in the free energy landscape, valid to determine second-order phase transitions such as charge-density wave of ferroelectric transitions.

3. [**Lead Telluride spectral properties**](http://sscha.eu/Tutorials/tutorial_spectral/)

    Here you will learn how to run the SSCHA minimization as an stand-alone program, and also how to calculate the phonon spectral functions, which are in the end what experiments probe. Several approaches to calculate the spectral function are exemplified.   

4. [**Tin Telluride with force fields**](http://sscha.eu/Tutorials/SnTe/)

    In this tutorial you will learn how to automatize a calculation with a python script using a force field. Also how to calculate the with a force field the free energy Hessian at different temperatures by scripting all the calculations. 

5. [**Sulfur hydride**](http://sscha.eu/Tutorials/Automatic_Calculations/)

    In this tutorial you can learn how to automatize a SSCHA minimization. The example works with the H$$_3$$S superconducting compound.

6. [**Lanthanum hydride**](http://sscha.eu/Tutorials/VariableCellRelaxation/)

    In this tutorial, which deals with LaH10, you can learn how to automatize a SSCHA relaxation also considering the lattice degrees of freedom. 
