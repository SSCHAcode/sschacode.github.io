---
layout: page
title: Tutorials with Jupyter
---

The following tutorials are on the form of jupyter notebooks. You can find the notebooks in the Tutorials directory of the [source code](https://github.com/SSCHAcode/python-sscha) from GitHub

We suggest that you start to run the tutorials in the following order:

1. [**Lead Telluride**](#Lead-Telluride)

    In this tutorial you will learn how you can set up a calculation from scratch starting from a CIF file. 

    It is illustrated how the ASE package can be used to prepare the input files to run the initial Born-Oppenheimer structural relaxation and harmonic phonon calculation. It is also illustrated how the output of these calculations can be used to start the SSCHA free energy minimization by creating first the ensemble and, second, calculating the energies, forces, and stress tensors for them. It is described how the latter calculations can be performed in three different ways: by using ASE to perform the calculations locally, by running the DFT calculations manually locally or in a cluster, and setting up an automatic submission to a cluster.

2. [**Lead Telluride structural instability**](#Lead-Telluride-structural-instability)

    In this tutorial you will understand how to calculate the free energy Hessian used to determine the stability of the system in the free energy landscape, valid to determine second-order phase transitions such as charge-density wave of ferroelectric transitions.

3. [**Lead Telluride spectral properties**](#Lead-Telluride-spectral-properties)

    Here you will learn how to run the SSCHA minimization as an stand-alone program, and also how to calculate the phonon spectral functions, which are in the end what experiments probe. Several approaches to calculate the spectral function are exemplified.   

4. [**Tin Telluride with force fields**](#Tin-Telluride-with-force-fields)

    In this tutorial you will learn how to automatize a calculation with a python script using a force field. Also how to calculate the with a force field the free energy Hessian at different temperatures by scripting all the calculations. 

5. [**Sulfur hydride**](#Sulfur-hydride)

    In this tutorial you can learn how to automatize a SSCHA minimization. The example works with the H$$_3$$S superconducting compound.

6. [**Lanthanum hydride**](#Lanthanum-hydride)

    In this tutorial, which deals with LaH$$_{10}$$, you can learn how to automatize a SSCHA relaxation also considering the lattice degrees of freedom.

<a name="Lead-Telluride"></a>
# 1. Lead Telluride

In this tutorial, we will set up, from scratch, a calculation of lead telluride (PbTe), a thermoelectric material with high thermoelectric efficiency. All the files needed for the calculation are in the directory Tutorials/PbTe of the *python-sscha* package.

Here, for carrying out the calculations, we will use Quantum ESPRESSO, ASE, and the SSCHA for thermodynamic properties. To setup ASE to work with espresso, please refer to the official [guide](https://wiki.fysik.dtu.dk/ase/ase/calculators/espresso.html).

We prepared the tutorial with the experimental structure in CIF format. You can download the starting CIF files from online databases or use ASE to build your structure.
In this case, we downloaded the structure "PbTe.cif"  from the [American Mineralogist Crystal Structure Database](http://rruff.geo.arizona.edu/AMS/amcsd.php).
The pseudopotentials used for the DFT calculations are in the pseudo_espresso folder. A full list of pseudopotentials is available on the [Quantum ESPRESSO website](https://www.quantum-espresso.org/pseudopotentials).

## Preparation
Both ASE and the SSCHA work in python, so we need to import them.
Moreover, since we want to use Quantum ESPRESSO, we also import the Quantum ESPRESSO calculator of ASE.
Indeed, you can replace the Quantum ESPRESSO calculator with any ASE calculator that can compute total energies, forces, and the stress tensors (the latter required only for cell relaxation). 

```python
# Load in the notebook all scientific libraries
# It can be replaced with:
# from numpy import *
# import numpy as np
# from matplotlib.pyplot import *
%pylab
```


```python
from __future__ import print_function
from __future__ import division
import sys,os

import ase
from ase.calculators.espresso import Espresso
from ase.visualize import view

# We import the basis modules for the SSCHA
import cellconstructor as CC
import cellconstructor.Structure
import cellconstructor.Phonons

# Import the SSCHA engine (we will use it later)
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
```

We will now set up the Quantum ESPRESSO calculator for ASE.  For now on, we will use it to relax the original structure within a simple DFT calculation (this is an experimental structure). Then, we will feed this calculator into the SSCHA code to perform the SSCHA minimization.

```python
# Lets define the pseudopotentials
pseudos = {"Pb": "Pb.upf",
          "Te": "Te.upf"}

# Now we define the parameters for the espresso calculations
input_params = {"ecutwfc" : 60, # The plane-wave wave-function cutoff
               "ecutrho": 240, # The density wave-function cutoff,
               "conv_thr": 1e-6, # The convergence for the DFT self-consistency
               "pseudo_dir" : "pseudo_espresso", # The directory of the pseudo potentials
               "tprnfor" : True, # Print the forces
               "tstress" : True # Print the stress tensor
               }

k_spacing = 0.2 #A^-1 The minimum distance in the Brillouin zone sampling

espresso_calc = Espresso(input_data = input_params, pseudopotentials = pseudos, kspacing = k_spacing)
```

We need to import the structure. 


```python
PbTe_atoms = ase.io.read("PbTe.cif")

# We can view the structure
view(PbTe_atoms)
```

As you may have noticed, this structure is in the conventional cell.
While it is useful for visualization purposes, it makes the SSCHA calculation harder, as more atoms are in the unit cell.
So we redefine the primitive cell with the *CellConstructor* package.
From an easy check on the structure, it is possible to recognize that the primitive vectors $$\vec{v'}$$ 
can be obtained from the conventional vectors $$\vec{v}$$ as follows:

$$
\vec {v'}_1 = \frac 12 \left(\vec v_1 + \vec v_2\right)
$$
$$
\vec {v'}_2 = \frac 12 \left(\vec v_1 - \vec v_2\right)
$$
$$
\vec {v'}_3 = \frac 12 \left(\vec v_1 + \vec v_3\right)
$$


```python
# Initialize a Cellconstructor Structure
struct = CC.Structure.Structure()
struct.generate_from_ase_atoms(PbTe_atoms)

# Define the new unit cell
new_cell = struct.unit_cell.copy()
new_cell[0,:] = .5 * struct.unit_cell[0,:] + .5*struct.unit_cell[1,:]
new_cell[1,:] = .5 * struct.unit_cell[0,:] - .5*struct.unit_cell[1,:]
new_cell[2,:] = .5 * struct.unit_cell[0,:] + .5*struct.unit_cell[2,:]

# Apply the new unit cell to the structure
# And remove duplicated atoms
struct.unit_cell = new_cell
struct.fix_coords_in_unit_cell()
PbTe_primitive = struct.get_ase_atoms()

view(PbTe_primitive)
```

You can see that the new unit cell is much smaller than the previous one, and we have only two atoms per cell, as expected from a rock salt structure.

We can relax the volume at ambient pressure with a variable cell relaxation in Quantum ESPRESSO. For this purpose, we will employ the ASE Espresso calculator to get familiar with it. Indeed the same operation can be performed using a standard espresso input or your favorite DFT program.

```python
# We override (or add) the key for the calculation type
input_params["calculation"] = "vc-relax"
# Generate once again the Espresso calculator
espresso_calc = Espresso(input_data = input_params, pseudopotentials = pseudos, kspacing = k_spacing)

# We attach the calculator to the cell
PbTe_primitive.set_calculator(espresso_calc)

# We run the relaxation
equilibrium_energy = PbTe_primitive.get_total_energy() # It should take few minutes
```


```python
# Ase should have created an input file espresso.pwi that you can check
# And redirect the output through espresso.pwo
# This relaxation should have changed the structure
# So we need to reload the structure from the output espresso file
PbTe_final= ase.io.read("espresso.pwo")

print("Volume before optimization: ", PbTe_primitive.get_volume(), " A^3")
print("Volume after optimization: ", PbTe_final.get_volume(), " A^3")
```

    Volume before optimization:  63.71002599999998  A^3
    Volume after optimization:  67.90591056326011  A^3


Now we have a structure in the primitive unit cell and relaxed.
We remark that the convergence parameters chosen for the minimization are too low to get accurate results, especially for the variable cell relaxation. Here the value of the pressure is overestimated by about 5 GPa; increase the wave-function cutoff to 70 Ry and the density cutoff to 280 Ry to get a more accurate result. We selected under converged parameters just to run very fast on a single processor as a demonstration. 

### Harmonic calculation

Now we need to perform a harmonic calculation to get the spectra.
This can be done with perturbation theory with Quantum ESPRESSO, directly with the ASE library, if we want to use the finite displacement approach, or with phonopy, if you want to correctly exploit the symmetries of the crystal.

In this example, we will show how to do this by exploiting Quantum ESPRESSO perturbation theory on DFT.
We will compute the dynamical matrix on a 2x2x2 q mesh (the equivalent of using a 2x2x2 supercell with the finite displacement approach)
and we will also compute effective charges as the system is an insulator.
We will then save the results in harmonic_dyn filenames.


```python
# We need an input file for phonon calculation in espresso.
ph_input = """
&inputph
    ! the final filename
    fildyn = "harmonic_dyn"
    
    ! the q mesh
    ldisp = .true.
    nq1 = 2 
    nq2 = 2
    nq3 = 2
    
    ! compute also the effective charges and the dielectric tensor
    epsil = .true.
&end
"""

# We write the input script and execute the phonon calculation program
with open("harmonic.phi", "w") as f:
    f.write(ph_input)

# Run the calculation
cmd = "mpirun -np 2 ph.x -npool 2 -i harmonic.phi > harmonic.pho"
res = os.system(cmd) # This command took with 4 i-7 processors about 5 minutes
```

You can see that espresso generated several files:
 * harmonic_dyn0
 * harmonic_dyn1 
 * ...
 
We import them with *CellConstructor*.

We need only to know the number of total q independent points by symmetries, that is,
the total number of dynamical matrix printed by Quantum ESPRESSO (harmonic_dyn0 does not count).
In this case, it should be 3.


```python
# Load the dynamical matrices just computed
harm_dyn = CC.Phonons.Phonons("harmonic_dyn", nqirr =3) # We load harmonic_dynX with X = 1,2,3

# Now we can print all the phonon frequencies
w_s, pols = harm_dyn.DiagonalizeSupercell()
# By default, the dynamical matrix is written in Ry/bohr^2, so the energy is in Ry
# To transform in cm-1 we can use the conversion factor provided by CellConstructor
print ("\n".join(["{:.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in  w_s]))
```

    -186.1459 cm-1
    -186.1459 cm-1
    -186.1459 cm-1
    -185.9474 cm-1
    -185.9474 cm-1
    -185.9474 cm-1
    -179.9445 cm-1
    -179.9445 cm-1
    -179.9445 cm-1
    -179.9445 cm-1
    -179.9445 cm-1
    -179.9445 cm-1
    -101.7456 cm-1
    -101.7456 cm-1
    -101.7456 cm-1
    -80.2691 cm-1
    -80.2691 cm-1
    -80.2691 cm-1
    -80.2691 cm-1
    -80.2691 cm-1
    -80.2691 cm-1
    48.0931 cm-1
    48.0931 cm-1
    48.0931 cm-1
    52.7799 cm-1
    52.7799 cm-1
    52.7799 cm-1
    52.7799 cm-1
    52.7799 cm-1
    52.7799 cm-1
    52.7799 cm-1
    52.7799 cm-1
    55.6362 cm-1
    55.6362 cm-1
    55.6362 cm-1
    55.6362 cm-1
    55.6362 cm-1
    55.6362 cm-1
    55.6362 cm-1
    55.6362 cm-1
    77.3157 cm-1
    77.3157 cm-1
    77.3157 cm-1
    77.3157 cm-1
    99.2432 cm-1
    99.2432 cm-1
    99.2432 cm-1
    99.2432 cm-1


What a mess?! many frequencies are imaginary (the negative ones).
We do not even have the acoustic frequencies at gamma equal to zero.

This is due to two reasons:

1. We have not converged parameters for the ab-initio simulation, so some of them may be an artifact.
2. PbTe may have a ferroelectric instability, so the system could be unstable.

If the system has a structural phase transition, the high-symmetry structure at $$T=0$$ K wants to distort the cell to break the symmetries, so no matter how well we converge the phonon calculations,
some imaginary frequencies will not be removed. 

In this case, the structure is in a saddle-point of the Born-Oppenheimer (BO) energy landscape. 
Many materials similar to PbTe are ferroelectric and have a structural phase transition (as SnSe, SnTe, etc.). PbTe is an exception, as it exhibits only an incipient ferroelectric transition but the structure is stable at the harmonic level (therefore, all the imaginary frequencies we found are just due to under converged parameters in the calculations).

Luckily, the SSCHA can deal with systems with imaginary frequencies, even if they are physically meaningful (instability).

Note that the presence of imaginary frequencies prevents the application of any harmonic approximation.

## The Self-Consistent Harmonic Approximation

Now we have all the ingredients to start the SSCHA calculation.
The SSCHA finds the optimal gaussian density matrix that minimizes the free energy of the system:
$$
\rho_{\mathcal R, \Upsilon}(\vec R) = \sqrt{\det(\Upsilon / 2\pi)} \exp \left[ -\frac 12 \sum_{\alpha\beta} (R_\alpha - \mathcal{R}_\alpha)\Upsilon_{\alpha\beta}(R_\beta - \mathcal{R}_\beta)\right]
$$
where $$\mathcal{R}$$ and $$\Upsilon$$ are, respectively, the average centroid position and the covariance matrix of the gaussian, while $$\alpha,\beta=1\cdots 3N$$ runs on both the atomic and the cartesian coordinates.

In the specific case, the $$\rho_{\mathcal R, \Upsilon}(\vec R)$$ density matrix can be represented by a positive definite dynamical matrix.
The equilibrium density matrix of any harmonic Hamiltonian is a Gaussian. Our density matrix is related one to one to
an auxiliary harmonic hamiltonian $$\mathcal H_{\mathcal R, \Upsilon}$$.
$$
\rho_{\mathcal R, \Upsilon} \longleftrightarrow \mathcal H_{\mathcal R, \Upsilon}
$$

### Preparation of the data

First, we need to obtain a good starting point for our density matrix $$\rho$$.
We can use the harmonic dynamical matrix we computed on the previous step (harmnic_dyn).

There is a problem: a harmonic hamiltonian must be positive definite to generate a density matrix, and the PbTe harmonic Hamiltonian is not.
The *CellConstructor* can force a dynamical matrix to be positive definite.

```python
# Load the dynamical matrix 
dyn = CC.Phonons.Phonons("harmonic_dyn", nqirr = 3)

# Apply the sum rule and symmetries
dyn.Symmetrize()

# Flip the imaginary frequencies into real ones
dyn.ForcePositiveDefinite()

# We can print the frequencies to show the magic:
w_s, pols = dyn.DiagonalizeSupercell()
print ("\n".join(["{:.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in  w_s]))
```

    0.0000 cm-1
    0.0000 cm-1
    0.0000 cm-1
    52.7799 cm-1
    52.7799 cm-1
    52.7799 cm-1
    52.7799 cm-1
    52.7799 cm-1
    52.7799 cm-1
    52.7799 cm-1
    52.7799 cm-1
    55.6362 cm-1
    55.6362 cm-1
    55.6362 cm-1
    55.6362 cm-1
    55.6362 cm-1
    55.6362 cm-1
    55.6362 cm-1
    55.6362 cm-1
    77.3157 cm-1
    77.3157 cm-1
    77.3157 cm-1
    77.3157 cm-1
    80.2691 cm-1
    80.2691 cm-1
    80.2691 cm-1
    80.2691 cm-1
    80.2691 cm-1
    80.2691 cm-1
    99.2432 cm-1
    99.2432 cm-1
    99.2432 cm-1
    99.2432 cm-1
    101.7456 cm-1
    101.7456 cm-1
    101.7456 cm-1
    179.9445 cm-1
    179.9445 cm-1
    179.9445 cm-1
    179.9445 cm-1
    179.9445 cm-1
    179.9445 cm-1
    186.1459 cm-1
    186.1459 cm-1
    186.1459 cm-1
    187.8042 cm-1
    187.8042 cm-1
    187.8042 cm-1


As you can see, we eliminated the imaginary frequencies. This is no more the harmonic dynamical matrix, however, it is a good starting point for defining the density matrix $$\rho_{\mathcal R, \Upsilon}(R)$$.

### The stochastic approach

To solve the self-consistent harmonic approximation we use a stochastic approach: we generate random ionic configurations distributed according to $$\rho_{\mathcal R, \Upsilon}$$ so that we can compute the free energy as
$$
F_{\mathcal R, \Upsilon} = F[{\mathcal H}_{\mathcal R, \Upsilon}] + \left<V(R) - \mathcal V(R)\right>_{\rho_{\mathcal R, \Upsilon}}
$$

Here $$F[\mathcal H_{\mathcal R, \Upsilon}]$$ is the free energy of the auxiliary harmonic hamiltonian $$\mathcal H_{\mathcal R, \Upsilon}$$, 
$$V(R)$$ is the BO energy landscape, $$\mathcal V(R)$$ is the harmonic potential of the auxiliary hamiltonian $$\mathcal H_{\mathcal R, \Upsilon}$$ and $$\left<\cdot\right>_{\rho_{\mathcal R, \Upsilon}}$$ is the average over the $$\rho_{\mathcal R, \Upsilon}$$ density matrix.

So we need to generate the ensemble randomly distributed according to $$\rho_{\mathcal R, \Upsilon}$$. This is done as follows:


```python
# We setup an ensemble for the SSCHA at T = 100 K using the density matrix  from the dyn dynamical matrix
ensemble = sscha.Ensemble.Ensemble(dyn, T0 = 100, supercell= dyn.GetSupercell())

# We generate 10 randomly displaced structures in the supercell
ensemble.generate(N = 10) 

# We can look at one configuration to have an idea on how they look like
# Hint: try to raise the temperature to 1000 K or 2000 K if you want to see a bigger distortion on the lattice
view(ensemble.structures[0].get_ase_atoms())
```

To minimize the free energy, we need to compute the quantity $$V(R)$$ (the energy) and its gradient (the atomic forces) for each ionic configuration generated within the ensemble. 

This can be done in many ways. You can do it manually:
You can save the ensemble as text files with the atomic coordinates and cell, copy that on a cluster, run your DFT calculations, parse the output files and tell the SSCHA code the values of energy and forces for each configuration manually. 
This is very general, in this way you can use whatever program to calculate energies and forces, but it requires a lot of interaction between the SSCHA code and the user.

A more interesting option is the possibility to run the calculation automatically with an ASE calculator. This is achieved by passing to the ensemble object the ASE calculator; it will do the calculation automatically. The drawback is that the calculations are executed on the same computer as the SSCHA code is installed. This can be an issue if you prefer running the SSCHA in your laptop, while you want to rely on an external cluster to perform the heavy DFT calculations.

It is also possible to configure a cluster for a remote calculation so that the code will automatically connect to the cluster, submit the calculations, and retrieve the results without any further interaction with the user. This latter option is the recommended one for large production runs. 

In the following subsections, We will follow all the possibilities, you can jump directly to the section you are more interested in.

### The ASE automatic calculation (local)
We start with the easiest possible thing. We already performed the calculation using ASE for the cell relaxation, we will set up it
now for a simple total energy DFT calculation, attach to our ensemble variable and compute the ensemble.

**NOTE**: Sometimes calculators fail randomly, the ensemble class will resubmit a single calculation for 5 times if it fails. If 5 consecutive fails are found, then the code raises an exception, showing the ASE error.



```python
# Lets setup the espresso calculator for ASE
# (You can substitute the following two lines with the calculator you prefer)
input_params["calculation"] = "scf" # Setup the simple DFT calculation

espresso_calc = Espresso(pseudopotentials = pseudos, input_data = input_params, kspacing = k_spacing, koffset = k_offset)

# Now we use the espresso calculator to compute all the configurations
ensemble.compute_ensemble(espresso_calc)
```

    conf 0 / 10
    conf 1 / 10
    conf 2 / 10
    conf 3 / 10
    conf 4 / 10
    conf 5 / 10
    conf 6 / 10
    conf 7 / 10
    conf 8 / 10
    conf 9 / 10



```python
# Now we can save the ensemble, to reload it in any later moment
ensemble.save("data_ensemble_ASE", population = 1)
```

### The manual calculation (local or on a cluster)
To perform a manual calculation we save each ionic configuration in a text file (atomic coordinates and the lattice vectors). 
You will need to parse them with a custom script to prepare appropriate input files for your favorite calculator and eventually send them in a cluster.
In this example, we will use the standard pw.x executable of quantum espresso for the DFT calculation, without relying on ASE. You can adapt it to your needs.




```python
# We save the ensemble as it is, before computing the forces
# The population flag allows us to save several ensembles inside the same directory
# As the SSCHA is an iterative algorithm, we may need more steps to converge,
# in this way we can save all the ensembles of the same SSCHA calculation inside the same directory
ensemble.save("data_ensemble_manual", population = 1)
```

Now let us have a look at the files created by the code in the ensemble directory 'data_ensemble_manual'.
We will notice files made as:

    u_population1_X.dat  
    scf_population1_X.dat   

Both files represent the ionic configurations. The Xs are the IDs (starting from 1) of the configurations. In the u_population1_X.dat file is stored the atomic displacement of each ion with respect to the centroid position. The scf_population1_X.dat is a file that contains the whole ionic displaced structure. In particular, in the first part, we have the atomic coordinates preceded by the atomic type (in Angstrom), while in the last part we have the lattice vectors (in Angstrom). This is the standard format to be used in Quantum ESPRESSO with ibrav=0, however, since it should not have any symmetry (is a randomly displaced configuration) it can be easily converted for input to any calculator.

So, we will use the scf_population1_X.dat to prepare our input scripts for Quantum ESPRESSO.

After we will get total energies, forces, and stresses for any ionic configurations, we must save the results in the following files:

    forces_population1_X.dat   
    pressures_population1_X.dat    
    energies_supercell_population1.dat


forces_population1_X.dat for each row must be filled with the force (in Ry/Bohr) on the corresponding atom of the same configuration, pressure_population1_X.dat must be filled with the 3x3 symmetric stress tensor (in Ry/Bohr^3) of the X configuration, and energies_supercell_population1.dat is a one-column file that contains, for each row, the total energy (in Ry) on the supercell of each configuration (in order of ID).

In the following code we will prepare the Quantum ESPRESSO pw.x input for each configuration.



```python
# Ok now we will use python to to parse the ensemble and generate input files for the ab-initio run

# First of all, we must prepare a generic header for the calculation
typical_espresso_header = """
&control
    calculation = "scf"
    tstress = .true.
    tprnfor = .true.
    disk_io = "none"
    pseudo_dir = "pseudo_espresso"
&end
&system
    nat = {}
    ntyp = 2
    ibrav = 0
    ecutwfc = 40
    ecutrho = 160
&end
&electrons
    conv_thr = 1d-6
    !diagonalization = "cg"
&end

ATOMIC_SPECIES
Pb 207.2 Pb.upf
Te 127.6 Te.upf

K_POINTS automatic
1 1 1 0 0 0
""".format(ensemble.structures[0].N_atoms) 
# We extract the number of atoms form the ensemble and the celldm(1) from the dynamical matrix (it is stored in Angstrom, but espresso wants it in Bohr)
# You can also read it on the fourth value of the first data line on the first dynamical matrix file (dyn_start_popilation1_1); In the latter case, it will be already in Bohr.

# Now we need to read the scf files
all_scf_files = [os.path.join("data_ensemble_manual", f) for f in os.listdir("data_ensemble_manual") if f.startswith("scf_")]

# In the previous line  I am reading all the files inside data_ensemble_manual os.listdir(data_ensemble_manual) and iterating over them (the f variable)
# I iterate only on the filenames that starts with scf_ 
# Then I join the directory name data_ensemble_manual to f. In unix it will be equal to data_ensemble_manual/scf_....
# (using os.path.join to concatenate path assure to have the correct behaviour independently on the operating system

# We will generate the input file in a new directory
if not os.path.exists("run_calculation"):
    os.mkdir("run_calculation")

for file in all_scf_files:
    # Now we are cycling on the scf_ files we found.
    # We must extract the number of the file
    # The file is the string "data_ensemble_manual/scf_population1_X.dat"
    # Therefore the X number is after the last "_" and before the "." character
    # We can split before the string file at each "_", isolate the last part "X.dat"
    # and then split it again on "." (obtaining ["X", "dat"]) and select the first element
    # then we convert the "X" string into an integer
    number = int(file.split("_")[-1].split(".")[0])
    
    # We decide the filename for the espresso input
    # We will call it run_calculation/espresso_run_X.pwi
    filename = os.path.join("run_calculation", "espresso_run_{}.pwi".format(number))
    
    # We start writing the file
    with open(filename, "w") as f:
        # We write the header
        f.write(typical_espresso_header)
        
        # Load the scf_population_X.dat file
        ff = open(file, "r")
        structure_lines = ff.readlines()
        ff.close()
        
        # Write the content on the espresso_run_X.pwi file
        # Note in the files we specify the units for both the cell and the structure [Angstrom]
        f.writelines(structure_lines)
        
```

The previous code will generate the input files inside a directory run_calculation for quantum espresso. Now we just need to run them.

You can copy them to your cluster, submit a job, and copy back the results.
For completeness, the next cell of code will perform the calculation locally, it could require some time.
If you pass them to a supercomputer, remember to copy also the pseudopotentials!


```python
directory = "run_calculation"
# Copy the pseudo
for file in os.listdir(directory):
    # Skip anything that is not an espresso input file
    if not file.endswith(".pwi"):
        continue 
        
    outputname = file.replace(".pwi", ".pwo")
    
    total_inputname = os.path.join(directory, file)
    total_outputname = os.path.join(directory, outputname)
    
    # Run the calculation (on 4 processors)
    cmd = "mpirun -np 4 pw.x -i {} > {}".format(total_inputname, total_outputname)
    print("Running: ", cmd) # On my laptop with 4 processors (i7) it takes about 1 minute for configuration
    os.system(cmd)
```

    Running:  /usr/bin/mpirun -np 4 pw.x -i run_calculation/espresso_run_3.pwi > run_calculation/espresso_run_3.pwo
    Running:  /usr/bin/mpirun -np 4 pw.x -i run_calculation/espresso_run_9.pwi > run_calculation/espresso_run_9.pwo
    Running:  /usr/bin/mpirun -np 4 pw.x -i run_calculation/espresso_run_4.pwi > run_calculation/espresso_run_4.pwo
    Running:  /usr/bin/mpirun -np 4 pw.x -i run_calculation/espresso_run_8.pwi > run_calculation/espresso_run_8.pwo
    Running:  /usr/bin/mpirun -np 4 pw.x -i run_calculation/espresso_run_10.pwi > run_calculation/espresso_run_10.pwo
    Running:  /usr/bin/mpirun -np 4 pw.x -i run_calculation/espresso_run_6.pwi > run_calculation/espresso_run_6.pwo
    Running:  /usr/bin/mpirun -np 4 pw.x -i run_calculation/espresso_run_7.pwi > run_calculation/espresso_run_7.pwo
    Running:  /usr/bin/mpirun -np 4 pw.x -i run_calculation/espresso_run_1.pwi > run_calculation/espresso_run_1.pwo
    Running:  /usr/bin/mpirun -np 4 pw.x -i run_calculation/espresso_run_2.pwi > run_calculation/espresso_run_2.pwo
    Running:  /usr/bin/mpirun -np 4 pw.x -i run_calculation/espresso_run_5.pwi > run_calculation/espresso_run_5.pwo


Either we run the calculation locally or we copied in a cluster and ran there, now we have the output files inside the directory run_calculations, we need to retrieve the energies, forces, and stress tensors (if any).

As written, we must convert the total energy of the supercell in Ry, the forces in Ry/Bohr, and the stress in Ry/Bohr^3.
Luckily quantum espresso already gives these quantities in the correct units, but be careful when using different calculators.
This problem does not arise when using automatic calculators, as the SSCHA and ASE will cooperate to convert the units to the correct one.
Now we will parse the Quantum ESPRESSO output looking for the energy, the forces, and the stress tensor.


```python
directory = "run_calculation"
output_filenames = [f for f in os.listdir(directory) if f.endswith(".pwo")] # We select only the output files
output_files = [os.path.join(directory, f) for f in output_filenames] # We add the directory/outpufilename to load them correctly

# We prepare the array of energies
energies = np.zeros(len(output_files)) 
for file in output_files:
    # Get the number of the configuration.
    id_number = int(file.split("_")[-1].split(".")[0]) # The same as before, we need the to extract the configuration number from the filename
    
    # Load the file
    ff = open(file, "r")
    lines = [l.strip() for l in ff.readlines()] # Read the whole file removing tailoring spaces
    ff.close()
    
    # Lets look for the energy (in espresso the first line that starts with !)
    # next is used to find only the first occurrence
    energy_line = next(l for l in lines if len(l) > 0 if l.split()[0] == "!")
    
    # Lets collect the energy (the actual number is the 5th item on the line, but python indexes start from 0)
    # note, also the id_number are saved starting from 1
    energies[id_number - 1] = float(energy_line.split()[4])
    
    # Now we can collect the force
    # We need the number of atoms
    nat_line = next( l for l in lines if len(l) > 0 if l.split()[0] == "number" and l.split()[2] == "atoms/cell" )
    nat = int(nat_line.split()[4])
    
    # Now allocate the forces and read them
    forces = np.zeros((nat, 3))
    forces_lines = [l for l in lines if len(l) > 0 if l.split()[0] == "atom"] # All the lines that starts with atom will contain a force
    for i in range(nat):
        forces[i, :] = [float(x) for x in forces_lines[i].split()[-3:]] # Get the last three number from the line containing the force
    
    # Now we can take the stress tensor
    stress = np.zeros((3,3))
    # We pick the index of the line that starts with the words total stress
    index_before_stress = next(i for i, l in enumerate(lines) if len(l) > 0 if l.split()[0] == "total" and l.split()[1] == "stress")
    # The stress tensor is located just after it
    for i in range(3):
        index = i + index_before_stress + 1
        stress[i, :] = [float(x) for x in lines[index].split()[:3]]

    # We can save the forces_population1_X.dat and pressures_population1_X.dat files
    force_file = os.path.join("data_ensemble_manual", "forces_population1_{}.dat".format(id_number))
    stress_file = os.path.join("data_ensemble_manual", "pressures_population1_{}.dat".format(id_number))
    np.savetxt(force_file, forces)
    np.savetxt(stress_file, stress)

# Now we read all the configurations, we can save the energy file
energy_file = os.path.join("data_ensemble_manual", "energies_supercell_population1.dat")
np.savetxt(energy_file, energies)
```

Now we have collected back the ensemble, and we can load it into the sscha once again.
This time, when we load the ensemble, we also need to specify the number of configurations.
If we made some mistakes in naming the files or some files have the wrong number of lines, the code would complain and stop.


```python
ensemble.load("data_ensemble_manual", population = 1, N = 10)
```

### The automatic submission to a cluster

The *python-sscha* code implements also the possibility to send the calculation to a remote cluster.
This option is useful, as installing the full sscha package on a cluster may be cumbersome. 

With this option, we will run all the SSCHA locally, and submit to the cluster only the DFT calculations.
At the current moment, the cluster automatic submission is compatible only with the SLURM queue system, and for now, only Quantum ESPRESSO input file generation has been tested.

First, we need to initialize the cluster, specifying all the relevant variables.
Here, we use the example of the MARCONI HPC of the CINECA, but can be edited to run on most clusters (tested in the Spanish MARE NOSTRUM, the French IRENE, and EKHI at the CFM)



```python
# Here we configure the cluster object MARCONI BROADWELL (we provide an example for GPU espresso running on MARCONI 100 on a separate directory)
import sscha.Cluster

my_hpc = sscha.Cluster.Cluster(pwd = None)

# We setup the connection info
my_hpc.hostname = "login.marconi.cineca.it@myuser" # The command to connect via ssh to the cluster
my_hpc.account_name = "IscrB_MYACRONYM" # The name of the project for the computation
my_hpc.workdir = "/marconi_work/IscrB_MYACRONYM/myuser/workdir" # the directory in which the calculations are performed

# Now we need to setup the espresso
# First we must tell the cluster where to find him:
my_hpc.binary = "pw.x -npool NPOOL -i  PREFIX.pwi > PREFIX.pwo"
# Then we need to specify if some modules must be loaded in the submission script
my_hpc.load_modules = """
# Here this is a bash script at the beginning of the submission
# We can load modules

module load espresso
"""

# All these information are independent from the calculation
# Now we need some more specific info, like the number of processors, pools and other stuff
my_hpc.n_cpu = 32 # We will use 32 processors
my_hpc.n_nodes = 1 #In 1 node
my_hpc.n_pool = 4 # This is an espresso specific tool, the parallel CPU are divided in 4 pools

# We can also choose in how many batch of jobs we want to submit simultaneously, and how many configurations for each job
my_hpc.batch_size = 40
my_hpc.job_number = 10
# In this way we submit 40 jobs, each one with 10 configurations

# We can specify the time limit for each job,
my_hpc.time = "00:30:00" # 30 minutes
```


```python
# Now we can compute the ensemble using this cluster configuration, in a similar way we made with ASE:
# We generate the ASE espresso calculator
input_params["calculation"] = "scf" # Setup the simple DFT calculation

espresso_calc = Espresso(pseudopotentials = pseudos, input_data = input_params, kspacing = k_spacing, koffset = k_offset)

# Now we use the espresso calculator with the cluster to send the calculation in the remote cluster.
ensemble.compute_ensemble(espresso_calc, cluster = my_hpc)
```

The last line will do the job. It will establish a connection to the cluster, prepare all the input files, send them, submit them through a SLURM queue system, and collect them back. 
The interesting feature of this procedure is that if some calculation fails, the automatic process will automatically resubmit it up to five times (this number can be changed as it is a property of the cluster, to avoid wasting CPU time if we setup a too low time limit. The variable is my_hpc.max_recalc.

You will need to connect without manually typing your password to use the cluster in this way. This is done by public/private key encryption with ssh, that is the safest way to go. However, some clusters do not allow connection with private/public keys (or the policy to obtain this kind of access requires a painful bureaucratic procedure) and impose the manual typing of the password each time. In this case, you can act in two ways: by setting up an ssh tunnel to connect (this is also refused by some clusters) or by explicitly telling python the password for the connection when you create the cluster object:
```python
my_hcp = sscha.Cluster.Cluster(pwd="mybeautifullpassword")
```

To run the code in this way, you need the ssh-pass package to be installed on your local computer.
This is not the best thing to do, as your password will be stored in clear text, and easily stolen from you. So before using this command check if it is impossible to establish a connection without the password request.


```python
# We can save the ensemble for further processing
ensemble.save("data_ensemble_cluster", population = 1)
```

## The SSCHA minimization
In the previous sections, we initialized the starting density matrix $$\rho_{\mathcal R,\Upsilon}$$, and we generated a randomly distributed ensemble. Then we used Quantum ESPRESSO to compute the energies, forces, and the stress tensors from it.

In this section we will use this ensemble to optimize the density matrix $$\rho_{\mathcal R, \Upsilon}$$ to minimize the free energy.
The minimization will continue until both we converge to a minimum or our ensemble does not describe sufficiently well anymore the updated $$\rho_{\mathcal R, \Upsilon}$$.
The ensemble degradation occurs when our updated parameters $$\mathcal R$$ and $$\Upsilon$$ change a lot with respect to the original value used to generate the ensemble.

This criterion is established by the Kong Liu ratio (see [Monacelli et al, PRB 98 (2), 024106](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.98.024106) for details on the definition).

In this case, $$\mathcal R$$ are fixed by symmetry, so our only degrees of freedom are in the $$\Upsilon$$ matrix. We force the system to ignore $$\mathcal R$$ for the minimization convergence.


```python
# Lets reset other calculation if you run this cell multiple times
ensemble.update_weights(dyn, 100) # Restore the original density matrix at T = 100 K
minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

# Ignore the structure minimization (is fixed by symmetry)
minimizer.minim_struct = False

# Setup the minimization parameter for the covariance matrix
minimizer.min_step_dyn = 0.05 # Values around 1 are good
#minimizer.precond_dyn = False
#minimizer.root_representation = "root2"

# Setup the threshold for the ensemble wasting
minimizer.kong_liu_ratio = 0.3 # Usually 0.5 is a good value

# Lest start the minimization
minimizer.init()
minimizer.run()
```

The previous calculation will print a lot of info on the standard output as the minimization goes.
However, after it finishes, you can plot the minimization features to have a hint of what happened:


```python
minimizer.plot_results()
```

Here you see several plots as a function of the minimization steps. 
The first one is the total Free energy (including anharmonic contributions) 
The second plot is the modulus of the free energy gradient with respect to the covariance $\Upsilon$ matrix. 
The third plot is the modulus of the free energy gradient with respect to the average atomic position $\mathcal R$.
The final plot is a measurement of the ensemble degrading.
You can notice that the gradient of the atomic position is always 0, this reflects the fact that the average position of the atoms must keep the
frequencies. 

You can also explicitly interrogate the code for quantities like Free energy and Stress tensor now that the system is ended.
To print some generic info you can use the finalize() method. Otherwise, you can access specific quantities:


```python
minimizer.finalize()
```

    
     * * * * * * * * 
     *             * 
     *   RESULTS   * 
     *             * 
     * * * * * * * * 
    
    
    Minimization ended after 10 steps
    
    Free energy = -4767586.25400589 +-       1.18792173 meV
    FC gradient modulus =   19709.10738392 +-       0.16384031 bohr^2
    Struct gradient modulus =       0.00000000 +-       0.00000000 meV/A
    Kong-Liu effective sample size =  3.204267119925592
    
    
     ==== STRESS TENSOR [GPa] ==== 
          4.37878027      0.00000000     -0.00000000                0.01248128      0.00000000      0.00000000
          0.00000000      4.37878027     -0.00000000    +-          0.00000000      0.01248128      0.00000000
         -0.00000000      0.00000000      4.37878027                0.00000000      0.00000000      0.01248128
    
     Ab initio average stress [GPa]:
          4.36834543      0.00000000     -0.00000000
          0.00000000      4.36834543     -0.00000000
         -0.00000000     -0.00000000      4.36834543
    



```python
print("The total free energy per unit cell is:", minimizer.get_free_energy(), " Ry")
print("The total stress tensor is [Ry/bohr^3]:")
print(minimizer.get_stress_tensor()[0])
print("And the stochastic error on the stress tensor is:")
print(minimizer.get_stress_tensor()[1])
print("The stocastic error of the free energy instead, was:", minimizer.get_free_energy(return_error = True)[1], " Ry")
```

    The total free energy per unit cell is: -350.4109991915717  Ry
    The total stress tensor is [Ry/bohr^3]:
    [[ 2.97663324e-04  2.71050543e-20 -1.35525272e-19]
     [ 2.71050543e-20  2.97663324e-04 -2.71050543e-20]
     [-2.16840434e-19  2.71050543e-20  2.97663324e-04]]
    And the stochastic error on the stress tensor is:
    [[8.48460188e-07 0.00000000e+00 0.00000000e+00]
     [0.00000000e+00 8.48460188e-07 0.00000000e+00]
     [0.00000000e+00 0.00000000e+00 8.48460188e-07]]
    The stocastic error of the free energy instead, was: 8.731060497873162e-05  Ry


## Saving the results

We completed the minimization, now we can save the final results.
In particular, we can save the density matrix $$\rho_{\mathcal R, \Upsilon}$$ (or rather the dynamical matrix and the centroid positions from which we can obtain the density matrix). 

These two quantities are both stored inside the auxiliary dynamical matrix. We can access it through:
```python
minimizer.dyn
```

We can, for example, show a new average structure (in this case will be equal to the beginning one, as it is fixed by symmetry) and print the frequencies of the effective matrix, to see how they changed.


```python
# Draw the 3D structure of the final average atomic positions
view(minimizer.dyn.structure.get_ase_atoms())

# We can save the dynamical matrix
minimizer.dyn.save_qe("dyn_pop1_")

# Print the frequencies before and after the minimization
w_old, p_old = ensemble.dyn_0.DiagonalizeSupercell() # This is the representation of the density matrix used to generate the ensemble
w_new, p_new = minimizer.dyn.DiagonalizeSupercell()

# We can now print them 
print(" Old frequencies |  New frequencies")
print("\n".join(["{:16.4f} | {:16.4f}  cm-1".format(w_old[i] * CC.Units.RY_TO_CM, w_new[i] * CC.Units.RY_TO_CM) for i in range(len(w_old))]))
```

     Old frequencies |  New frequencies
              0.0000 |          -0.0000  cm-1
              0.0000 |          -0.0000  cm-1
              0.0000 |           0.0000  cm-1
             52.7799 |          51.2496  cm-1
             52.7799 |          51.2496  cm-1
             52.7799 |          51.2496  cm-1
             52.7799 |          51.2496  cm-1
             52.7799 |          51.2496  cm-1
             52.7799 |          51.2496  cm-1
             52.7799 |          51.2496  cm-1
             52.7799 |          51.2496  cm-1
             55.6362 |          56.2187  cm-1
             55.6362 |          56.2187  cm-1
             55.6362 |          56.2187  cm-1
             55.6362 |          56.2187  cm-1
             55.6362 |          56.2187  cm-1
             55.6362 |          56.2187  cm-1
             55.6362 |          56.2187  cm-1
             55.6362 |          56.2187  cm-1
             77.3157 |          60.2082  cm-1
             77.3157 |          60.2082  cm-1
             77.3157 |          60.2082  cm-1
             77.3157 |          60.2082  cm-1
             80.2691 |          60.2082  cm-1
             80.2691 |          60.2082  cm-1
             80.2691 |          64.3616  cm-1
             80.2691 |          64.3616  cm-1
             80.2691 |          64.3616  cm-1
             80.2691 |          76.8430  cm-1
             99.2432 |          76.8430  cm-1
             99.2432 |          76.8430  cm-1
             99.2432 |          76.8430  cm-1
             99.2432 |          99.4508  cm-1
            101.7456 |          99.4508  cm-1
            101.7456 |          99.4508  cm-1
            101.7456 |          99.4508  cm-1
            179.9445 |         101.9649  cm-1
            179.9445 |         101.9649  cm-1
            179.9445 |         101.9649  cm-1
            179.9445 |         101.9649  cm-1
            179.9445 |         101.9649  cm-1
            179.9445 |         101.9649  cm-1
            186.1459 |         136.6898  cm-1
            186.1459 |         136.6898  cm-1
            186.1459 |         136.6898  cm-1
            187.8042 |         143.2382  cm-1
            187.8042 |         143.2382  cm-1
            187.8042 |         143.2382  cm-1


## Congratulations! 
You performed your first SSCHA relaxation! 
However, we still are not done. Most likely, the SSCHA stopped when the original ensemble does no longer describe the current density matrix.
We should repeat the procedure to extract a new ensemble with the last density matrix, compute energies and forces and then minimize again.
This procedure should be iterated until we converge.

There are many ways, you can manually iterate the ensemble generation, simply rerunning this notebook (but loading the new dynamical matrix we saved as dyn_pop1_).
However, the SSCHA code offers a module (Relax) that can automatize these iterations for you.
Indeed, in this case, you need to set up an automatic calculator with ASE, locally, or with the cluster option. 

<a name="Lead-Telluride-structural-instability"></a>
# 2. Lead Telluride structural instability

One of the main features of the SSCHA is that it provides a complete theoretical framework to study second-order phase transitions for structural instabilities. Examples of those are materials undergoing a charge-density wave, ferroelectric, or simply a structural transition. An example application for each one of this case with the SSCHA can be found in these papers:
1. [Bianco et. al. Nano Lett. 2019, 19, 5, 3098-3103](https://pubs.acs.org/doi/abs/10.1021/acs.nanolett.9b00504)
2. [Aseguinolaza et. al. Phys. Rev. Lett. 122, 075901](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.075901)
3. [Bianco et. al. Phys. Rev. B 97, 214101](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.214101)

According to Landau's theory of second-order phase transitions, a phase transition occurs when the free energy curvature around the high-symmetry structure on the direction of the order parameter becomes negative:
![](second_order.png)

For structural phase transitions, the order parameter is associated to phonon atomic displacements. So we just need to calculate the Free energy Hessian, as:

$$
\frac{\partial^2 F}{\partial R_a \partial R_b}.
$$

Here, $$a$$ and $$b$$ encode both atomic and Cartesian coordinates.
This quantity is very hard to compute with a finite difference approach, as it would require a SSCHA calculation for all possible atomic displacements (keeping atoms fixed). Also because finite difference approaches are hindered by the stochastic noise in the Free energy. Luckily, the SSCHA provides an analytical equation for the free energy Hessian, derived by Raffaello Bianco in the work [Bianco et. al. Phys. Rev. B 96, 014111](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.014111).
The free energy curvature can be written in matrix form as:

$$
\frac{\partial^2 F}{\partial {R_a}\partial {R_b}} = \Phi_{ab} + \sum_{cdefgh} \stackrel{(3)}{\Phi}_{acd}\Lambda_{cdef}[1 - \Lambda\stackrel{(4)}{\Phi}]^{-1}_{efgh} \stackrel{(3)}{\Phi}_{ghb}
$$

Here, $$\Phi$$ is the SCHA auxiliary force constant matrix obtained by the auxiliary harmonic Hamiltonian, $$\stackrel{(3,4)}{\Phi}$$ are the average of the 3rd and 4th derivative of the Born-Oppenheimer energy landscape on the SSCHA density matrix, while the $$\Lambda$$ tensor is a function of the frequencies of the auxiliary harmonic Hamiltonian.
Fortunately, this complex equation can be evaluated from the ensemble with a simple function call.

Lets see a practical example:

```python
%pylab
# Lets import all the sscha modules
import cellconstructor as CC
import cellconstructor.Phonons
import sscha, sscha.Ensemble

# We load the SSCHA dynamical matrix for the PbTe (the one after convergence)
dyn_sscha = CC.Phonons.Phonons("dyn_sscha", nqirr = 3)

# Now we load the ensemble
ensemble = sscha.Ensemble.Ensemble(dyn_sscha, T0 = 1000, supercell=dyn_sscha.GetSupercell())
ensemble.load("data_ensemble_final", N = 100, population = 5)

# If the SSCHA matrix was not the one used to compute the ensemble
# We must update the ensemble weights
# We can also use this function to simulate a different temperature.
ensemble.update_weights(dyn_sscha, T = 1000)

# ----------- COMPUTE THE FREE ENERGY HESSIAN -----------
dyn_hessian = ensemble.get_free_energy_hessian()
# -------------------------------------------------------

# We can save the free energy hessian as a dynamical matrix in quantum espresso format
dyn_hessian.save_qe("free_energy_hessian")
```

This code will do the trick. We can then print the frequencies of the hessian. If an imaginary frequency is present, then the system wants to spontaneously break the high symmetry phase. The frequencies in the free energy hessian are temperature dependent. Tracking the temperature at which an imaginary frequency appears the temperature at which the second-order phase transition occurs can be determined.

It is important to mention that by default the *bubble* approximation is assumed by the SSCHA code, meaning that in the equation above it is assumed that $$\stackrel{(4)}{\Phi}$$. This is an approximation that it is usually good, but needs to be checked. In order to include the $$\stackrel{(4)}{\Phi}$$ term the call to compute the Hessian needs to be modified as

```python
# ----------- COMPUTE THE FREE ENERGY HESSIAN -----------
dyn_hessian = ensemble.get_free_energy_hessian(include_v4 = True)
# -------------------------------------------------------
``` 

Including the $$\stackrel{(4)}{\Phi}$$ term for large supercells is time and memory consuming.

<a name="Lead-Telluride-spectral-properties"></a>
# 3. Lead Telluride spectral properties

In this tutorial we will learn how to perform static (free energy Hessian) and dynamic (spectral function) SSCHA phonon calculations. We will perform the calculations on PbTe in the rock-salt structure ([Ribeiro et al.,Phys. Rev. B 97, 014306](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.014306)) and, in order to speed-up the calculation, we will emply a force-field model based on the work of [Ai et al, Phys. Rev. B 90, 014308](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.014308). The force-field model can be downloaded and installed from [here](https://github.com/SSCHAcode/F3ToyModel). The material to run the example (dynamical matrices and path file) can be found in the folder python-sscha/Tutorials/Spectral_properties. Unless otherwise specified, the equations and the sections we will refer to are in the paper [Monacelli et al. arXiv:2103.03973](https://arxiv.org/abs/2103.03973) describing the SSCHA program.

We will work with a 4x4x4 supercell. To set up the interatomic force-field, we use the 8 (number of irreducible points in the 4x4x4 grid) matrices PbTe.ff.4x4x4.dyn#q to set up the interaction at harmonic level (obtained from first-principles), and we fix the parameters $$p_3$$, $$p_4$$ and $$p_{4\chi}$$ to set up the 3rd and 4th order anharmonic interaction. 

We perform the minimization of the free energy functional with respect to the auxiliary harmonic matrices (the atomic positions are fixed in the high-symmetry $$Fm\bar{3}m$$ rock-salt configuration) at $$300~\mathrm{K}$$. As we have learnt in previous tutorials, in order to obtain the final SSCHA dynamical matrices we need to perform the minimization in several steps, generating different populations. With this script we generate, from the temporary SSCHA dynamical matrices obtained at the end of a minimization, the next ensamble (a.k.a. population) to continue the minimization and we also compute the energy-forces for elements of this ensamble. The dynamical matrices to be used to generate the first ensamble can be the same dynamical matrices that define the toy-model harmonic force field. The parameters to be inserted in this script are: the population number that is going to be generated ("population", integer), and the number of elements in this ensamble (" n_random", integer).


```python
from __future__ import print_function

# Import the cellconstructor stuff
import cellconstructor as CC
import cellconstructor.Phonons

# Import the modules of the force field
import fforces as ff
import fforces.Calculator

# Import the modules to run the sscha
import sscha, sscha.Ensemble, sscha.SchaMinimizer
import sscha.Relax, sscha.Utilities


# Import Matplotlib to plot
import numpy as np
import matplotlib.pyplot as plt



# ========================= TOY MODEL DEFINITION ===========================
# Dynamical matrices that set up the harmonic part of the force-field
ff_dyn_name="PbTe.ff.4x4x4.dyn"
# Paramters that set up the anharmonic part of the force-field
p3 = -0.0140806002
p4 = -0.01089564958
p4x = 0.00254038964
# ====================================================================

# ======================== TO BE SET ================================
# population to be generated (integer number)
population=  
# dynamical matrices to be used to generate the population
dyn_sscha_name="PbTe.dyn"
# temperature
T=300.
# number of ensamble elements
n_random=
# dir where the population ensamble is stored
ens_savedir="ens_pop"+str(population)
# ====================================================================


INFO = """
We compute the population #{} for PbTe with the force-field defined by 
the harmonic matrices: {} 
and the parameters:
p3=  {}
p4=  {}
p4x= {}

The ensamble will be generated from dyn mat: {}
The temperature is {} K
The number of ensamble elements is {}.
The ensamble with energy/forces will be saved in: {}

""".format(population,ff_dyn_name,p3,p4,p4x,dyn_sscha_name,T,n_random,ens_savedir)


print(INFO)
print()
print(" ======== RUNNING ======== ")
print()


# Setup the harmonic part of the force-field 
ff_dynmat = CC.Phonons.Phonons(ff_dyn_name, 8)
ff_calculator = ff.Calculator.ToyModelCalculator(ff_dynmat)
# Setup the anharmonic part of the force-field
ff_calculator.type_cal = "pbtex"
ff_calculator.p3 = p3
ff_calculator.p4 = p4
ff_calculator.p4x = p4x
# Load matrices
dyn_sscha=CC.Phonons.Phonons( dyn_sscha_name,8)
# Generate the ensemble
supercell=dyn_sscha.GetSupercell()
ens = sscha.Ensemble.Ensemble(dyn_sscha, T, supercell)
ens.generate(n_random)
# Compute energy and forces for the ensemble elements
ens.get_energy_forces(ff_calculator , compute_stress = False)
# save population
ens.save(ens_savedir, population)
```

Once the population ensamble has been generated, the minimization has to be performed. The input file to perform the minimization as stand-alone code, without python scripting, is


```python
cat > input << EOF
&inputscha
    n_random =                                ! number of elements (integer)
    data_dir = "ens_pop#(population-number)"  ! population directory path (character)                     
    population = 1                            ! population number (integer)
    fildyn_prefix = "PbTe.dyn"                ! dyn mat that generated the population (character)
    nqirr = 8
    supercell_size =  4 4 4
    Tg = 300.0 
    T  = 300.0
    preconditioning = .true.
    root_representation = "normal"
    meaningful_factor = 1e-12
    gradi_op = "gc" 
    n_random_eff =                            ! Kong-Liu parameter minimum threshold
    minim_struc  =  .false. 
    print_stress = .false.
    eq_energy = 0.0
    lambda_a = 0.1
    lambda_w = 0.1
    max_ka = 10000
/
&utils
    save_freq_filename = "Frequencies.dat"
&end
EOF
```
and it is launched with:
```
sscha -i input --save-data data_saved > output
```

At the end of the minimizations, performed with several populations, we have:

The dynamical matrices that generated the last population #popnum

PbTe.dyn#q

The last population, stored in the folder

ens_pop#popnum

The SSCHA dynamical matrices obtained with the last minimization

PbTe.SSCHA.dyn#q

We have all the ingredients to perform the next step and calculate the free energy Hessian (including or not the 4th order SSCHA FCs terms).

At the end of the calculation we will have two main results

1. The free energy Hessian dynamical matrices:

   PbTe.Hessian.dyn#q


2. The third order SSCHA FCs (conveniently symmetrized), 
   crucial for the spectral analysis that we will perform:
   
   d3_realspace_sym.npy

This is the script to perform such a calculation, where the parameters to be inserted are: the ensamble data directory ("DATA_DIR",string), and the number of elements in the ensamble ("N_RANDOM",integer), whereas ("INCLUDE_V4", logical) can be left as False (we will use functionalities that use only the third order FCs, i.e. the ''bubble''):


```python
from __future__ import print_function
from __future__ import division

# Import the modules to read the dynamical matrix
import cellconstructor as CC
import cellconstructor.Phonons

# Import the SCHA modules
import sscha, sscha.Ensemble


# Here the input information
DATA_DIR = # path to the directory ens_pop#lastpop where the last population is stored
N_RANDOM = # number elements in the ensamble
DYN_PREFIX =  'PbTe.dyn'  # dyn mat that generated the last population
FINAL_DYN =   'PbTe.SSCHA.dyn'    # SSCHA dyn mat obtained with the last minimization 
SAVE_PREFIX = 'PbTe.Hessian.dyn'  # Free energy Hessian dynamical matrices
NQIRR = 8
Tg = 300
T =  300
POPULATION = # number of last population
INCLUDE_V4 = False # True to include the 4th-order SSCHA FC term to calculate the Hessian 

INFO = """
In this example we compute the free energy hessian.

The ensemble has been generated with the dynamical matrix at:
{}

And to compute the hessian we will use reweighting at:
{}

The original temperature was {} K, we use reweighting to {} K.
The ensemble, population {}, is located at: {}
The number of configuration is {}.
Do we include the v4 in the calculation? {}

The free energy Hessian will be saved in: {}

The (symmetrized) 3rd order FCs in d3_realspace_sym.npy

""".format(DYN_PREFIX, FINAL_DYN, Tg, T, POPULATION, DATA_DIR,
           N_RANDOM, INCLUDE_V4, SAVE_PREFIX)


print(INFO)
print()
print(" ======== RUNNING ======== ")
print()

print("Loading the original dynamical matrix...")
dyn = CC.Phonons.Phonons(DYN_PREFIX, NQIRR)
print("Loading the current dynamical matrix...")
final_dyn = CC.Phonons.Phonons(FINAL_DYN, NQIRR)

print("Loading the ensemble...")
ens = sscha.Ensemble.Ensemble(dyn, Tg, dyn.GetSupercell())
ens.load(DATA_DIR, POPULATION, N_RANDOM)
# If the ensemble was saved in binary format, load it with
# ens.load_bin(DATA_DIR, POPULATION)

print("Updating the importance sampling...")
ens.update_weights(final_dyn, T)

print("Computing the free energy hessian...")
# Set get_full_hessian to false to have only the odd correction
# Usefull if you want to study the convergence with the number of configuration
dyn_hessian = ens.get_free_energy_hessian(include_v4 = INCLUDE_V4,
                                          get_full_hessian = True,
                                          verbose = True)

print("Saving the hessian to {}...".format(SAVE_PREFIX))
dyn_hessian.save_qe(SAVE_PREFIX)
print("Done.")
```

The matrices PbTe.Hessian.dyn#q are a generalization of the standard harmonic dynamical matrices that include quantum and thermal effects on a static level, Eq.(61). As long as only the "bubble" term is included, Eq.(63), the SSCHA code is able to employ the Fourier interpolation technique to obtain the free energy Hessian on a generic __q__ , integrating on a generic **k** grid (Sec. IV-B, Eq. (66)). In order to perform such a calculations, we have to follow several preliminary steps (see Appendix E.2):

1. We "center" the 3rd order FCs (a step necessary to perform the Fourier interpolation)

2. We impose the acoustic sum rule (ASR)

With this script we load the 3rd oder FCs in d3_realspace_sym.npy, we center it and impose the ASR, and finally
we print it in the FC3 file.


```python
import cellconstructor as CC
import cellconstructor.ForceTensor
import cellconstructor.Structure
import cellconstructor.Spectral

import numpy as np


# Initialize the tensor3 object
# We need 2nd FCs of the used grid to configure the supercell.
# For example, we can use PbTe.dyn#q, or PbTe.SSCHA.dyn#q, or PbTe.Hessian.dyn#q
dyn = CC.Phonons.Phonons("PbTe.SSCHA.dyn",3) 
supercell = dyn.GetSupercell()
tensor3 = CC.ForceTensor.Tensor3(dyn.structure,
                                dyn.structure.generate_supercell(supercell),
                                supercell)

# Assign the tensor3 values
d3 = np.load("d3_realspace_sym.npy")*2.0
tensor3.SetupFromTensor(d3)

# Center and apply ASR
tensor3.Center()
tensor3.Apply_ASR()

# Print it
tensor3.WriteOnFile(fname="FC3",file_format='D3Q')

```

With PbTe.SSCHA.dyn#q and FC3 we have all the ingredients to perform the interpolated Hessian calculation.
As first thing, we can do a double check and verify that the centered (with imposed ASR) FC3 gives the same Hessian than the original calculation, for the same __q__ points and the integration **k**-grid commensurate with the supercell. We can use this input file:


```python
import cellconstructor as CC
import cellconstructor.ForceTensor
import cellconstructor.Structure
import cellconstructor.Spectral

dyn = CC.Phonons.Phonons("PbTe.SSCHA.dyn",8)

supercell = dyn.GetSupercell()
tensor3 = CC.ForceTensor.Tensor3(dyn.structure,
                                dyn.structure.generate_supercell(supercell),
                                supercell)

tensor3.SetupFromFile(fname="FC3",file_format='D3Q')


# integration grid
k_grid=[4,4,4]    

# q points in 2pi/Angstrom
list_of_q_points=[ [  0.0000000,  0.0000000,  0.0000000 ],
                   [ -0.0386763,  0.0386763, -0.0386763 ],
                   [  0.0773527, -0.0773527,  0.0773527 ],
                   [  0.0000000,  0.0773527,  0.0000000 ],
                   [  0.1160290, -0.0386763,  0.1160290 ],
                   [  0.0773527,  0.0000000,  0.0773527 ],
                   [  0.0000000, -0.1547054,  0.0000000 ],
                   [ -0.0773527, -0.1547054,  0.0000000 ]   ]


CC.Spectral.get_static_correction_along_path(dyn=dyn, 
                                             tensor3=tensor, 
                                             k_grid=k_grid, 
                                             q_path=list_of_q_points, 
                                             filename_st="v2_v2+d3static_freq.dat",
                                             T =300.0,
                                             print_dyn = False) # set true to print the Hessian dynamical matrices
                                                                # for each q point


```

This calculation (like all the calculations described below) 
can be performed in parallel on NPROC processors, using MPI. Called input.py the input file
we launch it with
```
mpirun -np NPROC python input.py > output
```

In the file *"v2_v2+d3static_freq.dat"* we have 8 rows (one for each q point).
The first column is the length of the path followed along these 8 points in $$2\pi$$/Å units.
After we have the SSCHA frequencies and the Hessian frequencies. They must coincide with the
frequencies that we have already calculated in PbTe.SSCHA.dyn#q and PbTe.Hessian.dyn#q. 

Up to now, we have not really used the Fourier interpolation ( the **q** and the __k__ grid points are commensurate with the supercell calculation) and, as a matter of fact, the centering+ASR imposition was not necessary. However, we can now compute the Hessian dynamical matrices and frequencies along a generic path, integrating on an arbitrary finer grid. We consider the path $$X-\Gamma-X$$ in the *"XGX_path.dat"* file
and we integrate on a $$20\times 20\times 20$$ grid (the path in *"XGX_path.dat"* is made of 1000 points. To speed up the calculations a path with less points can be used).


```python
import cellconstructor as CC
import cellconstructor.ForceTensor
import cellconstructor.Structure
import cellconstructor.Spectral

dyn = CC.Phonons.Phonons("PbTe.SSCHA.dyn",8)

supercell = dyn.GetSupercell()
tensor3 = CC.ForceTensor.Tensor3(dyn.structure,
                                dyn.structure.generate_supercell(supercell),
                                supercell)

tensor3.SetupFromFile(fname="FC3",file_format='D3Q')


# integration grid
k_grid=[20,20,20]    


CC.Spectral.get_static_correction_along_path(dyn=dyn, 
                                             tensor3=tensor, 
                                             k_grid=k_grid, 
                                             q_path_file="XGX_path.dat",
                                             filename_st="v2_v2+d3static_freq.dat",
                                             T =300.0,
                                             print_dyn = False) # set true to print the Hessian dynamical matrices
                                                                # for each q point
```

The result can be plotted to display the Hessian frequency dispersion along the path.

![](Hessian.png)

Notice in this plot the LO-TO splitting, which is calculated by the SSHCA code using the effective charges and the electronic permittivity tensor printed in the center zone dynamical matrix file (see Appendix E.3).

The static calculation is essentially devoted to identify the presence of structural instabilities (in this case, as we can see, we do not have instabilities along $$\Gamma-X$$).
To properly compute the phonon spectrum a dynamical calculation has to be performed. As a first thing, let us do a double-check calculation. Let us compute the spectral function $$\sigma(\Omega)$$ in $$X$$ and $$\Gamma$$ obtained in the static approximation Eqs.(72)-(74). The resulting spectral function must be composed of peaks around the Hessian frequencies.


```python
import cellconstructor as CC
import cellconstructor.ForceTensor
import cellconstructor.Structure
import cellconstructor.Spectral

dyn = CC.Phonons.Phonons("PbTe.SSCHA.dyn",8)

supercell = dyn.GetSupercell()
tensor3 = CC.ForceTensor.Tensor3(dyn.structure,
                                dyn.structure.generate_supercell(supercell),
                                supercell)

tensor3.SetupFromFile(fname="FC3",file_format='D3Q')


# integration grid
k_grid=[20,20,20]    

# X and G in 2pi/Angstrom
points=[[-0.1525326,  0.0,  0.0],
        [0.0       ,  0.0,  0.0]      ]

CC.Spectral.get_full_dynamic_correction_along_path(dyn=dyn, 
                                                   tensor3=tensor3, 
                                                   k_grid=k_grid,  
                                                   e1=100, de=0.1, e0=0,     # energy grid
                                                   sm1=1.0, sm0=1.0,  nsm=1, # smearing values
                                                   T=300,
                                                   q_path=points,                                                                                     
                                                   static_limit = True, #static approximation
                                                   notransl = True,  # projects out the acoustic zone center modes
                                                   filename_sp='static_spectral_func')
```

The input values e1, e0, de define the energy grid where the spectral function will be computed: initial, final and spacing value of the energy grid in cm$${}^{-1}$$, respectively. The input values sm0, sm1, nsm define the values used for $$\delta_{\text{se}}$$ of Eq. (75), the smearing used to compute the self-energy: initial, final and number of intermediate values between them, in cm$${}^{-1}$$ (however, as long as we consider the static approximation, the value of $$\delta_{\text{se}}$$ is immaterial). 

We decided to project out the part of the spectral function due to the pure translation modes (since they convey
a trivial information). The result is in the *"static_spectral_func_1.0.dat"* file. 
The first column is the distance of the followed reciprocal space path in $$2\pi$$/Å (an information that we are not going to use now). For each point we have the values $$\Omega$$ of the used energy grid (second column) and the spectral function $$\sigma(\Omega)$$ ( third column).

Here we plot $$\sigma(\Omega)$$ for the two points with black line. The colored vertical lines are the Hessian frequency values previously calculated in $$X$$ and $$\Gamma$$. Note that the higher peaks correspond to double degenerate frequencies.

![](specX.png)

![](specG.png)

Now we perform a full spectral calculation in $$\Gamma$$, Eq. (76), with the input


```python
import cellconstructor as CC
import cellconstructor.ForceTensor
import cellconstructor.Structure
import cellconstructor.Spectral

dyn = CC.Phonons.Phonons("PbTe.SSCHA.dyn",8)

supercell = dyn.GetSupercell()
tensor3 = CC.ForceTensor.Tensor3(dyn.structure,
                                dyn.structure.generate_supercell(supercell),
                                supercell)

tensor3.SetupFromFile(fname="FC3",file_format='D3Q')


# integration grid
k_grid=[20,20,20]    

# q point
G=[0.0,0.0,0.0]


CC.Spectral.get_full_dynamic_correction_along_path(dyn=dyn, 
                                           tensor3=tensor3, 
                                           k_grid=k_grid,  
                                           e1=145, de=0.1, e0=0,
                                           sm1=1, sm0=1,nsm=1,
                                           T=300,
                                           q_path=G,                                           
                                           notransl = True,
                                           filename_sp='full_spectral_func')
```

and plot the result (second and third column of the *"full_spectral_func_1.0.dat"* file). 

![](specG_full.png)

In this calculation we have calculated the full self-energy. We can also employ the no-mode-mixing approximation, 
Eqs. (78)-(80), discarding the off-diagonal elements of the  self-energy in the SSCHA mode basis set with the input 


```python
from __future__ import print_function
import cellconstructor as CC
import cellconstructor.ForceTensor
import cellconstructor.Structure
import cellconstructor.Spectral

dyn = CC.Phonons.Phonons("PbTe.SSCHA.dyn",8)

supercell = dyn.GetSupercell()
tensor3 = CC.ForceTensor.Tensor3(dyn.structure,
                                dyn.structure.generate_supercell(supercell),
                                supercell)

tensor3.SetupFromFile(fname="FC3",file_format='D3Q')


# integration grid
k_grid=[20,20,20]    

# 
G=[0.0,0.0,0.0]

CC.Spectral.get_diag_dynamic_correction_along_path(dyn=dyn, 
                                                   tensor3=tensor3,  
                                                   k_grid=k_grid, 
                                                   q_path=G,
                                                   T = 300.0, 
                                                   e1=145, de=0.1, e0=0,
                                                   sm1=1.0, nsm=1, sm0=1.0,
                                                   filename_sp = 'nomm_spectral_func')

```

The result is printed in the *"nomm_spectral_func_1.0.dat"* file. As before, the first column gives the distance along the path in reciprocal space (here inessential, we are doing a single point calculation), the second column the energy grid values, the third column the spectral function, and subsequently a column for each single mode contribution to the spectral function. Here we are interested in the modes 4, 5, 6 (the first three modes are the acoustic ones).

If we consider the sum of the spectral function of these three modes we essentially obtain the same result obtined with the full self-energy calculation. Therefore, we can safely use the no-mode-mixing approximation.
We can perform a spectral analysis mode by mode: 

![](specG_modes.png)

The modes 4 and 5 are degenerate (their sum gives the corresponding part of the spectral function showed in the previous figure).  As it is evident, while the mode 6 seems to have a Lorentzian character (with precise
center and linewidth), for the modes 4 and 5 there is a clear no-Lorentzian character. This is evident if we plot the Lorentzian spectral functions for these modes, using for example the "one-shot" approach, Eqs. (81), (84), (85), that we find in the file *"nomm_spectral_func_lorentz_one_shot_1.0.dat"*.
As we can see, the spectral function for the mode 6 is well described in the Lorentzian approximation, whereas modes 4 and 5 have a strong non-Lorentzian character.

![](specG_lorentzian.png)

To have a more complete picture we can compute the spectral function (in no-mode-mixing approximation)
along the $$X-\Gamma-X$$ path.


```python
from __future__ import print_function
import cellconstructor as CC
import cellconstructor.ForceTensor
import cellconstructor.Structure
import cellconstructor.Spectral

dyn = CC.Phonons.Phonons("PbTe.SSCHA.dyn",8)

supercell = dyn.GetSupercell()
tensor3 = CC.ForceTensor.Tensor3(dyn.structure,
                                dyn.structure.generate_supercell(supercell),
                                supercell)

tensor3.SetupFromFile(fname="FC3",file_format='D3Q')


# integration grid
k_grid=[20,20,20]    

CC.Spectral.get_diag_dynamic_correction_along_path(dyn=dyn, 
                                                   tensor3=tensor3,  
                                                   k_grid=k_grid, 
                                                   q_path_file="XGX_path.dat"
                                                   T = 300.0, 
                                                   e1=145, de=0.1, e0=0,
                                                   sm1=1.0, nsm=1, sm0=1.0,
                                                   filename_sp = 'nomm_spectral_func')
```

We plot the first three columns of *"nomm_spectral_func_1.0.dat"* with a colormap plot (the spectral function plotted as a color function)

![](spec_path.png)

We clearly see a satellite in $$\Gamma$$, which well reproduces what is seen in experiments ([Ribeiro et al.,Phys. Rev. B 97, 014306](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.014306)).

<a name="Tin-Telluride-with-force-fields"></a>
# 4. Tin Telluride with force fields

In this tutorial we are going to study the thermoelectric transition in SnTe.
To speedup the calculations, we will use a force-field that can mimic the physics of ferroelectric transitions in FCC lattices.

We will replicate the calculations performed in the paper by [Bianco et. al. Phys. Rev. B 96, 014111](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.014111).

The force field can be downloaded and installed from [here](https://github.com/SSCHAcode/F3ToyModel). 

## Initialization

As always, we need to initialize the working space.
This time we will initialize first the force field.
This force field needs the harmonic dynamical matrix to be initialized, and the higher order parameters.
We will initialize it in order to reproduce the results in [Bianco et. al. Phys. Rev. B 96, 014111](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.014111). The dynamical matrices needed for the harmonic part of the force field can be found in python-sscha/Tutorials/SnTe_ToyModel.

```python
# Import the cellconstructor stuff
import cellconstructor as CC
import cellconstructor.Phonons

# Import the modules of the force field
import fforces as ff
import fforces.Calculator

# Import the modules to run the sscha
import sscha, sscha.Ensemble, sscha.SchaMinimizer
import sscha.Relax, sscha.Utilities

# Load the dynamical matrix for the force field
ff_dyn = CC.Phonons.Phonons("ffield_dynq", 3)

# Setup the forcefield with the correct parameters
ff_calculator = ff.Calculator.ToyModelCalculator(ff_dyn)
ff_calculator.type_cal = "pbtex"
ff_calculator.p3 = 0.036475
ff_calculator.p4 = -0.022
ff_calculator.p4x = -0.014
```

We initialized a force field. For a detailed explanations of the parameters, refer to the [Bianco et. al. Phys. Rev. B 96, 014111](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.014111) paper.
Now ff_calculator behaves like any [ASE](https://wiki.fysik.dtu.dk/ase/) calculator, and can be used to compute forces and energies for our SSCHA minimization. Note: this force field is not able to compute stress, as it is defined only at fixed volume, so we cannot use it for a variable cell relaxation.

Now it is the time to initialize the SSCHA. 
We can start from the harmonic dynamical matrix we got for the force field. Remember, SSCHA dynamical matrices must be positive definite.
Since we are studying a system that has a spontaneous symmetry breaking at low temperature, the harmonic dynamical matrices will have imaginary phonons.
We must enforce phonons to be positive definite to start a SSCHA minimization.


```python
# Initialization of the SSCHA matrix
dyn_sscha = ff_dyn.Copy()
dyn_sscha.ForcePositiveDefinite()

# Apply also the ASR and the symmetry group
dyn_sscha.Symmetrize()
```

We must now prepare the sscha ensemble for the minimization. 
We will start with a $$T= 0K$$ simulation.


```python
ensemble = sscha.Ensemble.Ensemble(dyn_sscha, T0 = 0, supercell = dyn_sscha.GetSupercell())
```

We can now proceed with the sscha minimization.
Since we start from auxiliary dynamical matrices that are very far from the correct result, it is convenient to use a safer minimization scheme.
We will use the fourth root minimization, in which, instead of optimizing the auxiliary dynamical matrices themself, we will optimize their fourth root.
Then the dynamical matrix $$\Phi$$ will be obtained as:

$$
\Phi = \left(\sqrt[4]{\Phi}\right)^4.
$$

This constrains $$\Phi$$ to be positive definite during the minimization. Moreover, this minimization is more stable than the standard one. If you want further details, please look at [Monacelli et. al. Phys. Rev. B 98, 024106](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.98.024106).

We also change the Kong-Liu effective sample size threshold. This is a value that decrease during the minimization, as the parameters gets far away from the starting point. It estimates how the original ensemble is good to describe the new parameters. If this value goes below a given threshold, the minimization is stopped, and a new ensemble is extracted.


```python
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

# Lets setup the minimization on the fourth root
minim.root_representation = "root4" # Other possibilities are 'normal' and 'sqrt'

# To work correctly with the root4, we must deactivate the preconditioning on the dynamical matrix 
minim.precond_dyn = False

# Now we setup the minimization parameters
# Since we are quite far from the correct solution, we will use a small optimization step
minim.min_step_dyn = 1 # If the minimization ends with few steps (less than 10), decrease it, if it takes too much, increase it

# We decrease the Kong-Liu effective sample size below which the population is stopped
minim.kong_liu_ratio = 0.2 # Default 0.5 
```

We can setup the automatic relaxation, to avoid the need to restart the frequencies at each iteration.
We will also setup a custom function to save the frequencies at each iteration, to see how they evolves.
This is very usefull to understand if the algorithm is converged or not.

We will use ensembles of 1000 configurations for each population, and a maximum of 20 populations.


```python
relax = sscha.Relax.SSCHA(minim, 
                          ase_calculator = ff_calculator,
                          N_configs = 1000, 
                          max_pop = 20)

# Setup the custom function to print the frequencies at each step of the minimization
io_func = sscha.Utilities.IOInfo()
io_func.SetupSaving("frequencies.dat") # The file that will contain the frequencies is frequencies.dat

# Now tell relax to call the function to save the frequencies after each iteration
# CFP stands for Custom Function Post (Post = after the minimization step)
relax.setup_custom_functions(custom_function_post = io_func.CFP_SaveFrequencies)
```

We are ready to start the SSCHA minimization. This may take few minutes, depending on how powerfull is your PC. 
If you do not want to run this on your machine, skip it and pass to the following .


```python
relax.relax()

# Save the final dynamical matrix
relax.minim.dyn.save_qe("final_sscha_T0_")
```

## Plotting the results

We can plot the evolution of the frequencies, as well as the free energy and the gradient to see if the minimization ended correctly.



```python
# Import Matplotlib to plot
import numpy as np
import matplotlib.pyplot as plt

# Setup the interactive plotting mode
plt.ion()

# Lets plot the Free energy, gradient and the Kong-Liu effective sample size
relax.minim.plot_results()
```


![png](output_13_0.png)



![png](output_13_1.png)



![png](output_13_2.png)



![png](output_13_3.png)


As you can see, the free energy always decreases, and reaches a plateau. The gradient at the beginning increases, then it goes to zero very fast close to convergence.
To have an idea on how good is the ensemble, check the Kong-Liu effective sample size (the last plot). We can see that we went two times below the convergence threshold of 0.2 (200 out of 1000 configuration), this means that the we required three populations.
We have at the end 500 good configurations out of the original 1000. This means that the ensemble is still at its 50 % of efficiency.

Let's have a look on how the frequencies evolve now, by loading the frequencies.dat file that we created.


```python
frequencies = np.loadtxt("frequencies.dat")
N_steps, N_modes = frequencies.shape

#For each frequency, we plot it [we convert from Ry to cm-1]
plt.figure(dpi = 120)
for i_mode in range(N_modes):
    plt.plot(frequencies[:, i_mode] * CC.Units.RY_TO_CM)
plt.xlabel("Steps")
plt.ylabel("Frequencies [cm-1]")
plt.title("Evolution of the frequencies")
plt.tight_layout()
```


![png](output_15_0.png)


While at the beginning the frequencies were changing a lot, in the last population they changed more smoothly. 
This is a good sign of convergence. One frequency has a very small value. 
The cusps in the frequencies is the point in which we changed the ensemble. In the minimization we provided, this happened once slightly above step 100, in correspondence with the change of the ensemble into the last one.
We reached the final results with only 3 ensembles (3000 configurations energy/forces calculations).

# The instability

From the frequencies, we can see that we have one SSCHA frequency that is very low, below 10 cm$$^{-1}$$.
This is probably a sign of instability.
We remark: the sscha auxiliary frequencies are not the real frequencies observed in experiments, but rather are linked to the average displacements of atoms along that mode. In particular the average displacements of atoms can be computed from the SSCHA auxiliary frequencies as (includes both thermal and quantum fluctuations):

$$
\sigma_\mu = \sqrt{ \frac{1 + 2n_\mu}{2\omega_\mu}}
$$

where $$n_\mu$$ is the boson occupation number and $$\omega_\mu$$ is the SSCHA auxiliary frequency.

It is clear, that $$\omega_\mu$$ will always be positive, even if we have an instability. Since we have a very small mode in the SSCHA frequencies, it means that associated to that mode we have huge fluctuations. This can indicate an instability.
However, to test this we need to compute the free energy curvature along this mode. This can be obtained in one shot thanks to the theory developed in [Bianco et. al. Phys. Rev. B 96, 014111](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.014111). 

First of all, we generate a new ensemble with more configurations. To compute the hessian we will use an ensemble of 10000 configurations.


```python
# We reload the final result (no need to rerun the sscha minimization)
dyn_sscha_final = CC.Phonons.Phonons("final_sscha_T0_", 3)

# We reset the ensemble
ensemble = sscha.Ensemble.Ensemble(dyn_sscha_final, T0 = 0, supercell = dyn_sscha_final.GetSupercell())

# We need a bigger ensemble to properly compute the hessian
# Here we will use 10000 configurations
ensemble.generate(10000)

# We now compute forces and energies using the force field calculator
ensemble.get_energy_forces(ff_calculator, compute_stress = False)
```

Now we can compute the free energy hessian. 
We can choose if we neglect or not in the calculation the four phonon scattering process. Four phonon scattering processes require a huge memory allocation for big systems, that scales as $$(3 \cdot N)^4$$ with $$N$$ the number of atoms in the supercell. Moreover, it may require also more configurations to converge.

In almost all the systems we studied up to now, we found this four phonon scattering at high order to be negligible.
We remark, that the SSCHA minimization already includes four phonon scattering at the lowest order perturbation theory, thus neglecting this term only affects combinations of one or more four phonon scattering with two three phonon scatterings (high order diagrams).
For more details, see [Bianco et. al. Phys. Rev. B 96, 014111](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.014111). 



```python
dyn_hessian = ensemble.get_free_energy_hessian(include_v4 = False) # We neglect high-order four phonon scattering

# We can save it
dyn_hessian.save_qe("hessian")
```

We can print the eigenmodes of the free energy hessian to check if there is an instability:


```python
w_hessian, pols_hessian = dyn_hessian.DiagonalizeSupercell()

# Print all the frequency converting them into cm-1 (They are in Ry)
print("\n".join(["{:16.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in w_hessian]))
```

            -18.4862 cm-1
            -18.4862 cm-1
            -18.4862 cm-1
              0.0000 cm-1
              0.0000 cm-1
              0.0000 cm-1
             23.3091 cm-1
             23.3091 cm-1
             23.3091 cm-1
             23.3091 cm-1
             23.3091 cm-1
             23.3091 cm-1
             25.7898 cm-1
             25.7898 cm-1
             25.7898 cm-1
             48.9018 cm-1
             48.9018 cm-1
             48.9018 cm-1
             48.9018 cm-1
             48.9018 cm-1
             48.9018 cm-1
             48.9018 cm-1
             48.9018 cm-1
             52.5769 cm-1
             52.5769 cm-1
             52.5769 cm-1
             52.5769 cm-1
             52.5769 cm-1
             52.5769 cm-1
             64.6604 cm-1
             64.6604 cm-1
             64.6604 cm-1
             67.9293 cm-1
             67.9293 cm-1
             67.9293 cm-1
             67.9293 cm-1
             67.9293 cm-1
             67.9293 cm-1
             67.9293 cm-1
             67.9293 cm-1
             80.5268 cm-1
             80.5268 cm-1
             80.5268 cm-1
             80.5268 cm-1
            102.3613 cm-1
            102.3613 cm-1
            102.3613 cm-1
            102.3613 cm-1


Yes we have  imaginary phonons! We found an instability! You can check what happens if you include the fourth order.

# The phase transition

Up to now we studied the system at $$T = 0K$$ and we found that there is an instability. However, we can repeat the minimization at many temperatures, and track the phonon frequency to see which is the temperature at which the system becomes stable.

We can exploit the fact that our *python-sscha* package is a python library, and write a small script to automatize the calculation.

We will simulate the temperatures up to room temperature (300 K) with steps of 50 K.
Note, this will perform all the steps above 6 times, so it may take some minutes, depending on the PC (on a i3 from 2015, with one core, it took 2 hours).


```python
# Define the temperatures, from 50 to 300 K, 6 temperatures
temperatures = np.linspace(50, 300, 6)

lowest_hessian_mode = []
lowest_sscha_mode = []

# Perform a simulation at each temperature
t_old = 0
for T in temperatures:
    # Load the starting dynamical matrix
    dyn = CC.Phonons.Phonons("final_sscha_T{}_".format(int(t_old)), 3)
    
    # Prepare the ensemble
    ensemble = sscha.Ensemble.Ensemble(dyn, T0 = T, supercell = dyn.GetSupercell())
    
    # Prepare the minimizer 
    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
    minim.min_step_dyn = 0.002
    minim.kong_liu_ratio = 0.5
    #minim.root_representation = "root4"
    #minim.precond_dyn = False
    
    # Prepare the relaxer (through many population)
    relax = sscha.Relax.SSCHA(minim, ase_calculator = ff_calculator, N_configs=1000, max_pop=5)
    
    # Relax
    relax.relax()
    
    # Save the dynamical matrix
    relax.minim.dyn.save_qe("final_sscha_T{}_".format(int(T)))
    
    # Recompute the ensemble for the hessian calculation
    ensemble = sscha.Ensemble.Ensemble(relax.minim.dyn, T0 = T, supercell = dyn.GetSupercell())
    ensemble.generate(5000)
    ensemble.get_energy_forces(ff_calculator, compute_stress = False)
    
    # Get the free energy hessian
    dyn_hessian = ensemble.get_free_energy_hessian(include_v4 = False)
    dyn_hessian.save_qe("hessian_T{}_".format(int(T)))

    # Get the lowest frequencies for the sscha and the free energy hessian
    w_sscha, pols_sscha = relax.minim.dyn.DiagonalizeSupercell()
    # Get the structure in the supercell
    superstructure = relax.minim.dyn.structure.generate_supercell(relax.minim.dyn.GetSupercell()) #
    
    # Discard the acoustic modes
    acoustic_modes = CC.Methods.get_translations(pols_sscha, superstructure.get_masses_array())
    w_sscha = w_sscha[~acoustic_modes]
    
    lowest_sscha_mode.append(np.min(w_sscha) * CC.Units.RY_TO_CM) # Convert from Ry to cm-1
    
    w_hessian, pols_hessian = dyn_hessian.DiagonalizeSupercell()
    # Discard the acoustic modes
    acoustic_modes = CC.Methods.get_translations(pols_hessian, superstructure.get_masses_array())
    w_hessian = w_hessian[~acoustic_modes]
    lowest_hessian_mode.append(np.min(w_hessian) * CC.Units.RY_TO_CM) # Convert from Ry to cm-1
    
    t_old = T

# We prepare now the file to save the results
freq_data = np.zeros( (len(temperatures), 3))
freq_data[:, 0] = temperatures
freq_data[:, 1] = lowest_sscha_mode
freq_data[:, 2] = lowest_hessian_mode

# Save results on file
np.savetxt("hessian_vs_temperature.dat", freq_data, header = "T [K]; SSCHA mode [cm-1]; Free energy hessian [cm-1]")
```

We can now load and plot the results.



```python
hessian_data = np.loadtxt("hessian_vs_temperature.dat")
plt.figure(dpi = 120)
plt.plot(hessian_data[:,0], hessian_data[:,1], label = "Min SCHA freq", marker = ">")
plt.plot(hessian_data[:,0], hessian_data[:,2], label = "Free energy curvature", marker = "o")
plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
plt.xlabel("Temperature [K]")
plt.ylabel("Frequency [cm-1]")
plt.legend()
plt.tight_layout()
```


![png](output_26_0.png)


From the previous plot we can easily see that the phase transition occurs at about 170 K.
It is worth noting, as pointed out, that the SSCHA frequency is always positive definite, and no divergency is present in correspondance of the transition.
Moreover, the Free energy curvature is more noisy than the SSCHA one. Mainly because the error is bigger on small frequencies for the fact that
$$
\omega \sim \sqrt{\Phi}.
$$
But this is also due to the computation of the free energy curvature itself, which requires the third order force constant tensor that requires more configurations to converge.

Be aware, if you study phase transition in charge density waves, like [NbS2](https://pubs.acs.org/doi/abs/10.1021/acs.nanolett.9b00504) or [TiSe2](https://arxiv.org/abs/1910.12709),
or thermoelectric materials like [SnSe](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.075901) or [SnS](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.214307) usually the transition temperature depends strongly on the supercell size.

For the Landau theory of phase transition, since the SSCHA is a mean-field approach, we expect that around the transition the critical exponent of the temperature goes as
$$
\omega \sim T^\frac 1 2.
$$

Thus it is usually better to plot the temperature versus the square of the frequency:


```python
hessian_data = np.loadtxt("hessian_vs_temperature.dat")
plt.figure(dpi = 120)
plt.plot(hessian_data[:,0], np.sign(hessian_data[:,2]) * hessian_data[:,2]**2, label = "Free energy curvature", marker = "o")
plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
plt.xlabel("Temperature [K]")
plt.ylabel("$\omega^2$ [cm-2]")
plt.legend()
plt.tight_layout()
```


![png](output_28_0.png)


As seen, the linear interpolation between 100 and 250 K is much better in this plot, allowing for a more precise determination of the transition temperature.

<a name="Sulfur-hydride"></a>
# 5. Sulfur hydride

Sulfur hydride is one of the most interesting system of physics, it has one of the highest superconducting $$T_c$$'s ever measured (203 K).
Thanks to the very light mass of hydrogen atoms, its physical properties are dominated by anharmonicity and quantum fluctuations, making it one of the best systems to be tackled with SSCHA.

In this simple tutorial we will learn how to run an automatic quantum relaxation of the structure. 
We will follow the lines of these works: [Errea et al. Nature volume 532, pages81–84(2016)](https://www.nature.com/articles/nature17175) and [Bianco et al Phys. Rev. B 97, 214101](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.214101).

We provide the H3S.scf file, that contains the lattice and the atomic position of the high symmetry phase of H3S. 
The harmonic dynamical matrix are given as computed by quantum espresso, however, they can be recomputed live by running the get_phonons.py script (needs ASE and quantum-espresso phonon package). All these files are provided in python-sscha/Tutorials/H3S.

You will see that also this case has an instability at gamma, with a phonon mode that is imaginary.

# How to automatize the relaxation

In the first [PbTe tutorial](http://sscha.eu/Tutorials/Tutorial_PbTe/) we saw how to prepare the SSCHA calculation from scratch and how to manipulate ensembles and the
minimization.
However, all the procedure could be a bit clumsy. In this section we will show how to automatize everything.
We will complete the PbTe simulation to get some useful results.

The automatic relaxation is performed through the Relax module, that implements two classes, one for the static cell relaxation and the other for the variable cell relaxation


```python
# Import ASE to setup the automatic calculator
import ase
from ase.calculators.espresso import Espresso

# Lets import cellconstructor
import cellconstructor as CC
import cellconstructor.Phonons

# Import the sscha
import sscha, sscha.Ensemble, sscha.Relax, sscha.SchaMinimizer

from numpy import *
import numpy as np
import matplotlib.pyplot as plt

plt.ion() #Setup interactive plot mode
```

To make everything clear, we will start again from the harmonic dynamical matrix for our structure: matdynX.    
We need to impose the sum rule and to force it to be positive definite before starting a SSCHA calculation


```python
pseudo = {"H":"H.pbe-rrkjus_psl.1.0.0.UPF",
         "S" : "S.pbe-nl-rrkjus_psl.1.0.0.UPF"}

input_data = {"ecutwfc" : 35,
              "ecutrho" : 350,
              "occupations" : "smearing",
              "input_dft" : "blyp",
              "mixing_beta" : 0.2,
              "conv_thr" : 1e-9,
              "degauss" : 0.02,
              "smearing" : "mp",
              "pseudo_dir" : ".",
             "tprnfor" : True}


espresso_calc = Espresso(pseudopotentials =pseudo,
                        input_data = input_data,
                        kspacing = 0.06,
                        koffset = (1,1,1))
```


```python
# Let us load the starting dynamical matrix
dyn = CC.Phonons.Phonons("matdyn", nqirr = 3)

# Apply the sum rule and delete the imaginary modes
dyn.Symmetrize()
dyn.ForcePositiveDefinite()

# Generate the ensemble and the minimizer objects
ensemble = sscha.Ensemble.Ensemble(dyn, 1000, supercell = dyn.GetSupercell())
minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

# We setup all the minimization parameters
minimizer.min_step_dyn = 0.001
minimizer.kong_liu_ratio = 0.5
```

You can automatize also the cluster calculation by setting up a cluster object 
like we did in the previous tutorials.
If you keep it as None (as done in the following cell) the calculation will be runned locally.
Remember, if you use the cluster, you need to copy first the pseudopotential in the working directory of the cluster.


```python
# Here we prepare a cluster
# Here we configure the cluster object MARCONI
import sscha.Cluster
my_hpc = sscha.Cluster.Cluster(pwd = None)

# We setup the connection info
my_hpc.hostname = "ekhi" # The command to connect via ssh to the cluster
#my_hpc.account_name = "IscrB_COMRED" # The name of the project for the computation
my_hpc.workdir = "/scratch/lorenzo/my_h3s_dir" # the directory in which the calculations are performed

# Now we need to setup the espresso
# First we must tell the cluster where to find him:
my_hpc.binary = "pw.x -npool NPOOL -i  PREFIX.pwi > PREFIX.pwo"
# Then we need to specify if some modules must be loaded in the submission script
my_hpc.load_modules = """
# Here this is a bash script at the beginning of the submission
# We can load modules

module load QuantumESPRESSO
export OMP_NUM_THREADS=1
"""

# All these information are independent from the calculation
# Now we need some more specific info, like the number of processors, pools and other stuff
my_hpc.n_cpu = 40 # We will use 32 processors
my_hpc.n_nodes = 1 #In 1 node
my_hpc.n_pool = 20 # This is an espresso specific tool, the parallel CPU are divided in 4 pools

# We can also choose in how many batch of jobs we want to submit simultaneously, and how many configurations for each job
my_hpc.batch_size = 10
my_hpc.job_number = 20
# In this way we submit 10 jobs, each one with 10 configurations (overall 100 configuration at time)

# We give 25 seconds of timeout
my_hpc.set_timeout(30)

# We can specify the time limit for each job,
my_hpc.time = "00:20:00" # 5 minutes

my_hpc.setup_workdir()
```


```python
# We prepare the Automatic relaxation
relax = sscha.Relax.SSCHA(minimizer, ase_calculator = espresso_calc,
                         N_configs = 200,
                         max_pop = 6,
                         save_ensemble = True,
                         cluster = my_hpc)
```

Lets see the parameters:
    
    minimizer     : it is the SSCHA_Minimizer, containing the settings of each minimization.
    N_configs     : the number of configurations to generate at each run (we call them populations)
    max_pop       : The maximum number of populations after wich the calculation is stopped even if not converged.
    save_ensemble : If True, after each energy and force calculation, the ensemble will be saved.
    cluster       : The cluster object to be used for submitting the configurations.
    

If no cluster is provided (like in this case) the calculation is performed locally.
The save_ensemble keyword will save the ensemble inside the directory specified in the running command
```python
relax.relax(ensemble_loc = "directory_of_the_ensemble")
```

The calculation can be run just by calling
```python
relax.relax()
```
And the code will proceed with the automatic SSCHA minimization.
However, before doing so, we will setup a custom function.
These are usefull to manipulate a bit the minimization, printing or saving extra info
during the minimization.
In particular, we will setup a function to be called after each minimization step of each population, that
will save the current frequencies of the auxiliary dynamical matrix in an array, so that we will be able
to plot the whole frequency evolution after the minimization.

The following cell will start the computation. Consider that it may take long time.
Running espresso in parallel with 4 processors on a Intel(R) Core(TM) i7-4790K CPU at 4.00GHz the single ab-initio run takes a bit more than 1 minute. In this example we are running 20 configurations per population and up to 6 populations. 
So the overall time is about 2 hours and half.

For heavier relaxations it is strongly suggested to use the cluster module, or to directly run the SSCHA from a more powerful computer.

Thanks to the line inside our custom function
```python
np.savetxt("all_frequencies.dat", all_frequencies)
```
we will save at each minimization step the evolution of the frequencies.
So we can load and plot them even during the minimization, to get a clue on how frequencies of the SSCHA dynamical matrix are evolving.

**NOTE**: These frequencies are not the physical frequencies of the phonon quasiparticles, but just the frequencies of the auxiliary dynamical matrix. To obtain the real physical frequencies, look at the [structural instability](http://sscha.eu/Tutorials/StructuralInstability/) and the (spectral properties)[http://sscha.eu/Tutorials/tutorial_spectral/] tutorials.


```python
# We reset the frequency array
all_frequencies = []

# We define a function that will be called by the code after each minimization step, 
# passing to us the minimizer at that point
# We will store the frequencies, to plot after the minimization their evolution.
def add_current_frequencies(minimizer):
    # Get the frequencies
    w, p = minimizer.dyn.DiagonalizeSupercell()
    
    all_frequencies.append(w)
    
    # In this way the file will be updated at each step
    np.savetxt("all_frequencies.dat", all_frequencies)

# We add this function to the relax ojbect
relax.setup_custom_functions(custom_function_post = add_current_frequencies)

# Now we are ready to
# ***************** RUN *****************
relax.relax(ensemble_loc = "data_ensemble_autorelax")
# ***************************************
```


```python
# We can save the frequency list for furture analisys
# NOTE: This will overwrite the one given in the examples!
np.savetxt("all_frequencies.dat", all_frequencies)
```




```python
# Lets plot the minimization
relax.minim.plot_results()
```


```python
# Lets plot all the frequencies in their evolution
all_frequencies = np.loadtxt("all_frequencies.dat")
freqs = np.array(all_frequencies).T 
n_freqs = shape(freqs)[0]

plt.figure(dpi = 100)
plt.title("Frequencies")
plt.ylabel("Energy [cm-1]")
plt.xlabel("Steps")
for i in range(n_freqs):
    plt.plot(freqs[i, :] * CC.Units.RY_TO_CM)

```




<img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAyAAAAJYCAYAAACadoJwAAAgAElEQVR4Xuy9CZglV1n//6397r3N9OyTyQZJIAlJEEKC5AcoARcii48KUQJEwoMgIqsgAsom4IKALIKAgiIooqAPuOFfSQBJIAkhIftk9unp/a61/5/3VNXt2z09M919b9/pnv5WUnPqVp06p+pT595+v/We8x4tjuMYXEiABEiABEiABEiABEiABEigDwQ0CpA+UGYVJEACJEACJEACJEACJEACigAFCBsCCZAACZAACZAACZAACZBA3whQgPQNNSsiARIgARIgARIgARIgARKgAGEbIAESIAESIAESIAESIAES6BsBCpC+oWZFJEACJEACJEACJEACJEACFCBsAyRAAiRAAiRAAiRAAiRAAn0jQAHSN9SsiARIgARIgARIgARIgARIgAKEbYAESIAESIAESIAESIAESKBvBChA+oaaFZEACZAACZAACZAACZAACVCAsA2QAAmQAAmQAAmQAAmQAAn0jQAFSN9QsyISIAESIAESIAESIAESIAEKELYBEiABEiABEiABEiABEiCBvhGgAOkbalZEAiRAAiRAAiRAAiRAAiRAAcI2QAIkQAIkQAIkQAIkQAIk0DcCFCB9Q82KSIAESIAESIAESIAESIAEKEDYBkiABEiABEiABEiABEiABPpGgAKkb6hZEQmQAAmQAAmQAAmQAAmQAAUI2wAJkAAJkAAJkAAJkAAJkEDfCFCA9A01KyIBEiABEiABEiABEiABEqAAYRsgARIgARIgARIgARIgARLoGwEKkL6hZkUkQAIkQAIkQAIkQAIkQAIUIGwDJEACJEACJEACJEACJEACfSNAAdI31KyIBEiABEiABEiABEiABEiAAoRtgARIgARIgARIgARIgARIoG8EKED6hpoVkQAJkAAJkAAJkAAJkAAJUICwDZAACZAACZAACZAACZAACfSNAAVI31CzIhIgARIgARIgARIgARIgAQoQtgESIAESIAESIAESIAESIIG+EaAA6RtqVkQCJEACJEACJEACJEACJEABwjZAAiRAAiRAAiRAAiRAAiTQNwIUIH1DzYpIgARIgARIgARIgARIgAQoQNgGSIAESIAESIAESIAESIAE+kaAAqRvqFkRCZAACZAACZAACZAACZAABQjbAAmQAAmQAAmQAAmQAAmQQN8IUID0DTUrIgESIAESIAESIAESIAESoABhGyABEiABEiABEiABEiABEugbAQqQvqFmRSRAAiRAAiRAAiRAAiRAAhQgbAMkQAIkQAIkQAIkQAIkQAJ9I0AB0jfUrIgESIAESIAESIAESIAESIAChG2ABEiABEiABEiABEiABEigbwQoQPqGmhWRAAmQAAmQAAmQAAmQAAlQgLANkAAJkAAJkAAJkAAJkAAJ9I0ABUjfULMiEiABEiABEiABEiABEiABChC2ARIgARIgARIgARIgARIggb4RoADpG2pWRAIkQAIkQAIkQAIkQAIkQAHCNkACJEACfSSgadpJa7vmmmvw3//93328oo1X1fXXX4/Pf/7z+N///V88+clP3ngAeMckQAIkcJoJUICc5gfA6kmABDYWgUyAvOhFL1r0xi+44AK86U1v2lhQ+ny3FCB9Bs7qSIAESGABAQoQNgkSIAES6COBTIDEcdzHWllVJ4HDhw9jZmYGZ511FvL5POGQAAmQAAn0mQAFSJ+BszoSIIGNTYACZGM/f949CZAACZAAQAHCVkACJEACfSSwHAHywAMP4Pzzz8fTn/50/MM//APe/va34x//8R9x8OBBvPrVr8YHPvCB9pX/67/+Kz784Q/j//7v/1CtVrFjxw485znPwVve8hYMDw8fd4cTExN461vfiq985SuYmprCOeecg5e//OVqtW0b5557LqT+bPnd3/1dvOtd78Jf//VfQ7owLVx27tyJI0eOIAiC447dfffdeO9734tvfvObGBsbw+DgIJ72tKfh937v93DhhRfOy//JT34Sv/7rv44/+IM/UPVId7T/+I//QL1ex8UXX4y3ve1t+Nmf/dlFn9iPfvQjxUTqES/HwMAAHvWoRykOv/VbvwXDMNR5J+uCJfV88IMfxBe/+EXcf//90HUdj33sY/GKV7wCv/qrv3pcvQ8//HD73vbv3688Ktu2bcNP/uRP4rWvfa16flxIgARIgATmE6AAYYsgARIggT4SWIkAedKTnoRaraaEhwxSl+Xyyy+HiAJZXve61+GP/uiP4DgOfuInfgJbtmzBHXfcoQSEGMA333wzNm/e3L5LER9XX3017r33Xmzfvl0NxJ6cnFSD30WAiJDplQAR4fSCF7wAnufhsssuU+Xu27dPCaVSqYSvf/3r6lqyJRMgL37xi/HV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OHbsmBpzISF6n/KUp6i38SIYOhfJ89a3vhX/9E//hKmpKTXwWroyydwcElVqMQHiuq4a3yAiRDwB4gURb4pMDihduSTcbmcXrKw+ETrSbUs8ERK9SsoXw15mGpd5NZ797GerfbKsVIDIuSIcpPuYCCm5PxEREgpY6vjN3/xNJQpkWawLVnat4pWR8Lp/93d/p0Ibi3CSEMLCQwSHRL0SrhlzERnSfUvEiQgyOfbEJz5RCRW5Py4kQAIkQALHE+AYELYKEiABEiCBNgERENJVaTEBQkwkQAIkQAIk0AsCFCC9oMgySIAESOAMIUABcoY8SN4GCZAACaxhAhQga/jh8NJIgARIoN8EKED6TZz1kQAJkMDGI0ABsvGeOe+YBEiABE5IgAKEjYMESIAESGC1CVCArDZhlk8CJEACJEACJEACJEACJNAmQAHCxkACJEACJEACJEACJEACJNA3AhQgfUPNikiABEiABEiABEiABEiABChA2AZIgARIgARIgARIgARIgAT6RoACpG+oWREJkAAJkAAJkAAJkAAJkAAFyBpqA+Pj4/jGN76BPXv2IJ/Pr6Er46WQAAmQAAmQAAmQAAkIgWazib179+Laa6/Fpk2bCGUFBChAVgBttU75/Oc/j+uvv361ime5JEACJEACJEACJEACPSLwuc99Di984Qt7VNrGKoYCZA0975tvvhlPfvKTIQ36wgsvXENXxkshARIgARIgARIgARIQAvfcc496Yfytb30LV199NaGsgAAFyAqgrdYp3//+93HFFVfgtttuw+WXX75a1bBcEiABEiABEiABEiCBFRKgvbZCcB2nUYB0z7BnJbBB9wwlCyIBEiABEiABEiCBVSFAe617rBQg3TPsWQls0D1DyYJIgARIgARIgARIYFUI0F7rHisFSPcMe1YCG3TPULIgEiABEiABEiABElgVArTXusdKAdI9w56VwAbdM5QsiARIgARIgARIgARWhQDtte6xUoB0z7BnJbBB9wwlCyIBEiABEiABEiCBVSFAe617rBQg3TPsWQls0D1DyYJIgARIgARIgARIYFUI0F7rHisFSPcMe1YCG3TPULIgEiABEiABEiABElgVArTXusdKAdI9w56VwAbdM5QsiARIgARIgARIgARWhQDtte6xUoB0z7BnJbBB9wwlCyIBEiABEiABEiCBVSFAe617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wASYABNgAm1LgO21+FXbNAGyZs0arFixAnfddZc48nZsAtygj82IYzABJsAEmAATYAJMYD4JsL0Wn37TBEg2m8VNN92Ea6+9Nn4ul0gK3KCXSEVzMZkAE2ACTIAJMIFFS4DttfhV1zQBQqudn3vuubjhhhvi53KJpMANeolUNBeTCTABJsAEmAATWLQE2F6LX3VNEyBUOa961atwxx134OUvf3n8nC6BFLhBL4FK5iIyASbABJgAE2ACi5oA22vxq69pAoSy9sMf/hBvfvObsXbtWqxbt04MSK/faBasn/3sZ/FL0SYpcINuk4rkYjABJsAEmAATYAJtS4DttfhV2zQBQuLj8ssvh+d5yGQyM86CtXfv3vilaJMUuEG3SUVyMZgAE2ACTIAJMIG2JcD2WvyqbZoAOeWUU0Tu/vmf/xmbN2+On9MlkAI36CVQyVxEJsAEmAATYAJMYFETYHstfvU1TYCk02nceOON+PM///P4uZwmhUKhgM997nPYunUrHn74YRw6dAhXXXUVbrvttimxyQtDM3LdeuutII8LdQm7+uqrcd11103pFtasuLOBwA16NpQ4DhNgAkyACTABJsAE5o8A22vx2TdNgLzgBS/Aa17zGnziE5+In8tpUti1axdOOOEErFy5Elu2bBHrjcwkQN773vfiq1/9Kt7xjnfgvPPOwwMPPCCECoV/+ctfnpB6s+LOBgI36NlQ4jhMgAkwASbABJgAE5g/AmyvxWffNAFy77334q1vfSt++tOfNqULlmVZGBwcxOrVq+G6LgzDmFaAPP744zjjjDPwvve9DzfffHOV2Ac+8AF86UtfwqOPPorTTz9dhDcr7myriRv0bElxPCbABJgAE2ACTIAJzA8Bttfic2+aAHnb296GRx55BE899RRe+MIXYv369dPOgnX77bfHLsXRBMjHP/5xfOYzn8GOHTuExyTadu7ciY0bN4Luf/rTnxbBzYo72wJyg54tKY7HBJgAE2ACTIAJMIH5IcD2WnzuTRMgqqoeM3c0DS+NuYi7HU2AvOIVrxBeDhojMnlbvnw5zjzzTNx9993iVrPizrZ83KBnS4rjMQEmwASYABNgAkxgfgiwvRafe9MESPyszT6FowkQ6l5lmqYYrD55O+uss+A4juh6RVuz4k5XkoMHD4L2+o28RW95y1tEXilvvDEBJsAEmAATYAJM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JE888YR85zvfkdNOO81QjAXIUUcdJXpO47TbbrvJ/vvv37Tv9NNPl2uvvVZOPPFEOfDAA2XJkiVy8803i+6/5pprJlUyOPSksHESBCAAAQhAAAIQaBsB4jV71DNegKgwePe73y35fL5Gq1AoyD777CMrVqyQ5cuXSyaTqQmQiy66SC699NINkn388cdl7733ljPOOEOuuuqqWtozzzxTrr76ann00Udlzz33bLl0cOiWkXECBCAAAQhAAAIQaCsB4jV73DNegIyH6JxzzpErr7zSNJ/abrvtmgSIihCdent7xzxdj19++eXy3HPPNTXb0iZeu+yyi0xExIxlGIe2d2gsQAACEIAABCAAgTQJEK/Z0+1aAfLpT39afvzjH8sbb7wh/f39NQGi69pMSydtenXWWWeZZlWN07x580wtx+uvvz6qBLbeemvZd9995a677mq5dHDolpFxAgQgAAEIQAACEGgrAeI1e9xdKUCefPJJ0wTriCOOkJ/85CeG4osvvignnXSSHH300bLDDjvIq6++KjfccIOok5177rmmw3k8afOqXC4njzzyyKgS0H4m5XJZtJnWhqbXXntNdG6cNF/HH3+8sat2mCAAAQhAAAIQgAAEOosAAsS+PLpOgKxdu1YOOOAAIzC0FmPHHXccl6KOdHXIIYfIQw89JE8//bTsuuuuJq0utaZj6dKlo87VDunar+SZZ57ZYOlcfPHFcskll4yZBgFi79hYgAAEIAABCEAAAmkQQIDYU+0qAaKdz7X51O9+9zu588475dBDD90owdtvv12OPPJIuf766+Xkk0826akB2Sg2EkAAAhCAAAQgAIEZSQABYl+sXSNASqWSERL33Xef6fuhw+1OZHrsscfMiFeXXXaZXHjhheYU+oBMhBxpIAABCEAAAhCAwMwjgACxL9OuECCVSsV820NrM7773e/KZz7zmQmTu+2228w3RPSDg/pBQp1UiFxxxRXjjoKlx1WwtDrh0K0SIz0EIAABCEAAAhBoLwHiNXveM16A+L4vxx13nPzoRz+S6667Tk455ZQxqWm/ja222qrpmDbZmjt3rvlo4bPPPmu+dq6T9h3Rka7G+w7IsmXLZK+99mq5dHDolpFxAgQgAAEIQAACEGgrAeI1e9wzXoCcffbZ8o1vfEPe//73yxe+8IVRxD70oQ+ZDuU6+tWqVavksMMO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AAEIQAACEJhCAgiQKYTPpSEAAQh0EoGnn35aLr30Ulm6dKm8/PLLMnv2bNlhhx3k4IMPlgULFsi2224r//u//yv/9V//JZ/73Odkp5126qTskxcIQAACEJgmBBAg06SgyCYEIACBNAn87ne/kw984AOy2WabGXGx8847y8qVK+WJJ56QX/ziF3LbbbeZ44sXL5YvfOEL8sADD5htJghAAAIQgECrBBAgrRIjPQQgAIEZSOCII46QX/3qV/LUU0/JnDlzmu6wUChIuVyWTTbZBAEyA8ueW4IABCDQbgIIkHYT53oQgAAEOpDA29/+dunp6ZFly5aNm7uLL75YLrnkklHHb7rpJlNrotOLL75o0tx5552mBkXFzHHHHSf/9E//JPl8vnau1p4888wzRvR86Utfkt/85jfm+n/7t38rixYtkv7+/lraF154Qf7xH//RpF2xYoVsuumm8s53vlO+8pWvyGGHHdaBNMkSBCAAAQhsiAACBP+AAAQgAAH56Ec/appV3X///XLggQeOSeSxxx6Tb33rW3LjjTfKhRdeKHvssYdJp+l32WUXee655+SAAw6QbDYr8+fPN31Gfv/734sKlHnz5skdd9whjuOYc1SAPP7440ZoHHTQQTJ37lz57W9/K9/97nflQx/6kPzyl7806bTmRcXGunXr5NRTTzV9UlSEaJMxvda5555L6UEAAhCAwDQjgACZZgVGdiEAAQikQeDXv/61fPCDH5RKpSL77LOPEQX777+/EQNbbbVV7ZIb6gPysY99zPQZ+cMf/iCbb7557Zxvf/vbcsYZZ8hdd91lhEgsQB588EE5++yz5etf/3ot7Ze//GVZuHCh6ej+8Y9/XB599FGTn3//9383tSNMEIAABCAw/QkgQKZ/GXIHEIAABBIh8PDDD5vg/+677zY1DjplMhk57bTTjEjQmo3xBMiaNWvkzW9+s/zDP/yDXHDBBU35Wb16tey+++6mtkKbVzUKkFdeecXUlMTT8uXLZeuttzYd3W+44Qb5y1/+Ykbb0iZeV111lemHwgQBCEAAAtObAAJkepcfuYcABCCQOAHf9+XPf/6zaY6lwuPZZ581/Tq0H8d4AkSbRL3vfe/bYF4++9nPyi233FITIDqkbyx0Gk/Ukbje7XrhpAAAA+NJREFU/e53yz333GN2X3TRRXLFFVcYMfTe977X1KJ86lOfkre+9a2J3zsGIQABCEAgfQIIkPQZcwUIQAAC05bAqlWrZNdddzW1GypExhMg2n9Dm2xpP41jjjlmzPvdZpttTH8OnbQPyEQFiKbXDuu33367aLOte++9V0qlksmLihomCEAAAhCYXgQQINOrvMgtBCAAgbYTeNe73iX/93//J8PDw6YDunYwH/kdEB3xSvuKnHzyyXLddddtNI8qQFRMbKwJ1liG3njjDXnPe95jRIiOusUEAQhAAALTiwACZHqVF7mFAAQgkAoBrVU49NBDxfO8Jvs6stXf/M3fmD4c2rn8hz/8oXz60582HyY8+uijm9Jq06j//u//lkceeaQ2QlacoFgsis5xH45YgIzXCf1nP/uZHHnkkbJ27Vrp6+sz/U8ap8MPP9yMmrV+/fpUeGAUAhCAAATSI4AASY8tliEAAQhMGwIqMrQj+VFHHWUEh/a3ePrpp02fDe1EroJAR7nSviFve9vbTA2Edk7v7e01fT/0y+kqVnQ4Xe3XcdJJJxk7g4OD5uOGP/7xj+VHP/qRqHDQSQWIDus7e/ZsOfjgg81QvvEwvDoal/b/0CF7f/rTn9aadakI0m+FaM3J97//fbP/2muvnTaMySgEIAABCIQEECB4AgQgAAEImJGvtFZjyZIlplnUwMCAaVKl39o455xzzDKerrnmGrnyyivNCFXVatV85yP+EOFrr70ml19+ufnmh9rRGg8VJype9IODW2yxRU2AaL8ObcqlQ/Tqhwj1Q4XHHnusfO1rXzPCRKfnn3/edEDXYYJffvllcV3XfHPkxBNPlC9+8YtGKDFBAAIQgMD0IoAAmV7lRW4hAAEIzAgC8ZfQVVQwQQACEIBAdxFAgHRXeXO3EIAABDqCAAKkI4qBTEAAAhCYEgIIkCnBzkUhAAEIdDcBBEh3lz93DwEIdDcBBEh3lz93DwEIQGBKCCBApgQ7F4UABCDQEQQQIB1RDGQCAhCAAAQgAAEIQAAC3UEAAdId5cxdQgACEIAABCAAAQhAoCMIIEA6ohjIBAQgAAEIQAACEIAABLqDAAKkO8qZu4QABCAAAQhAAAIQgEBHEECAdEQxkAkIQAACEIAABCAAAQh0BwEESHeUM3cJAQhAAAIQgAAEIACBjiCAAOmIYiATEIAABCAAAQhAAAIQ6A4CCJDuKGfuEgIQgAAEIAABCEAAAh1BAAHSEcVAJiAAAQhAAAIQgAAEINAdBBAg3VHO3CUEIAABCEAAAhCAAAQ6ggACpCOKgUxAAAIQgAAEIAABCECgOwj8/8EpamnAeHDwAAAAAElFTkSuQmCC" width="640">


# The Structural phase transition

We got the minimum of the free energy, now we can compute the free energy curvature around this high symmetry structure, to discover if there are some imaginary phonons. 
In order to do this, we need the  free energy Hessian:

$$
\frac{\partial^2 F}{\partial R_a \partial R_b} \approx \Phi_{ab} + \sum_{pqrs}\stackrel{(3)}{\Phi}_{apq} \Lambda_{pqrs} \stackrel{(3)}{\Phi}_{rsb}
$$

This is obtained neglecting 4 phonon scattering (that may be relevant for highly anharmonic systems).

You can compute this quantity using the get_free_energy_hessian subroutine of the ensemble class.
However, you need a bigger ensemble to get reliable results. We will use the last ensemble, but for production calculations it is advisable to generate a new, bigger, ensemble.

Here, the flag include_v4 tells if you need to include the V4


```python
# Compute the free energy hessian
free_energy_hessian = relax.minim.ensemble.get_free_energy_hessian(include_v4 = False)
```


```python
# Now we can save the free energy hessian like it was a dynamical matrix
free_energy_hessian.save_qe("hessian")

# We can print its eigenvalues to see if there are imaginary phonons
w, p = free_energy_hessian.DiagonalizeSupercell()
print(w * CC.Units.RY_TO_CM)
```

    [-3.86241653e-05 -2.54094408e-05  1.93744849e-05  3.13241198e+02
      3.13241198e+02  3.13241198e+02  3.13241198e+02  3.13241198e+02
      3.13241198e+02  3.30348694e+02  3.30348694e+02  3.30348694e+02
      4.15342661e+02  4.15342661e+02  4.15342661e+02  4.15342661e+02
      4.15342661e+02  4.15342661e+02  5.45615319e+02  5.45615319e+02
      5.45615319e+02  5.45615319e+02  5.45615319e+02  5.45615319e+02
      7.42041373e+02  7.42041373e+02  7.42041373e+02  7.52648326e+02
      7.52648326e+02  7.52648326e+02  8.43506806e+02  8.43506806e+02
      8.43506806e+02  8.43506806e+02  8.43506806e+02  8.43506806e+02
      1.08903503e+03  1.08903503e+03  1.08903503e+03  1.08903503e+03
      1.08903503e+03  1.08903503e+03  1.20924095e+03  1.20924095e+03
      1.20924095e+03  1.21606945e+03  1.21606945e+03  1.21606945e+03
      1.21606945e+03  1.21606945e+03  1.21606945e+03  1.22775643e+03
      1.22775643e+03  1.22775643e+03  1.22775643e+03  1.22775643e+03
      1.22775643e+03  1.23314055e+03  1.23314055e+03  1.23314055e+03
      1.23314055e+03  1.23314055e+03  1.23314055e+03  1.25031872e+03
      1.25031872e+03  1.25031872e+03  1.25031872e+03  1.25031872e+03
      1.25031872e+03  1.34998026e+03  1.34998026e+03  1.34998026e+03
      1.43074091e+03  1.43074091e+03  1.43074091e+03  1.43074091e+03
      1.43074091e+03  1.43074091e+03  1.53256354e+03  1.53256354e+03
      1.53256354e+03  1.57033447e+03  1.57033447e+03  1.57033447e+03
      1.57033447e+03  1.57033447e+03  1.57033447e+03  1.63331416e+03
      1.63331416e+03  1.63331416e+03  1.63331416e+03  1.63331416e+03
      1.63331416e+03  1.67069280e+03  1.67069280e+03  1.67069280e+03]


No imaginary phonon! The structure is stable at this pressure.

Remember, to get converged result, you should study this result as a function of the number of configurations and as a function of the supercell size!


<a name="Lanthanum-hydride"></a>
# 6. Lanthanum hydride

In this tutorial we will replicate the results of the work [Errea et. al, Nature 578, 66–69(2020)](https://www.nature.com/articles/s41586-020-1955-z), where they show how the rhombohedral structure $$R-3m$$ of LaH10, local minimum of the Born-Oppenheimer (BO) energy landscape,  collapses into the higher symmetry $$Fm-3m$$ at low pressure, where this phase is unstable within the harmonic approximation. 

**Note**: We will use underconverged parameters to be able run the calculation on a local computer, for production runs you must study the convergence, especially with the k point sampling and the supercell size.

For this tutorial, we provide a dynamical matrix, obtained by a sscha relaxation at fixed cell (see tutorial on [H3S](http://sscha.eu/Tutorials/Automatic_Calculations/) for example). Also the pseudopotentials used are provided. These input files can be found in python-sscha/Tutorials/LaH10.  


```python
%pylab
from __future__ import print_function
```

## The ab-initio parameters

Here we setup the calculator to compute the energy, forces, and pressures of the BO energy landscape. Please, refer to [Espresso pw.x guide](https://www.quantum-espresso.org/Doc/INPUT_PW.html) for the detailed description on the input parameters, and [ASE Espresso calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/espresso.html) for a description on how to properly setup an Espresso calculator suited for your application.


```python
import ase
from ase.calculators.espresso import Espresso

import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax

pseudo = {"H": "H.pbe-rrkjus_psl.1.0.0.UPF",
         "La" : "La.pbe-spfn-rrkjus_psl.1.0.0.UPF"}
input_params = {"tstress" : True, # Print the stress in the output
                "tprnfor" : True, # Print the forces in the output
                "ecutwfc" : 35,  #The wavefunction energy cutoff for plane-waves (Ry)
                "ecutrho" : 350, # The density energy cutoff (Ry)
                "mixing_beta" : 0.2,  # The mixing parameter in the self-consistent calculation
                "conv_thr" : 1e-9,    # The energy convergence threshold (Ry)
                "degauss" : 0.02,  # Smearing temperature (Ry)
                "smearing" : "mp",
                "pseudo_dir" : ".",
                "occupations" : "smearing",
               "disk_io" : "none"}

k_points = (8,8,8) # The k points grid (you can alternatively specify a kspacing)
k_offset = (1,1,1) # The offset of the grid (can increase convergence)

espresso_calc = Espresso(pseudopotentials = pseudo, input_data = input_params, 
                        kpts = k_points, koffset = k_offset)
```

## The preparation of the minimization
In the following cell we prepare the minimization parameters and the ensemble.


```python
# We now load the dynamical matrix
dyn = CC.Phonons.Phonons("dyn")
dyn.Symmetrize() #Enforce the sum rule

# We prepare the ensemble
ensemble = sscha.Ensemble.Ensemble(dyn, T0 = 0, supercell = dyn.GetSupercell())

# We prepare the sscha minimizer
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

# We set up the minimization parameters
minim.min_step_dyn = 0.05     # The minimization step on the dynamical matrix
minim.min_step_struc = 0.05   # The minimization step on the structure
minim.kong_liu_ratio = 0.5     # The parameter that estimates whether the ensemble is still good
minim.gradi_op = "all" # Check the stopping condition on both gradients
minim.meaningful_factor = 0.2 # How much small the gradient should be before I stop?
```

We setup the standard sscha minimiztion. Now we must prepare the calculator for the automatic relaxation exactly like we did for the H3S example.
Remember you can always specify a cluster for the automatic calculation:

*In case you want to use a cluster, please remember to upload the pseudos on the cluster working directory!*

Skip the following cell if you do not want to setup a cluster, and run the calculation locally.


```python
# Here we prepare a cluster
# Here we configure the cluster object MARCONI
import sscha.Cluster
my_hpc = sscha.Cluster.Cluster(pwd = None)

# We setup the connection info
my_hpc.hostname = "ekhi" # The command to connect via ssh to the cluster
#my_hpc.account_name = "IscrB_COMRED" # The name of the project for the computation
my_hpc.workdir = "/scratch/lorenzo/my_calculation" # the directory in which the calculations are performed

# Now we need to setup the espresso
# First we must tell the cluster where to find him:
my_hpc.binary = "pw.x -npool NPOOL -i  PREFIX.pwi > PREFIX.pwo"
# Then we need to specify if some modules must be loaded in the submission script
my_hpc.load_modules = """
# Here this is a bash script at the beginning of the submission
# We can load modules

module load QuantumESPRESSO
export OMP_NUM_THREADS=1
"""

# All these information are independent from the calculation
# Now we need some more specific info, like the number of processors, pools and other stuff
my_hpc.n_cpu = 32 # We will use 32 processors
my_hpc.n_nodes = 1 #In 1 node
my_hpc.n_pool = 16 # This is an espresso specific tool, the parallel CPU are divided in 4 pools

# We can also choose in how many batch of jobs we want to submit simultaneously, and how many configurations for each job
my_hpc.batch_size = 20
my_hpc.job_number = 20
# In this way we submit 10 jobs, each one with 10 configurations (overall 100 configuration at time)

# We give 25 seconds of timeout
my_hpc.set_timeout(25)

# We can specify the time limit for each job,
my_hpc.time = "00:10:00" # 5 minutes

# Create the working directory if none on the cluster 
# And check the connection
my_hpc.setup_workdir()
```

## Prepare the automatic relaxation
Now we prepare the relaxation object. 
We use the same object (SSCHA) inside the Relax module as we did for the H3S example.
set my_hpc to None if you want to run the calculation locally


```python
# Decomment the following line if you did not set up the cluster
#my_hpc = None

relax = sscha.Relax.SSCHA(minim, ase_calculator = espresso_calc,
                         N_configs = 400,
                         max_pop = 20,
                         save_ensemble = True,
                         cluster = my_hpc)
```

This time we are interested in plotting the symmetry analisys, as we evolve the minimization


```python
import spglib
print ("The original spacegroup is:", spglib.get_spacegroup(dyn.structure.get_ase_atoms(), 0.05))
```

    The original spacegroup is: R-3m (166)


We create a custom function to print the spacegroup after each iteration of the minimization.
In this way we can follow the evolution of the dynamical matrix as it evolves.
We use a threshold for symmetries of 0.05 Å. Remember that the SSCHA is a stochastic method, the atomic position is affected by stochastic noise.
If you want to increase your accuracy in the identification of the space group, you should accordingly increase the number of configurations, remember that the stochastic noise scales as $$1/\sqrt{N_{configs}}$$.


```python
# we define a function that prints the space group during the optimization
space_groups = []
def print_spacegroup(minim):
    spgroup = spglib.get_spacegroup(minim.dyn.structure.get_ase_atoms(), 0.05)
    space_groups.append(spgroup)
    
    # We can save them in the output at each minimization step
    f = open("space_group.dat", "w")
    f.writelines(["{}) {}\n".format(i+1, x) for i,x in enumerate(space_groups)])
    f.close()
    
relax.setup_custom_functions(custom_function_post = print_spacegroup)
```

## Perform the variable cell relaxation
We are ready to start a variable cell relaxation. There are two different variable cell relaxations implemented in the SSCHA code: target pressure or fixed volume.
In the target pressure calculation, the SSCHA adapts the cell until the stress stress tensor is uniform and reproduces the pressure we want. In this calculation the volume changes. In the fixed volume, instead, the SSCHA optimizes the cell parameters keeping the overall volume unchanged.

We will use the latter in this example, but feel free to experiment by selecting fix_volume to false, and manually change the target_pressure argument. We prepared another command (commented) to perform a variable cell relaxation with target pressure.

The advantage of keeping the volume fixed is that we reduce the number of degrees of freedom in the cell optimization.

The static_bulk_modulus is a flag that allow the program to estimate, given the stress tensor, how to change the unit cell to reach the optimal value in the lowest number of steps. A good value is the static bulk modulus, that is the derivative pressure with respect to the volume times the volume. The code expects it in GPa.
Usually for high pressure materials, the bulk modulus is around hundreds of GPa, while for ice at ambient pressure is about 10 GPa. A high value of the bulk modulus will mean a slower change in the unit cell, so if you have the fealing that the unit cell is not changing a lot between sequent steps, try to reduce it.  


```python
# Now we can run the calculation!!!
# In this case we fix the volume (we optimize lattice parameters)
# But you can also fixe the target pressure (as done in the commented line)
import os
if not os.path.exists("ensembles"):
    os.mkdir("ensembles")
relax.vc_relax(fix_volume = True, static_bulk_modulus = 120, ensemble_loc = "ensembles")
#relax.vc_relax(target_press = 120, static_bulk_modulus = 200, ensemble_loc = "ensembles")
```

## Results
The minimization is done, now we can study the evolution of the rhombohedral angle in subsequent populations
![Rombohedral angle](romb_ang.png)

As you can see, the angle is getting close to 60 degrees. 
We can see if it recognize a closer high symmetry structure.



```python
relax.minim.finalize()
relax.minim.plot_results()
```

    
     * * * * * * * * 
     *             * 
     *   RESULTS   * 
     *             * 
     * * * * * * * * 
    
    
    Minimization ended after 432 steps
    
    Free energy = -1853233.86250700 +-       6.13461465 meV
    FC gradient modulus =       3.00827542 +-       3.10405062 bohr^2
    Struct gradient modulus =      40.59413639 +-     200.71832520 meV/A
    Kong-Liu effective sample size =  396.12419553802914
    
    
     ==== STRESS TENSOR [GPa] ==== 
        161.29075385      0.00000000      0.00000000                0.34837928      0.00000000      0.00000000
          0.00000000    161.29075385     -0.00000000    +-          0.00000000      0.34837928      0.00000000
          0.00000000     -0.00000000    161.29043593                0.00000000      0.00000000      0.55826675
    
     Ab initio average stress [GPa]:
        157.36742136      0.00000000      0.00000000
          0.00000000    157.36742136     -0.00000000
          0.00000000     -0.00000000    157.43887873
    


The previous command, as always, plots the results of the minimizations and prints the information on the last minimization, as the free energy and the final stress tensor.

It is printing two kinds of stress tensors: the real one, and the average over the ab-initio stresses.

The latter is just the stochastic average over the ab-initio stresses, however, it does not include the kinetic contribution, so it is not the real stress (see [Monacelli et. al. Phys. Rev. B 98, 024106](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.98.024106)).


```python
spglib.get_spacegroup(relax.minim.dyn.structure.get_ase_atoms(), 0.05)
```




    u'Fm-3m (225)'




```python
spglib.get_spacegroup(dyn.structure.get_ase_atoms(), 0.1)
```




    u'R-3m (166)'




```python
view(relax.minim.dyn.structure.get_ase_atoms())
```


```python
from ase.visualize import view
```


```python
relax.minim.dyn.structure.unit_cell
```




    array([[ 1.79783817e+00, -1.03798235e+00,  2.91372421e+00],
           [-2.43331324e-15,  2.07596470e+00,  2.91372421e+00],
           [-1.79783817e+00, -1.03798235e+00,  2.91372421e+00]])
