---
layout: page
title: Old Tutorials
---

In order to get familiar with the code we strongly suggest that before attempting a real calculation you run the tutorials provided here.

We suggest that you start to run the tutorials in the following order:

1. [**Lead Telluride**](#Lead-Telluride)

    In this tutorial you will learn how you can set up a calculation from scratch starting from a CIF file. 

    It is illustrated how the ASE package can be used to prepare the input files to run the initial Born-Oppenheimer structural relaxation and harmonic phonon calculation. It is also illustrated how the output of these calculations can be used to start the SSCHA free energy minimization by creating first the ensemble and, second, calculating the energies, forces, and stress tensors for them. It is described how the latter calculations can be performed in three different ways: by using ASE to perform the calculations locally, by running the DFT calculations manually locally or in a cluster, and setting up an automatic submission to a cluster.

2. [**Lead Telluride structural instability**](#Lead-Telluride-structural-instability)

    In this tutorial you will understand how to calculate the free energy Hessian used to determine the stability of the system in the free energy landscape, valid to determine second-order phase transitions such as charge-density wave of ferroelectric transitions.

3. [**Lead Telluride spectral properties**](http://sscha.eu/Tutorials/tutorial_spectral/)

    Here you will learn how to run the SSCHA minimization as an stand-alone program, and also how to calculate the phonon spectral functions, which are in the end what experiments probe. Several approaches to calculate the spectral function are exemplified.   

4. [**Tin Telluride with force fields**](http://sscha.eu/Tutorials/SnTe/)

    In this tutorial you will learn how to automatize a calculation with a python script using a force field. Also how to calculate the with a force field the free energy Hessian at different temperatures by scripting all the calculations. 

5. [**Sulfur hydride**](http://sscha.eu/Tutorials/Automatic_Calculations/)

    In this tutorial you can learn how to automatize a SSCHA minimization. The example works with the H$$_3$$S superconducting compound.

6. [**Lanthanum hydride**](http://sscha.eu/Tutorials/VariableCellRelaxation/)

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


