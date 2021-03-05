# Quantum variable cell relaxation of LaH10

In this tutorial we will replicate the results of the work [Errea et. al, Nature 578, 66–69(2020)](https://www.nature.com/articles/s41586-020-1955-z), where they show how the rombohedral structure R-3m of LaH10, global minimum of the Born-Oppenheimer (BO) energy landscape, to collapse into the higher symmetry Fm-3m at low pressure, where this phase is unstable within the harmonic approximation. 

*Note*: We will use underconverged parameters to be able run the calculation on a local computer, for production runs you must study the convergence, especially with K point sampling and the supercell size.

For this tutorial, we provide a dynamical matrix, obtained by a sscha relaxation at fixed cell (see tutorial on H3S for example). 


```python
%pylab
from __future__ import print_function
```

    Using matplotlib backend: Qt5Agg
    Populating the interactive namespace from numpy and matplotlib


## The ab-initio parameters
Here we setup the calculator to compute the energy, forces and pressures of the BO energy landscape. Please, refer to [Espresso pw.x guide](https://www.quantum-espresso.org/Doc/INPUT_PW.html) for the detailed description on the input parameters, and [ASE Espresso calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/espresso.html) for a description on how to properly setup an Espresso calculator suited for your application.


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

*In the case you want to use this cluster, please remember to upload the pseudos on the cluster working directory!*

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
We use a threshold for symmetries of 0.05 A. Remember that the SSCHA is a stochastic method, the atomic position is affected by stochastic noise.
If you want to increase your accuracy in the identification of the spacegroup, you should accordingly increase the number of configurations, remember that the stochastic noise scales as $1/\sqrt{N_{configs}}$.


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

We will use the latter in this example, but feel free to experiment by selecting fix_volume to false, and manually change the target_pressure argument. I prepared another command (commented) to perform a variable cell relaxation with target pressure.

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
The minimization is done, now we can study the evolution of the Rombohedral angle in subsequent populations
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




```python

```
