---
layout: page
title: Frequently Asked Questions (FAQs)
---

1. [How do I start a calculation if the dynamical matrices have imaginary frequencies?](#How-do-I-start-a-calculation-if-the-dynamical-matrices-have-imaginary-frequencies?)
2. [What are the reasonable values for the steps (lambda_a, lambda_w, min_step_dyn and min_step_struc)?](#What-are-the-reasonable-values-for-the-steps-(lambda_a,-lambda_w,-min_step_dyn-and-min_step_struc)?)
3. [In a variable cell optimization, what is a reasonable value for the bulk modulus?](#In-a-variable-cell-optimization,-what-is-a-reasonable-value-for-the-bulk-modulus?)
4. [The code stops saying it has found imaginary frequencies, how do I fix it?](#The-code-stops-saying-it-has-found-imaginary-frequencies,-how-do-I-fix-it?)
5. [Why the gradient sometimes increases during a minimization?](#Why-the-gradient-sometimes-increases-during-a-minimization?)
6. [The gradients on my simulations are increasing a lot, why is this happening?](#The-gradients-on-my-simulations-are-increasing-a-lot,-why-is-this-happening?)
7. [How do I check if my calculations are well converged?](#How-do-I-check-if-my-calculations-are-well-converged?)
8. [How do I plot the frequencies of the dynamical matrix during the minimization?](#How-do-I-plot-the-frequencies-of-the-dynamical-matrix?)
9. [What is the final error on the structure or the dynamical matrix of a SCHA minimization?](#What-is-the-final-error-on-the-structure-or-the-dynamical-matrix-of-a-SCHA-minimization?)
10. [How does the error over the gradients scale with the number of configurations?](#How-does-the-error-over-the-gradients-scale-with-the-number-of-configurations?)
11. [When I relax the cell, is it necessary for the gradients to reach zero before making a step with the new cell?](#When-I-relax-the-cell,-is-it-necessary-for-the-gradients-to-reach-zero-before-making-a-step-with-the-new-cell?)
12. [I cannot remove the pressure anisotropy after relaxing the cell, what is happening?](#I-cannot-remove-the-pressure-anisotropy-after-relaxing-the-cell,-what-is-happening?)
13. [How may I run a calculation neglecting symmetries?](#How-may-I-run-a-calculation-neglecting-symmetries?)
14. [In which units are the lattice vectors, the atomic positions, and the mass of the atoms in the dynamical matrix file?](#In-which-units-are-the-lattice-vectors,-the-atomic-positions,-and-the-mass-of-the-atoms-in-the-dynamical-matrix-file?)
15. [What is the difference between the different kind of minimization (preconditioning and root_representation)?](#What-is-the-difference-between-the-different-kind-of-minimization-(preconditioning-and-root_representation)?)
16. [How do I lock modes from m to n in the minimization?](#How-do-I-lock-modes-from-m-to-n-in-the-minimization?)
17. [How do I lock a special atom in the minimization?](#How-do-I-lock-a-special-atom-in-the-minimization?)
18. [How do I understand if I have to generate a new population or the minimization converged?](#How-do-I-understand-if-I-have-to-generate-a-new-population-or-the-minimization-converged?)
19. [How do I choose the appropriate value of Kong-Liu effective sample size or ratio?](#How-do-I-choose-the-appropriate-value-of-Kong-Liu-effective-sample-size-or-ratio?)
20. [How do I understand if the free energy hessian calculation is converged?](#How-do-I-understand-if-the-free-energy-hessian-calculation-is-converged?)
21. [How can I add more configurations to an existing ensembe?](#How-can-I-add-more-configurations-to-an-existing-ensembe?)
22. [How do I fix the random number generator seed to make a calculation reproducible?](#How-do-I-fix-the-random-number-generator-seed-to-make-a-calculation-reproducible?) 
23. [When I try to start a sscha calculation or import cellconstructor I get ModuleNotFoundError: No module named 'symph'
](#symph-not-found) 

<a name="How-do-I-start-a-calculation-if-the-dynamical-matrices-have-imaginary-frequencies?"></a>
# How do I start a calculation if the dynamical matrices have imaginary frequencies? 

A good starting point for the SSCHA minimization is usually provided by the dynamical matrices obtained in a standard harmonic calculation. However, they can have imaginary frequencies. This may be related to both instabilities (the structure is a saddle-point of the Born-Oppenheimer energy landscape) or to a not well converged choice of the parameters for computing the harmonic frequencies. In both cases, it is very easy to get a new dynamical matrix that is positive definite and can be used as a starting point. An example is provided in the [Tutorial on H3S](http://sscha.eu/Tutorials/Automatic_Calculations/). Assuming your not positive definite dynamical matrix is in Quantum Espresso format “harm1” … “harmN” (with N the number of irreducible q points), you can generate a positive definite dynamical matrix “positive1” … “positiveN” with the following python script that uses CellConstructor.

```
# Load the cellconstructor library
import cellconstructor as CC
import cellconstructor.Phonons
   
# Load the harmonic not-positive definite dynamical matrix
# We are reading 6 dynamical matrices
harm = CC.Phonons.Phonons("harm", nqirr = 6)

# Apply the acoustic sum rule and the symmetries
harm.Symmetrize()

# Force the frequencies to be positive definite
harm.ForcePositiveDefinite()

# Save the final dynamical matrix, ready to be used in a sscha run
harm.save_qe("positive")
```

The previous script (that we can save into *script.py*) will generate the positive definite matrix ready for the sscha run. It may be executed with

```
python script.py
```

<a name="What-are-the-reasonable-values-for-the-steps-(lambda_a,-lambda_w,-min_step_dyn-and-min_step_struc)?"></a>
# What are the reasonable values for the steps (lambda_a, lambda_w, min_step_dyn and min_step_struc)?

The code minimizes using a Newton method: preconditioned gradient descent. Thanks to an analytical evaluation of the hessian matrix, the step is rescaled so that the theoretical best step is close to 1. In other words: **one is theoretically the  best (and the default) choice for the steps**. However, the SSCHA is a stochastic algorithm, therefore, if the ensemble is too small, or the gradient is very big, this step could bring you outside the region in which the ensemble is describing well the physics very soon.

Since SSCHA can exploit the reweighting technique, and the most computational expensive part of the algorithm is the computation of forces and energies if an *ab initio* approach is taken, it is often much better to use a small step (smaller than the optimal one). **Good values of the steps are usually around 0.01 and 0.1**. Rule of thumb: the minimization should not end because it went outside the stochastic regime before at least 10 steps have been made. This will depend on the system, the number of configurations, and how far from the correct solution you are.

**lambda_w** is the step in the atomic positions (stand-alone program input).

**lambda_a** is the step in the dynamical matrix (stand-alone program input).

If you are using a python script, the equivalent variables are the attributes of the sscha.SchaMinimizer.SSCHA_Minimizer class.

**min_step_struc** is the step in the atomic positions.

**min_step_dyn** is the step in the dynamical matrix.

<a name="In-a-variable-cell-optimization,-what-is-a-reasonable-value-for-the-bulk-modulus?"></a>
# In a variable cell optimization, what is a reasonable value for the bulk modulus?

The bulk modulus is just an indicative parameter used to guess the optimal step of the lattice parameters in order to converge as quickly as possible. It is expressed in GPa. You can find online the bulk modulus for many materials. Find a material similar to the one you are studying and look if the bulk modulus is available in the literature.

Usual values are between 10 GPa and 100 GPa for systems at ambient conditions. Diamond has a bulk modulus of about 500 GPa. High pressure hydrides have a bulk modulus of around 500 GPa as well.

If you have no idea on the bulk modulus, you can easily compute it by doing two static *ab initio* calculations at very close volumes (by varying the cell size), and then computing the differences between the pressures:
\\( B = - \Omega \frac{dP}{d\Omega} \\) where \\(\Omega\\) is the unit-cell volume and \\(P\\) is the pressure (in GPa).

<a name="The-code-stops-saying-it-has-found-imaginary-frequencies,-how-do-I-fix-it?"></a>
# The code stops saying it has found imaginary frequencies, how do I fix it?

This means that your step is too large. You can reduce the step of the minimization. An alternative (often more efficient) is to switch to the root representation. In this way the square root of the dynamical matrix is minimized, and the total dynamical matrix is positive definite in the whole minimization by construction.

In the namelist input you activate this minimization with the following keywords inside the &inputscha namelist
```
preconditioning = .false.
root_representation = "root4"
```
if you are running the SSCHA as a stand-alone program

If you are using a python script instead, you should setup the following attributes of the sscha.SchaMinimizer.SSCHA_Minimizer class
```
minim.preconditioning = False
minim.root_representation = "root4"
```

It is possible that the optimal step size for the root_representation is different than the other one.

<a name="Why-the-gradient-sometimes-increases-during-a-minimization?"></a>
# Why the gradient sometimes increases during a minimization?

Nothing in principle assures that a gradient should always go down. It is possible that at the beginning of the calculation, when we are far from the solution, one of the gradients increases. However, when we get closer to the solution, indeed the gradients must decrease. If this does not happen it could be due to the ensemble that has fewer configurations than necessary. In this case, a good choice is to increase the number of the effective sample size (the kong-liu ratio), in order to stop the minimization when the gradient starts increasing, or to increase the number of configurations in the ensemble.

In any case, what must decrease is the free energy. If you see that the gradient is increasing but the free energy decreases, then the minimization is correct. However, if both the gradient and the free energy are increasing, something is wrong. This could be due to a too big step size. Then try to reduce the value of **lambda_a** and **lambda_w** (in the input file) or **min_step_dyn** and **min_step_struc** (in the python script). It could also be due to a wasted ensemble, in this case, check the value of the Kong-Liu effective sample size. If it is below or around 0.5, then try to increase the threshold at which the calculation is stopped, **kong_liu_ratio** (in the python script) or **N_random_eff** (in the input file), or increase the number of configurations for the next population.

<a name="The-gradients-on-my-simulations-are-increasing-a-lot,-why-is-this-happening?"></a>
# The gradients on my simulations are increasing a lot, why is this happening?

See the previous question.

<a name="How-do-I-check-if-my-calculations-are-well-converged?"></a>
# How do I check if my calculations are well converged?

In general, if the gradient goes to zero and the Kong Liu ratio is above 0.5 probably your calculation converged very well.
There are some cases (especially in systems with many atoms) in which it is difficult to have an ensemble sufficiently big to reach this condition.
In these cases, you can look at the history of the frequencies in the last populations and see whether they are changing or not. 

To plot the frequencies look at [this answer](#How-do-I-plot-the-frequencies-of-the-dynamical-matrix?).




<a name="How-do-I-plot-the-frequencies-of-the-dynamical-matrix?"></a>
# How do I plot the frequencies of the dynamical matrix during the optimization?
To check if the SSCHA is converging, you should plot the dynamical matrix's frequencies during the minimization.
In particular, you should look if, between different populations, the evolution of each frequency is consistent. If it seems that frequencies are evolving randomly from a population to the next one, you should increase the number of configurations, otherwise, you can keep the number fixed.

The code can print the frequencies at each step.
If you run the code with an input script, you should provide in the &utils tag the filename for the frequencies:

```
       &utils
           save_frequencies = "freqs.dat"
       &utils
```

You can use the same function from the python script by calling a custom function that saves the frequencies after each optimization step. The Utilities module of the SSCHA offers this function:

```
IO_freq = sscha.Utilities.IOInfo()
IO_freq.SetupSaving("freqs.dat")

# Initialize the minimizer as minim [...]
minim.run(custom_function_post = IO_freq.CFP_SaveFrequencies)
```

The code here is providing the SSCHA code a function (IO_freq.CFP_SaveFrequencies) that is called after each minimization step. This function saves all the frequencies of the current dynamical matrix in the file specified by IO_freq.SetupSaving("freqs.dat").

To plot the results, the SSCHA offers an executable script, installed together with the code. Just run:

```
plot_frequencies.py freqs.dat
```

And the code will plot all the frequencies. You can also pass more than one file. In this case, the frequencies are concatenated.Plotting the frequencies of the dynamical matrix is a very good way to check if the algorithm is converging correctly.


<a name="What-is-the-final-error-on-the-structure-or-the-dynamical-matrix-of-a-SCHA-minimization?"></a>
# What is the final error on the structure or the dynamical matrix of a SCHA minimization?

This is a difficult question. The best way to estimate the error is to generate a new ensemble with the same number of configuration at the end of the minimization and check how the final optimized solution changes with this new ensemble. This is also a good way to test if the solution is actually converged to the correct solution. The magnitude of the changes in the auxiliary dynamical matrix’s frequencies and structure is an accurate estimation on the stochastic error.

You can always split the ensemble in two and run two minimizations with the two half of the ensembe to get a hint on the error on the structure or on the dynamical matrix.
To split the ensemble, refer to the *FAQ* about the error on the hessian matrix.

<a name="How-does-the-error-over-the-gradients-scale-with-the-number-of-configurations?"></a>
# How does the error over the gradients scale with the number of configurations?

The error scales as any stochastic method, with the inverse of the square root of the number of configurations. So to double the accuracy, the number of configurations must be multiplied by 4.

<a name="When-I-relax-the-cell,-is-it-necessary-for-the-gradients-to-reach-zero-before-making-a-step-with-the-new-cell?"></a>
# When I relax the cell, is it necessary for the gradients to reach zero before making a step with the new cell?

In general it is good to have a reasonable dynamical matrix before starting with a variable cell relaxation. The best strategy is to perform a fixed cell relaxation with few configurations until you are close to the final solution (the gradients are comparable with their errors). Then you can start a variable cell relaxation and submit new populations in the suggested new cell even if the previous one was not perfectly converged.

<a name="I-cannot-remove-the-pressure-anisotropy-after-relaxing-the-cell,-what-is-happening?"></a>
# I cannot remove the pressure anisotropy after relaxing the cell, what is happening?

Variable cell calculation is a tricky algorithm. It could be that your bulk modulus is strongly anisotropic, so the algorithm has difficulties in optimizing well. In general the stress tensor is also affected by stochastic error, so it is impossible to completely remove anisotropy. However, a converged result is one in which the residual anisotropy in the stress tensor is comparable to the stochastic error on the stress tensor. If you are not able to converge, you can either increase the number of configurations, modify the bulk_modulus parameter (increase it if the stress change too much between two populations, decrease it if it does not changes enough) or fix the overall volume (by using the fix_volume flag in the &relax namespace or in the vc_relax method if you are using the python script).
    Fixing the volume can improve the convergence of the variable cell algorithm.

<a name="How-may-I-run-a-calculation-neglecting-symmetries?"></a>
# How may I run a calculation neglecting symmetries?

You can tell the code to neglect symmetries with the `neglect_symmetries = .true.` flag.
In the python script this is done setting the attribute *neglect_symmetries* of sscha.SchaMinimizer.SSCHA_Minimizer to False.

<a name="In-which-units-are-the-lattice-vectors,-the-atomic-positions,-and-the-mass-of-the-atoms-in-the-dynamical-matrix-file?"></a>
# In which units are the lattice vectors, the atomic positions, and the mass of the atoms in the dynamical matrix file?

The dynamical matrix follows the Quantum ESPRESSO units. They are Rydberg atomic units (unit of mass is 1/2  the electron mass, energy is Ry, positions are in Bohr). However, Quantum ESPRESSO may have an *ibrav* not equal to zero (the third number in the header of the dynamical matrix). In this case, please, refer to the espresso *ibrav* guide in the PW.x input description: <https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm199>

<a name="What-is-the-difference-between-the-different-kind-of-minimization-(preconditioning-and-root_representation)?"></a>
# What is the difference between the different kind of minimizations (preconditioning and root_representation)?

 The target of a SSCHA minimization is to find the ionic density matrix \\(\rho(\Phi, \vec {\mathcal R})\\) that minimizes the total
    free energy. It may happen, if we are using a too big step for the auxiliar dynamical matrix \\(\Phi\\) that it becomes not positive definite.
    This may be due to the stochastic noise during the minimization.
  In order to avoid this, you may set **root_representation** to either **sqrt** or **root4** (inside &inputscha namespace or the SSCHA_Minimizer object)
    In this way, instead of minimizing the \\(\Phi\\) matrix, we minimize with respect to \\(\sqrt{\Phi}\\) or \\(\sqrt[4]{\Phi}\\).
    Therefore the new dynamical matrix are constrained in a space that is positive definite. Moreover, it has been proved that \\(\sqrt[4]{\Phi}\\)
    minimization is better conditioned than the original, and thus should reach the minimum faster.

Alternatively, a similar effect to the speedup in the minimization obtained with **root4** is possible using preconditioning (by setting **preconditioning** or **precond_dyn** to True in the input file or the python script, respectively). This way also the single minimization step runs faster, as it avoids passing in the root space of the dynamical matrix (but indeed, you can have imaginary frequencies).

Since the gradient computation is much slower (especially for system with more than 80 atoms in the supercell) without  preconditioning,
    it is possible to combine the preconditioning with the root representation to have a faster gradient computation and to guarantee that
    the dynamical matrix is positive definite by construction at each step.
    However, in this way the good condition number obtained by the preconditioning (or the root4 representation) is spoiled. For this reason, when using the preconditioning, avoid using **root4**, and chose instead **sqrt** as root_representation.

The default values are:
```
&inputscha
        root_representation = "normal"
        preconditioning = .true.
&end
```
or in python     
```
# The ensemble has been loaded as ens
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ens)
minim.root_representation = "normal"
minim.precond_dyn = True
```

<a name="How-do-I-lock-modes-from-m-to-n-in-the-minimization?"></a>
# How do I lock modes from m to n in the minimization?

Constrains to the minimization within the mode space may be added both in the input script and directly by the python version.
    
In the input script, inside the namespace **&utils**, you should add, for instance, **mu_free_start = 30** and **mu_free_end = 36** to optimize only between mode 30 and 36 (for each q point). You can also use the keywords **mu_lock_start** and **mu_lock_end** to freeze only a subset of modes. You can also choose if you want to freeze only the dynamical matrix or also the structure relaxation along those directions, by picking: **project_dyn = .true.** and **project_structure = .false.**. In this way, you freeze only the dynamical matrix along the specified modes, but not the structure.

Modes may be also locked within the python scripting. Look at the LockModes example in the Examples directory.

<a name="How-do-I-lock-a-special-atom-in-the-minimization?"></a>
# How do I lock a special atom in the minimization?

More complex constrains than mode locking may be activated in the minimization, but their use is limited within the python scripting.
    You can write your own constraining function that will be applied to the structure gradient or to the dynamical matrix gradient.
    This function should take as input the two gradients (dynamical matrix and structure) and operate directly on them.
    Then it can be passed to the minimization engine as *custom_function_gradient*.

```
LIST_OF_ATOMS_TO_FIX = [0, 2, 3]
def fix_atoms(gradient_dyn, gradient_struct):
# Fix the atoms in the list
gradient_struct[LIST_OF_ATOMS_TO_FIX, :] = 0
minim.run( custom_function_gradient = fix_atoms )
```

Here, `minim` is the `SSCHA_Minimizer` class. In this case we only fix the structure gradient. However, in this way the overall gradient will have a translation (acoustic sum rule is violated). Be very careful when doing this kind of constrains, and check if it is really what you want.

A more detailed and working example that fixes also the degrees of freedom of the dynamical matrix is reported in the FixAtoms example.

<a name="How-do-I-understand-if-I-have-to-generate-a-new-population-or-the-minimization-converged?"></a>
# How do I understand if I have to generate a new population or the minimization converged?

In general, if the code stops because the gradient is much below the error (less then 1%), then it is converged (with a Kong-Liu threshold ratio of at least 0.5). If the code ends the minimization because it went outside the stochastic criteria, a new population is required.
    There are cases in which you use too few configurations to reach a small gradient before wasting the ensemble. If this is the case, print the frequencies during the minimizations (using the &utils card with `save_freq_filename` attribute). You may compare subsequent minimizations, if the frequencies are randomly moving between different minimization (and you cannot identify a trend in none of them), then you reach the limit of accuracy of the ensemble. Frequencies are a much better parameter to control for convergence than free energy, as the free energy close to the minimum is quadratic.

<a name="How-do-I-choose-the-appropriate-value-of-Kong-Liu-effective-sample-size-o-ratio?"></a>
# How do I choose the appropriate value of Kong-Liu effective sample size or ratio?

The Kong-Liu (KL) effective sample size is an estimation on how good is the extracted set of configurations to describe the BO landscape around the current values of dynamical matrix and the centroid position. After the ensemble is generated, the KL sample size matches with the actual number of configurations. However, as the minimization goes on, the KL sample size is reduced. The code stops when the KL sample size is below a certain threshold.

The default value of Kong-Liu threshold ratio is 0.5 (effective sample size = 0.5 the original number of configurations). This is a good and safe value for most situations. However, if you are very far from the minimum and the gradient is big, you can trust it even if it is very noisy. For this reason you can lower the Kong-Liu ratio to 0.2 or 0.1. However, notice that by construction the KL effective sample size is always bigger than 2. For this reason if you use 10 configurations, and you set a threshold ratio below 0.2, you will never reach the threshold, and your minimization will continue forever (going into a very bad regime where you are minimizing something that is completely random). On the other side, on some very complex system close to the minimum, it could be safe to increase the KL ratio even at 0.6.

<a name="How-do-I-understand-if-the-free-energy-hessian-calculation-is-converged?"></a>
# How do I understand if the free energy hessian calculation is converged?


The free energy hessian requires normally many more configurations than the SCHA minimization itself. First of all, to run the free energy Hessian, the SSCHA minimization must end with a gradient that can be decreased indefinitely without decreasing the KL below 0.7 /0.8.
    Then you can estimate the error by repeating the hessian calculation with half of the ensemble and check how the frequencies of the hessian changes. This is also a good check for the final error on the frequencies.

You can split your ensemble in two by using the split function.
```
    import sscha, sscha.Ensemble

    # Load the dynamical matrix as dyn
    # [...]

    # ens is the Ensemble() class correctly initialized.
    # We can for example load it
    # Assuming it is stored in 'data_dir' with population 1 and 1000 configurations
    # We assume to have loaded the original dynamical matrix dyn and to know the generating temperature T
    ens = sscha.Ensemble.Ensemble(dyn, T, dyn.GetSupercell())
    ens.load("data_dir", population = 1, N = 1000)

    # We create a mask that selects which configurations to take
    first_half_mask = np.zeros(ens.N, dtype = bool)
    first_half_mask[: ens.N // 2] = True

    # We create also the mask for the second half
    # by taking the not operation on the first_half_mask
    second_half_mask = ~first_half_mask

    # Now we split the ensemble
    ens_first_half = ens.split(first_half_mask)
    ens_second_half = ens.split(second_half_mask)

    # We can save the two half ensembles as population 2 and 3.
    ens_first_half.save("data_dir", population = 2)
    ens_second_half.save("data_dir", population = 3)
```
This simple script will generate two ensembles inside `data_dir` directory with population 2 and 3, each one containing the first
    and the second half of the ensemble with population 1 respectively. You can perform then your calculation of the free energy hessian
    with both the ensemble to estimate the error on the frequencies and the polarization vectors.

<a name="How-can-I-add-more-configurations-to-an-existing-ensembe?"></a>
# How can I add more configurations to an existing ensembe?

You can use the split and merge functions of the Ensemble class.
    First of all you generate a new ensemble, you compute the energy and force for that ensemble,
    then you merge it inside another one.
```
    # Load the original ensemble (first population with 1000 configurations)
    ens = sscha.Ensemble.Ensemble(dynmat, T, dynmat.GetSupercell())
    ens.load("data_dir", population = 1, N = 1000)

    # Generate a new ensemble with other 1000 configurations
    new_ensemble = sscha.Ensemble.Ensemble(dynmat, T, dynmat.GetSupercell())
    new_ensemble.generate(1000)

    # Compute the energy and forces for the new ensemble
    # For example in this case we assume to have initialized 'calc' as an ASE calculator.
    # But you can also save it with a different population,
    # manually compute energy and forces, and then load again the ensemble.
    new_ensemble.get_energy_forces(calc)

    # Merge the two ensembles
    ens.merge(new_ensemble)

    # Now ens contains the two ensembles. You can save it or directly use it for a SSCHA calculation
    ens.save("data_dir", population = 2)
```

Indeed, to avoid mistakes, when merging the ensemble you must be careful that the dynamical matrix and the temperature
    used to generate both ensembles are the same.

<a name="How-do-I-fix-the-random-number-generator-seed-to-make-a-calculation-reproducible?"></a>
# How do I fix the random number generator seed to make a calculation reproducible?

As for version 1.0, this can be achieved only by using the python script.
    Since python uses numpy as random number generator, you can, at the beginning of the script that generates the ensemble, use the following:
```
    import numpy as np

    X = 0
    np.random.seed(seed = X)
```
where `X` is the integer used as a seed. By default, if not specified, it is initialized with None that it is equivalent of initializing with the current local time.


<a name="symph-not-found"></a>
# When I try to start a sscha calculation or import cellconstructor I get ModuleNotFoundError: No module named 'symph'


This error is spotting that the symph module is missing. Symph is an extention built and installed toghether with cellconstructor. Try to reinstall cellconstructor and look for errors in the output. Note that symph is not compiled correctly if you do not have BLAS and LAPACK libraries. Follow the installation guide for further details on how to properly install CellConstructor and python-sscha. 


