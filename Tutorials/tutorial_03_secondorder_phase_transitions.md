::: {.document}
::: {.documentwrapper}
::: {.bodywrapper}
::: {.body role="main"}
::: {#hands-on-session-3-calculations-of-second-order-phase-transitions-with-the-sscha .section}
Hands-on-session 3 - Calculations of second-order phase transitions with the SSCHA[¶](#hands-on-session-3-calculations-of-second-order-phase-transitions-with-the-sscha "Permalink to this headline"){.headerlink}
==================================================================================================================================================================================================================

In this hands-on, we learn how to calculate second-order phase
transitions within the SSCHA.

::: {#structural-instability-calculation-of-the-hessian .section}
Structural instability: calculation of the Hessian[¶](#structural-instability-calculation-of-the-hessian "Permalink to this headline"){.headerlink}
---------------------------------------------------------------------------------------------------------------------------------------------------

According to Landau's theory, a second-order phase transition occurs
when the free energy curvature around the high-symmetry structure on the
direction of the order parameter becomes negative:

![[Fig. 7 ]{.caption-number}[Landau's theory of second-order phase
transitions.]{.caption-text}[¶](#id3 "Permalink to this image"){.headerlink}](_images/second_order.png)

For structural *displacive* phase transitions, the order parameter is
associated to phonon atomic displacements:

::: {.math .notranslate .nohighlight}
\\\[\\frac{\\partial\^2 F}{\\partial R_a \\partial R_b}\\\]
:::

Thus, the Free energy Hessian is the central quantity to study
second-order phase transitions. The SSCHA provides an analytical
equation for the free energy Hessian, derived by Raffaello Bianco in the
work [Bianco et. al. Phys. Rev. B 96,
014111](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.014111){.reference
.external}. The free energy curvature can be written as:

::: {.math .notranslate .nohighlight}
\\\[\\frac{\\partial\^2 F}{\\partial {R_a}\\partial {R_b}} = \\Phi\_{ab}
+ \\sum\_{cdef} \\stackrel{(3)}{\\Phi}\_{acd}\[1 -
\\Lambda\\stackrel{(4)}{\\Phi}\]\^{-1}\_{cdef}
\\stackrel{(3)}{\\Phi}\_{efb}\\\]
:::

Fortunately, this complex equation can be evaluated from the ensemble
with a simple function call:

::: {.highlight-python .notranslate}
::: {.highlight}
    ensemble.get_free_energy_hessian()
:::
:::

Lets see a practical example, first we calculate the SSCHA dynamical
matrix for the SnTe:

To speedup the calculations, we will use a force-field that can mimic
the physics of ferroelectric transitions in FCC lattices.

We begin importing some libraries:

::: {.highlight-python .notranslate}
::: {.highlight}
    #!/usr/bin/env python
    # -*- coding: utf-8 -*-
    #
    #  SSCHA_exercise_Calculus.py
    #
    # Import the cellconstructor stuff
    import cellconstructor as CC
    import cellconstructor.Phonons
    import cellconstructor.ForceTensor
    import cellconstructor.Structure
    import cellconstructor.Spectral

    # Import the modules of the force field
    import fforces as ff
    import fforces.Calculator

    # Import the modules to run the sscha
    import sscha, sscha.Ensemble, sscha.SchaMinimizer
    import sscha.Relax, sscha.Utilities

    import spglib
    from ase.visualize import view

    # Import Matplotlib to plot
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import timeit
:::
:::

Next we set some variables for the calculation:

::: {.highlight-python .notranslate}
::: {.highlight}
    #Setting the variables:
    #Setting the temperature in Kelvin:
    Temperature = 0
    #Setting the number of configurations:
    configurations = 50
    #Setting the names and location of the files:
    Files_dyn_SnTe = "ffield_dynq"
    #Set the number of irreducible q (reated to the supercell size):
    nqirr = 3
    #Setting the frequencies output file:
    File_frequencies = "frequencies.dat"
    #Setting the dynamical matrix output filename:
    File_final_dyn = "final_sscha_T{}_".format(int(Temperature))
:::
:::

Now we need to calculate the SSCHA dynamical matrix. For that we follow
some steps:

1.  First we prepare the Toy model force field that substitutes the
    usual *ab-initio* for this tutorial. This force field needs the
    harmonic dynamical matrix to be initialized, and the higher order
    parameters. Finally, the dynamical matrix for the minimization is
    loaded and readied. Since we are studying a system that has a
    spontaneous symmetry breaking at low temperature, the harmonic
    dynamical matrices will have imaginary phonons. We must enforce
    phonons to be positive definite to start a SSCHA minimization.

::: {.highlight-python .notranslate}
::: {.highlight}
    # Load the dynamical matrix for the force field
    ff_dyn = CC.Phonons.Phonons("ffield_dynq", 3)

    # Setup the forcefield with the correct parameters
    ff_calculator = ff.Calculator.ToyModelCalculator(ff_dyn)
    ff_calculator.type_cal = "pbtex"
    ff_calculator.p3 = 0.036475
    ff_calculator.p4 = -0.022
    ff_calculator.p4x = -0.014

    # Initialization of the SSCHA matrix
    dyn_sscha = CC.Phonons.Phonons(Files_dyn_SnTe, nqirr)
    # Flip the imaginary frequencies into real ones
    dyn_sscha.ForcePositiveDefinite()
    # Apply the ASR and the symmetry group
    dyn_sscha.Symmetrize()
:::
:::

2.  The next step is to create the ensembles for the specified
    temperature. As an extra, we also look for the space group of the
    structure.

::: {.highlight-python .notranslate}
::: {.highlight}
    ensemble = sscha.Ensemble.Ensemble(dyn_sscha,
            T0 = Temperature, supercell = dyn_sscha.GetSupercell())
    # Detect space group
    symm=spglib.get_spacegroup(dyn_sscha.structure.get_ase_atoms(),
            0.005)
    print('Initial SG = ', symm)
:::
:::

3.  Next comes the minimization step. Here we can set the fourth root
    minimization, in which, instead of optimizing the auxiliary
    dynamical matrices themselves, we will optimize their fourth root.

    ::: {.math .notranslate .nohighlight}
    \\\[\\Phi = \\left({\\sqrt\[4\]{\\Phi}}\\right)\^4\\\]
    :::

    This constrains the dynamical matrix to be positive definite during
    the minimization. Next the automatic relaxation is set with the
    option here to use the Sobol sequence for the ensemble generation.

    We also set a custom function to save the frequencies at each
    iteration, to see how they evolves. This is very useful to
    understand if the algorithm is converged or not.

    Then the dynamical matrix of the converged minimization is saved in
    a file, and finally we take a look at the space group and the
    structure.

::: {.highlight-python .notranslate}
::: {.highlight}
    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

    # Now we setup the minimization parameters
    # Since we are quite far from the correct solution,
    # we will use a small optimization step
    minim.set_minimization_step(0.25)

    # Reduce the threshold for the gradient convergence
    minim.meaningful_factor = 0.01

    # If the minimization ends with few steps (less than 10),
    # decrease it, if it takes too much, increase it

    # We decrease the Kong-Liu effective sample size below
    # which the population is stopped
    minim.kong_liu_ratio = 0.5 # Default 0.5
    # We relax the structure
    relax = sscha.Relax.SSCHA(minim,
                      ase_calculator = ff_calculator,
                      N_configs = configurations,
                      max_pop = 50)

    # Setup the custom function to print the frequencies
    # at each step of the minimization
    io_func = sscha.Utilities.IOInfo()
    io_func.SetupSaving(File_frequencies)
    # The file that will contain the frequencies is frequencies.dat

    # Now tell relax to call the function to save the frequencies
    # after each iteration
    # CFP stands for Custom Function Post (Post = after the minimization step)
    #relax.setup_custom_functions(custom_function_post = io_func.CFP_SaveFrequencies)
    relax.setup_custom_functions(custom_function_post = io_func.CFP_SaveAll)
    # Finally we do all the free energy calculations.
    relax.relax()
    #relax.vc_relax(static_bulk_modulus=40, fix_volume = False)

    # Save the final dynamical matrix
    relax.minim.dyn.save_qe(File_final_dyn)
    # Detect space group
    symm=spglib.get_spacegroup(relax.minim.dyn.structure.get_ase_atoms(),
                0.005)
    print('New SG = ', symm)
    view(relax.minim.dyn.structure.get_ase_atoms())
:::
:::

This code will calculate the SSCHA minimization with the
*ff_calculator*. We cat use **sscha-plot-data.py** to take a look at the
minimization.

::: {.highlight-bash .notranslate}
::: {.highlight}
    python sscha-plot-data.py frequencies.dat
:::
:::

![](_images/Figure_1_N.png)

![](_images/Figure_2_N.png)

Note: this force field model is not able to compute stress, as it is
defined only at fixed volume, so we cannot use it for a variable cell
relaxation.

**Now we can search for instabilities.**

If we have a very small mode in the SSCHA frequencies, it means that
associated to that mode we have huge fluctuations. This can indicate an
instability. However, to test this we need to compute the free energy
curvature along this mode. This can be obtained in one shot thanks to
the theory developed in [Bianco et. al. Phys. Rev. B 96,
014111.](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.014111){.reference
.external}

For that we create another program to do the job.

As before, we begin importing some libraries and setting variables:

::: {.highlight-python .notranslate}
::: {.highlight}
    #!/usr/bin/env python
    # -*- coding: utf-8 -*-
    #
    #  SSCHA_exercise_Unstable.py
    #
    # Import the cellconstructor stuff
    import cellconstructor as CC
    import cellconstructor.Phonons
    import cellconstructor.ForceTensor
    import cellconstructor.Structure
    import cellconstructor.Spectral

    # Import the modules of the force field
    import fforces as ff
    import fforces.Calculator

    # Import the modules to run the sscha
    import sscha, sscha.Ensemble, sscha.SchaMinimizer
    import sscha.Relax, sscha.Utilities

    import spglib
    from ase.visualize import view

    # Import Matplotlib to plot
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import timeit

    #Setting the variables:
    #Setting the temperature in Kelvin:
    Temperature = 0
    #Setting the number of configurations:
    configurations = 50
    #Setting the names and location of the files:
    Files_dyn_SnTe = "ffield_dynq"
    #Set the number of irreducible q (reated to the supercell size):
    nqirr = 3
    #Setting the frequencies output file:
    File_frequencies = "frequencies.dat"
    #Setting the dynamical matrix output filename:
    File_final_dyn = "final_sscha_T{}_".format(int(Temperature))
:::
:::

Now we look for that instability:

1.  The *ff_calculator* toy potential is defined as we have seen in the
    previous program.

::: {.highlight-python .notranslate}
::: {.highlight}
    # Load the dynamical matrix for the force field
    ff_dyn = CC.Phonons.Phonons("ffield_dynq", 3)

    # Setup the forcefield with the correct parameters
    ff_calculator = ff.Calculator.ToyModelCalculator(ff_dyn)
    ff_calculator.type_cal = "pbtex"
    ff_calculator.p3 = 0.036475
    ff_calculator.p4 = -0.022
    ff_calculator.p4x = -0.014

    # Initialization of the SSCHA matrix
    dyn_sscha = CC.Phonons.Phonons(Files_dyn_SnTe, nqirr)
    dyn_sscha.ForcePositiveDefinite()

    # Apply also the ASR and the symmetry group
    dyn_sscha.Symmetrize()
:::
:::

2.  Next, we will load the dynamical matrix calculated previously with
    the *ff_calculator* toy potential, so there is no need to calculate
    it again.

::: {.highlight-python .notranslate}
::: {.highlight}
    # The SSCHA dynamical matrix is needed (the one after convergence)
    # We reload the final result (no need to rerun the sscha minimization)
    dyn_sscha_final = CC.Phonons.Phonons(File_final_dyn, nqirr)
:::
:::

3.  Then, as the Hessian calculation is more sensible, we generate a new
    ensemble with more configurations. To compute the hessian we will
    use an ensemble of 10000 configurations. Note here that we can use
    less if we use Sobol sequence or we can load a previously generated
    ensemble.

::: {.highlight-python .notranslate}
::: {.highlight}
    # We reset the ensemble
    ensemble = sscha.Ensemble.Ensemble(dyn_sscha_final, T0 = Temperature,
                        supercell = dyn_sscha_final.GetSupercell())

    # We need a bigger ensemble to properly compute the hessian
    # Here we will use 10000 configurations
    ensemble.generate(5000, sobol = True)
    #ensemble.generate(10000, sobol = False)
    #We could also load the ensemble with
    # ensemble.load("data_ensemble_final", N = 100, population = 5)
:::
:::

4.  We now compute forces and energies using the force field calculator.

::: {.highlight-python .notranslate}
::: {.highlight}
    # We now compute forces and energies using the force field calculator
    ensemble.get_energy_forces(ff_calculator, compute_stress = False)
:::
:::

5.  Finally the free energy hessian is calculated in the *hessian*
    function. We can choose if we neglect or not in the calculation the
    four phonon scattering process. Four phonon scattering processes
    require a huge memory allocation for big systems, that scales as
    (3N)\^4 with N the number of atoms in the supercell. Moreover, it
    may require also more configurations to converge.

    In almost all the systems we studied up to now, we found this four
    phonon scattering at high order to be negligible. We remark, that
    the SSCHA minimization already includes four phonon scattering at
    the lowest order perturbation theory, thus neglecting this term only
    affects combinations of one or more four phonon scattering with two
    three phonon scatterings (high order diagrams). For more details,
    see [Bianco et. al. Phys. Rev. B
    96, 014111.](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.014111){.reference
    .external}

    We can then print the frequencies of the hessian. If an imaginary
    frequency is present, then the system wants to spontaneously break
    the high symmetry phase.

::: {.highlight-python .notranslate}
::: {.highlight}
    print("Updating the importance sampling...")
    # If the sscha matrix was not the one used to compute the ensemble
    # We must update the ensemble weights
    # We can also use this function to simulate a different temperature.
    ensemble.update_weights(dyn_sscha_final, Temperature)
    # ----------- COMPUTE THE FREE ENERGY HESSIAN -----------
    print("Computing the free energy hessian...")
    dyn_hessian = ensemble.get_free_energy_hessian(include_v4 = False)
    # We neglect high-order four phonon scattering
    #dyn_hessian = ensemble.get_free_energy_hessian(include_v4 = True,
    #              get_full_hessian = True,verbose = True) # Full calculus
    # We can save the free energy hessian as a dynamical matrix
    # in quantum espresso format
    dyn_hessian.save_qe("hessian")
    # -------------------------------------------------------
    # We calculate the frequencies of the hessian:
    w_hessian, pols_hessian = dyn_hessian.DiagonalizeSupercell()

    # Print all the frequency converting them into cm-1 (They are in Ry)
    print("\n".join(["{:16.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in w_hessian]))
:::
:::

The frequencies in the free energy hessian are temperature dependent.

We can look at the eigenmodes of the free energy hessian to check if we
have imaginary phonons. If there are negative frequencies then we found
an instability. You can check what happens if you include the fourth
order.

::: {.topic}
Exercise

The Sobol sequences reduces the number of configurations by doing a
better mapping of the gaussian than a random distribution. By uniformity
spreading the samplings with a low discrepancy sequence like Sobol it is
possible to reduce the number of configurations needed. Low discrepancy
sequences tend to sample space "more uniformly" than random numbers.
Algorithms that use such sequences may have superior convergence. You
can test this in the calculation of the hessian by changing the number
of configurations and the mapping scheme in the *ensemble.generate()*
function.
:::
:::

::: {#second-order-phase-transition .section}
Second order phase transition[¶](#second-order-phase-transition "Permalink to this headline"){.headerlink}
----------------------------------------------------------------------------------------------------------

Up to now we studied the system at T=0K and we found that there is an
instability. However, we can repeat the minimization at many
temperatures, and track the phonon frequency to see which is the
temperature at which the system becomes stable.

Again we load and set the variables. Now the we have several
temperatures so we store them in an array:

::: {.highlight-python .notranslate}
::: {.highlight}
    #!/usr/bin/env python
    # -*- coding: utf-8 -*-
    #
    #  SSCHA_exercise_Unstable.py
    #
    # Import the cellconstructor stuff
    import cellconstructor as CC
    import cellconstructor.Phonons
    import cellconstructor.ForceTensor
    import cellconstructor.Structure
    import cellconstructor.Spectral

    # Import the modules of the force field
    import fforces as ff
    import fforces.Calculator

    # Import the modules to run the sscha
    import sscha, sscha.Ensemble, sscha.SchaMinimizer
    import sscha.Relax, sscha.Utilities

    import spglib
    from ase.visualize import view

    # Import Matplotlib to plot
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import timeit

    #Setting the variables:
    #Setting the temperature in Kelvin:
    Temperature = 0
    #Setting the number of configurations:
    configurations = 50
    #Setting the names and location of the files:
    Files_dyn_SnTe = "ffield_dynq"
    #Set the number of irreducible q (reated to the supercell size):
    nqirr = 3
    #Setting the frequencies output file:
    File_frequencies = "frequencies.dat"
    #Setting the dynamical matrix output filename:
    File_final_dyn = "final_sscha_T{}_".format(int(Temperature))
    sobol = False
    sobol_scatter = False
:::
:::

1.  Like in the previous program, first we prepare the Toy model force
    field

::: {.highlight-python .notranslate}
::: {.highlight}
    # Load the dynamical matrix for the force field
    ff_dyn = CC.Phonons.Phonons("ffield_dynq", 3)

    # Setup the forcefield with the correct parameters
    ff_calculator = ff.Calculator.ToyModelCalculator(ff_dyn)
    ff_calculator.type_cal = "pbtex"
    ff_calculator.p3 = 0.036475
    ff_calculator.p4 = -0.022
    ff_calculator.p4x = -0.014
:::
:::

2.  We are going to need a range of temperatures for this calculation:

::: {.highlight-python .notranslate}
::: {.highlight}
    # Define the temperatures, from 50 to 300 K, 6 temperatures
    temperatures = np.linspace(50, 300, 6)

    lowest_hessian_mode = []
    lowest_sscha_mode = []

    # Perform a simulation at each temperature
    t_old = Temperature
:::
:::

3.  In the next part we condense the calculation of the hessians in a
    loop for different temperatures. In the end, it searches for the
    lowest non acoustic frequency to save with the correspondent
    auxiliar sscha frequency.

::: {.highlight-python .notranslate}
::: {.highlight}
    for Temperature in temperatures:
        # Load the starting dynamical matrix
        dyn = CC.Phonons.Phonons(File_final_dyn.format(int(t_old)), nqirr)

        # Prepare the ensemble
        ensemble = sscha.Ensemble.Ensemble(dyn, T0 = Temperature,
                                      supercell = dyn.GetSupercell())

        # Prepare the minimizer
        minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
        minim.min_step_struc = 0.05
        minim.min_step_dyn = 0.002
        minim.kong_liu_ratio = 0.5
        minim.meaningful_factor = 0.000001
        #minim.root_representation = "root4"
        #minim.precond_dyn = False
        #minim.minim_struct = True
        #minim.neglect_symmetries = True
        minim.enforce_sum_rule = True  # Lorenzo's solution to the error

        # Prepare the relaxer (through many population)
        relax = sscha.Relax.SSCHA(minim, ase_calculator = ff_calculator,
                      N_configs=configurations, max_pop=20)

        # Relax
        relax.relax(sobol = sobol, sobol_scramble = sobol_scatter)
        #relax.relax()

        # Save the dynamical matrix
        relax.minim.dyn.save_qe(File_final_dyn.format(int(Temperature)))

        # Detect space group
        symm=spglib.get_spacegroup(relax.minim.dyn.structure.get_ase_atoms(),
                                          0.005)
        print('Current SG = ', symm,' at T=',int(Temperature))

        # Recompute the ensemble for the hessian calculation
        ensemble = sscha.Ensemble.Ensemble(relax.minim.dyn, T0 = Temperature,
                                supercell = dyn.GetSupercell())
        ensemble.generate(configurations, sobol = sobol,
                                          sobol_scramble = sobol_scatter)
        ensemble.get_energy_forces(ff_calculator, compute_stress = False)
        #gets the energies and forces from ff_calculator

        #update weights!!!
        ensemble.update_weights(relax.minim.dyn, Temperature)
        # Get the free energy hessian
        dyn_hessian = ensemble.get_free_energy_hessian(include_v4 = False)
        #free energy hessian as in Bianco paper 2017
        dyn_hessian.save_qe("hessian_T{}_".format(int(Temperature)))

        # Get the lowest frequencies for the sscha and the free energy hessian
        w_sscha, pols_sscha = relax.minim.dyn.DiagonalizeSupercell() #dynamical matrix
        # Get the structure in the supercell
        superstructure = relax.minim.dyn.structure.generate_supercell(relax.minim.dyn.GetSupercell())

        # Discard the acoustic modes
        acoustic_modes = CC.Methods.get_translations(pols_sscha,
                                      superstructure.get_masses_array())
        w_sscha = w_sscha[~acoustic_modes]

        lowest_sscha_mode.append(np.min(w_sscha) * CC.Units.RY_TO_CM) # Convert from Ry to cm-1

        w_hessian, pols_hessian = dyn_hessian.DiagonalizeSupercell() #recomputed dyn for hessian
        # Discard the acoustic modes
        acoustic_modes = CC.Methods.get_translations(pols_hessian,
                                        superstructure.get_masses_array())
        w_hessian = w_hessian[~acoustic_modes]
        lowest_hessian_mode.append(np.min(w_hessian) * CC.Units.RY_TO_CM) # Convert from Ry to cm-1
        #print ("\n".join(["{:.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in pols_hessian]))
        #exit()

        t_old = Temperature
    # We prepare now the file to save the results
    freq_data = np.zeros( (len(temperatures), 3))
    freq_data[:, 0] = temperatures
    freq_data[:, 1] = lowest_sscha_mode
    freq_data[:, 2] = lowest_hessian_mode

    # Save results on file
    np.savetxt("{}_hessian_vs_temperature.dat".format(configurations),
                freq_data, header = "T [K]; SSCHA mode [cm-1]; Free energy hessian [cm-1]")
:::
:::

4.  Finally we make a graphic output of the data.

::: {.highlight-python .notranslate}
::: {.highlight}
    hessian_data = np.loadtxt("{}_hessian_vs_temperature.dat".format(configurations))

    plt.figure(dpi = 120)
    plt.plot(hessian_data[:,0], hessian_data[:,1],
                              label = "Min SCHA freq", marker = ">")
    plt.plot(hessian_data[:,0], hessian_data[:,2],
                              label = "Free energy curvature", marker = "o")
    plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
    plt.xlabel("Temperature [K]")
    plt.ylabel("Frequency [cm-1]")
    plt.legend()
    plt.tight_layout()
    plt.savefig('{}_Temp_Freq.png'.format(configurations))
    #plt.show()

    plt.figure(dpi = 120)
    plt.plot(hessian_data[:,0], np.sign(hessian_data[:,2]) * hessian_data[:,2]**2,
                        label = "Free energy curvature", marker = "o")
    plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
    plt.xlabel("Temperature [K]")
    plt.ylabel("$\omega^2$ [cm-2]")
    plt.legend()
    plt.tight_layout()
    plt.savefig('{}_Temp_Omeg.png'.format(configurations))
    #plt.show()
:::
:::

We will simulate the temperatures up to room temperature (300 K) with
steps of 50 K. Note, this will perform all the steps above 6 times, so
it may take some minutes, depending on the PC (on a i3 from 2015, with
one core, it took 2 hours). If it takes too long you can reduce the
number of steps by changing the temperature array in *Temperature_i =
np.linspace(50, 300, 6)*.

![[Fig. 8 ]{.caption-number}[Frequencies versus
Temperatures]{.caption-text}[¶](#id4 "Permalink to this image"){.headerlink}](_images/5000_Temp_Freq.png)

In [[Frequencies versus Temperatures]{.std
.std-ref}](#fig-results1){.reference .internal} we can see that the
phase transition is between 100K and 150K. We see that the data points
do not drawn a linear figure. We can increase the number of Temperature
points to locate the exact transition temperature, but there is another
better way to find it.

![[Fig. 9 ]{.caption-number}[squared Frequencies versus
Temperatures.]{.caption-text}[¶](#id5 "Permalink to this image"){.headerlink}](_images/5000_Temp_Omeg.png)

For the Landau theory of phase transition, since the SSCHA is a
mean-field approach, we expect that around the transition the critical
exponent of the temperature goes as

::: {.math .notranslate .nohighlight}
\\\[\\omega \\sim \\sqrt{\\Phi}\\\]
:::

For this reason is better to plot the temperature versus the square of
the frequency as in [[squared Frequencies versus Temperatures.]{.std
.std-ref}](#fig-results2){.reference .internal} This makes the graph
lineal and so we can easily estimate the critic temperature by linear
interpolation.

We are using only 50 configurations in the ensemble. Note that this
makes a fast calculation but is a low number for this calculations
because the free energy calculations are more noisy than the SSCHA
frequencies. This is due to the fact that the computation of the free
energy requires the third order force constant tensor, and that requires
more configurations to converge.

::: {.topic}
Exercise

How the calculation of the free energy changes with the number of
configurations?
:::

![[Fig. 10 ]{.caption-number}[Evolution of the lowest *soft* frequency
in relation to the number of configurations in the ensemble with a
stable configuration. The line is the media and the shade is the
standard
deviation.]{.caption-text}[¶](#id6 "Permalink to this image"){.headerlink}](_images/Conf_Freq.png)

::: {.topic}
Exercise

Plot the Hessian phonon dispersion

![](_images/dispersion.png)
:::

![[Fig. 11 ]{.caption-number}[Workflow of the SSCHA objects for: A) Free
energy minimization; B) Structural instabilities search; C) Temperature
loop for second order phase transition code. Dotted lines are functions
within objects. The dotted and dashed lines indicate the relationship of
the dynamic matrix to the
objects.]{.caption-text}[¶](#id7 "Permalink to this image"){.headerlink}](_images/Hands-on3.png)
:::
:::
:::
:::
:::

::: {.sphinxsidebar role="navigation" aria-label="main navigation"}
::: {.sphinxsidebarwrapper}
::: {.relations}
### Related Topics

-   [Documentation overview](index.html)
    -   Previous: [Hands-on-session 2 - Advanced free energy
        minimization](tutorial_02_advanced_submission.html "previous chapter")
    -   Next: [Hands-on-session 4 - Calculation of spectral properties
        with the Self Consistent Harmonic
        Approximation](tutorial_04_calculation_spectrum.html "next chapter")
:::

::: {#searchbox style="display: none" role="search"}
### Quick search {#searchlabel}

::: {.searchformwrapper}
:::
:::
:::
:::

::: {.clearer}
:::
:::

::: {.footer}
©2023, Lorenzo Monacelli. \| Powered by [Sphinx
4.2.0](http://sphinx-doc.org/) & [Alabaster
0.7.12](https://github.com/bitprophet/alabaster) \| [Page
source](_sources/tutorial_03_secondorder_phase_transitions.rst.txt)
:::
