---
layout: page
title: PbTe structural instability tutorial
---

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

