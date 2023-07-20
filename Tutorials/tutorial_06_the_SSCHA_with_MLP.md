::: {.document} ::: {.documentwrapper} ::: {.bodywrapper} ::: {.body
role=\"main\"} :::
{\#hands-on-session-6-the-sscha-with-machine-learning-potentials
.section} Hands-on-session 6 - The SSCHA with machine learning
potentials\[¶\](\#hands-on-session-6-the-sscha-with-machine-learning-potentials
\"Permalink to this headline\"){.headerlink}
============================================================================================================================================================================
The minimization of the variational free energy demands a large number
of single-point density functional theory (DFT) calculations. These
calculations are performed on supercells, repetitions of the primitive
unit cell. DFT calculations can become very costly very fast if we need
to increase the size of the supercell. This can happen in case we have
very slowly decaying second-order force constants, large primitive unit
cells, or simply very low symmetry. In some of these cases, DFT is
prohibitive due to the large number of atoms per calculation or we
simply need a very large number of configurations to converge our
results (for example when we need to compute free energy hessian to
check the dynamical stability of the system). In the last ten years,
there has been a large amount of research put into developing
machine-learned (ML) interatomic potentials. Contrary to the traditional
interatomic potentials, they do not have a fixed analytical form and
thus are much more flexible and transferable. They are usually trained
on a very large number of DFT data and have very good accuracy. ML
potentials are considerably slower than the traditional interatomic
potentials, however still orders of magnitude faster than DFT, with a
considerably better scaling with a number of atoms. The synergy between
SSCHA and machine learning interatomic potentials is obvious. If we can
use the machine learning interatomic potentials as a calculator for
forces, stresses, and energies we can go to much larger supercells and
numbers of configurations. The stochastic sampling employed by SSCHA
gives a very good method for obtaining training sets needed to train
machine learning interatomic potentials. The force, energy and stresses
errors produced by ML interatomic potentials will influence SSCHA
results less due to the averaging effects (in case the errors are not
biased). There are a number of freely available implementations of ML
interatomic potentials (\[Gaussian Approximation
Potentials\](https://libatoms.github.io/GAP/){.reference .external},
\[NequIP\](https://nequip.readthedocs.io/en/latest/){.reference
.external},
\[pacemaker\](https://pacemaker.readthedocs.io/en/latest/){.reference
.external}, etc.), and at this point, they can be used without a large
prior knowledge of the theory behind ML potentials. :::
{\#hands-on-exercise .section} Hands-on
exercise\[¶\](\#hands-on-exercise \"Permalink to this
headline\"){.headerlink}
\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\--
For this exercise we will be using \[Gaussian Approximation
Potentials\](https://libatoms.github.io/GAP/){.reference .external},
however, the framework can be applied to any other type of ML
interatomic potential. In the exercise, we will obtain the training data
from the Tersoff interatomic potential, instead of the DFT. We have
provided starting dynamical matrices calculated for the structure at 0
K. Now we will calculate training and test ensemble with Tersoff
potential: ::: {.highlight-python .notranslate} ::: {.highlight} from
quippy.potential import Potential import cellconstructor as CC import
cellconstructor.Phonons import sscha, sscha.Ensemble temperature = 0.0
\# Temperature at which we generate SSCHA configurations nconf_train =
1000 \# Number of configurations in the training set nconf_test = 500 \#
Number of configurations in the test set \# Load the Tersoff potential
that we want to fit with ML GAP pot = Potential(\'IP Tersoff\',
param_filename=\'./06_the_SSCHA_with_machine_learning_potentials/ip.parms.Tersoff.xml\')
\# Load dynamical matrices dyn_prefix =
\'./06_the_SSCHA_with_machine_learning_potentials/start_dyn\' nqirr = 3
dyn = CC.Phonons.Phonons(dyn_prefix, nqirr) \# Generate training
ensemble ensemble_train = sscha.Ensemble.Ensemble(dyn, T0=temperature,
supercell = dyn.GetSupercell()) ensemble_train.generate(N = nconf_train)
ensemble_train.compute_ensemble(pot, compute_stress = True,
stress_numerical = False, cluster = None, verbose = True) \# This line
will save ensemble in correct format
ensemble_train.save_enhanced_xyz(\'train.xyz\', append_mode = False,
stress_key = \"stress\", forces_key = \"forces\", energy_key =
\"energy\") \# Generate test ensemble ensemble_test =
sscha.Ensemble.Ensemble(dyn, T0=temperature, supercell =
dyn.GetSupercell()) ensemble_test.generate(N = nconf_test)
ensemble_test.compute_ensemble(pot, compute_stress = True,
stress_numerical = False, cluster = None, verbose = True)
ensemble_test.save_enhanced_xyz(\'test.xyz\', append_mode = False,
stress_key = \"stress\", forces_key = \"forces\", energy_key =
\"energy\") ::: ::: The training of the ML interatomic potential can be
done with a command gap_fit which should be available after installing
quippy-ase. This command takes a large number of arguments so it is
easier to make a bash script. We will name it train.sh: :::
{.highlight-bash .notranslate} ::: {.highlight} \#!/bin/bash gap_fit
energy_parameter_name=energy force_parameter_name=forces \\
stress_parameter_name=stress virial_parameter_name=virial \\
do_copy_at_file=F sparse_separate_file=T gp_file=GAP.xml \\
at_file=train.xyz e0_method=\"average\" \\ default_sigma={0.001 0.03
0.03 0} sparse_jitter=1.0e-8 \\ gap={soap cutoff=4.2 n_sparse=200
covariance_type=dot_product \\ sparse_method=cur_points delta=0.205
zeta=4 l_max=4 \\ n_max=8 atom_sigma=0.5 cutoff_transition_width=0.8 \\
add_species } ::: ::: The meaning of each argument is not important
right now, but can be easily looked up on the official website . We run
the training command: ::: {.highlight-bash .notranslate} :::
{.highlight} bash train.sh ::: ::: Note, the training is memory
intensive, so you may need to allocate extra memory on your virtual
machine if you are employing Quantum Mobile. 4Gb of Ram are required.
You may need to restart the virtual machine. This should take a minute
or so. Once it is finished, if the memory was enough and the command
typed correctly, one should obtain the \`GAP.xml\`{.docutils .literal
.notranslate} file in the working directory containing the GAP ML
interatomic potential. We can use test.xyz file to check how well our ML
potential reproduces data with this simple script: :::
{.highlight-python .notranslate} ::: {.highlight} import numpy as np
import ase from ase import Atoms from quippy.potential import Potential
import matplotlib import matplotlib.pyplot as plt from
matplotlib.gridspec import GridSpec fpaths =
matplotlib.font_manager.findSystemFonts() infile = \'test.xyz\' \# test
datasets \# Read in .xyz files using ase method atoms =
ase.io.read(infile, \':\', format=\'extxyz\') nconf = len(atoms)
print(\'Number of configurations in the dataset: \' + str(nconf)) natoms
= \[len(at.symbols) for at in atoms\] \# Load in newly trained GAP
potential gap_file = \'./GAP.xml\' pot = Potential(\'IP GAP\',
param_filename=gap_file) \# Read in potential \# Collect previously
calculated (with Tersoff) atomic properties dft_energies =
\[atom.get_potential_energy() for atom in atoms\] dft_forces =
\[atom.get_forces() for atom in atoms\] dft_stress =
\[atom.get_stress()\[0:3\] for atom in atoms\] \# Now recalculate them
with GAP en_gap = \[\] forces_gap = \[\] stress_gap = \[\] for i in
range(nconf): if(i%100 == 0): print(\'Configuration: \', i + 1) \# Make
ase Atoms object atoms_gap = Atoms(symbols = atoms\[i\].symbols, cell =
atoms\[i\].cell,\\ scaled_positions =
atoms\[i\].get_scaled_positions(),\\ calculator = pot, pbc = True) \#
Calculate total energies of structures with GAP
en_gap.append(atoms_gap.get_potential_energy()) \# Calculate forces on
atoms forces_gap.append(atoms_gap.get_forces()) \# Calculate stress and
only take diagonal elements
stress_gap.append(atoms_gap.get_stress()\[0:3\]) \# Calculate errors
energy_errors = np.zeros_like(en_gap) forces_errors =
np.zeros_like(forces_gap) GPa = 1.60217733e-19\*1.0e21 stress_errors =
np.zeros_like(stress_gap) for i in range(nconf): \# Calculate energy
errors energy_errors\[i\] = (atoms\[i\].get_potential_energy() -\\
en_gap\[i\])/natoms\[i\] \# Calculate errors on forces
forces_errors\[i\] = atoms\[i\].get_forces() - forces_gap\[i\] \#
Calculate errors on stress stress_errors\[i\] = dft_stress\[i\] -
stress_gap\[i\] \# Function to plot Tersoff vs GAP results def
plot_comparison(ax, data1, data2, data3, \\ xlabel = \'Original energy
(eV)\', ylabel = \'ML energy (eV)\'): sizes =
np.array(data3/np.amax(data3))\*2.0 ax.scatter(data1, data2, marker =
\'o\', s = sizes, c = \'red\') lims = \[np.min(\[ax.get_xlim(),
ax.get_ylim()\]),\\ np.max(\[ax.get_xlim(), ax.get_ylim()\])\] \# now
plot both limits against eachother ax.plot(lims, lims, \'k-\',
alpha=0.75, zorder=0) ax.set_xlabel(xlabel) ax.set_ylabel(ylabel) \#
Function to plot histogram of errors def plot_error_histogram(ax, x,
nbins, xlabel): import scipy.stats as st from scipy.stats import norm
ax.hist(x, density=True, bins=nbins) mu, std = norm.fit(x) xmin, xmax =
ax.get_xlim() x1 = np.linspace(xmin, xmax, 100) p = norm.pdf(x1, mu,
std) ax.plot(x1, p, \'k\', linewidth=2) ax.set_ylabel(\"Probability\")
ax.set_xlabel(xlabel) \# Sometimes forces arrays can be ragged list,
this will flatten them def flatten(xs): res = \[\] def loop(ys): for i
in ys: if isinstance(i, list): loop(i) elif(isinstance(i, np.ndarray)):
loop(i.tolist()) else: res.append(i) loop(xs) return res \# Plot stuff
plt.rcParams\[\"font.family\"\] = \"Times New Roman\"
plt.rcParams\[\'mathtext.fontset\'\] = \"stix\"
plt.rcParams.update({\'font.size\': 16}) fig =
plt.figure(figsize=(6.4\*3.0, 4.8\*2.0)) gs1 = GridSpec(2, 3) ax00 =
fig.add_subplot(gs1\[0, 0\]) plot_comparison(ax00, dft_energies, en_gap,
energy_errors, \\ xlabel = \'Original energy (eV)\', ylabel = \'ML
energy (eV)\') ax01 = fig.add_subplot(gs1\[0, 1\]) plot_comparison(ax01,
dft_forces, forces_gap, forces_errors, \\ xlabel = r\'Original force
(eV/\$\\AA\$)\', ylabel = r\'ML force (eV/\$\\AA\$)\') ax02 =
fig.add_subplot(gs1\[0, 2\]) plot_comparison(ax02,
np.array(flatten(dft_stress))\*GPa, \\
np.array(flatten(stress_gap))\*GPa,
np.array(flatten(stress_errors))\*GPa, \\ xlabel = \'Original stress
(GPa)\', ylabel = \'ML stress (GPa)\') ax10 = fig.add_subplot(gs1\[1,
0\]) plot_error_histogram(ax10, energy_errors, 100, \'Energy error
(eV/atom)\') ax11 = fig.add_subplot(gs1\[1, 1\]) flattened_forces =
flatten(forces_errors) plot_error_histogram(ax11, flattened_forces, 100,
\'Force error (eV/\$\\AA\$)\') ax12 = fig.add_subplot(gs1\[1, 2\])
plot_error_histogram(ax12, np.array(\[item for sublist in stress_errors
for item in sublist\])\*GPa,\\ 100, \'Stress error (GPa)\')
fig.savefig(\'test.pdf\') plt.show() ::: ::: In the upper panel figures,
ideally, we would like points to lie on the diagonal. When fitting
interatomic potential we aim for normal distribution of errors centered
at 0 (without bias) with as small as possible standard deviation. We
should have very nice results for energies and forces. Now that we are
happy with the potential let us use it to relax SSCHA at 2000 K: :::
{.highlight-python .notranslate} ::: {.highlight} from quippy.potential
import Potential import cellconstructor as CC import
cellconstructor.Phonons \# Import the SSCHA engine (we will use it
later) import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax \#
Declare SSCHA variables temperature = 2000.0 nconf = 1000 max_pop =
10000 \# Load in the GAP potential gap_file = \'./GAP.xml\' pot =
Potential(\'IP GAP\', param_filename=gap_file) \# Load in the SSCHA
dynamical matrices dyn_prefix =
\'./06_the_SSCHA_with_machine_learning_potentials/start_dyn\' nqirr = 3
dyn = CC.Phonons.Phonons(dyn_prefix, nqirr) \# Relax the structure at
2000 K ensemble = sscha.Ensemble.Ensemble(dyn, T0=temperature, supercell
= dyn.GetSupercell()) ensemble.generate(N = nconf) minimizer =
sscha.SchaMinimizer.SSCHA_Minimizer(ensemble) minimizer.min_step_dyn =
0.1 minimizer.kong_liu_ratio = 0.5 minimizer.meaningful_factor = 0.001
minimizer.max_ka = 100000 relax = sscha.Relax.SSCHA(minimizer,
ase_calculator = pot, N_configs = nconf, max_pop = max_pop,
save_ensemble = True) relax.vc_relax(ensemble_loc=\'Ensemble_location\')
relax.minim.dyn.save_qe(\'final_dyn\') \# We can check minimization
procedure relax.minim.plot_results(save_filename = \'sscha\', plot =
False) ::: ::: We have relaxed SSCHA at 2000 K. We can check that
everything went well in \"sscha\" file. However, we do not know whether
this is correct. We need to check how our ML potential performs at 2000
K. ::: {.topic} Exercise: Let\'s create a dataset of SSCHA-generated
configuration at 2000 K using GAP relaxed dynamical matrices and compute
it using Tersoff potential. Next, check the performance of the GAP ML
potential against this new dataset. ::: ::: {.topic} Excercise: We
should see GAP performing quite worse compared to the test.xyz case. How
can we improve GAP potential? Let\'s do it. ::: ::: {.topic} Excercise:
How do the Tersoff phonons compare to GAP phonons? ::: ::: {.topic}
Excercise: Does this translate to larger supercells? ::: ::: ::: ::: :::
::: ::: {.sphinxsidebar role=\"navigation\" aria-label=\"main
navigation\"} ::: {.sphinxsidebarwrapper} ::: {.relations} \#\#\#
Related Topics - \[Documentation overview\](index.html) - Previous:
\[Hands-on-session 5 - Raman and Infrared spectra with the
Time-Dependent Self Consistent Harmonic
Approximation\](tutorial_05_ramanir.html \"previous chapter\") - Next:
\[Hands-on-session 7: Calculation of the electron-phonon interaction and
superconducting properties with the
SSCHA\](tutorial_07_simple_electron_phonon.html \"next chapter\") :::
::: {\#searchbox style=\"display: none\" role=\"search\"} \#\#\# Quick
search {\#searchlabel} ::: {.searchformwrapper} ::: ::: ::: ::: :::
{.clearer} ::: ::: ::: {.footer} ©2023, Lorenzo Monacelli. \\\| Powered
by \[Sphinx 4.2.0\](http://sphinx-doc.org/) & \[Alabaster
0.7.12\](https://github.com/bitprophet/alabaster) \\\| \[Page
source\](\_sources/tutorial_06_the_SSCHA_with_MLP.rst.txt) :::
