::: {.document}
::: {.documentwrapper}
::: {.bodywrapper}
::: {.body role="main"}
::: {#hands-on-session-7-calculation-of-the-electron-phonon-interaction-and-superconducting-properties-with-the-sscha .section}
Hands-on-session 7: Calculation of the electron-phonon interaction and superconducting properties with the SSCHA[¶](#hands-on-session-7-calculation-of-the-electron-phonon-interaction-and-superconducting-properties-with-the-sscha "Permalink to this headline"){.headerlink}
===============================================================================================================================================================================================================================================================================

In this hands-on-session we will learn to calculate the electron-phonon
interaction and superconductivity properties in strongly anharmonic
systems combining electron-phonon matrix elements calculated within
density-functional perturbation theory (DFPT), as implemented in Quantum
Espresso (QE), and the anharmonic phonon frequencies and polarization
vectors obtained with the SSCHA. For that purpose we will use a slightly
modified version of Quantum Espresso version 5.1.0, which includes some
extra features developed by us that can be used to combine the
electron-phonon matrix elements with real space SSCHA force constants,
calculate the [\\(\\alpha\^2F(\\omega)\\)]{.math .notranslate
.nohighlight}, calculate [\\(T_c\\)]{.math .notranslate .nohighlight}
with empirical equations, and solve isotropic Migdal-Eliashberg
equations.

::: {#calculation-of-the-electron-phonon-matrix-elements .section}
Calculation of the electron-phonon matrix elements[¶](#calculation-of-the-electron-phonon-matrix-elements "Permalink to this headline"){.headerlink}
----------------------------------------------------------------------------------------------------------------------------------------------------

As an example we will use a non-converged calculation on PdH, a strongly
anharmonic superconductor, in which the Pd atoms form a fcc lattice and
H atoms sit at the octahedral interstitial sites. The crystal structure
is the rock-salt one, whose space group is [\\(Fm\\bar{3}m\\)]{.math
.notranslate .nohighlight}.

We will perform first a harmonic phonon calculation and the calculation
of the electron-phonon coupling constant for the irreducible q points in
a 2x2x2 grid. In order to know how many irreducible q points are there
for a given crystal one can use the kpoints.x program of QE (for
instance, in the QE version distributed with the SSCHA, one can find it
in qe-5.1.0_elph/PW/tools). In this case there are only three q points
in the irreducible grid. For a particular q point, the calculation of
the electron-phonon matrix elements is performed in three steps:

1.  Perform a standard DFT calculation of the crystal structure

2.  Perform a standard DFPT to calculate the harmonic dynamical matrix
    at a particular point in the grid as well as the derivative of the
    Kohn-Sham (KS) potential for this particular q point. The latter
    will be needed for the electron-phonon calculation.

3.  Perform a non-self-consistent calculation of the band structure in a
    finer electronic k point grid, read the derivative of the KS
    potential, read the harmonic dynamical matrix to obtain the phonon
    frequencies and polarization vectors, and calculate the
    electron-phonon matrix elements.

The input we will use for the standard DFT calculation that we will use
is the following:

::: {.highlight-fortran .notranslate}
::: {.highlight}
    &control
        ! Type of calculation
        calculation     =       'scf'
        ! Show more details in the output
        verbosity       =       'high'
        ! Calculate stress tensor
        tstress         =       .true.
        ! Calculate forces
        tprnfor         =       .true.
        ! Prefix for tmp files
        prefix          =       'pdh'
        ! Location of the pseudopotentials
        pseudo_dir      =       './'
        ! Folder for the tmp files
        outdir          =       './tmp'
    /
    &system
        ! Type of unit cell
        ibrav           =       2
        ! Lattice parameter in a0, Bohr length
        celldm(1)       =       7.80
        ! Number of atoms
        nat             =       2
        ! Number of atom types
        ntyp            =       2
        ! Plane-wave cutoff in Ry
        ecutwfc         =       30.0
        ! Density cutoff in Ry
        ecutrho         =       300.0
        ! It is a metal so use smearing
        occupations     =       'smearing'
        ! Type of smearing
        smearing        =       'mp'
        ! Broadening of the smearing in Ry
        degauss         =       0.020
    /
    &electrons
        ! Parameter for the DFT scf cycle
        mixing_beta     =       0.7
        ! Energy threshold to stop the scf cycle
        conv_thr        =       1.0d-8
    /
    ATOMIC_SPECIES
        Pd  106.42   Pd.pz-nd-rrkjus.UPF
        H   1.00794  H.pz-rrkjus.UPF
    ATOMIC_POSITIONS {crystal}
        Pd       0.0000        0.0000        0.0000
        H        0.5000        0.5000        0.5000
    K_POINTS {automatic}
        10 10 10 1 1 1
:::
:::

The reader is referred to the official [QE
guide](https://www.quantum-espresso.org/Doc/INPUT_PW.html){.reference
.external} to check the all the details about the input files used, even
if a short description is provided here. The pseudopotentials can be
found in the folder 07_simple_electron_phonon/pseudos/ .

If this input file was named as pw.in, we would run QE as follows:

::: {.highlight-console .notranslate}
::: {.highlight}
    qe-5.1.0_elph/bin/pw.x < pw.in > pw.out
:::
:::

This will calculate the Kohn-Sham potential, needed for the phonon and
electron-phonon calculations, apart from the usual total energy and
forces of the structure. Note that this is not a converged calculation.
One should carefully check the convergence with respect to the cutoffs,
smearing, k-point grids, etc.

The second step is to calculate the harmonic dynamical matrix within
DFPT. This is a model input file to calculate it in the first (the
[\\(\\Gamma\\)]{.math .notranslate .nohighlight} point) q point of the
2x2x2 grid.

::: {.highlight-fortran .notranslate}
::: {.highlight}
    Phonon calculation on the 1st point of a 2x2x2 q point grid
    &inputph
        ! Prefix for tmp files (same as for pw.x)
        prefix           = 'pdh'
        ! Folder for tmp files (same as for pw.x)
        outdir           = './tmp/'
        ! Name of dynamical matrices calculated
        fildyn           = "harmonic_dyn"
        ! Mass of 1st atom type in m.a.u
        amass(1)         = 106.42
        ! Mass of 2nd atom type in m.a.u
        amass(2)         = 1.00794
        ! File where the derivative of the KS potential will be stored
        fildvscf         = 'pdh_dv'
        ! Calculate the phonons in a grid nq1 x nq2 x nq3 grid of q points
        ldisp            = .true.
        nq1              = 2
        nq2              = 2
        nq3              = 2
        ! Threshold for the self-consistent loop
        tr2_ph           = 1.0d-16
        ! First q point to calculate
        start_q          = 1
        ! Last q point to calculate
        last_q           = 1
    &end
:::
:::

The reader is referred to the official [QE
guide](https://www.quantum-espresso.org/Doc/INPUT_PH.html){.reference
.external} for more details on the input parameters. If this input file
was named as ph.in, we would run QE as follows:

::: {.highlight-console .notranslate}
::: {.highlight}
    qe-5.1.0_elph/bin/ph.x < ph.in > ph.out
:::
:::

As an output we will obtain the dynamical matrix at the
[\\(\\Gamma\\)]{.math .notranslate .nohighlight} point stored in the
file harmonic_dyn1.

Once we have the dynamical matrix calculated and the derivative of the
KS potential stored we can calculate the electron-phonon matrix elements
using the modified version of QE. The input file is the following:

::: {.highlight-fortran .notranslate}
::: {.highlight}
    Phonon calculation on the 1st point of a 2x2x2 q point grid
    &inputph
        ! Prefix for tmp files (same as for pw.x)
        prefix           = 'pdh'
        ! Folder for tmp files (same as for pw.x)
        outdir           = './tmp/'
        ! Name of dynamical matrices calculated
        fildyn           = "harmonic_dyn"
        ! Mass of 1st atom type in m.a.u
        amass(1)         = 106.42
        ! Mass of 2nd atom type in m.a.u
        amass(2)         = 1.00794
        ! File where the derivative of the KS potential will be stored
        fildvscf         = 'pdh_dv'
        ! Calculate the phonons in a grid nq1 x nq2 x nq3 grid of q points
        ldisp            = .true.
        nq1              = 2
        nq2              = 2
        nq3              = 2
        ! Threshold for the self-consistent loop
        tr2_ph           = 1.0d-16
        ! First q point to calculate
        start_q          = 1
        ! Last q point to calculate
        last_q           = 1
        ! Do not calculate dynamical matrix
        trans            = .false.
        ! Type of electron-phonon interaction
        electron_phonon  = 'simple'
        ! Minimum Gaussian broadening in Ry for the double Dirac delta
        el_ph_sigma      = 0.004
        ! The number of Gaussian broadenings that will be studied
        el_ph_nsigma     = 25
        ! nk1 x nk2 x nk3 is the grid for the non-scf calculation
        ! used in the electron-phonon calculations
        nk1              = 20
        nk2              = 20
        nk3              = 20
        ! k1, k2, k3 determine whether the grid is shifted from Gamma
        k1               = 1
        k2               = 1
        k3               = 1
    &end
:::
:::

The reader is referred to the official [QE
guide](https://www.quantum-espresso.org/Doc/INPUT_PH.html){.reference
.external} for more details on the input parameters. If this input file
was named as elph.in, we would run QE as follows:

::: {.highlight-console .notranslate}
::: {.highlight}
    qe-5.1.0_elph/bin/ph.x < elph.in > elph.out
:::
:::

As an output we will obtain several files giving information on the
phonon linewidth coming from the electron-phonon interaction and so on.
Most of them can be obtained with the standard version of QE, but the
modified version we are providing here prints also the
'fildyn'.elph.d.mat.'q point number' files, in this case
harmonic_dyn1.elph.d.mat.1. This file is important for our purpose as it
is necessary to combine the SSCHA dynamical matrices with the obtained
electron-phonon matrix elements. What it contains explicitly is the
following:

::: {.math .notranslate .nohighlight}
\\\[\\Delta \^{ab}(\\mathbf{q}) =
\\frac{1}{N\_{\\mathrm{F}}N\_{\\mathbf{k}}}\\sum\_{n,n\',\\mathbf{k}}d\^{a}\_{n\\mathbf{k},n\'\\mathbf{k}+\\mathbf{q}}d\^{b}\_{n\'\\mathbf{k}+\\mathbf{q},n\\mathbf{k}}\\delta(\\epsilon
\_{n\\mathbf{k}} - \\epsilon \_{\\mathrm{F}})\\delta(\\epsilon
\_{n\'\\mathbf{k}+\\mathbf{q}} - \\epsilon \_{\\mathrm{F}}),\\\]
:::

where

::: {.math .notranslate .nohighlight}
\\\[d\^{a}\_{n\\mathbf{k}, n\'\\mathbf{k}+\\mathbf{q}} = \\langle
n\\mathbf{k} \\vert \\frac{\\delta V\_{KS}}{\\delta u
\^{a}(\\mathbf{q})}\\vert n\'\\mathbf{k} + \\mathbf{q} \\rangle\\\]
:::

are the electron-phonon matrix elements between differen KS states
[\\(\\vert n \\mathbf{k} \\rangle\\)]{.math .notranslate .nohighlight}
([\\(n\\)]{.math .notranslate .nohighlight} is aband index and
[\\(\\mathbf{k}\\)]{.math .notranslate .nohighlight} the wave number) of
the derivative of the KS potential with respect to the Fourier
transformed displacement in Cartesian basis. In the above equations
lower case latin indexes ([\\(a\\)]{.math .notranslate .nohighlight}
...) denote both atoms in the unit cell as well as Cartesian indexes.
Above [\\(N_F\\)]{.math .notranslate .nohighlight} is the density of
states (DOS) at the Fermi level per spin, [\\(\\epsilon
\_{\\mathrm{F}}\\)]{.math .notranslate .nohighlight} is the Fermi
energy, [\\(N\_{\\mathbf{k}}\\)]{.math .notranslate .nohighlight} the
number of k points in the sum. The file contains first the Gaussian
broadening used in the calculation (DOS and double Dirac delta in the
equation), the DOS at the Fermi level calculated with that broadening,
and later prints the elements of the
[\\(\\Delta\^{ab}(\\mathbf{q})\\)]{.math .notranslate .nohighlight}
matrix. Then it continues with the same data for the next broadening
calculated.

::: {.topic}
Exercise

-   Calculate the electron-phonon matrix elements for the other q points
    in the 2x2x2 grid.
:::
:::

::: {#the-sscha-calculation .section}
The SSCHA calculation[¶](#the-sscha-calculation "Permalink to this headline"){.headerlink}
------------------------------------------------------------------------------------------

This system, even with this unconverged parameters, is extremely
anharmonic and the SSCHA strongly renormalizes the phonon spectrum.

::: {.topic}
Exercise

-   Perform a SSCHA calculation on a 2x2x2 supercell to obtain
    renormalized phonon frequencies.
:::

As performing the SSCHA even with this unconverged parameters may take a
considerable time, we provide auxiliary SSCHA dynamical matrices
obtained with few configurations in 07_simple_electron_phonon/sscha/
with the name sscha_T0.0_dyn\*.
:::

::: {#combine-the-sscha-dynamical-matrices-with-the-electron-phonon-matrix-elements .section}
Combine the SSCHA dynamical matrices with the electron-phonon matrix elements[¶](#combine-the-sscha-dynamical-matrices-with-the-electron-phonon-matrix-elements "Permalink to this headline"){.headerlink}
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

We will now combine the SSCHA dynamical matrices with the
electron-phonon matrix elements calculates previously. In order to do
that (the reason will be apparent later) we will first Fourier transform
the SSCHA dynamical matrices and create the real space SSCHA force
constants. We can do that with the q2r.x code of QE. Let's first copy
the harmonic_dyn0 obtained in the phonon calculations, file that
contains the list of q points in the 2x2x2 grid, to the folder where the
SSCHA dynamical matrices are and let's name it following the new
notation:

::: {.highlight-console .notranslate}
::: {.highlight}
    cp {PATH_TO_HARMONIC_CALCULATION}/harmonic_dyn0 {PATH_TO_SSCHA_RESULTS}/sscha_T0.0_dyn0
:::
:::

Now we are ready to perform the Fourier transform. Let's prepare the
input for q2r.x:

::: {.highlight-fortran .notranslate}
::: {.highlight}
    &input
      ! Name of dynamical matrices
      fildyn  =  'sscha_T0.0_dyn'
      ! Type of ASR imposed
      zasr    =  'crystal'
      ! Name of obtained force constants
      flfrc   =  'sscha_T0.0.fc'
    /
:::
:::

More details about the input for q2r.x can be found
[here](https://www.quantum-espresso.org/Doc/INPUT_Q2R.html){.reference
.external}. If the input file was named as q2r.in, we would run it as

::: {.highlight-console .notranslate}
::: {.highlight}
    qe-5.1.0_elph/bin/q2r.x < q2r.in > q2r.out
:::
:::

This will create the sscha_T0.0.fc file with the real space force
constants.

Now we are ready to combine these SSCHA real space force constants with
the electron-phonon matrix elements calculated with DFPT. For that the
easiest thing to do is to copy the SSCHA force constants file
(sscha_T0.0.fc) and the files where the [\\(\\Delta
\^{ab}(\\mathbf{q})\\)]{.math .notranslate .nohighlight} matrices were
stored (harmonic_dyn\*.elph.d.mat.\*) to a new folder. We will use the
elph_fc.x, which is written by us based on matdyn.x of QE and is not
present in the QE distribution, to perform this calculation. The input
for this code looks as follows:

::: {.highlight-fortran .notranslate}
::: {.highlight}
    Calculation of superconducting properties with elph_fc.x
    &input
        ! ASR type imposed
        asr          =   'crystal'
        ! Mass of 1st atom in m.a.u
        amass(1)     =    106.42
        ! Mass of 2nd atom in m.a.u
        amass(2)     =    1.00794
        ! File with the SSCHA FCs
        flfrc        =   'sscha_T0.0.fc'
        ! Prefix for the files with the elph
        fildyn       =   'harmonic_dyn'
        ! Number of broadenings for the phonon Dirac Delta
        nbroad       =   1
        ! Gaussian broadening for the Dirac delta in cm-1
        minbroad     =   10
    /
    .004 ! Broadening for the electrons chosen
    3    ! Number of q points in IBZ followed by the q list (in 2pi/a) and multiplicity
    0.000000000000000E+00   0.000000000000000E+00   0.000000000000000E+00 1
    0.500000000000000E+00  -0.500000000000000E+00   0.500000000000000E+00 4
    0.000000000000000E+00  -0.100000000000000E+01   0.000000000000000E+00 3
:::
:::

If the input file was named as elph_fc.in, we would run the code as
follows:

::: {.highlight-console .notranslate}
::: {.highlight}
    qe-5.1.0_elph/bin/elph_fc.x < elph_fc.in > elph_fc.out
:::
:::

The code Fourier transforms the SSCHA real space force constants to the
list of q points provided in order to obtain the dynamica matrices at
these points, and combines them with the
[\\(\\Delta\^{ab}(\\mathbf{q})\\)]{.math .notranslate .nohighlight}
matrices to calculate the Eliashberg function
[\\(\\alpha\^2F(\\omega)\\)]{.math .notranslate .nohighlight} as

::: {.math .notranslate .nohighlight}
\\\[\\alpha\^2F(\\omega) = \\frac{1}{N\_{\\mathbf{q}}}\\sum
\_{ab\\mathbf{q}\\mu} \\frac{e\^a\_{\\mu}(\\mathbf{q})\\Delta
\^{ab}(\\mathbf{q})e\^b\_{\\mu}(\\mathbf{q})\^\*}{2\\omega\_{\\mu}(\\mathbf{q})\\sqrt{m_am_b}}
\\delta(\\omega-\\omega\_{\\mu}(\\mathbf{q})),\\\]
:::

where [\\(\\omega\_{\\mu}(\\mathbf{q})\\)]{.math .notranslate
.nohighlight} and [\\(e\^a\_{\\mu}(\\mathbf{q})\\)]{.math .notranslate
.nohighlight} are, respectively, the phonon frequencies and polarization
vectors obtained diagonalizing the Fourier interpolated SSCHA force
constants at point q, and [\\(m_a\\)]{.math .notranslate .nohighlight}
the masses of the atoms. The Dirac delta on the equation is approximated
with a Gaussian of 10 cm-1 broadening, following the input parameter. In
the output elph_fc.out the code gives the SSCHA dynamical matrix at each
q point, the phonon linewidth (HWHM), the contribution to the
electron-phonon coupling of each mode, etc. Note that the code skips the
Gamma point and does not consider it in the calculation. The reason is
that the equation used at this point is divergent (see discussion in
Appendix C in
[arXiv:2303.02621](https://arxiv.org/abs/2303.02621){.reference
.external}). The code also calculates the total electron-phonon coupling
constant [\\(\\lambda\\)]{.math .notranslate .nohighlight} and
[\\(\\omega\_{log}\\)]{.math .notranslate .nohighlight}. It also gives
the value of the the superconducting critical temperature calculated for
different values of the Coulomb pseudopotential [\\(\\mu\^\*\\)]{.math
.notranslate .nohighlight} within the semiempirical McMillan and
Allen-Dynes formulas.

The code also prints the calculated Eliashberg function, phonon density
of states, partial electron-phonon coupling constant, as well as the
projection of both the phonon DOS and Eliashberg function on different
atoms. We can do for example plots like this one with this data:

![[Fig. 15 ]{.caption-number}[Figure with the
[\\(\\alpha\^2F(\\omega)\\)]{.math .notranslate .nohighlight} Eliashberg
function, its partial contributions from Pd and H atoms. The partial
contribution to the electron-phonon coupling constant is also plotted,
[\\(\\lambda(\\omega)\\)]{.math .notranslate
.nohighlight}.]{.caption-text}[¶](#id3 "Permalink to this image"){.headerlink}](_images/a2F_PdH.pdf)

::: {.topic}
Exercise

-   Make a figure as the one above but with a different smearing for the
    double Dirac delta on the electronic states. Note that the
    electron-phonon matrix elements were calculated for smearings
    proportional to 0.04 Ry.
:::

The fact that the elph_fc.x code works with real space SSCHA force
constants allows us to combine the calculation of the electron-phonon
matrix elements in a small supercell with electron-phonon matrix
elements in a finer q point grid trivially.

::: {.topic}
Exercise

-   Calculate the electron-phonon coupling constant using the
    electron-phonon matrix elements calculated in a 4x4x4 q point grid
    combining it with the SSCHA force constants on a 2x2x2 supercell.
    For that use the data in
    07_simple_electron_phonon/elph_matrix_elements_444 .
:::
:::

::: {#solution-of-isotropic-migdal-eliashberg-equations .section}
Solution of isotropic Migdal-Eliashberg equations[¶](#solution-of-isotropic-migdal-eliashberg-equations "Permalink to this headline"){.headerlink}
--------------------------------------------------------------------------------------------------------------------------------------------------

Once the Eliashberg function [\\(\\alpha\^2F(\\omega)\\)]{.math
.notranslate .nohighlight} has been calculated combining the SSCHA and
the electron-phonon matrix elements, one can easily solve isotropic
Migdal-Eliashberg equations, which are not semiempirical. We provide a
utility to perform this calculation as well, ME.x. This is the input
file that needs to be prepared:

::: {.highlight-fortran .notranslate}
::: {.highlight}
    &inputme
         ! First temperature for the calculation
         t_first      =  1.0
         ! Last temperature for the calculation
         t_last       =  80.00
         ! Total number of temperatures
         t_number     =  80
         ! Cutoff for Matsubara frequencies in Ha
         wc_cutoff    =  0.05
         ! Initial guess for the gap in meV
         gap_guess    =  10.
         ! File with the Eliashberg function
         a2f_filename =  "a2F.10.dat"
         ! mu*
         mu_star      =  0.10
    /
:::
:::

Note that the cutoff for the Matsubara frequencies should be around 10
times the highest phonon frequency, not bigger. If this input file was
named as ME.in, we would run the code as follows:

::: {.highlight-console .notranslate}
::: {.highlight}
    qe-5.1.0_elph/bin/ME.x < ME.in > ME.out
:::
:::

In the output the code will calculate the superconducting gap as a
function of temperature. In the provided file t_gap.dat the gap as a
function of temperature is provided. This can be used to generate plots
like this one to estimate the critical temperature from the temperature
at which the gap closes.

![[Fig. 16 ]{.caption-number}[Superconducting gap as a function of
temperature.]{.caption-text}[¶](#id4 "Permalink to this image"){.headerlink}](_images/gap_PdH.pdf)
:::

::: {#important-remarks .section}
Important remarks[¶](#important-remarks "Permalink to this headline"){.headerlink}
----------------------------------------------------------------------------------

The calculations above are not converged and are just meant to
illustrate the use of the several codes. In a proper calculation one
should take into account the following points:

1.  The electron-phonon coupling constant needs to be converged with the
    number of k points for the electrons and the smearing used for the
    Double delta. These parameters enter into the equation of the
    [\\(\\Delta\^{ab}(\\mathbf{q})\\)]{.math .notranslate .nohighlight}
    matrices. The typical thing is to calculate [\\(\\lambda\\)]{.math
    .notranslate .nohighlight} for different k point grids as a function
    of the smearing and see for which low value of the smearing
    [\\(\\lambda\\)]{.math .notranslate .nohighlight} plateaus. Note
    that the physical limit is the one with infinite number of k points
    and 0 smearing.

2.  In this example we have used auxiliary SSCHA dynamical matrices to
    incorporate anharmonic effects into the calculation of
    electron-phonon properties. However, we could have also used those
    dynamical matrices that come from the Hessian of the SSCHA free
    energy. Sometimes the differences are minor, but in other cases it
    may be important. In this cases the best approach is to calculate
    the Eliashberg function with the full spectral function as described
    recently in
    [arXiv:2303.07962](https://arxiv.org/abs/2303.07962){.reference
    .external}.
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
    -   Previous: [Hands-on-session 6 - The SSCHA with machine learning
        potentials](tutorial_06_the_SSCHA_with_MLP.html "previous chapter")
    -   Next: [Hands-on-session 8: EPIq - Anharmonicity in
        electron-phonon coupling related
        properties](tutorial_08_epiq.html "next chapter")
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
source](_sources/tutorial_07_simple_electron_phonon.rst.txt)
:::
