::: {.document}
::: {.documentwrapper}
::: {.bodywrapper}
::: {.body role="main"}
::: {#hands-on-session-8-epiq-anharmonicity-in-electron-phonon-coupling-related-properties .section}
Hands-on-session 8: EPIq - Anharmonicity in electron-phonon coupling related properties[¶](#hands-on-session-8-epiq-anharmonicity-in-electron-phonon-coupling-related-properties "Permalink to this headline"){.headerlink}
===========================================================================================================================================================================================================================

::: {#introduction .section}
Introduction[¶](#introduction "Permalink to this headline"){.headerlink}
------------------------------------------------------------------------

In this hands-on session we learn how to include anharmonic effects
calculated within the SSCHA in the calculation of electron-phonon
coupling related properties using
[EPIq](https://the-epiq-team.gitlab.io/epiq-site/){.reference
.external}.

[![logo_epiq](_images/SSCHA_epiq_Logo.png)](_images/SSCHA_epiq_Logo.png){.reference
.internal .image-reference}

In some systems the first principles calculation of electron-phonon
coupling matrix elements can be demanding. EPIq (Electron-Phonon wannier
Interpolation over k and q-points) is an open-source software that
allows to speed up the calculation of electron--phonon coupling related
properties using the Wannier interpolation technique. Details on the
interpolation scheme can be found
[here](https://the-epiq-team.gitlab.io/epiq-site/docs/th_found){.reference
.external}. Within EPIq, it is possible to include anharmonic
corrections to the dynamical matrices as calculated within the SSCHA.
:::

::: {#requirements .section}
Requirements[¶](#requirements "Permalink to this headline"){.headerlink}
------------------------------------------------------------------------

In the interest of time, in this hands-on session the following starting
data are at your disposal:

1.  Electron-phonon matrix elements
    [\\(g\^{\\nu}\_{m,n}(\\mathbf{k},\\mathbf{q})\\)]{.math .notranslate
    .nohighlight} computed from first principles.

2.  Wannier interpolation files which encode the trasformation to the
    optimally smooth subspace, [\\(U\_{mn}\\)]{.math .notranslate
    .nohighlight} : [\\(\\ket{\\mathbf{R}n} = \\frac{1}{\\sqrt{N_k\^w}}
    \\sum\_{\\bf k=1}\^{N_k\^w}\\sum\_{m=1}\^{N\_{\\rm
    w}}e\^{-i\\mathbf{k}\\cdot
    \\mathbf{R}}U\_{mn}(\\mathbf{k})\|\\psi\_{{\\bf
    k}m}\\rangle\\)]{.math .notranslate .nohighlight}

3.  Anharmonic dynamical matrices
    [\\(D\^{SCHA}\_{\\mu,\\nu}(\\mathbf{k},\\mathbf{q})\\)]{.math
    .notranslate .nohighlight}.

4.  Harmonic dynamical matrices (as a reference)
    [\\(D\^{HARM}\_{\\mu,\\nu}(\\mathbf{k},\\mathbf{q})\\)]{.math
    .notranslate .nohighlight}.

::: {.admonition .attention}
Attention

A folder prepared for you with these data for the present tutorial can
be downloaded
[08_EPIq](https://ehubox.ehu.eus/s/Y48Wc8iX9Z76jqN?path=%2F08_EPIq){.reference
.external} folder in the shared cloud. Right click on
`tutorial_data`{.docutils .literal .notranslate} on the navigation bar
and download the whole folder. It contains:

> <div>
>
> 1.  The electron-phonon coupling matrix elements can be found in the
>     `mat_elem`{.docutils .literal .notranslate} folder. Each file
>     corresponds to a different [\\(\\mathbf{q}\\)]{.math .notranslate
>     .nohighlight}-point in the first Brillouin zone.
>
> 2.  The `Wannier`{.docutils .literal .notranslate} folder contains the
>     files `(.eig, .chk)`{.docutils .literal .notranslate}.
>
> 3.  The dynamical matrices are stored in the `dyn_mat`{.docutils
>     .literal .notranslate} directory. `dynq*`{.docutils .literal
>     .notranslate} files are harmonic
>     ([\\(D\^{SCHA}\_{\\mu,\\nu}(\\mathbf{k},\\mathbf{q})\\)]{.math
>     .notranslate .nohighlight}) while `MoS2.Hessian.dyn*`{.docutils
>     .literal .notranslate} are anharmonic dynamical matrices computed
>     with the SSCHA code
>     ([\\(D\^{SCHA}\_{\\mu,\\nu}(\\mathbf{k},\\mathbf{q})\\)]{.math
>     .notranslate .nohighlight}).
>
> Place all the downloaded material where you intend to run the
> tutorial. A suggested structure is for example:
>
> </div>

::: {.highlight-console .notranslate}
::: {.highlight}
    handson_8/
        |
        + ----- epiq/
        |   |
        |   + ----- bin/
        |   |
        |   + ----- src/
        |
        |
        + ----- MoS2/
            |
            + ----- MoS2.eig
            |
            + ----- MoS2.chk
            |
            + ----- mat_elem/
            |   |
            |   + ----- MoS2_elph.mat.1_q*
            |
            + ----- dyn_mat/
                |
                + ----- dynq*
                |
                + ----- MoS2.Hessian.dyn*
:::
:::
:::
:::

::: {#about-epiq .section}
About EPIq[¶](#about-epiq "Permalink to this headline"){.headerlink}
--------------------------------------------------------------------

![[Fig. 17 ]{.caption-number}[Epiq site:
<https://the-epiq-team.gitlab.io/epiq-site/>]{.caption-text}[¶](#id3 "Permalink to this image"){.headerlink}
¶ Epiq paper:
<http://arxiv.org/abs/2306.15462>](_images/qr-code-epiq-both.png)

The Electron-Phonon Intepolation over q package exploits Wannier
interpolation to obtain many proprieties in solids. Bloch theorem allows
to describe an infinite system in real space with a continuum of
Hamiltonian for different k-points in the Brillouin zone.

The properties of a material refers to a certain observable averaged
over the sample.

Thanks to the Bloch theorem it is often convenient to perform the
average in the reciprocal k-space. In other words as a sum over the
whole Brillouin zone of the quantities defined at each k-point.

::: {.math .notranslate .nohighlight}
\\\[\\langle O \\rangle=\\frac{1}{N_k}\\sum\_{\\mathbf{k}}
F_O(\\mathbf{k})\\\]
:::

The quality of the averaging depends on the finesse of the sampling
[\\(N_k\\)]{.math .notranslate .nohighlight} and of the smoothness of
the integrand function [\\(F_O(\\mathbf{k})\\)]{.math .notranslate
.nohighlight} (a constant function is totally described with just one
k-point).

[![sampling](_images/sampling.png){.align-center}](_images/sampling.png){.reference
.internal .image-reference}

Wannier interpolation is an efficient way to refined the Brillouin zone
sampling.

[![wannier_interpolation](_images/wannier_interpolation.png){.align-center}](_images/wannier_interpolation.png){.reference
.internal .image-reference}

::: {#calculations-available-in-epiq .section}
### Calculations available in epiq[¶](#calculations-available-in-epiq "Permalink to this headline"){.headerlink}

1.  Adiabatic (static) and non-adiabatic (dynamic) force constant
    matrices.

2.  **Electron-phonon contribution to the phonon linewidth** and related
    quantities.

3.  Isotropic and **anisotropic Eliashberg equations**.

4.  Double Resonant Raman scattering.

5.  Electron lifetime and relaxation time.
:::

::: {#epiq-workflow .section}
### EPIq workflow[¶](#epiq-workflow "Permalink to this headline"){.headerlink}

[![wannier](_images/workflow.png){.align-center}](_images/workflow.png){.reference
.internal .image-reference}

The core steps of any calculation employing EPIq are performed using the
main executable `epiq.x`{.docutils .literal .notranslate} and consists
in two main stages.

1.  A preliminary step where electron-phonon coupling matrix elements
    and the Hamiltonian are Fourier-transformed to real space and
    written to file.

2.  The electron-phonon coupling matrix elements and the Hamiltonian are
    transformed back to reciprocal space to compute the property of
    interest at arbitrary k- and q- values.
:::

::: {#epiq-input-file .section}
### EPIq input file[¶](#epiq-input-file "Permalink to this headline"){.headerlink}

The input file is divided in three namelists:

1.  &control, specifying what calculation the code will perform

2.  &electrons, specifying the electronic parameters of the property to
    be computed

3.  &phonons, specifying the phonons parameters for the property to be
    computed

Finally, the last lines of the input indicates the electron momentum
(k-) and phonon momentum (q-) meshes on which the matrix elements are
interpolated to.
:::
:::

::: {#let-s-practice-calculation-of-electron-phonon-coupling-related-properties-for-doped-monolayer-rm-mos-2 .section}
Let's practice: calculation of electron-phonon coupling related properties for doped monolayer [\\({\\rm MoS}\_2\\)]{.math .notranslate .nohighlight}[¶](#let-s-practice-calculation-of-electron-phonon-coupling-related-properties-for-doped-monolayer-rm-mos-2 "Permalink to this headline"){.headerlink}
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

[![wannier](_images/mos2.png){.align-center}](_images/mos2.png){.reference
.internal .image-reference}

In this tutorial we calculate electron-phonon coupling related
properties for doped monolayer [\\({\\rm MoS}\_2\\)]{.math .notranslate
.nohighlight}, and evaluate the effect of anharmonicity on them.

::: {#wannier-interpolation .section}
### Wannier interpolation[¶](#wannier-interpolation "Permalink to this headline"){.headerlink}

As explained in the previous section, any calculation within EPIq starts
with a preliminary stepreferred to as `dump`{.docutils .literal
.notranslate}. During this step, the Hamiltonian and the electron-phonon
coupling matrix elements in real space are computed and written to file,
to be used in any subsequent calculation. This step is performed only
once. The input file is the following:

::: {.highlight-fortran .notranslate}
::: {.highlight}
    &control
        dump_gR=.true.,
        prefix='MoS2',
        elphmat_dir="./mat_elem/"
      /
      &electrons
      /
      &phonons
        nq1=8,
        nq2=8,
        nq3=1,
      /
:::
:::

Notice, in particular, the following paramters:

-   `dump_gR`{.docutils .literal .notranslate} in the namelist
    `control`{.docutils .literal .notranslate} tell the program to save
    the auxiliary file containing the Hamiltonian and the
    electron-phonon coupling matrix elements in real space. Since we are
    not calculating any property of the system the namelist
    `electron`{.docutils .literal .notranslate} and `phonons`{.docutils
    .literal .notranslate} are essentially empty. Only
    `nq1,nq2,nq3`{.docutils .literal .notranslate} have to be supplied
    in order to specifie the q-points mesh where the input
    electron-phonon coupling matrix elements were computed (in this
    case, a 8x8x1 q-grid).

::: {.admonition .note}
Note

In order to keep everything tidy you can use keep the el-ph matrix
elements in a separate folder and use the variable:

&control \> elphmat_dir="\<path-to-folder\>"
:::

Run this preliminary calculation using:

::: {.highlight-console .notranslate}
::: {.highlight}
    mpirun -n <NPROC> {$path_to_epiq}/bin/epiq.x -inp dump.in > dump.out
:::
:::

After the dump has ended, some output files have been produced.

-   The output file is binary and named `G_and_H.bin`{.docutils .literal
    .notranslate} contains the Hamiltonian and the electron-phonon
    coupling matrix elements in the Wannier representation. If, however,
    `ascii_G_and_H=.true.`{.docutils .literal .notranslate} is added in
    the `&control`{.docutils .literal .notranslate} namelist then the
    produced output file is readable and is named
    `G_and_H.asc`{.docutils .literal .notranslate}.

-   .dat files containing the real space localization of the Hamiltonian
    and the electron-phonon coupling matrix elements.
:::

::: {#check-real-space-localization .section}
### Check real space localization[¶](#check-real-space-localization "Permalink to this headline"){.headerlink}

If the Wannier transformation is well converged, the matrix elements are
optimally localized in real space. Always check their localization.

::: {.topic}
Excercise:

Using the two-columns files:

> <div>
>
> "MoS2_gw_R\_ph.mu\*.dat.pe_1"
>
> :   [\\(\|R\| \\qquad \\sum\_{m,n} \|
>     g\^{\\nu}\_{m,n}(0,\\mathbf{R})\|\^2\\)]{.math .notranslate
>     .nohighlight}
>
> "MoS2_gw_R\_el.mu\*.dat.pe_1"
>
> :   [\\(\|r\| \\qquad \\sum\_{m,n} \|
>     g\^{\\nu}\_{m,n}(\\mathbf{r},0)\|\^2\\)]{.math .notranslate
>     .nohighlight}
>
> </div>

plot the averaged modulus of the electron-phonon matrix elements as a
function of the distance in real space [\\(\|R\|\\)]{.math .notranslate
.nohighlight}. Are they localized?
:::
:::

::: {#phonon-linewidth-calculation .section}
### Phonon linewidth calculation[¶](#phonon-linewidth-calculation "Permalink to this headline"){.headerlink}

We would now focus on one of the system properties that can be
calculated with EPIq: the phonon linewidth
[\\(\\gamma\_{\\mathbf{q},\\nu}\\)]{.math .notranslate .nohighlight}. We
will consider the Allen phonon linewidth, which is defined by the
following equation:

::: {.math .notranslate .nohighlight}
\\\[\\gamma\_{\\mathbf{q},\\nu} = \\frac{4 \\pi
\\omega\_{\\mathbf{q},\\nu}}{N_k}\\sum\_{m,n}\\sum\_{\\mathbf{k}}\|g\^{\\nu}\_{m,n}(\\mathbf{k},\\mathbf{q})\|\^2\\delta(\\epsilon\_{\\mathbf{k}+\\mathbf{q},m}-\\epsilon\_{F})\\delta(\\epsilon\_{\\mathbf{k},n}-\\epsilon\_{F})\\\]
:::

where the electron-phonon coupling is defined from the deformation
potential as:

::: {.math .notranslate .nohighlight}
\\\[g\^{\\nu}\_{m,n}(\\mathbf{k},\\mathbf{q}) =
\\sum_s\\mathbf{e}\^s\_{\\mathbf{q},\\nu}\\cdot\\mathbf{d}\^s\_{m,n}(\\mathbf{k},\\mathbf{q})/\\sqrt{2M_s\\omega\_{\\mathbf{q},\\nu}}\\\]
:::

::: {#example-of-input-file-for-linewidth-calculation-of-monolayer-rm-mos-2 .section}
#### Example of input file for linewidth calculation of monolayer [\\({\\rm MoS}\_2\\)]{.math .notranslate .nohighlight}[¶](#example-of-input-file-for-linewidth-calculation-of-monolayer-rm-mos-2 "Permalink to this headline"){.headerlink}

We first calculate the mode-resolved [\\(\\gamma\\)]{.math .notranslate
.nohighlight} at the M-point of the Brillouin zone. In the following
example we perform the calculation for phonon of momentum
[\\(\\mathbf{q}=\\mathbf{M}\\)]{.math .notranslate .nohighlight} for two
values of the electronic smearing. The input parameters are explained in
detail in the [EPIq
manual](https://the-epiq-team.gitlab.io/epiq-site/docs/manual/){.reference
.external} . Here is the input file:

::: {.highlight-fortran .notranslate}
::: {.highlight}
    &control
      prefix='MoS2',
      calculation='ph_linewidth',
      read_dumped_gr=.true.,
      dump_gR=.false.,
      elphmat_dir="./mat_elem/"
      out2json=.true.
     /

     &electrons
          ngauss=0,
          sigma_min=0.01,
          sigma_max=0.05
          nsigma=5,
      /

      &phonons
          use_alternative_dyn=.true.
          prefix_alt_dyn='./dyn_mat/dynq'
          Fourier_interp_dyn=.true.,
          nq1=8,nq2=8,nq3=1,
      /

      k-points
      automatic
      4 4 1 0 0 0

      q-points
      crystal
      1
      0.5 0 0 1 ! M
:::
:::

The linewidth calculation is then started by:

::: {.highlight-console .notranslate}
::: {.highlight}
    mpirun -n <NPROC> {$path_to_epiq}/bin/epiq.x -inp lw.in > lw.out
:::
:::

Note that this calculation is done within the harmonic approximation, as
we are using the dynamical matrices indicated by the variable
`prefix_alt_dyn='./dyn_mat/dynq'`{.docutils .literal .notranslate}.
:::

::: {#parameters .section}
#### Parameters[¶](#parameters "Permalink to this headline"){.headerlink}

> <div>
>
> -   The linewidth calculation is selected by setting the
>     `calculation`{.docutils .literal .notranslate} parameter equal to
>     `'ph_linewidth'`{.docutils .literal .notranslate} in the
>     `&control`{.docutils .literal .notranslate} namelist.
>
> -   In the namelist `control`{.docutils .literal .notranslate},
>     `read_dumped_gr`{.docutils .literal .notranslate}, which is set to
>     `.true.`{.docutils .literal .notranslate}, indicating that
>     electron-phonon coupling matrix elements can be read from
>     `G_and_H.bin`{.docutils .literal .notranslate} )
>
> -   In the namelist `electrons`{.docutils .literal .notranslate},
>     `sigma_min`{.docutils .literal .notranslate},
>     `sigma_max`{.docutils .literal .notranslate}, `ngauss`{.docutils
>     .literal .notranslate} and `nsigma`{.docutils .literal
>     .notranslate} specify maximum, minimum, type and number of
>     electronic smearing values to use. `efermi`{.docutils .literal
>     .notranslate} and `ef_from_input`{.docutils .literal .notranslate}
>     specify the initial guess for the Fermi level ( the Fermi level
>     calculated by Quantum ESPRESSO is usually a good guess) and
>     whether the Fermi level should be re-calculated by epiq (
>     `ef_from_input`{.docutils .literal .notranslate} equal to
>     `.false.`{.docutils .literal .notranslate}) or set from input (
>     `ef_from_input`{.docutils .literal .notranslate} equal to
>     `.true.`{.docutils .literal .notranslate}). The variable
>     `thr_compute_k`{.docutils .literal .notranslate} is used to
>     restrict the calculation only to k-points possessing at least one
>     eigenvalue in the specified energy region near the Fermi level.
>
> -   In the namelist `phonons`{.docutils .literal .notranslate},
>     `Fourier_interp_dyn`{.docutils .literal .notranslate} equal to
>     `.true.`{.docutils .literal .notranslate} asks to interpolate the
>     dynamical matrices for the q-points that do not belong to the
>     Wannier grid. Alternatively, *EPIq* gives the opportunity to read
>     eigenvalues and eigenvectors produced by `matdyn.x`{.docutils
>     .literal .notranslate} of the Quantum ESPRESSO package
>     (`matdyn.eig`{.docutils .literal .notranslate} file), putting
>     `Fourier_interp_dyn=.false.`{.docutils .literal .notranslate} and
>     `read_modes=.true.`{.docutils .literal .notranslate} in input, or
>     even to directly read dynamical matrices from a
>     `matdyn.dyn`{.docutils .literal .notranslate} file, putting
>     `Fourier_interp_dyn=.false.`{.docutils .literal .notranslate} and
>     `read_modes=.false.`{.docutils .literal .notranslate} .
>
> </div>
:::

::: {#output-files .section}
#### Output files[¶](#output-files "Permalink to this headline"){.headerlink}

If the `out2json`{.docutils .literal .notranslate} flat is set to
`.true.`{.docutils .literal .notranslate}, the file
`MoS2_lambda.json`{.docutils .literal .notranslate} will be produced. It
can be automatically parsed using python as in the following lines where
the variable `q_pts`{.docutils .literal .notranslate} is a "dict" whose
entries are the results of the calculation for each q-point.

::: {.highlight-python .notranslate}
::: {.highlight}
    import json
    ff = './MoS2_lambda.json'
    with open(ff,'r') as f:
        q_pts = json.load(f)

    first_q = q_pts["1"]

    print("Fraction coordinate of the q point:", first_q["xq_frac"])
    print("Electronic temperaturs used:", first_q["T"])
    print("Results for the first mode:", first_q["1"])
    print("Frequency of the second mode in meV:", first_q["2"]["freq"])
:::
:::

Here, a simple python script to plot the linewidth esteemed with the
Allen formula:

::: {.highlight-python .notranslate}
::: {.highlight}
    import matplotlib.pyplot as plt
    import numpy as np
    for q in list(q_pts.values())[:]:
        for mod in range(1,10):
            plt.plot( q["T"],q[f'{mod}']["gamma_allen"],label=r"$\omega_"+f"{mod}$={q[f'{mod}']['freq']:.0f} (meV)",linestyle='--' )
        plt.xlabel('T  (eV)')
        plt.ylabel(r'$\gamma(T)$  (meV)')
        plt.legend(title=f"q={ np.round(np.array(q['xq_frac']),2) }")
        plt.show()
:::
:::

The other output file, `MoS2_lambda.d`{.docutils .literal .notranslate},
contains all the properties calculated by EPIq. The way this file is
formatted is specified in output file, `lw.out`{.docutils .literal
.notranslate}. Notice in particular that the first column contains the
electronic smearing, while the second column contains the Allen
linewidth.

::: {.topic}
Exercise:

Perfom a convergence study of the linewidth at [\\(\\mathbf{M}\\)]{.math
.notranslate .nohighlight} as function of the smearing and the k-mesh
density. What is the minimum temperature at which the linewidth of the
eighth mode is to be considered at convergence with a k-mesh of 8x8x1 ?
And of 12x12x1?
:::

::: {.topic}
Exercise:

Which mode shows larger smearing dependence?
:::
:::

::: {#inclusion-of-anharmonicity .section}
#### Inclusion of anharmonicity[¶](#inclusion-of-anharmonicity "Permalink to this headline"){.headerlink}

Now, we want to observe what changes with the inclusion of **anharmonic
effects**. To this aim, we need to correctly specify the prefix of the
Free-energy Hessian matrices calculated within the SSCHA using
`prefix_alt_dyn='./dyn_mat/MoS2.Hessian.dyn'`{.docutils .literal
.notranslate}.

::: {.topic}
Exercise:

Compute the linewidth at [\\(\\mathbf{q}=\\mathbf{M}\\)]{.math
.notranslate .nohighlight} using anharmonic dynamical matrix computed
thanks to the SSCHA code. Why it seems like the first and the second
mode are exchanged with respect to the harmonic dynamical matrices?
Which are the modes presenting a larger anharmonic correction?
:::
:::

::: {#dispersion-along-a-line .section}
#### Dispersion along a line[¶](#dispersion-along-a-line "Permalink to this headline"){.headerlink}

Now we want to perform the calculation of phonon linewidth along a
certain crystalline direction and produce a plot like this one:

[![wannier](_images/lw.png){.align-center}](_images/lw.png){.reference
.internal .image-reference}

where the thickness of the lines is proportional to the calculated
phonon linewidth. In order to to this, we repeat the phonon linewidth
calculation, this time considering the whole [\\(\\Gamma\\)]{.math
.notranslate .nohighlight} -M path:

> <div>
>
> ::: {.highlight-fortran .notranslate}
> ::: {.highlight}
>     &control
>        dump_gR=.false.,
>        read_dumped_gr=.true.,
>        prefix='MoS2',
>        calculation='ph_linewidth',
>        elphmat_dir="./mat_elem/"
>     &end
>     &electrons
>        ngauss=0,
>        sigma_min=0.01,
>        sigma_max=0.02
>        nsigma=2,
>     &end
>     &phonons
>        use_alternative_dyn=.true.
>        prefix_alt_dyn='./dyn_mat/MoS2.Hessian.dyn'
>        Fourier_interp_dyn=.true.,
>        nq1=8,nq2=8,nq3=1,
>     &end
>     k-points
>     automatic
>     8 8 1  0  0  0
>     q-points
>     crystal
>     11
>     0.001 0 0 1
>     0.05 0 0 1
>     0.10 0 0 1
>     ...
>     ...
>     0.50 0 0 1
> :::
> :::
>
> </div>

Once the linewidth calculation has finished, we can obtain a plottable
file using the `linewidth_path.x`{.docutils .literal .notranslate}
post-processing tool. Here is an example of input file:

> <div>
>
> ::: {.highlight-fortran .notranslate}
> ::: {.highlight}
>     &input_lambda
>       prefix='MoS2'
>       lkp_sequential=.true.
>       sigma_min=0.01,
>       sigma_max=0.02,
>       nsigma=2,
>       chosen_sigma=0.01
>     &end
>      3.159998   0.000000   0.000000
>     -1.579999   2.736639   0.000000
>      0.000000   0.000000  19.141895
>      crystal
>      11
>      0.001 0 0 1
>      0.05 0 0 1
>      0.10 0 0 1
>      ...
>      ...
>      0.50 0 0 1
> :::
> :::
>
> </div>

Notice the parameter `chosen_sigma`{.docutils .literal .notranslate},
which specifies what smearing will be used to produce the plottable
file, and the lattice parameters at the end of the namelist. The
post-processing is executed as follows:

::: {.highlight-console .notranslate}
::: {.highlight}
    {$path_to_epiq}/bin/linewidth_path.x < path.in > path.out
:::
:::

Finally, use the following gnuplot script plots the q-resolved linewidth
for the acoustic modes:

> <div>
>
> ::: {.highlight-none .notranslate}
> ::: {.highlight}
>     set ylabel '{/Symbol w}(meV)'
>     set xlabel 'Gamma - M'
>     set style fill transparent solid 0.25
>     pl for [i=0:9] 'MoS2_lw_path.d' every 9::i u 1:2 w l lt rgb 'black' title ''
>     repl  for [i=0:9] 'MoS2_lw_path.d' every 9::i u ($1):($2-$3/2):($2+$3/2)\
>     w filledc lt rgb 'red' title ''
> :::
> :::
>
> </div>

::: {.topic}
Exercise:

Try to produce two plots of the whole phonon spectrum, comparing the
harmonic and the anharmonic result. Do you observe any differences?
:::
:::

::: {#advanced-tutorial-migdal-eliashberg-calculation-using-sscha-hessian-matrices .section}
#### Advanced tutorial: Migdal-Eliashberg calculation using SSCHA Hessian matrices[¶](#advanced-tutorial-migdal-eliashberg-calculation-using-sscha-hessian-matrices "Permalink to this headline"){.headerlink}

*EPIq* also allows to solve the anisotropic Eliashberg equations on the
imaginary axis in order to calculate the [\\(\\mathbf{k}\\)]{.math
.notranslate .nohighlight}-resolved superconducting gap. We give an
example for the superconducting gap of doped monolayer
[\\(\\textrm{MoS}\_2\\)]{.math .notranslate .nohighlight} at T= 1 K.

> <div>
>
> ::: {.highlight-fortran .notranslate}
> ::: {.highlight}
>     &control
>         dump_gR=.false.,
>         read_dumped_gr=.true.,
>         prefix='MoS2',
>         calculation='migdal_eliashberg',
>         elphmat_dir="./mat_elem/"
>         &end
>         &electrons
>         theta=1,
>         efermi=-2.0,
>         ef_from_input=.false.,
>     &end
>     &phonons
>         use_alternative_dyn=.false.
>         prefix_alt_dyn='./dyn_mat/MoS2.Hessian.dyn'
>         read_modes=.false.,
>         nq1=8,
>         nq2=8,
>         nq3=1,
>     &end
>     &input_migdal
>         initialize='step',
>         gap_threshold=0.025,
>         sigma_me=0.1,
>         ME_Fermi_thickness=0.4,
>         mustar=0.1,
>         nmatsu=64,
>         nitermax=100,
>         alpha_mix_me=0.5,
>         gap_init=4.0,
>     &end
>     nkfs
>     64 64 1
>     40
> :::
> :::
>
> </div>

Run the *EPIq* calculation as follows:

::: {.highlight-console .notranslate}
::: {.highlight}
    mpirun -n 4 {$path_to_epiq}/bin/epiq.x <input_ME > out_ME
:::
:::

Note that for the `migdal_eliashberg`{.docutils .literal .notranslate}
calculation, the number of [\\(\\mathbf{k}\\)]{.math .notranslate
.nohighlight}-points (40 in this case) must be a multiple of the mpi
processes (4 in this case).

Notice the following parameters in the input file:

-   In the `&electrons`{.docutils .literal .notranslate} namelist,
    `theta=1.0`{.docutils .literal .notranslate} specifies the
    temperature of the calculation in Kelvin.

-   In the `&input_migdal`{.docutils .literal .notranslate} namelist,
    `sigma_me`{.docutils .literal .notranslate} specifies the broadening
    (in eV) to be used in the calculation and
    `ME_Fermi_thickness`{.docutils .literal .notranslate} the energy
    range around the Fermi level where electron eigenvalues are
    considered in the calculation. `nmatsu`{.docutils .literal
    .notranslate} indicates the number of Matsubara frequencies to be
    employed in the sum ( see
    [here](https://the-epiq-team.gitlab.io/epiq-site/docs/calc/eliash/){.reference
    .external} for further details: )

-   The code solves the equation generating a certain number of random
    [\\(\$\$\\mathbf{k}\\)]{.math .notranslate .nohighlight}-points,
    defined by the `64 64 1`{.docutils .literal .notranslate} grid in
    input and having an eigenvalue on the Fermi surface.

EPIq also calculates the Fermi surface using the grid specified after
`nkfs`{.docutils .literal .notranslate} (64 64 1 here). After EPIq is
done, two output files are produced:

-   MgB2.bxsf contains the Fermi surface in the XCrysDen .bxsf format.

-   MgB2_ME.d contains the Fermi-surface resolved Eliashberg gap.

You can produce a plot of the Fermi-surface resolved superconducting gap
using the `plot_ME_fs.x`{.docutils .literal .notranslate} post
processing. Execute it as:

::: {.highlight-console .notranslate}
::: {.highlight}
    {$path_to_epiq}/bin/plot_ME_fs.x
:::
:::

This is, for example, what you get for [\\({\\rm MgB}\_2\\)]{.math
.notranslate .nohighlight}, with its famous double gap:

[![logo_epiq](_images/fs_mgb2.png)](_images/fs_mgb2.png){.reference
.internal .image-reference}

Try to obtain the same kind of plot for [\\({\\rm MoS}\_2\\)]{.math
.notranslate .nohighlight} using
[fermisurfer](https://mitsuaki1987.github.io/fermisurfer/){.reference
.external} ([github
repository](https://github.com/mitsuaki1987/fermisurfer.git/){.reference
.external}) typing:

::: {.highlight-console .notranslate}
::: {.highlight}
    fermisurfer MoS2.frmsf
:::
:::

::: {.topic}
Exercise:

Perform a convergence study of the average superconducting gap as
function of the number of Matsubara frequencies and k-points. What is
the required value of nmatsu and k-points to have a precision better
than 0.1 meV?
:::

::: {.admonition .note}
Note

In order to perform a complete calculation of anharmonic electron phonon
coupling related properties within SSCHA+EPIq requires the following
steps:

1.  [The SSCHA code](http://sscha.eu/){.reference .external}
    [\\(\\rightarrow\\)]{.math .notranslate .nohighlight} First, compute
    the free energy Hessian within the stochastic self-consistent
    harmonic approximation (SSCHA).

2.  [Quantum ESPRESSO
    code](https://www.quantum-espresso.org/){.reference .external}
    [\\(\\rightarrow\\)]{.math .notranslate .nohighlight} Then, compute
    electron-phonon coupling matrix elements following the instructions
    reported in the EPIq site
    <https://the-epiq-team.gitlab.io/epiq-site/docs/tutorials/step1/>.

3.  [Wannier90 code](http://www.wannier.org/){.reference .external}
    [\\(\\rightarrow\\)]{.math .notranslate .nohighlight} Identify the
    unitary transformation connecting the Bloch and the maximally-smooth
    gauge (required to interpolate the electron-phonon coupling matrix
    elements in the Wannier basis).

4.  [EPIq code](https://the-epiq-team.gitlab.io/epiq-site/){.reference
    .external} [\\(\\rightarrow\\)]{.math .notranslate .nohighlight}
    Perform the electron-phonon copuling interpolation in the Wannier
    basis and calculate physical properties including the anharmonic
    correction from SSCHA.

All these open-source software can be downloaded following the
instruction in each website.

If you want to know more about the full procedure, please take a look at
the tutorials on the [EPIq
site.](https://the-epiq-team.gitlab.io/epiq-site/docs/tutorials/tutorials/){.reference
.external}
:::
:::
:::
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
    -   Previous: [Hands-on-session 7: Calculation of the
        electron-phonon interaction and superconducting properties with
        the
        SSCHA](tutorial_07_simple_electron_phonon.html "previous chapter")
    -   Next: [Hands-on-session 9 - Thermal conductivity calculations
        with the SSCHA](tutorial_09_TC.html "next chapter")
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
source](_sources/tutorial_08_epiq.rst.txt)
:::
