---
layout: page
title: Calculation of the electron-phonon interaction and superconducting properties with the SSCHA
---

This tutorial was prepared for the [2023 SSCHA School](http://sscha.eu/Schools/2023/home/) by Ion Errea. You can see here the video os the hands-on session:

<iframe width="560" height="315" src="https://www.youtube.com/embed/oSRAx6eNFtc" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>

The material needed for this tutorial can be downloaded [here](https://ehubox.ehu.eus/s/Y48Wc8iX9Z76jqN).

<p>In this hands-on-session we will learn to calculate the electron-phonon interaction and superconductivity properties in strongly anharmonic systems combining electron-phonon matrix elements calculated within density-functional perturbation theory (DFPT), as implemented in Quantum Espresso (QE), and the anharmonic phonon frequencies and polarization vectors obtained with the SSCHA. For that purpose we will use a slightly modified version of Quantum Espresso version 5.1.0, which includes some extra features developed by us that can be used to combine the electron-phonon matrix elements with real space SSCHA force constants, calculate the <span class="math notranslate nohighlight">\(\alpha^2F(\omega)\)</span>, calculate <span class="math notranslate nohighlight">\(T_c\)</span>  with empirical equations, and solve isotropic Migdal-Eliashberg equations.</p>
<section id="calculation-of-the-electron-phonon-matrix-elements">
<h2>Calculation of the electron-phonon matrix elements<a class="headerlink" href="#calculation-of-the-electron-phonon-matrix-elements" title="Permalink to this headline"> </a></h2>
<p>As an example we will use a non-converged calculation on PdH, a strongly anharmonic superconductor, in which the Pd atoms form a fcc lattice and H atoms sit at the octahedral interstitial sites. The crystal structure is the rock-salt one, whose space group is <span class="math notranslate nohighlight">\(Fm\bar{3}m\)</span>.</p>
<p>We will perform first a harmonic phonon calculation and the calculation of the electron-phonon coupling constant for the irreducible q points in a 2x2x2 grid. In order to know how many irreducible q points are there for a given crystal one can use the kpoints.x program of QE (for instance, in the QE version distributed with the SSCHA, one can find it in qe-5.1.0_elph/PW/tools). In this case there are only three q points in the irreducible grid. For a particular q point, the calculation of the electron-phonon matrix elements is performed in three steps:</p>
<ol class="arabic simple">
<li><p>Perform a standard DFT calculation of the crystal structure</p></li>
<li><p>Perform a standard DFPT to calculate the harmonic dynamical matrix at a particular point in the grid as well as the derivative of the Kohn-Sham (KS) potential for this particular q point. The latter will be needed for the electron-phonon calculation.</p></li>
<li><p>Perform a non-self-consistent calculation of the band structure in a finer electronic k point grid, read the derivative of the KS potential, read the harmonic dynamical matrix to obtain the phonon frequencies and polarization vectors, and calculate the electron-phonon matrix elements.</p></li>
</ol>
<p>The input we will use for the standard DFT calculation that we will use is the following:</p>
<div class="highlight-fortran notranslate"><div class="highlight"><pre><span></span>&amp;control
    ! Type of calculation
    calculation     =       &#39;scf&#39;
    ! Show more details in the output
    verbosity       =       &#39;high&#39;
    ! Calculate stress tensor
    tstress         =       .true.
    ! Calculate forces
    tprnfor         =       .true.
    ! Prefix for tmp files
    prefix          =       &#39;pdh&#39;
    ! Location of the pseudopotentials
    pseudo_dir      =       &#39;./&#39;
    ! Folder for the tmp files
    outdir          =       &#39;./tmp&#39;
/
&amp;system
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
    occupations     =       &#39;smearing&#39;
    ! Type of smearing
    smearing        =       &#39;mp&#39;
    ! Broadening of the smearing in Ry
    degauss         =       0.020
/
&amp;electrons
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
</pre></div>
</div>
<p>The reader is referred to the official  <a class="reference external" href="https://www.quantum-espresso.org/Doc/INPUT_PW.html">QE guide</a> to check the all the details about the input files used, even if a short description is provided here. The pseudopotentials can be found in the folder 07_simple_electron_phonon/pseudos/ .</p>
<p>If this input file was named as pw.in, we would run QE as follows:</p>
<div class="highlight-console notranslate"><div class="highlight"><pre><span></span><span class="go">qe-5.1.0_elph/bin/pw.x &lt; pw.in &gt; pw.out</span>
</pre></div>
</div>
<p>This will calculate the Kohn-Sham potential, needed for the phonon and electron-phonon calculations, apart from the usual total energy and forces of the structure. Note that this is not a converged calculation. One should carefully check the convergence with respect to the cutoffs, smearing, k-point grids, etc.</p>
<p>The second step is to calculate the harmonic dynamical matrix within DFPT. This is a model input file to calculate it in the first (the <span class="math notranslate nohighlight">\(\Gamma\)</span> point) q point of the 2x2x2 grid.</p>
<div class="highlight-fortran notranslate"><div class="highlight"><pre><span></span><span class="n">Phonon</span> <span class="n">calculation</span> <span class="n">on</span> <span class="n">the</span> <span class="mi">1</span><span class="n">st</span> <span class="n">point</span> <span class="n">of</span> <span class="n">a</span> <span class="mi">2</span><span class="n">x2x2</span> <span class="n">q</span> <span class="n">point</span> <span class="n">grid</span>
<span class="p">&amp;</span><span class="n">inputph</span>
    <span class="c">! Prefix for tmp files (same as for pw.x)</span>
    <span class="n">prefix</span>           <span class="o">=</span> <span class="s1">&#39;pdh&#39;</span>
    <span class="c">! Folder for tmp files (same as for pw.x)</span>
    <span class="n">outdir</span>           <span class="o">=</span> <span class="s1">&#39;./tmp/&#39;</span>
    <span class="c">! Name of dynamical matrices calculated</span>
    <span class="n">fildyn</span>           <span class="o">=</span> <span class="s2">&quot;harmonic_dyn&quot;</span>
    <span class="c">! Mass of 1st atom type in m.a.u</span>
    <span class="n">amass</span><span class="p">(</span><span class="mi">1</span><span class="p">)</span>         <span class="o">=</span> <span class="mi">10</span><span class="mf">6.42</span>
    <span class="c">! Mass of 2nd atom type in m.a.u</span>
    <span class="n">amass</span><span class="p">(</span><span class="mi">2</span><span class="p">)</span>         <span class="o">=</span> <span class="mf">1.00794</span>
    <span class="c">! File where the derivative of the KS potential will be stored</span>
    <span class="n">fildvscf</span>         <span class="o">=</span> <span class="s1">&#39;pdh_dv&#39;</span>
    <span class="c">! Calculate the phonons in a grid nq1 x nq2 x nq3 grid of q points</span>
    <span class="n">ldisp</span>            <span class="o">=</span> <span class="p">.</span><span class="n">true</span><span class="p">.</span>
    <span class="n">nq1</span>              <span class="o">=</span> <span class="mi">2</span>
    <span class="n">nq2</span>              <span class="o">=</span> <span class="mi">2</span>
    <span class="n">nq3</span>              <span class="o">=</span> <span class="mi">2</span>
    <span class="c">! Threshold for the self-consistent loop</span>
    <span class="n">tr2_ph</span>           <span class="o">=</span> <span class="mf">1.0d-16</span>
    <span class="c">! First q point to calculate</span>
    <span class="n">start_q</span>          <span class="o">=</span> <span class="mi">1</span>
    <span class="c">! Last q point to calculate</span>
    <span class="n">last_q</span>           <span class="o">=</span> <span class="mi">1</span>
<span class="p">&amp;</span><span class="k">end</span>
</pre></div>
</div>
<p>The reader is referred to the official <a class="reference external" href="https://www.quantum-espresso.org/Doc/INPUT_PH.html">QE guide</a> for more details on the input parameters. If this input file was named as ph.in, we would run QE as follows:</p>
<div class="highlight-console notranslate"><div class="highlight"><pre><span></span><span class="go">qe-5.1.0_elph/bin/ph.x &lt; ph.in &gt; ph.out</span>
</pre></div>
</div>
<p>As an output we will obtain the dynamical matrix at the <span class="math notranslate nohighlight">\(\Gamma\)</span> point stored in the file harmonic_dyn1.</p>
<p>Once we have the dynamical matrix calculated and the derivative of the KS potential stored we can calculate the electron-phonon matrix elements using the modified version of QE. The input file is the following:</p>
<div class="highlight-fortran notranslate"><div class="highlight"><pre><span></span><span class="n">Phonon</span> <span class="n">calculation</span> <span class="n">on</span> <span class="n">the</span> <span class="mi">1</span><span class="n">st</span> <span class="n">point</span> <span class="n">of</span> <span class="n">a</span> <span class="mi">2</span><span class="n">x2x2</span> <span class="n">q</span> <span class="n">point</span> <span class="n">grid</span>
<span class="p">&amp;</span><span class="n">inputph</span>
    <span class="c">! Prefix for tmp files (same as for pw.x)</span>
    <span class="n">prefix</span>           <span class="o">=</span> <span class="s1">&#39;pdh&#39;</span>
    <span class="c">! Folder for tmp files (same as for pw.x)</span>
    <span class="n">outdir</span>           <span class="o">=</span> <span class="s1">&#39;./tmp/&#39;</span>
    <span class="c">! Name of dynamical matrices calculated</span>
    <span class="n">fildyn</span>           <span class="o">=</span> <span class="s2">&quot;harmonic_dyn&quot;</span>
    <span class="c">! Mass of 1st atom type in m.a.u</span>
    <span class="n">amass</span><span class="p">(</span><span class="mi">1</span><span class="p">)</span>         <span class="o">=</span> <span class="mi">10</span><span class="mf">6.42</span>
    <span class="c">! Mass of 2nd atom type in m.a.u</span>
    <span class="n">amass</span><span class="p">(</span><span class="mi">2</span><span class="p">)</span>         <span class="o">=</span> <span class="mf">1.00794</span>
    <span class="c">! File where the derivative of the KS potential will be stored</span>
    <span class="n">fildvscf</span>         <span class="o">=</span> <span class="s1">&#39;pdh_dv&#39;</span>
    <span class="c">! Calculate the phonons in a grid nq1 x nq2 x nq3 grid of q points</span>
    <span class="n">ldisp</span>            <span class="o">=</span> <span class="p">.</span><span class="n">true</span><span class="p">.</span>
    <span class="n">nq1</span>              <span class="o">=</span> <span class="mi">2</span>
    <span class="n">nq2</span>              <span class="o">=</span> <span class="mi">2</span>
    <span class="n">nq3</span>              <span class="o">=</span> <span class="mi">2</span>
    <span class="c">! Threshold for the self-consistent loop</span>
    <span class="n">tr2_ph</span>           <span class="o">=</span> <span class="mf">1.0d-16</span>
    <span class="c">! First q point to calculate</span>
    <span class="n">start_q</span>          <span class="o">=</span> <span class="mi">1</span>
    <span class="c">! Last q point to calculate</span>
    <span class="n">last_q</span>           <span class="o">=</span> <span class="mi">1</span>
    <span class="c">! Do not calculate dynamical matrix</span>
    <span class="n">trans</span>            <span class="o">=</span> <span class="p">.</span><span class="n">false</span><span class="p">.</span>
    <span class="c">! Type of electron-phonon interaction</span>
    <span class="n">electron_phonon</span>  <span class="o">=</span> <span class="s1">&#39;simple&#39;</span>
    <span class="c">! Minimum Gaussian broadening in Ry for the double Dirac delta</span>
    <span class="n">el_ph_sigma</span>      <span class="o">=</span> <span class="mf">0.004</span>
    <span class="c">! The number of Gaussian broadenings that will be studied</span>
    <span class="n">el_ph_nsigma</span>     <span class="o">=</span> <span class="mi">25</span>
    <span class="c">! nk1 x nk2 x nk3 is the grid for the non-scf calculation</span>
    <span class="c">! used in the electron-phonon calculations</span>
    <span class="n">nk1</span>              <span class="o">=</span> <span class="mi">20</span>
    <span class="n">nk2</span>              <span class="o">=</span> <span class="mi">20</span>
    <span class="n">nk3</span>              <span class="o">=</span> <span class="mi">20</span>
    <span class="c">! k1, k2, k3 determine whether the grid is shifted from Gamma</span>
    <span class="n">k1</span>               <span class="o">=</span> <span class="mi">1</span>
    <span class="n">k2</span>               <span class="o">=</span> <span class="mi">1</span>
    <span class="n">k3</span>               <span class="o">=</span> <span class="mi">1</span>
<span class="p">&amp;</span><span class="k">end</span>
</pre></div>
</div>
<p>The reader is referred to the official <a class="reference external" href="https://www.quantum-espresso.org/Doc/INPUT_PH.html">QE guide</a> for more details on the input parameters. If this input file was named as elph.in, we would run QE as follows:</p>
<div class="highlight-console notranslate"><div class="highlight"><pre><span></span><span class="go">qe-5.1.0_elph/bin/ph.x &lt; elph.in &gt; elph.out</span>
</pre></div>
</div>
<p>As an output we will obtain several files giving information on the phonon linewidth coming from the electron-phonon interaction and so on. Most of them can be obtained with the standard version of QE, but the modified version we are providing here prints also the ‘fildyn’.elph.d.mat.’q point number’ files, in this case harmonic_dyn1.elph.d.mat.1. This file is important for our purpose as it is necessary to combine the SSCHA dynamical matrices with the obtained electron-phonon matrix elements. What it contains explicitly is the following:</p>
<div class="math notranslate nohighlight">
\[\Delta ^{ab}(\mathbf{q}) = \frac{1}{N_{\mathrm{F}}N_{\mathbf{k}}}\sum_{n,n',\mathbf{k}}d^{a}_{n\mathbf{k},n'\mathbf{k}+\mathbf{q}}d^{b}_{n'\mathbf{k}+\mathbf{q},n\mathbf{k}}\delta(\epsilon _{n\mathbf{k}} - \epsilon _{\mathrm{F}})\delta(\epsilon _{n'\mathbf{k}+\mathbf{q}} - \epsilon _{\mathrm{F}}),\]</div>
<p>where</p>
<div class="math notranslate nohighlight">
\[d^{a}_{n\mathbf{k}, n'\mathbf{k}+\mathbf{q}} = \langle n\mathbf{k} \vert \frac{\delta V_{KS}}{\delta u ^{a}(\mathbf{q})}\vert n'\mathbf{k} + \mathbf{q} \rangle\]</div>
<p>are the electron-phonon matrix elements between differen KS states <span class="math notranslate nohighlight">\(\vert n \mathbf{k}  \rangle\)</span> (<span class="math notranslate nohighlight">\(n\)</span> is aband index and <span class="math notranslate nohighlight">\(\mathbf{k}\)</span> the wave number) of the derivative of the KS potential with respect to the Fourier transformed displacement in Cartesian basis. In the above equations lower case latin indexes (<span class="math notranslate nohighlight">\(a\)</span> …) denote both atoms in the unit cell as well as Cartesian indexes. Above <span class="math notranslate nohighlight">\(N_F\)</span> is the density of states (DOS) at the Fermi level per spin, <span class="math notranslate nohighlight">\(\epsilon _{\mathrm{F}}\)</span> is the Fermi energy, <span class="math notranslate nohighlight">\(N_{\mathbf{k}}\)</span> the number of k points in the sum. The file contains first the Gaussian broadening used in the calculation (DOS and double Dirac delta in the equation), the DOS at the Fermi level calculated with that broadening, and later prints the elements of the <span class="math notranslate nohighlight">\(\Delta^{ab}(\mathbf{q})\)</span> matrix. Then it continues with the same data for the next broadening calculated.</p>
<div class="topic">
<p class="topic-title">Exercise</p>
<ul class="simple">
<li><p>Calculate the electron-phonon matrix elements for the other q points in the 2x2x2 grid.</p></li>
</ul>
</div>
</section>
<section id="the-sscha-calculation">
<h2>The SSCHA calculation<a class="headerlink" href="#the-sscha-calculation" title="Permalink to this headline"> </a></h2>
<p>This system, even with this unconverged parameters, is extremely anharmonic and the SSCHA strongly renormalizes the phonon spectrum.</p>
<div class="topic">
<p class="topic-title">Exercise</p>
<ul class="simple">
<li><p>Perform a SSCHA calculation on a 2x2x2 supercell to obtain renormalized phonon frequencies.</p></li>
</ul>
</div>
<p>As performing the SSCHA even with this unconverged parameters may take a considerable time, we provide auxiliary SSCHA dynamical matrices obtained with few configurations in 07_simple_electron_phonon/sscha/ with the name sscha_T0.0_dyn*.</p>
</section>
<section id="combine-the-sscha-dynamical-matrices-with-the-electron-phonon-matrix-elements">
<h2>Combine the SSCHA dynamical matrices with the electron-phonon matrix elements<a class="headerlink" href="#combine-the-sscha-dynamical-matrices-with-the-electron-phonon-matrix-elements" title="Permalink to this headline"> </a></h2>
<p>We will now combine the SSCHA dynamical matrices with the electron-phonon matrix elements calculates previously. In order to do that (the reason will be apparent later) we will first Fourier transform the SSCHA dynamical matrices and create the real space SSCHA force constants. We can do that with the q2r.x code of QE. Let’s first copy the harmonic_dyn0 obtained in the phonon calculations, file that contains the list of q points in the 2x2x2 grid, to the folder where the SSCHA dynamical matrices are and let’s name it following the new notation:</p>
<div class="highlight-console notranslate"><div class="highlight"><pre><span></span><span class="go">cp {PATH_TO_HARMONIC_CALCULATION}/harmonic_dyn0 {PATH_TO_SSCHA_RESULTS}/sscha_T0.0_dyn0</span>
</pre></div>
</div>
<p>Now we are ready to perform the Fourier transform. Let’s prepare the input for q2r.x:</p>
<div class="highlight-fortran notranslate"><div class="highlight"><pre><span></span><span class="p">&amp;</span><span class="n">input</span>
  <span class="c">! Name of dynamical matrices</span>
  <span class="n">fildyn</span>  <span class="o">=</span>  <span class="s1">&#39;sscha_T0.0_dyn&#39;</span>
  <span class="c">! Type of ASR imposed</span>
  <span class="n">zasr</span>    <span class="o">=</span>  <span class="s1">&#39;crystal&#39;</span>
  <span class="c">! Name of obtained force constants</span>
  <span class="n">flfrc</span>   <span class="o">=</span>  <span class="s1">&#39;sscha_T0.0.fc&#39;</span>
<span class="o">/</span>
</pre></div>
</div>
<p>More details about the input for q2r.x can be found <a class="reference external" href="https://www.quantum-espresso.org/Doc/INPUT_Q2R.html">here</a>. If the input file was named as q2r.in, we would run it as</p>
<div class="highlight-console notranslate"><div class="highlight"><pre><span></span><span class="go">qe-5.1.0_elph/bin/q2r.x &lt; q2r.in &gt; q2r.out</span>
</pre></div>
</div>
<p>This will create the sscha_T0.0.fc file with the real space force constants.</p>
<p>Now we are ready to combine these SSCHA real space force constants with the electron-phonon matrix elements calculated with DFPT. For that the easiest thing to do is to copy the SSCHA force constants file (sscha_T0.0.fc) and the files where the <span class="math notranslate nohighlight">\(\Delta ^{ab}(\mathbf{q})\)</span> matrices were stored (harmonic_dyn*.elph.d.mat.*) to a new folder. We will use the elph_fc.x, which is written by us based on matdyn.x of QE and is not present in the QE distribution, to perform this calculation. The input for this code looks as follows:</p>
<div class="highlight-fortran notranslate"><div class="highlight"><pre><span></span><span class="n">Calculation</span> <span class="n">of</span> <span class="n">superconducting</span> <span class="n">properties</span> <span class="n">with</span> <span class="n">elph_fc</span><span class="p">.</span><span class="n">x</span>
<span class="p">&amp;</span><span class="n">input</span>
    <span class="c">! ASR type imposed</span>
    <span class="n">asr</span>          <span class="o">=</span>   <span class="s1">&#39;crystal&#39;</span>
    <span class="c">! Mass of 1st atom in m.a.u</span>
    <span class="n">amass</span><span class="p">(</span><span class="mi">1</span><span class="p">)</span>     <span class="o">=</span>    <span class="mi">10</span><span class="mf">6.42</span>
    <span class="c">! Mass of 2nd atom in m.a.u</span>
    <span class="n">amass</span><span class="p">(</span><span class="mi">2</span><span class="p">)</span>     <span class="o">=</span>    <span class="mf">1.00794</span>
    <span class="c">! File with the SSCHA FCs</span>
    <span class="n">flfrc</span>        <span class="o">=</span>   <span class="s1">&#39;sscha_T0.0.fc&#39;</span>
    <span class="c">! Prefix for the files with the elph</span>
    <span class="n">fildyn</span>       <span class="o">=</span>   <span class="s1">&#39;harmonic_dyn&#39;</span>
    <span class="c">! Number of broadenings for the phonon Dirac Delta</span>
    <span class="n">nbroad</span>       <span class="o">=</span>   <span class="mi">1</span>
    <span class="c">! Gaussian broadening for the Dirac delta in cm-1</span>
    <span class="n">minbroad</span>     <span class="o">=</span>   <span class="mi">10</span>
<span class="o">/</span>
<span class="p">.</span><span class="mi">004</span> <span class="c">! Broadening for the electrons chosen</span>
<span class="mi">3</span>    <span class="c">! Number of q points in IBZ followed by the q list (in 2pi/a) and multiplicity</span>
<span class="mf">0.000000000000000E+00</span>   <span class="mf">0.000000000000000E+00</span>   <span class="mf">0.000000000000000E+00</span> <span class="mi">1</span>
<span class="mf">0.500000000000000E+00</span>  <span class="o">-</span><span class="mf">0.500000000000000E+00</span>   <span class="mf">0.500000000000000E+00</span> <span class="mi">4</span>
<span class="mf">0.000000000000000E+00</span>  <span class="o">-</span><span class="mf">0.100000000000000E+01</span>   <span class="mf">0.000000000000000E+00</span> <span class="mi">3</span>
</pre></div>
</div>
<p>If the input file was named as elph_fc.in, we would run the code as follows:</p>
<div class="highlight-console notranslate"><div class="highlight"><pre><span></span><span class="go">qe-5.1.0_elph/bin/elph_fc.x &lt; elph_fc.in &gt; elph_fc.out</span>
</pre></div>
</div>
<p>The code Fourier transforms the SSCHA real space force constants to the list of q points provided in order to obtain the dynamica matrices at these points, and combines them with the <span class="math notranslate nohighlight">\(\Delta^{ab}(\mathbf{q})\)</span> matrices to calculate the Eliashberg function <span class="math notranslate nohighlight">\(\alpha^2F(\omega)\)</span> as</p>
<div class="math notranslate nohighlight">
\[\alpha^2F(\omega) = \frac{1}{N_{\mathbf{q}}}\sum _{ab\mathbf{q}\mu} \frac{e^a_{\mu}(\mathbf{q})\Delta ^{ab}(\mathbf{q})e^b_{\mu}(\mathbf{q})^*}{2\omega_{\mu}(\mathbf{q})\sqrt{m_am_b}} \delta(\omega-\omega_{\mu}(\mathbf{q})),\]</div>
<p>where <span class="math notranslate nohighlight">\(\omega_{\mu}(\mathbf{q})\)</span> and <span class="math notranslate nohighlight">\(e^a_{\mu}(\mathbf{q})\)</span> are, respectively, the phonon frequencies and polarization vectors obtained diagonalizing the Fourier interpolated SSCHA force constants at point q, and <span class="math notranslate nohighlight">\(m_a\)</span> the masses of the atoms. The Dirac delta on the equation is approximated with a Gaussian of 10 cm-1 broadening, following the input parameter. In the output elph_fc.out the code gives the SSCHA dynamical matrix at each q point, the phonon linewidth (HWHM), the contribution to the electron-phonon coupling of each mode, etc. Note that the code skips the Gamma point and does not consider it in the calculation. The reason is that the equation used at this point is divergent (see discussion in Appendix C in <a class="reference external" href="https://arxiv.org/abs/2303.02621">arXiv:2303.02621</a>). The code also calculates the total electron-phonon coupling constant <span class="math notranslate nohighlight">\(\lambda\)</span> and <span class="math notranslate nohighlight">\(\omega_{log}\)</span>. It also gives the value of the the superconducting critical temperature calculated for different values of the Coulomb pseudopotential <span class="math notranslate nohighlight">\(\mu^*\)</span> within the semiempirical McMillan and Allen-Dynes formulas.</p>
<p>The code also prints the calculated Eliashberg function, phonon density of states, partial electron-phonon coupling constant, as well as the projection of both the phonon DOS and Eliashberg function on different atoms. We can do for example plots like this one with this data:</p>
<figure class="align-default" id="id3">
<img alt="../figures_07/a2F_PdH.pdf" src="../figures_07/a2F_PdH.pdf" />
<figcaption>
<p><span class="caption-number">Fig. 15 </span><span class="caption-text">Figure with the <span class="math notranslate nohighlight">\(\alpha^2F(\omega)\)</span> Eliashberg function, its partial contributions from Pd and H atoms. The partial contribution to the electron-phonon coupling constant is also plotted, <span class="math notranslate nohighlight">\(\lambda(\omega)\)</span>.</span><a class="headerlink" href="#id3" title="Permalink to this image"> </a></p>
</figcaption>
</figure>
<div class="topic">
<p class="topic-title">Exercise</p>
<ul class="simple">
<li><p>Make a figure as the one above but with a different smearing for the double Dirac delta on the electronic states. Note that the electron-phonon matrix elements were calculated for smearings proportional to 0.04 Ry.</p></li>
</ul>
</div>
<p>The fact that the elph_fc.x code works with real space SSCHA force constants allows us to combine the calculation of the electron-phonon matrix elements in a small supercell with electron-phonon matrix elements in a finer q point grid trivially.</p>
<div class="topic">
<p class="topic-title">Exercise</p>
<ul class="simple">
<li><p>Calculate the electron-phonon coupling constant using the electron-phonon matrix elements calculated in a 4x4x4 q point grid combining it with the SSCHA force constants on a 2x2x2 supercell. For that use the data in 07_simple_electron_phonon/elph_matrix_elements_444 .</p></li>
</ul>
</div>
</section>
<section id="solution-of-isotropic-migdal-eliashberg-equations">
<h2>Solution of isotropic Migdal-Eliashberg equations<a class="headerlink" href="#solution-of-isotropic-migdal-eliashberg-equations" title="Permalink to this headline"> </a></h2>
<p>Once the Eliashberg function <span class="math notranslate nohighlight">\(\alpha^2F(\omega)\)</span> has been calculated combining the SSCHA and the electron-phonon matrix elements, one can easily solve isotropic Migdal-Eliashberg equations, which are not semiempirical. We provide a utility to perform this calculation as well, ME.x. This is the input file that needs to be prepared:</p>
<div class="highlight-fortran notranslate"><div class="highlight"><pre><span></span><span class="p">&amp;</span><span class="n">inputme</span>
     <span class="c">! First temperature for the calculation</span>
     <span class="n">t_first</span>      <span class="o">=</span>  <span class="mf">1.0</span>
     <span class="c">! Last temperature for the calculation</span>
     <span class="n">t_last</span>       <span class="o">=</span>  <span class="mi">8</span><span class="mf">0.00</span>
     <span class="c">! Total number of temperatures</span>
     <span class="n">t_number</span>     <span class="o">=</span>  <span class="mi">80</span>
     <span class="c">! Cutoff for Matsubara frequencies in Ha</span>
     <span class="n">wc_cutoff</span>    <span class="o">=</span>  <span class="mf">0.05</span>
     <span class="c">! Initial guess for the gap in meV</span>
     <span class="n">gap_guess</span>    <span class="o">=</span>  <span class="mi">1</span><span class="mf">0.</span>
     <span class="c">! File with the Eliashberg function</span>
     <span class="n">a2f_filename</span> <span class="o">=</span>  <span class="s2">&quot;a2F.10.dat&quot;</span>
     <span class="c">! mu*</span>
     <span class="n">mu_star</span>      <span class="o">=</span>  <span class="mf">0.10</span>
<span class="o">/</span>
</pre></div>
</div>
<p>Note that the cutoff for the Matsubara frequencies should be around 10 times the highest phonon frequency, not bigger. If this input file was named as ME.in, we would run the code as follows:</p>
<div class="highlight-console notranslate"><div class="highlight"><pre><span></span><span class="go">qe-5.1.0_elph/bin/ME.x &lt; ME.in &gt; ME.out</span>
</pre></div>
</div>
<p>In the output the code will calculate the superconducting gap as a function of temperature. In the provided file t_gap.dat the gap as a function of temperature is provided. This can be used to generate plots like this one to estimate the critical temperature from the temperature at which the gap closes.</p>
<figure class="align-default" id="id4">
<img alt="../figures_07/gap_PdH.pdf" src="../figures_07/gap_PdH.pdf" />
<figcaption>
<p><span class="caption-number">Fig. 16 </span><span class="caption-text">Superconducting gap as a function of temperature.</span><a class="headerlink" href="#id4" title="Permalink to this image"> </a></p>
</figcaption>
</figure>
</section>
<section id="important-remarks">
<h2>Important remarks<a class="headerlink" href="#important-remarks" title="Permalink to this headline"> </a></h2>
<p>The calculations above are not converged and are just meant to illustrate the use of the several codes. In a proper calculation one should take into account the following points:</p>
<ol class="arabic simple">
<li><p>The electron-phonon coupling constant needs to be converged with the number of k points for the electrons and the smearing used for the Double delta. These parameters enter into the equation of the <span class="math notranslate nohighlight">\(\Delta^{ab}(\mathbf{q})\)</span> matrices. The typical thing is to calculate <span class="math notranslate nohighlight">\(\lambda\)</span> for different k point grids as a function of the smearing and see for which low value of the smearing <span class="math notranslate nohighlight">\(\lambda\)</span> plateaus. Note that the physical limit is the one with infinite number of k points and 0 smearing.</p></li>
<li><p>In this example we have used auxiliary SSCHA dynamical matrices to incorporate anharmonic effects into the calculation of electron-phonon properties. However, we could have also used those dynamical matrices that come from the Hessian of the SSCHA free energy. Sometimes the differences are minor, but in other cases it may be important. In this cases the best approach is to calculate the Eliashberg function with the full spectral function as described recently in <a class="reference external" href="https://arxiv.org/abs/2303.07962">arXiv:2303.07962</a>.</p></li>
</ol>
</section>

