---
layout: page
title: Setting up automatic SSCHA simulations on clusters
---

This tutorial will cover more advanced code features, like the SSCHA code’s interoperability with a high-performance computer (HPC). The tutorial is divided into two sections. In the first section, we will perform a free energy minimization manually; then we will learn how to automatize the interaction with a cluster to run ab initio calculations automatically.  

This tutorial was prepared for the [2023 SSCHA School](http://sscha.eu/Schools/2023/home/) by Lorenzo Monacelli. You can see here the video os the hands-on session:

<iframe width="560" height="315" src="https://www.youtube.com/embed/Z0K_9zrGcJI" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>

The material needed for this tutorial can be downloaded [here](https://ehubox.ehu.eus/s/Y48Wc8iX9Z76jqN).

<section id="manual-submission">
<h2>Manual submission<a class="headerlink" href="#manual-submission" title="Permalink to this headline"> </a></h2>
<dl class="simple">
<dt>The SSCHA calculation comprises three main steps iterated until convergence:</dt><dd><ol class="arabic simple">
<li><p>The generation of a random ensemble of ionic configurations</p></li>
<li><p>Calculations of energies and forces on the ensemble</p></li>
<li><p>The SSCHA free energy minimization</p></li>
</ol>
</dd>
</dl>
<p>In the first hands-on session, you configured the code to do these iterations automatically.
Thanks to the ASE EMT force field, the code can automatically compute energies, forces, and stress tensors without user interaction.</p>
<p>However, if you need to compute energies and forces from an <em>ab initio</em> calculation like DFT, you may want to
run the DFT code on a different machine, like a cluster.</p>
<p>You can use the manual submission if you want more control over the procedure.</p>
<p>We will compute the sulfur hydride (superconductor with <span class="math notranslate nohighlight">\(T_c = 203\)</span> K), using a DFT code like quantum Espresso to calculate energy and forces.</p>
<p>The harmonic phonons (computed using quantum Espresso) is provided in the directory <strong>02_manual_submission</strong>, where you can find the input and output files of the quantum espresso calculation to calculate the harmonic phonons, and the dynamical matrices, named
dyn_h3s_harmonic_1, dyn_h3s_harmonic_2 and dyn_h3s_harmonic_3.</p>
<p>They respect the naming convention so that each file contains a different q point: since we are using a 2x2x2 mesh to sample the Brillouin zone of phonons, the different q points are ordered in three separate files, each one grouping the <em>star</em> of q (the q points related by symmetry operations).</p>
<p>We start by plotting the dispersion of the harmonic dynamical matrix.
Please write in a file the following script and run it.</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span><span class="o">,</span> <span class="nn">cellconstructor.Phonons</span>
<span class="kn">import</span> <span class="nn">cellconstructor.ForceTensor</span>
<span class="kn">import</span> <span class="nn">ase</span><span class="o">,</span> <span class="nn">ase.dft</span>

<span class="kn">import</span> <span class="nn">matplotlib.pyplot</span> <span class="k">as</span> <span class="nn">plt</span>
<span class="kn">import</span> <span class="nn">numpy</span> <span class="k">as</span> <span class="nn">np</span>

<span class="n">dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="s2">&quot;dyn_h3s_harmonic_&quot;</span><span class="p">,</span> <span class="mi">3</span><span class="p">)</span> <span class="c1"># Load 3 files</span>

<span class="n">PATH</span> <span class="o">=</span> <span class="s2">&quot;GHNPGN&quot;</span>
<span class="n">N_POINTS</span> <span class="o">=</span> <span class="mi">1000</span>

<span class="c1"># Use ASE to get the q points from the path</span>
<span class="n">band_path</span> <span class="o">=</span> <span class="n">ase</span><span class="o">.</span><span class="n">dft</span><span class="o">.</span><span class="n">kpoints</span><span class="o">.</span><span class="n">bandpath</span><span class="p">(</span><span class="n">PATH</span><span class="p">,</span>
    <span class="n">dyn</span><span class="o">.</span><span class="n">structure</span><span class="o">.</span><span class="n">unit_cell</span><span class="p">,</span>
    <span class="n">N_POINTS</span><span class="p">)</span>

<span class="c1"># Get the q points in cartesian coordinates</span>
<span class="n">q_path</span> <span class="o">=</span> <span class="n">band_path</span><span class="o">.</span><span class="n">cartesian_kpts</span><span class="p">()</span>

<span class="c1"># Get the values of x axis and labels for plotting the band path</span>
<span class="n">x_axis</span><span class="p">,</span> <span class="n">xticks</span><span class="p">,</span> <span class="n">xlabels</span> <span class="o">=</span> <span class="n">band_path</span><span class="o">.</span><span class="n">get_linear_kpoint_axis</span><span class="p">()</span>

<span class="c1"># Perform the interpolation of the dynamical matrix along the q_path</span>
<span class="n">frequencies</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">ForceTensor</span><span class="o">.</span><span class="n">get_phonons_in_qpath</span><span class="p">(</span><span class="n">dyn</span><span class="p">,</span> <span class="n">q_path</span><span class="p">)</span>

<span class="c1"># Plot the dispersion</span>
<span class="n">fig</span> <span class="o">=</span> <span class="n">plt</span><span class="o">.</span><span class="n">figure</span><span class="p">()</span>
<span class="n">ax</span> <span class="o">=</span> <span class="n">plt</span><span class="o">.</span><span class="n">gca</span><span class="p">()</span>
<span class="n">ax</span><span class="o">.</span><span class="n">set_title</span><span class="p">(</span><span class="s2">&quot;Harmonic H3S Phonon dispersion&quot;</span><span class="p">)</span>
<span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">frequencies</span><span class="o">.</span><span class="n">shape</span><span class="p">[</span><span class="o">-</span><span class="mi">1</span><span class="p">]):</span>
   <span class="n">ax</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="n">x_axis</span><span class="p">,</span> <span class="n">frequencies</span><span class="p">[:,</span> <span class="n">i</span><span class="p">],</span> <span class="n">color</span> <span class="o">=</span> <span class="s1">&#39;r&#39;</span><span class="p">)</span>


<span class="k">for</span> <span class="n">x</span> <span class="ow">in</span> <span class="n">xticks</span><span class="p">:</span>
   <span class="n">ax</span><span class="o">.</span><span class="n">axvline</span><span class="p">(</span><span class="n">x</span><span class="p">,</span> <span class="mi">0</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="n">color</span><span class="o">=</span><span class="s1">&#39;k&#39;</span><span class="p">,</span> <span class="n">lw</span><span class="o">=</span><span class="mf">0.4</span><span class="p">)</span> <span class="c1"># Plot vertical lines for each high-symmetry point</span>

<span class="c1"># Set the labels of the axis as the Brilluin zone letters</span>
<span class="n">ax</span><span class="o">.</span><span class="n">set_xticks</span><span class="p">(</span><span class="n">xticks</span><span class="p">)</span>
<span class="n">ax</span><span class="o">.</span><span class="n">set_xticklabels</span><span class="p">(</span><span class="n">xlabels</span><span class="p">)</span>

<span class="n">ax</span><span class="o">.</span><span class="n">set_ylabel</span><span class="p">(</span><span class="s2">&quot;Frequency [cm-1]&quot;</span><span class="p">)</span>
<span class="n">ax</span><span class="o">.</span><span class="n">set_xlabel</span><span class="p">(</span><span class="s2">&quot;q-path&quot;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">tight_layout</span><span class="p">()</span>
<span class="n">plt</span><span class="o">.</span><span class="n">savefig</span><span class="p">(</span><span class="s2">&quot;harmonic_h3s_dispersion.png&quot;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">show</span><span class="p">()</span>
</pre></div>
</div>
<p>You should see the figure <a class="reference internal" href="#harmonic-disp"><span class="std std-ref">Dispersion of the harmonic phonons of H3S</span></a>.</p>
<figure class="align-center" id="id1">
<span id="harmonic-disp"></span><a class="reference internal image-reference" href="../figures_02/harmonic_h3s_dispersion.png"><img alt="Phonon dispersion with imaginary modes." src="../figures_02/harmonic_h3s_dispersion.png" style="width: 50%;" /></a>
<figcaption>
<p><span class="caption-number">Fig. 6 </span><span class="caption-text">Dispersion of the harmonic phonons of H3S</span><a class="headerlink" href="#id1" title="Permalink to this image"> </a></p>
</figcaption>
</figure>
<p>The dispersion presents imaginary phonons throughout most of the Brillouin zone.
To start the SSCHA, we need a <strong>positive</strong> definite dynamical matrix.
Since the starting point for the SSCHA does not matter, we may flip the phonons to be positive:</p>
<div class="math notranslate nohighlight">
\[\Phi_{ab} = \sum_\mu \sqrt{m_am_b} \left|\omega_\mu\right|^2 e_\mu^a e_\mu^b\]</div>
<p>where <span class="math notranslate nohighlight">\(m_a\)</span> is the mass of the a-th atom, <span class="math notranslate nohighlight">\(\omega_\mu\)</span> is the frequency of the dynamical matrix, and <span class="math notranslate nohighlight">\(e_\mu\)</span> is the corresponding eigenvector. This operation can be performed with the command</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="n">dyn</span><span class="o">.</span><span class="n">ForcePositiveDefinite</span><span class="p">()</span>
</pre></div>
</div>
<p>and save the results into <code class="docutils literal notranslate"><span class="pre">start_sscha1</span></code>, <code class="docutils literal notranslate"><span class="pre">start_sscha2</span></code>, and <code class="docutils literal notranslate"><span class="pre">start_sscha3</span></code> with</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="n">dyn</span><span class="o">.</span><span class="n">save_qe</span><span class="p">(</span><span class="s2">&quot;start_sscha&quot;</span><span class="p">)</span>
</pre></div>
</div>
<div class="topic">
<p class="topic-title">Exercise</p>
<p>Plot the phonon dispersion of the positive definite dynamical matrix obtained in this way.
Save the resulting dynamical matrix as ‘start_sscha’ to continue with the following section.</p>
</div>
<section id="ensemble-generation">
<h3>Ensemble generation<a class="headerlink" href="#ensemble-generation" title="Permalink to this headline"> </a></h3>
<p>Now that we have a good starting point for the dynamical matrix, we are ready to
generate the first ensemble to start the free energy optimization.
Here is a script to generate the ensemble.</p>
<p>The following script supposes that you saved the dynamical matrix after enforcing them to be positive definite as “start_sscha”. However, you can edit the script to read the harmonic dynamical matrices and impose the positiveness within the same script.</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span><span class="o">,</span> <span class="nn">cellconstructor.Phonons</span>
<span class="kn">import</span> <span class="nn">sscha</span><span class="o">,</span> <span class="nn">sscha.Ensemble</span>
<span class="kn">import</span> <span class="nn">numpy</span> <span class="k">as</span> <span class="nn">np</span>

<span class="c1"># Fix the seed so that we all generate the same ensemble</span>
<span class="n">np</span><span class="o">.</span><span class="n">random</span><span class="o">.</span><span class="n">seed</span><span class="p">(</span><span class="mi">0</span><span class="p">)</span>

<span class="c1"># Load the dynamical matrix</span>
<span class="n">dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="s2">&quot;start_sscha&quot;</span><span class="p">,</span> <span class="n">nqirr</span><span class="o">=</span><span class="mi">3</span><span class="p">)</span>

<span class="c1">#[ apply here the needed changes to dyn ]</span>

<span class="c1"># Prepare the ensemble</span>
<span class="n">temperature</span> <span class="o">=</span> <span class="mi">300</span> <span class="c1"># 300 K</span>
<span class="n">ensemble</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Ensemble</span><span class="o">.</span><span class="n">Ensemble</span><span class="p">(</span><span class="n">dyn</span><span class="p">,</span> <span class="n">temperature</span><span class="p">)</span>

<span class="c1"># Generate the ensemble</span>
<span class="n">number_of_configurations</span> <span class="o">=</span> <span class="mi">10</span>
<span class="n">ensemble</span><span class="o">.</span><span class="n">generate</span><span class="p">(</span><span class="n">number_of_configurations</span><span class="p">)</span>

<span class="c1"># Save the ensemble into a directory</span>
<span class="n">save_directory</span> <span class="o">=</span> <span class="s2">&quot;data&quot;</span>
<span class="n">population_id</span> <span class="o">=</span> <span class="mi">1</span>
<span class="n">ensemble</span><span class="o">.</span><span class="n">save</span><span class="p">(</span><span class="n">save_directory</span><span class="p">,</span> <span class="n">population_id</span><span class="p">)</span>
</pre></div>
</div>
<p>If you try to run the code, you can face an error telling you that the dynamical matrix does not satisfy the acoustic sum rule (ASR).
This occurs because quantum Espresso does not impose the ASR by default. However, we can enforce the acoustic sum rule with the following:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="n">dyn</span><span class="o">.</span><span class="n">Symmetrize</span><span class="p">()</span>
</pre></div>
</div>
<p>Besides the ASR, this function will also impose all the symmetries on the dynamical matrix, ensuring it is correct.</p>
<div class="topic">
<p class="topic-title">Exercise</p>
<p>Impose the acoustic sum rule and the symmetries and generate the ensemble.
Either add this after loading the dynamical matrix or do it once overriding the ‘start_sscha’ files.</p>
</div>
</section>
<section id="calculation-of-energies-and-forces">
<h3>Calculation of energies and forces<a class="headerlink" href="#calculation-of-energies-and-forces" title="Permalink to this headline"> </a></h3>
<p>Very good; if you imposed the sum rule correctly, the ensemble should have been correctly generated.
The script should have created the <em>data</em> directory and two sets of dynamical matrices:</p>
<ol class="arabic simple">
<li><p>dyn_start_population1_x</p></li>
<li><p>dyn_end_population1_x</p></li>
</ol>
<p>where x goes from 1 to 3. These are the same dynamical matrix as the original one.
In particular, dyn_start is the dynamical matrix used to generate the ensemble, and dyn_end is the final dynamical matrix after the free energy optimization. Since we did not run the sscha, they are the same.</p>
<p>If we look inside the <em>data</em> directory, we find:</p>
<ol class="arabic simple">
<li><p>energies_supercell_population1.dat</p></li>
<li><p>scf_population1_x.dat</p></li>
<li><p>u_population1_x.dat</p></li>
</ol>
<p>where x counts from 1 to the total number of configurations, the energies_supercell file contains any structure’s total DFT energy (in Ry).
Since we have not yet performed DFT calculations, it is full of 0s.</p>
<p><em>u_population1_x.dat</em> files contain the cartesian displacements of each atom in the supercell with respect to the average position.
We will not touch this file, but the sscha uses it to load the ensemble much faster when we have many configurations and big systems.</p>
<p>The last files are the <em>scf_population1_x.dat</em>, containing the ionic positions, including the atomic type, in Cartesian coordinates.</p>
<p>This file contains the structure in the supercell; it is already in the standard quantum espresso format, so you can attach this text to the header file of the quantum espresso input to have a complete input file for this structure.
However, you can easily manipulate this file to adapt it to your favorite programs, like VASP, ABINIT, SIESTA, CP2K, CASTEP, or any other.</p>
<p>You can visualize a structure using ASE and Cellconstructor:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">ase</span><span class="o">,</span> <span class="nn">ase.visualize</span>
<span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span><span class="o">,</span> <span class="nn">cellconstructor.Structure</span>

<span class="n">struct</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Structure</span><span class="o">.</span><span class="n">Structure</span><span class="p">()</span>
<span class="n">struct</span><span class="o">.</span><span class="n">read_scf</span><span class="p">(</span><span class="s2">&quot;data/scf_population1_1.dat&quot;</span><span class="p">)</span>
<span class="n">ase_struct</span> <span class="o">=</span> <span class="n">struct</span><span class="o">.</span><span class="n">get_ase_atoms</span><span class="p">()</span>
<span class="n">ase</span><span class="o">.</span><span class="n">visualize</span><span class="o">.</span><span class="n">view</span><span class="p">(</span><span class="n">ase_struct</span><span class="p">)</span>
</pre></div>
</div>
<p>Indeed, using the same trick, you can export the structure in any file format that ASE support, including input files for different programs mentioned above.</p>
<p>Here, we will use quantum Espresso. The header file for the quantum espresso calculation is in <em>espresso_header.pwi</em>.
Remember that the configurations are in the supercell, so the number of atoms (here 32 instead of 4) and any extensive parameter like the k-point mesh should be rescaled accordingly. Here we employ an 8x8x8 k-mesh for the electronic calculation, while to compute the harmonic phonons with a unit cell calculation, we use a 16x16x16 k-mesh since the sscha configurations are 2x2x2 bigger than the original one, and thus the Brillouin zone is a factor 0.5x0.5x0.5 smaller.</p>
<p>You can append each scf file to this header to get the espresso input.</p>
<p>We have only ten configurations; in production runs, using at least hundreds of configurations per ensemble is appropriate.
Therefore, it is impractical to create the input file for each of them manually.</p>
<div class="highlight-bash notranslate"><div class="highlight"><pre><span></span><span class="ch">#!/bin/bash</span>

<span class="nv">HEADER_FILE</span><span class="o">=</span>espresso_header.pwi
<span class="nv">DATA_DIR</span><span class="o">=</span>data
<span class="nv">POPULATION</span><span class="o">=</span><span class="m">1</span>

<span class="c1"># Define a directory in which to save all the input files</span>
<span class="nv">TARGET_DIRECTORY</span><span class="o">=</span><span class="nv">$DATA_DIR</span>/input_files_population<span class="nv">$POPULATION</span>

mkdir -p <span class="nv">$TARGET_DIRECTORY</span>

<span class="k">for</span> file <span class="k">in</span> <span class="sb">`</span>ls <span class="nv">$DATA_DIR</span>/scf_population<span class="si">${</span><span class="nv">POPULATION</span><span class="si">}</span>*.dat<span class="sb">`</span>
<span class="k">do</span>
    <span class="c1"># Extract the configuration index</span>
    <span class="c1"># (the grep command returns only the expression</span>
    <span class="c1">#  that matches the regular expression from the file name)</span>
    <span class="nv">index</span><span class="o">=</span><span class="sb">`</span><span class="nb">echo</span> <span class="nv">$file</span> <span class="p">|</span> grep -oP <span class="s1">&#39;(?&lt;=population1_).*(?=\.dat)&#39;</span><span class="sb">`</span>

    <span class="nv">target_input_file</span><span class="o">=</span><span class="nv">$TARGET_DIRECTORY</span>/structure_<span class="si">${</span><span class="nv">index</span><span class="si">}</span>.pwi

    <span class="c1"># Copy the template header file</span>
    cp <span class="nv">$HEADER_FILE</span> <span class="nv">$target_input_file</span>

    <span class="c1"># Attach after the header the structure</span>
    cat <span class="nv">$file</span> &gt;&gt; <span class="nv">$target_input_file</span>
<span class="k">done</span>
</pre></div>
</div>
<p>Executing this script, you have created a directory inside the data dir called <em>input_files_population1</em>
which contains all the input files for quantum Espresso.</p>
<p>You can run these with your own laptop if you have a good computer.
However, the calculation is computationally demanding: each configuration contains plenty of atoms and no symmetries at all,
as they are snapshots of the quantum/thermal motion of the nuclei.
The alternative is to copy these files on a cluster and submit a calculation there.</p>
<p>The espresso files are run with the command</p>
<div class="highlight-bash notranslate"><div class="highlight"><pre><span></span>mpirun -np NPROC pw.x -i input_file.pwi &gt; output_file.pwo
</pre></div>
</div>
<p>where NPROC is the number of processors in which we want to run.
<strong>Remember to copy the pseudopotential in the same directory where you run the pw.x executable.</strong></p>
<p>However, we skip this part now (try it yourself later!)</p>
<p>We provide the output files in the folder <em>output_espresso</em></p>
<p>Once we have the output files from Espresso, we need to save the energies, forces, and stress tensors in the ensemble directory.</p>
<div class="highlight-bash notranslate"><div class="highlight"><pre><span></span><span class="ch">#!/bin/bash</span>

<span class="nv">N_CONFIGS</span><span class="o">=</span><span class="m">10</span>
<span class="nv">POPULATION</span><span class="o">=</span><span class="m">1</span>
<span class="nv">PATH_TO_DIR</span><span class="o">=</span><span class="s2">&quot;data/output_espresso&quot;</span>
<span class="nv">N_ATOMS</span><span class="o">=</span><span class="m">32</span>


<span class="nv">ENERGY_FILE</span><span class="o">=</span><span class="s2">&quot;data/energies_supercell_population</span><span class="si">${</span><span class="nv">POPULATION</span><span class="si">}</span><span class="s2">.dat&quot;</span>

<span class="c1"># Clear the energy file</span>
rm -rf <span class="nv">$ENERGY_FILE</span>

<span class="k">for</span> i <span class="k">in</span> <span class="sb">`</span>seq <span class="m">1</span> <span class="m">10</span><span class="sb">`</span>
<span class="k">do</span>
    <span class="nv">filename</span><span class="o">=</span><span class="si">${</span><span class="nv">PATH_TO_DIR</span><span class="si">}</span>/structure_<span class="nv">$i</span>.pwo
    <span class="nv">force_file</span><span class="o">=</span>data/forces_population<span class="si">${</span><span class="nv">POPULATION</span><span class="si">}</span>_<span class="nv">$i</span>.dat
    <span class="nv">stress_file</span><span class="o">=</span>data/pressures_population<span class="si">${</span><span class="nv">POPULATION</span><span class="si">}</span>_<span class="nv">$i</span>.dat

    <span class="c1"># Get the total energy</span>
    grep ! <span class="nv">$filename</span> <span class="p">|</span> awk <span class="s1">&#39;{print $5}&#39;</span> &gt;&gt; <span class="nv">$ENERGY_FILE</span>
    grep force <span class="nv">$filename</span> <span class="p">|</span> grep atom <span class="p">|</span> awk <span class="s1">&#39;{print $7, $8, $9}&#39;</span> &gt; <span class="nv">$force_file</span>
    grep <span class="s2">&quot;total   stress&quot;</span> <span class="nv">$filename</span> -A3 <span class="p">|</span> tail -n +2 <span class="p">|</span> awk <span class="s1">&#39;{print $1, $2, $3}&#39;</span> &gt; <span class="nv">$stress_file</span>
<span class="k">done</span>
</pre></div>
</div>
<p>This script works specifically for quantum Espresso. It extracts energy, forces, and the stress tensor and fills the files <em>data/energies_supercell_population1.dat</em>, <em>forces_population1_X.dat</em>, and <em>pressures_population1_X.dat</em> with the results obtained from the output file of quantum Espresso.</p>
<p>The units of measurement are</p>
<ol class="arabic simple">
<li><p>Ry for the energy (in the supercell)</p></li>
<li><p>Ry/Bohr for the forces</p></li>
<li><p>Ry/Bohr^3 for the stress tensor</p></li>
</ol>
<p>Here, we do not need conversion, as these are the default units quantum Espresso gives. However, remember to convert correctly to these units if you use a different program, like VASP.</p>
</section>
<section id="free-energy-minimization">
<h3>Free energy minimization<a class="headerlink" href="#free-energy-minimization" title="Permalink to this headline"> </a></h3>
<p>We have the ensemble ready to be loaded back into the Python script and start a minimization.
This is done with the following scripts</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">sscha</span><span class="o">,</span> <span class="nn">sscha.Ensemble</span><span class="o">,</span> <span class="nn">sscha.SchaMinimizer</span>
<span class="kn">import</span> <span class="nn">sscha.Utilities</span>
<span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span><span class="o">,</span> <span class="nn">cellconstructor.Phonons</span>

<span class="n">POPULATION</span> <span class="o">=</span> <span class="mi">1</span>

<span class="n">dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="s2">&quot;start_sscha&quot;</span><span class="p">,</span> <span class="mi">3</span><span class="p">)</span>
<span class="n">ensemble</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Ensemble</span><span class="o">.</span><span class="n">Ensemble</span><span class="p">(</span><span class="n">dyn</span><span class="p">,</span> <span class="mi">0</span><span class="p">)</span>
<span class="n">ensemble</span><span class="o">.</span><span class="n">load</span><span class="p">(</span><span class="s2">&quot;data&quot;</span><span class="p">,</span> <span class="n">population</span> <span class="o">=</span> <span class="n">POPULATION</span><span class="p">,</span> <span class="n">N</span> <span class="o">=</span> <span class="mi">10</span><span class="p">)</span>

<span class="n">minim</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">SchaMinimizer</span><span class="o">.</span><span class="n">SSCHA_Minimizer</span><span class="p">(</span><span class="n">ensemble</span><span class="p">)</span>
<span class="n">minim</span><span class="o">.</span><span class="n">init</span><span class="p">()</span>

<span class="c1"># Save the minimization details</span>
<span class="n">ioinfo</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Utilities</span><span class="o">.</span><span class="n">IOInfo</span><span class="p">()</span>
<span class="n">ioinfo</span><span class="o">.</span><span class="n">SetupSaving</span><span class="p">(</span><span class="s2">&quot;minim_</span><span class="si">{}</span><span class="s2">&quot;</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="n">POPULATION</span><span class="p">))</span>


<span class="n">minim</span><span class="o">.</span><span class="n">run</span><span class="p">(</span><span class="n">custom_function_post</span> <span class="o">=</span> <span class="n">ioinfo</span><span class="o">.</span><span class="n">CFP_SaveAll</span><span class="p">)</span>
<span class="n">minim</span><span class="o">.</span><span class="n">finalize</span><span class="p">()</span>
<span class="n">minim</span><span class="o">.</span><span class="n">dyn</span><span class="o">.</span><span class="n">save_qe</span><span class="p">(</span><span class="s2">&quot;final_sscha_dyn_population</span><span class="si">{}</span><span class="s2">_&quot;</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="n">POPULATION</span><span class="p">))</span>
</pre></div>
</div>
<p>You can plot the results of the minimization with</p>
<div class="highlight-console notranslate"><div class="highlight"><pre><span></span><span class="go">sscha-plot-data.py minim_1</span>
</pre></div>
</div>
<p><strong>Congratulations!</strong> You run your first completely manual SSCHA run.</p>
<p>The output file informs us that minimization ended because the ensemble is <em>out of the stochastic criteria</em>. This means that the dynamical matrix
changed a sufficient amount that the original ensemble was not good enough anymore to describe the free energy of the new dynamical matrix; therefore, a new ensemble should be extracted.</p>
<p>In the early days of the SSCHA, this procedure should have been iterated repeatedly until convergence.
Nowadays, we have a fully automatic procedure that can automatize all these steps configuring the ssh connection to
a cluster.</p>
</section>
</section>
<section id="automatic-submission-with-a-cluster">
<h2>Automatic submission with a cluster<a class="headerlink" href="#automatic-submission-with-a-cluster" title="Permalink to this headline"> </a></h2>
<p>In the previous section, you made all the steps to run a sscha calculation manually.
This consists of iterating through the following steps:</p>
<ol class="arabic simple">
<li><p>generating the input files for Espresso,</p></li>
<li><p>transferring them to a cluster,</p></li>
<li><p>submitting the calculations,</p></li>
<li><p>retrieving the outputs,</p></li>
<li><p>reload the ensemble</p></li>
<li><p>run the free energy minimization</p></li>
</ol>
<p>In this section, we learn how to automatize these passages.
We must set up the interaction between the SSCHA library and the HPC cluster running the DFT calculations.
As of June 2023, this automatic interaction is only supported for quantum Espresso and SLURM-based clusters. However, writing plugins to support different DFT codes and cluster schedulers should be easy.</p>
<p>The configuration of the DFT parameter has been introduced in the previous hands-on session; thus, we skip and provide a file called <em>espresso_calculator.py</em>, which defines a function <em>get_h3s_calculator</em> returning the
calculator object for quantum Espresso with the input parameters for H3S.</p>
<p>We focus instead on the configuration of the cluster.
Create a new file called <em>cluster.py</em>. The following script provides an example to connect to a cluster with username <code class="docutils literal notranslate"><span class="pre">sschauser</span></code> and login node <code class="docutils literal notranslate"><span class="pre">my.beautiful.cluster.eu</span></code>:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span><span class="o">,</span> <span class="nn">cellconstructor.Phonons</span>
<span class="kn">import</span> <span class="nn">sscha</span>
<span class="kn">import</span> <span class="nn">sscha.Cluster</span>

<span class="kn">import</span> <span class="nn">sys</span><span class="o">,</span> <span class="nn">os</span>


<span class="k">def</span> <span class="nf">configure_cluster</span><span class="p">(</span><span class="n">cluster_workdir</span> <span class="o">=</span> <span class="s2">&quot;H3S&quot;</span><span class="p">):</span>
    <span class="n">cluster</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Cluster</span><span class="o">.</span><span class="n">Cluster</span><span class="p">(</span><span class="n">hostname</span> <span class="o">=</span> <span class="s2">&quot;sschauser@my.beautiful.cluster.eu&quot;</span><span class="p">)</span>

    <span class="n">cluster</span><span class="o">.</span><span class="n">use_memory</span> <span class="o">=</span> <span class="kc">True</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">ram</span> <span class="o">=</span> <span class="mi">180000</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">use_partition</span> <span class="o">=</span> <span class="kc">True</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">partition_name</span> <span class="o">=</span> <span class="s2">&quot;workstations&quot;</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">account_name</span> <span class="o">=</span> <span class="s2">&quot;my_allocation_resources&quot;</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">n_nodes</span> <span class="o">=</span> <span class="mi">1</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">use_cpu</span> <span class="o">=</span> <span class="kc">False</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">custom_params</span><span class="p">[</span><span class="s2">&quot;get-user-env&quot;</span><span class="p">]</span> <span class="o">=</span> <span class="kc">None</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">custom_params</span><span class="p">[</span><span class="s2">&quot;cpus-per-task&quot;</span><span class="p">]</span> <span class="o">=</span> <span class="mi">2</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">custom_params</span><span class="p">[</span><span class="s2">&quot;ntasks-per-node&quot;</span><span class="p">]</span> <span class="o">=</span> <span class="mi">48</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">time</span> <span class="o">=</span> <span class="s2">&quot;12:00:00&quot;</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">n_cpu</span> <span class="o">=</span> <span class="mi">48</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">n_pool</span> <span class="o">=</span> <span class="mi">48</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">job_number</span> <span class="o">=</span> <span class="mi">12</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">batch_size</span> <span class="o">=</span> <span class="mi">2</span>

    <span class="n">home_workdir</span><span class="o">=</span><span class="n">os</span><span class="o">.</span><span class="n">path</span><span class="o">.</span><span class="n">join</span><span class="p">(</span><span class="s2">&quot;$HOME&quot;</span><span class="p">,</span> <span class="n">cluster_workdir</span><span class="p">)</span>
    <span class="n">scratch_workdir</span> <span class="o">=</span> <span class="n">os</span><span class="o">.</span><span class="n">path</span><span class="o">.</span><span class="n">join</span><span class="p">(</span><span class="s2">&quot;/scratch/$USER/&quot;</span><span class="p">,</span> <span class="n">cluster_workdir</span><span class="p">)</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">workdir</span> <span class="o">=</span> <span class="n">home_workdir</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">add_set_minus_x</span> <span class="o">=</span> <span class="kc">True</span>  <span class="c1"># Avoid the set -x</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">load_modules</span> <span class="o">=</span> <span class="sa">f</span><span class="s2">&quot;&quot;&quot;</span>

<span class="s2">module purge</span>
<span class="s2">module load intel</span>
<span class="s2">module load intel-mpi</span>
<span class="s2">module load intel-mkl</span>
<span class="s2">module load quantum-espresso/6.8.0-mpi</span>

<span class="s2">export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK</span>

<span class="s2">mkdir -p </span><span class="si">{</span><span class="n">scratch_workdir</span><span class="si">}</span><span class="s2"></span>
<span class="s2">cp $HOME/espresso/pseudo/* </span><span class="si">{</span><span class="n">scratch_workdir</span><span class="si">}</span><span class="s2">/</span>
<span class="s2">&quot;&quot;&quot;</span>

    <span class="k">def</span> <span class="nf">cp_files</span><span class="p">(</span><span class="n">lbls</span><span class="p">):</span>
        <span class="n">extrain</span> <span class="o">=</span> <span class="sa">f</span><span class="s2">&quot;cd </span><span class="si">{</span><span class="n">scratch_workdir</span><span class="si">}</span><span class="se">\n</span><span class="s2">&quot;</span>
        <span class="n">extraout</span> <span class="o">=</span> <span class="s2">&quot;sleep 1</span><span class="se">\n</span><span class="s2">&quot;</span>
        <span class="k">for</span> <span class="n">lbl</span> <span class="ow">in</span> <span class="n">lbls</span><span class="p">:</span>
            <span class="n">extrain</span> <span class="o">+=</span> <span class="sa">f</span><span class="s2">&quot;cp </span><span class="si">{</span><span class="n">home_workdir</span><span class="si">}</span><span class="s2">/</span><span class="si">{</span><span class="n">lbl</span><span class="si">}</span><span class="s2">.pwi </span><span class="si">{</span><span class="n">scratch_workdir</span><span class="si">}</span><span class="s2">/</span><span class="se">\n</span><span class="s2">&quot;</span>
            <span class="n">extraout</span> <span class="o">+=</span> <span class="sa">f</span><span class="s2">&quot;mv </span><span class="si">{</span><span class="n">scratch_workdir</span><span class="si">}</span><span class="s2">/</span><span class="si">{</span><span class="n">lbl</span><span class="si">}</span><span class="s2">.pwo </span><span class="si">{</span><span class="n">home_workdir</span><span class="si">}</span><span class="s2">/</span><span class="se">\n</span><span class="s2">&quot;</span>

        <span class="k">return</span> <span class="n">extrain</span><span class="p">,</span> <span class="n">extraout</span>

    <span class="c1"># Add the possibility to copy the input files</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">additional_script_parameters</span> <span class="o">=</span> <span class="n">cp_files</span>

    <span class="c1"># Force to open a shell when executing ssh commands</span>
    <span class="c1"># (Otherwise the cluster will not load the module environment)</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">use_active_shell</span> <span class="o">=</span> <span class="kc">True</span>
    <span class="n">cluster</span><span class="o">.</span><span class="n">setup_workdir</span><span class="p">()</span>

    <span class="c1"># Check the communication</span>
    <span class="k">if</span> <span class="ow">not</span> <span class="n">cluster</span><span class="o">.</span><span class="n">CheckCommunication</span><span class="p">():</span>
        <span class="k">raise</span> <span class="ne">ValueError</span><span class="p">(</span><span class="s2">&quot;Impossible to connect to the cluster.&quot;</span><span class="p">)</span>

    <span class="k">return</span> <span class="n">cluster</span>
</pre></div>
</div>
<p>This file contains all the information to connect with the cluster that you can customize to adapt to your HPC center.</p>
<p>Let us dive a bit into the options.</p>
<p>The first thing to know how to configure is the ssh host connection. For example, if I connect to a cluster using the command</p>
<div class="highlight-bash notranslate"><div class="highlight"><pre><span></span>ssh sschauser@my.beautiful.cluster.eu
</pre></div>
</div>
<p>You have to specify the entire string <em>sschauser&#64;my.beautiful.cluster.eu</em> inside the <code class="docutils literal notranslate"><span class="pre">hostname</span></code> key
at the first definition of the cluster. If you have an ssh config file enabled, you can substitute the hostname with the name in the configuration file corresponding to a HostName inside <code class="docutils literal notranslate"><span class="pre">.ssh/config</span></code> located in your home directory.</p>
<p>The best procedure is to enable a public-private key <strong>without</strong> encryption. You can activate the encryption if you have a wallet system in your PC that keeps the password saved, but in this case, the user must log in with the screen unlocked to work.</p>
<p>If the HPC does not allow you to configure a pair of ssh keys for the connection and requires the standard username/password connection,
you can add the <code class="docutils literal notranslate"><span class="pre">pwd</span></code> keyword in the definition of the cluster. This is not encouraged, as you will store your password in clear text inside the script (so if you are in a shared workstation, remember to limit the read access to your scripts to other users, and <em>do not</em> send the script accidentally to other people with your password).</p>
<p>For example:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="n">cluster</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Cluster</span><span class="o">.</span><span class="n">Cluster</span><span class="p">(</span><span class="n">hostname</span><span class="o">=</span><span class="s2">&quot;sschauser@my.beautiful.cluster.eu&quot;</span><span class="p">,</span> <span class="n">pwd</span><span class="o">=</span><span class="s2">&quot;mybeautifulpassword&quot;</span><span class="p">)</span>
</pre></div>
</div>
<p>The other options are all standard SLURM configurations, as the amount of ram, name of partition, and account for the submission, number of nodes, total time, and custom parameters specific for each cluster
These parameters are transformed into the submission script for slurm as</p>
<div class="highlight-bash notranslate"><div class="highlight"><pre><span></span><span class="c1">#SLURM --time=12:00:00</span>
<span class="c1">#SLURM --get-user-env</span>
<span class="c1">#SLURM --cpus-per-task=2</span>
<span class="c1"># [...]</span>
</pre></div>
</div>
<p>Most variables have the <em>use_xxx</em> attribute; if set to False, the corresponding option is not printed.
In the last version of SSCHA, if you manually edit a variable, it should automatically set the corresponding <em>use_xxx</em> to true.</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="n">cluster</span><span class="o">.</span><span class="n">use_partition</span> <span class="o">=</span> <span class="kc">True</span>
<span class="n">cluster</span><span class="o">.</span><span class="n">partition_name</span> <span class="o">=</span> <span class="s2">&quot;workstations&quot;</span>
<span class="n">cluster</span><span class="o">.</span><span class="n">account_name</span> <span class="o">=</span> <span class="s2">&quot;my_allocation_resources&quot;</span>
</pre></div>
</div>
<p>Most clusters must run on specific partitions; in this case, activate the partition flag with the <code class="docutils literal notranslate"><span class="pre">use_partition</span></code> variable and specify the appropriate <code class="docutils literal notranslate"><span class="pre">partition_name</span></code>. Also, most of the time, the computational resources are related to specific accounts indicated with <code class="docutils literal notranslate"><span class="pre">account_name</span></code>.</p>
<p>Particular attention needs to be taken to the following parameters</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="n">cluster</span><span class="o">.</span><span class="n">n_nodes</span> <span class="o">=</span> <span class="mi">1</span>
<span class="n">cluster</span><span class="o">.</span><span class="n">time</span> <span class="o">=</span> <span class="s2">&quot;12:00:00&quot;</span>
<span class="n">cluster</span><span class="o">.</span><span class="n">n_cpu</span> <span class="o">=</span> <span class="mi">48</span>
<span class="n">cluster</span><span class="o">.</span><span class="n">n_pool</span> <span class="o">=</span> <span class="mi">48</span>
<span class="n">cluster</span><span class="o">.</span><span class="n">job_number</span> <span class="o">=</span> <span class="mi">12</span>
<span class="n">cluster</span><span class="o">.</span><span class="n">batch_size</span> <span class="o">=</span> <span class="mi">2</span>
</pre></div>
</div>
<p>These parameters are specific for the kind of calculation.</p>
<blockquote>
<div><ol class="arabic simple">
<li><p>n_nodes specifies the number of nodes</p></li>
<li><p>time specifies the total time</p></li>
<li><p>n_cpu specify how many processors to call quantum Espresso with</p></li>
<li><p>n_pool is the number of pools for the quantum espresso parallelization; it should be the greatest common divisor between the number of CPUs and K points.</p></li>
<li><p>batch_size how many pw.x calculations to group in the same job (executed one after the other without queue time).</p></li>
<li><p>job_number how many jobs will be submitted simultaneously (executed in parallel, but with queue time).</p></li>
</ol>
</div></blockquote>
<p>The total time requested must be roughly the time expected for a single calculation multiplied by the batch_size.
It is convenient to overshoot the requested time, as some configurations may take a bit more time.</p>
<p>The workdir is the directory in which all the input files are copied inside the cluster.
This cluster uses a local scratch for the submission (the job must copy all the input on a local scratch of the node and then copy back the results to the shared filesystem).
If no local scratch is requried, then we can set the working directory (usually a shared scratch) with the command</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="n">cluster</span><span class="o">.</span><span class="n">workdir</span> <span class="o">=</span> <span class="s2">&quot;/scratch/myuser/&quot;</span>
</pre></div>
</div>
<p>However, this submission script (as ekhi) must work on a shared workdir, which is inside the home directory.
Therefore, we must tell the cluster to copy the files from the workdir to the local scratch before and after each calculation.
This is done by setting a custom function, executed for each calculation</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="k">def</span> <span class="nf">cp_files</span><span class="p">(</span><span class="n">lbls</span><span class="p">):</span>
    <span class="n">extrain</span> <span class="o">=</span> <span class="sa">f</span><span class="s2">&quot;cd </span><span class="si">{</span><span class="n">scratch_workdir</span><span class="si">}</span><span class="se">\n</span><span class="s2">&quot;</span>
    <span class="n">extraout</span> <span class="o">=</span> <span class="s2">&quot;sleep 1</span><span class="se">\n</span><span class="s2">&quot;</span>
    <span class="k">for</span> <span class="n">lbl</span> <span class="ow">in</span> <span class="n">lbls</span><span class="p">:</span>
        <span class="n">extrain</span> <span class="o">+=</span> <span class="sa">f</span><span class="s2">&quot;cp </span><span class="si">{</span><span class="n">home_workdir</span><span class="si">}</span><span class="s2">/</span><span class="si">{</span><span class="n">lbl</span><span class="si">}</span><span class="s2">.pwi </span><span class="si">{</span><span class="n">scratch_workdir</span><span class="si">}</span><span class="s2">/</span><span class="se">\n</span><span class="s2">&quot;</span>
        <span class="n">extraout</span> <span class="o">+=</span> <span class="sa">f</span><span class="s2">&quot;mv </span><span class="si">{</span><span class="n">scratch_workdir</span><span class="si">}</span><span class="s2">/</span><span class="si">{</span><span class="n">lbl</span><span class="si">}</span><span class="s2">.pwo </span><span class="si">{</span><span class="n">home_workdir</span><span class="si">}</span><span class="s2">/</span><span class="se">\n</span><span class="s2">&quot;</span>

    <span class="k">return</span> <span class="n">extrain</span><span class="p">,</span> <span class="n">extraout</span>

<span class="c1"># Add the possibility to copy the input files</span>
<span class="n">cluster</span><span class="o">.</span><span class="n">additional_script_parameters</span> <span class="o">=</span> <span class="n">cp_files</span>
</pre></div>
</div>
<p>Each cluster must load modules to run a calculation. All the modules and other commands to run before the calculations are stored in the text variable <em>load_modules</em></p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="n">cluster</span><span class="o">.</span><span class="n">load_modules</span> <span class="o">=</span> <span class="sa">f</span><span class="s2">&quot;&quot;&quot;</span>
<span class="s2">module purge</span>
<span class="s2">module load intel</span>
<span class="s2">module load intel-mpi</span>
<span class="s2">module load intel-mkl</span>
<span class="s2">module load quantum-espresso/6.8.0-mpi</span>

<span class="s2">export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK</span>

<span class="s2">mkdir -p </span><span class="si">{</span><span class="n">scratch_workdir</span><span class="si">}</span><span class="s2"></span>
<span class="s2">cp $HOME/espresso/pseudo/* </span><span class="si">{</span><span class="n">scratch_workdir</span><span class="si">}</span><span class="s2">/</span>
<span class="s2">&quot;&quot;&quot;</span>
</pre></div>
</div>
<p>The specific of the modules to load depends on the cluster,
in this case, we also create the local scratch directory and copy the pseudopotential.</p>
<p>To check the connection and set up the working directory (create it on the cluster if it does not exist) use the</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="n">cluster</span><span class="o">.</span><span class="n">setup_workdir</span><span class="p">()</span>

<span class="k">if</span> <span class="ow">not</span> <span class="n">cluster</span><span class="o">.</span><span class="n">CheckCommunication</span><span class="p">():</span>
    <span class="k">raise</span> <span class="ne">ValueError</span><span class="p">(</span><span class="s2">&quot;Cluster connection failed!&quot;</span><span class="p">)</span>
</pre></div>
</div>
<div class="topic">
<p class="topic-title">Exercise</p>
<p>Customize the cluster.py file to connect to the ekhi server, following the instructions provided in the <a class="reference internal" href="Ekhi.html#ekhiconfig"><span class="std std-ref">ekhi guide</span></a>.</p>
</div>
<section id="how-to-submit-a-calculation-with-a-cluster-automatically">
<h3>How to submit a calculation with a cluster automatically<a class="headerlink" href="#how-to-submit-a-calculation-with-a-cluster-automatically" title="Permalink to this headline"> </a></h3>
<p>Now that we have seen how to configure the cluster, it is time to start an actual calculation.
We can use this option to directly evaluate the ensemble generated manually before in the following way:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span><span class="o">,</span> <span class="nn">cellconstructor.Phonons</span>
<span class="kn">import</span> <span class="nn">sscha</span><span class="o">,</span> <span class="nn">sscha.Ensemble</span>

<span class="c1"># Import the two python scripts for the cluster and espresso configurations</span>
<span class="kn">import</span> <span class="nn">espresso_calculator</span>
<span class="kn">import</span> <span class="nn">cluster</span>

<span class="c1"># Generate an ensemble with 10 configurations</span>
<span class="n">dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="s2">&quot;start_sscha&quot;</span><span class="p">,</span> <span class="mi">3</span><span class="p">)</span>
<span class="n">ensemble</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Ensemble</span><span class="o">.</span><span class="n">Ensemble</span><span class="p">(</span><span class="n">dyn</span><span class="p">,</span> <span class="mi">300</span><span class="p">)</span>
<span class="n">ensemble</span><span class="o">.</span><span class="n">generate</span><span class="p">(</span><span class="mi">10</span><span class="p">)</span>

<span class="c1"># Get the espresso and cluster configurations</span>
<span class="n">espresso_config</span> <span class="o">=</span> <span class="n">espresso_calculator</span><span class="o">.</span><span class="n">get_calculator</span><span class="p">()</span>
<span class="n">cluster_config</span> <span class="o">=</span> <span class="n">cluster</span><span class="o">.</span><span class="n">configure_cluster</span><span class="p">()</span>

<span class="c1"># Compute the ensemble</span>
<span class="n">ensemble</span><span class="o">.</span><span class="n">compute_ensemble</span><span class="p">(</span><span class="n">espresso_config</span><span class="p">,</span> <span class="n">cluster</span><span class="o">=</span><span class="n">cluster_config</span><span class="p">)</span>

<span class="c1"># Save the ensemble (using population 2 to avoid overwriting the other one)</span>
<span class="n">ensemble</span><span class="o">.</span><span class="n">save</span><span class="p">(</span><span class="s2">&quot;data&quot;</span><span class="p">,</span> <span class="mi">2</span><span class="p">)</span>
</pre></div>
</div>
<p>As seen here, once the cluster is configured (but this needs to be done only once),
it is straightforward to compute the ensemble’s energy, forces, and stresses.</p>
<p>While the calculation is running, the temporary files copied from/to the cluster are stored in a directory
that is local_workdir. This is, by default, called cluster_work.
They are called ESP_x.pwi EXP_x.pwo, the input and output files, and with ESP_x.sh you have the SLURM submission script.</p>
<p>Indeed, as you have seen in the previous hands-on session, it is possible to use the cluster keyword also in the <em>SSCHA</em> object of the <em>Relax</em> module to automatize all the procedures.</p>
<p>The following script runs the complete automatic relaxation of the SSCHA.</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span><span class="o">,</span> <span class="nn">cellconstructor.Phonons</span>
<span class="kn">import</span> <span class="nn">sscha</span><span class="o">,</span> <span class="nn">sscha.Ensemble</span>
<span class="kn">import</span> <span class="nn">sscha.SchaMinimizer</span><span class="o">,</span> <span class="nn">sscha.Relax</span>

<span class="c1"># Import the two python scripts for the cluster and espresso configurations</span>
<span class="kn">import</span> <span class="nn">espresso_calculator</span>
<span class="kn">import</span> <span class="nn">cluster</span>

<span class="c1"># Generate an ensemble with 10 configurations</span>
<span class="n">dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="s2">&quot;start_sscha&quot;</span><span class="p">,</span> <span class="mi">3</span><span class="p">)</span>
<span class="n">ensemble</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Ensemble</span><span class="o">.</span><span class="n">Ensemble</span><span class="p">(</span><span class="n">dyn</span><span class="p">,</span> <span class="mi">300</span><span class="p">)</span>

<span class="c1"># Get the espresso and cluster configurations</span>
<span class="n">espresso_config</span> <span class="o">=</span> <span class="n">espresso_calculator</span><span class="o">.</span><span class="n">get_calculator</span><span class="p">()</span>
<span class="n">cluster_config</span> <span class="o">=</span> <span class="n">cluster</span><span class="o">.</span><span class="n">configure_cluster</span><span class="p">()</span>

<span class="c1"># Setup the minimizer</span>
<span class="n">minimizer</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">SchaMinimizer</span><span class="o">.</span><span class="n">SSCHA_Minimizer</span><span class="p">(</span><span class="n">ensemble</span><span class="p">)</span>

<span class="c1"># Setup the automatic relaxation</span>
<span class="n">relax</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Relax</span><span class="o">.</span><span class="n">SSCHA</span><span class="p">(</span><span class="n">minimizer</span><span class="p">,</span> <span class="n">espresso_config</span><span class="p">,</span>
        <span class="n">N_configs</span><span class="o">=</span><span class="mi">10</span><span class="p">,</span>
        <span class="n">max_pop</span><span class="o">=</span><span class="mi">3</span><span class="p">,</span>
        <span class="n">save_ensemble</span><span class="o">=</span><span class="kc">True</span><span class="p">,</span>
        <span class="n">cluster</span><span class="o">=</span><span class="n">cluster_config</span><span class="p">)</span>

<span class="c1"># Setup the IO to save the minimization data and the frequencies</span>
<span class="n">ioinfo</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Utilities</span><span class="o">.</span><span class="n">IOInfo</span><span class="p">()</span>
<span class="n">ioinfo</span><span class="o">.</span><span class="n">SetupSaving</span><span class="p">(</span><span class="s2">&quot;minimization_data&quot;</span><span class="p">)</span>

<span class="c1"># Activate the data saving in the minimization</span>
<span class="n">relax</span><span class="o">.</span><span class="n">setup_custom_functions</span><span class="p">(</span><span class="n">custom_function_post</span><span class="o">=</span><span class="n">ioinfo</span><span class="o">.</span><span class="n">CFP_SaveAll</span><span class="p">)</span>

<span class="c1"># Perform the NVT simulation</span>
<span class="n">relax</span><span class="o">.</span><span class="n">relax</span><span class="p">(</span><span class="n">get_stress</span><span class="o">=</span><span class="kc">True</span><span class="p">)</span>

<span class="c1"># Save the data</span>
<span class="n">relax</span><span class="o">.</span><span class="n">minim</span><span class="o">.</span><span class="n">finalize</span><span class="p">()</span>
<span class="n">relax</span><span class="o">.</span><span class="n">minim</span><span class="o">.</span><span class="n">dyn</span><span class="o">.</span><span class="n">save_qe</span><span class="p">(</span><span class="s2">&quot;final_dyn&quot;</span><span class="p">)</span>
</pre></div>
</div>
<p>As for this NVT, you can also use <em>vc_relax</em> for the NPT simulation or the NVT with variable cell shape.</p>
</section>
</section>


