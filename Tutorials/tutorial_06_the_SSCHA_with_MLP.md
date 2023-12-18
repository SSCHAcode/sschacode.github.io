---
layout: page
title: The SSCHA with machine learning potentials
---

This tutorial was prepared for the [2023 SSCHA School](http://sscha.eu/Schools/2023/home/) by Đorđe Dangić. You can see here the video os the hands-on session:

<iframe width="560" height="315" src="https://www.youtube.com/embed/gaIT7gRECho" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>

The material needed for this tutorial can be downloaded [here](https://ehubox.ehu.eus/s/Y48Wc8iX9Z76jqN).

<p>The minimization of the variational free energy demands a large number of single-point density functional theory (DFT) calculations. These calculations are performed on supercells, repetitions of the primitive unit cell. DFT calculations can become very costly very fast if we need to increase the size of the supercell. This can happen in case we have very slowly decaying second-order force constants, large primitive unit cells, or simply very low symmetry. In some of these cases, DFT is prohibitive due to the large number of atoms per calculation or we simply need a very large number of configurations to converge our results (for example when we need to compute free energy hessian to check the dynamical stability of the system).</p>
<p>In the last ten years, there has been a large amount of research put into developing machine-learned (ML) interatomic potentials. Contrary to the traditional interatomic potentials, they do not have a fixed analytical form and thus are much more flexible and transferable. They are usually trained on a very large number of DFT data and have very good accuracy. ML potentials are considerably slower than the traditional interatomic potentials, however still orders of magnitude faster than DFT, with a considerably better scaling with a number of atoms.</p>
<p>The synergy between SSCHA and machine learning interatomic potentials is obvious. If we can use the machine learning interatomic potentials as a calculator for forces, stresses, and energies we can go to much larger supercells and numbers of configurations. The stochastic sampling employed by SSCHA gives a very good method for obtaining training sets needed to train machine learning interatomic potentials. The force, energy and stresses errors produced by ML interatomic potentials will influence SSCHA results less due to the averaging effects (in case the errors are not biased).</p>
<p>There are a number of freely available implementations of ML interatomic potentials (<a class="reference external" href="https://libatoms.github.io/GAP/">Gaussian Approximation Potentials</a>, <a class="reference external" href="https://nequip.readthedocs.io/en/latest/">NequIP</a>, <a class="reference external" href="https://pacemaker.readthedocs.io/en/latest/">pacemaker</a>, etc.), and at this point, they can be used without a large prior knowledge of the theory behind ML potentials.</p>
<section id="hands-on-exercise">
<h2>Hands-on exercise<a class="headerlink" href="#hands-on-exercise" title="Permalink to this headline"> </a></h2>
<p>For this exercise we will be using <a class="reference external" href="https://libatoms.github.io/GAP/">Gaussian Approximation Potentials</a>, however, the framework can be applied to any other type of ML interatomic potential. In the exercise, we will obtain the training data from the Tersoff interatomic potential, instead of the DFT.</p>
<p>We have provided starting dynamical matrices calculated for the structure at 0 K. Now we will calculate training and test ensemble with Tersoff potential:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">from</span> <span class="nn">quippy.potential</span> <span class="kn">import</span> <span class="n">Potential</span>
<span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span>
<span class="kn">import</span> <span class="nn">cellconstructor.Phonons</span>
<span class="kn">import</span> <span class="nn">sscha</span><span class="o">,</span> <span class="nn">sscha.Ensemble</span>

<span class="n">temperature</span> <span class="o">=</span> <span class="mf">0.0</span>   <span class="c1"># Temperature at which we generate SSCHA configurations</span>
<span class="n">nconf_train</span> <span class="o">=</span> <span class="mi">1000</span>  <span class="c1"># Number of configurations in the training set</span>
<span class="n">nconf_test</span> <span class="o">=</span> <span class="mi">500</span>    <span class="c1"># Number of configurations in the test set</span>

<span class="c1"># Load the Tersoff potential that we want to fit with ML GAP</span>
<span class="n">pot</span> <span class="o">=</span> <span class="n">Potential</span><span class="p">(</span><span class="s1">&#39;IP Tersoff&#39;</span><span class="p">,</span>
        <span class="n">param_filename</span><span class="o">=</span><span class="s1">&#39;./06_the_SSCHA_with_machine_learning_potentials/ip.parms.Tersoff.xml&#39;</span><span class="p">)</span>
<span class="c1"># Load dynamical matrices</span>
<span class="n">dyn_prefix</span> <span class="o">=</span> <span class="s1">&#39;./06_the_SSCHA_with_machine_learning_potentials/start_dyn&#39;</span>
<span class="n">nqirr</span> <span class="o">=</span> <span class="mi">3</span>
<span class="n">dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="n">dyn_prefix</span><span class="p">,</span> <span class="n">nqirr</span><span class="p">)</span>

<span class="c1"># Generate training ensemble</span>
<span class="n">ensemble_train</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Ensemble</span><span class="o">.</span><span class="n">Ensemble</span><span class="p">(</span><span class="n">dyn</span><span class="p">,</span> <span class="n">T0</span><span class="o">=</span><span class="n">temperature</span><span class="p">,</span>
        <span class="n">supercell</span> <span class="o">=</span> <span class="n">dyn</span><span class="o">.</span><span class="n">GetSupercell</span><span class="p">())</span>
<span class="n">ensemble_train</span><span class="o">.</span><span class="n">generate</span><span class="p">(</span><span class="n">N</span> <span class="o">=</span> <span class="n">nconf_train</span><span class="p">)</span>
<span class="n">ensemble_train</span><span class="o">.</span><span class="n">compute_ensemble</span><span class="p">(</span><span class="n">pot</span><span class="p">,</span> <span class="n">compute_stress</span> <span class="o">=</span> <span class="kc">True</span><span class="p">,</span>
        <span class="n">stress_numerical</span> <span class="o">=</span> <span class="kc">False</span><span class="p">,</span> <span class="n">cluster</span> <span class="o">=</span> <span class="kc">None</span><span class="p">,</span> <span class="n">verbose</span> <span class="o">=</span> <span class="kc">True</span><span class="p">)</span>

<span class="c1"># This line will save ensemble in correct format</span>
<span class="n">ensemble_train</span><span class="o">.</span><span class="n">save_enhanced_xyz</span><span class="p">(</span><span class="s1">&#39;train.xyz&#39;</span><span class="p">,</span> <span class="n">append_mode</span> <span class="o">=</span> <span class="kc">False</span><span class="p">,</span>
        <span class="n">stress_key</span> <span class="o">=</span> <span class="s2">&quot;stress&quot;</span><span class="p">,</span> <span class="n">forces_key</span> <span class="o">=</span> <span class="s2">&quot;forces&quot;</span><span class="p">,</span>
        <span class="n">energy_key</span> <span class="o">=</span> <span class="s2">&quot;energy&quot;</span><span class="p">)</span>

<span class="c1"># Generate test ensemble</span>
<span class="n">ensemble_test</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Ensemble</span><span class="o">.</span><span class="n">Ensemble</span><span class="p">(</span><span class="n">dyn</span><span class="p">,</span> <span class="n">T0</span><span class="o">=</span><span class="n">temperature</span><span class="p">,</span>
        <span class="n">supercell</span> <span class="o">=</span> <span class="n">dyn</span><span class="o">.</span><span class="n">GetSupercell</span><span class="p">())</span>
<span class="n">ensemble_test</span><span class="o">.</span><span class="n">generate</span><span class="p">(</span><span class="n">N</span> <span class="o">=</span> <span class="n">nconf_test</span><span class="p">)</span>
<span class="n">ensemble_test</span><span class="o">.</span><span class="n">compute_ensemble</span><span class="p">(</span><span class="n">pot</span><span class="p">,</span> <span class="n">compute_stress</span> <span class="o">=</span> <span class="kc">True</span><span class="p">,</span>
        <span class="n">stress_numerical</span> <span class="o">=</span> <span class="kc">False</span><span class="p">,</span> <span class="n">cluster</span> <span class="o">=</span> <span class="kc">None</span><span class="p">,</span> <span class="n">verbose</span> <span class="o">=</span> <span class="kc">True</span><span class="p">)</span>

<span class="n">ensemble_test</span><span class="o">.</span><span class="n">save_enhanced_xyz</span><span class="p">(</span><span class="s1">&#39;test.xyz&#39;</span><span class="p">,</span> <span class="n">append_mode</span> <span class="o">=</span> <span class="kc">False</span><span class="p">,</span>
        <span class="n">stress_key</span> <span class="o">=</span> <span class="s2">&quot;stress&quot;</span><span class="p">,</span> <span class="n">forces_key</span> <span class="o">=</span> <span class="s2">&quot;forces&quot;</span><span class="p">,</span>
        <span class="n">energy_key</span> <span class="o">=</span> <span class="s2">&quot;energy&quot;</span><span class="p">)</span>
</pre></div>
</div>
<p>The training of the ML interatomic potential can be done with a command gap_fit which should be available after installing quippy-ase. This command takes a large number of arguments so it is easier to make a bash script. We will name it train.sh:</p>
<div class="highlight-bash notranslate"><div class="highlight"><pre><span></span><span class="ch">#!/bin/bash</span>

gap_fit <span class="nv">energy_parameter_name</span><span class="o">=</span>energy <span class="nv">force_parameter_name</span><span class="o">=</span>forces <span class="se">\</span>
    <span class="nv">stress_parameter_name</span><span class="o">=</span>stress <span class="nv">virial_parameter_name</span><span class="o">=</span>virial <span class="se">\</span>
    <span class="nv">do_copy_at_file</span><span class="o">=</span>F <span class="nv">sparse_separate_file</span><span class="o">=</span>T <span class="nv">gp_file</span><span class="o">=</span>GAP.xml <span class="se">\</span>
    <span class="nv">at_file</span><span class="o">=</span>train.xyz <span class="nv">e0_method</span><span class="o">=</span><span class="s2">&quot;average&quot;</span> <span class="se">\</span>
    <span class="nv">default_sigma</span><span class="o">={</span><span class="m">0</span>.001 <span class="m">0</span>.03 <span class="m">0</span>.03 <span class="m">0</span><span class="o">}</span> <span class="nv">sparse_jitter</span><span class="o">=</span><span class="m">1</span>.0e-8 <span class="se">\</span>
    <span class="nv">gap</span><span class="o">={</span>soap <span class="nv">cutoff</span><span class="o">=</span><span class="m">4</span>.2 <span class="nv">n_sparse</span><span class="o">=</span><span class="m">200</span> <span class="nv">covariance_type</span><span class="o">=</span>dot_product <span class="se">\</span>
        <span class="nv">sparse_method</span><span class="o">=</span>cur_points <span class="nv">delta</span><span class="o">=</span><span class="m">0</span>.205 <span class="nv">zeta</span><span class="o">=</span><span class="m">4</span> <span class="nv">l_max</span><span class="o">=</span><span class="m">4</span> <span class="se">\</span>
        <span class="nv">n_max</span><span class="o">=</span><span class="m">8</span> <span class="nv">atom_sigma</span><span class="o">=</span><span class="m">0</span>.5 <span class="nv">cutoff_transition_width</span><span class="o">=</span><span class="m">0</span>.8 <span class="se">\</span>
        add_species <span class="o">}</span>
</pre></div>
</div>
<p>The meaning of each argument is not important right now, but can be easily looked up on the official website <a class="reference external" href="https://libatoms.github.io/GAP/gap_fit.html">https://libatoms.github.io/GAP/gap_fit.html</a>. We run the training command:</p>
<div class="highlight-bash notranslate"><div class="highlight"><pre><span></span>bash train.sh
</pre></div>
</div>
<p>Note, the training is memory intensive, so you may need to allocate extra memory on your virtual machine if you are employing Quantum Mobile. 4Gb of Ram are required.
You may need to restart the virtual machine.</p>
<p>This should take a minute or so. Once it is finished, if the memory was enough and the command typed correctly, one should obtain the <code class="docutils literal notranslate"><span class="pre">GAP.xml</span></code> file in the working directory containing the GAP ML interatomic potential. We can use test.xyz file to check how well our ML potential reproduces data with this simple script:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">numpy</span> <span class="k">as</span> <span class="nn">np</span>
<span class="kn">import</span> <span class="nn">ase</span>
<span class="kn">from</span> <span class="nn">ase</span> <span class="kn">import</span> <span class="n">Atoms</span>
<span class="kn">from</span> <span class="nn">quippy.potential</span> <span class="kn">import</span> <span class="n">Potential</span>
<span class="kn">import</span> <span class="nn">matplotlib</span>
<span class="kn">import</span> <span class="nn">matplotlib.pyplot</span> <span class="k">as</span> <span class="nn">plt</span>
<span class="kn">from</span> <span class="nn">matplotlib.gridspec</span> <span class="kn">import</span> <span class="n">GridSpec</span>
<span class="n">fpaths</span> <span class="o">=</span> <span class="n">matplotlib</span><span class="o">.</span><span class="n">font_manager</span><span class="o">.</span><span class="n">findSystemFonts</span><span class="p">()</span>

<span class="n">infile</span> <span class="o">=</span> <span class="s1">&#39;test.xyz&#39;</span> <span class="c1"># test datasets</span>

<span class="c1"># Read in .xyz files using ase method</span>
<span class="n">atoms</span> <span class="o">=</span> <span class="n">ase</span><span class="o">.</span><span class="n">io</span><span class="o">.</span><span class="n">read</span><span class="p">(</span><span class="n">infile</span><span class="p">,</span> <span class="s1">&#39;:&#39;</span><span class="p">,</span> <span class="nb">format</span><span class="o">=</span><span class="s1">&#39;extxyz&#39;</span><span class="p">)</span>
<span class="n">nconf</span> <span class="o">=</span> <span class="nb">len</span><span class="p">(</span><span class="n">atoms</span><span class="p">)</span>
<span class="nb">print</span><span class="p">(</span><span class="s1">&#39;Number of configurations in the dataset: &#39;</span> <span class="o">+</span> <span class="nb">str</span><span class="p">(</span><span class="n">nconf</span><span class="p">))</span>
<span class="n">natoms</span> <span class="o">=</span> <span class="p">[</span><span class="nb">len</span><span class="p">(</span><span class="n">at</span><span class="o">.</span><span class="n">symbols</span><span class="p">)</span> <span class="k">for</span> <span class="n">at</span> <span class="ow">in</span> <span class="n">atoms</span><span class="p">]</span>

<span class="c1"># Load in newly trained GAP potential</span>
<span class="n">gap_file</span> <span class="o">=</span> <span class="s1">&#39;./GAP.xml&#39;</span>
<span class="n">pot</span> <span class="o">=</span> <span class="n">Potential</span><span class="p">(</span><span class="s1">&#39;IP GAP&#39;</span><span class="p">,</span> <span class="n">param_filename</span><span class="o">=</span><span class="n">gap_file</span><span class="p">)</span> <span class="c1"># Read in potential</span>

<span class="c1"># Collect previously calculated (with Tersoff) atomic properties</span>
<span class="n">dft_energies</span> <span class="o">=</span> <span class="p">[</span><span class="n">atom</span><span class="o">.</span><span class="n">get_potential_energy</span><span class="p">()</span> <span class="k">for</span> <span class="n">atom</span> <span class="ow">in</span> <span class="n">atoms</span><span class="p">]</span>
<span class="n">dft_forces</span> <span class="o">=</span> <span class="p">[</span><span class="n">atom</span><span class="o">.</span><span class="n">get_forces</span><span class="p">()</span> <span class="k">for</span> <span class="n">atom</span> <span class="ow">in</span> <span class="n">atoms</span><span class="p">]</span>
<span class="n">dft_stress</span> <span class="o">=</span> <span class="p">[</span><span class="n">atom</span><span class="o">.</span><span class="n">get_stress</span><span class="p">()[</span><span class="mi">0</span><span class="p">:</span><span class="mi">3</span><span class="p">]</span> <span class="k">for</span> <span class="n">atom</span> <span class="ow">in</span> <span class="n">atoms</span><span class="p">]</span>

<span class="c1"># Now recalculate them with GAP</span>
<span class="n">en_gap</span> <span class="o">=</span> <span class="p">[]</span>
<span class="n">forces_gap</span> <span class="o">=</span> <span class="p">[]</span>
<span class="n">stress_gap</span> <span class="o">=</span> <span class="p">[]</span>
<span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">nconf</span><span class="p">):</span>
    <span class="k">if</span><span class="p">(</span><span class="n">i</span><span class="o">%</span><span class="mi">100</span> <span class="o">==</span> <span class="mi">0</span><span class="p">):</span>
        <span class="nb">print</span><span class="p">(</span><span class="s1">&#39;Configuration: &#39;</span><span class="p">,</span> <span class="n">i</span> <span class="o">+</span> <span class="mi">1</span><span class="p">)</span>
    <span class="c1"># Make ase Atoms object</span>
    <span class="n">atoms_gap</span> <span class="o">=</span> <span class="n">Atoms</span><span class="p">(</span><span class="n">symbols</span> <span class="o">=</span> <span class="n">atoms</span><span class="p">[</span><span class="n">i</span><span class="p">]</span><span class="o">.</span><span class="n">symbols</span><span class="p">,</span> <span class="n">cell</span> <span class="o">=</span> <span class="n">atoms</span><span class="p">[</span><span class="n">i</span><span class="p">]</span><span class="o">.</span><span class="n">cell</span><span class="p">,</span>\
    <span class="n">scaled_positions</span> <span class="o">=</span> <span class="n">atoms</span><span class="p">[</span><span class="n">i</span><span class="p">]</span><span class="o">.</span><span class="n">get_scaled_positions</span><span class="p">(),</span>\
    <span class="n">calculator</span> <span class="o">=</span> <span class="n">pot</span><span class="p">,</span> <span class="n">pbc</span> <span class="o">=</span> <span class="kc">True</span><span class="p">)</span>
    <span class="c1"># Calculate total energies of structures with GAP</span>
    <span class="n">en_gap</span><span class="o">.</span><span class="n">append</span><span class="p">(</span><span class="n">atoms_gap</span><span class="o">.</span><span class="n">get_potential_energy</span><span class="p">())</span>
    <span class="c1"># Calculate forces on atoms</span>
    <span class="n">forces_gap</span><span class="o">.</span><span class="n">append</span><span class="p">(</span><span class="n">atoms_gap</span><span class="o">.</span><span class="n">get_forces</span><span class="p">())</span>
    <span class="c1"># Calculate stress and only take diagonal elements</span>
    <span class="n">stress_gap</span><span class="o">.</span><span class="n">append</span><span class="p">(</span><span class="n">atoms_gap</span><span class="o">.</span><span class="n">get_stress</span><span class="p">()[</span><span class="mi">0</span><span class="p">:</span><span class="mi">3</span><span class="p">])</span>

<span class="c1"># Calculate errors</span>
<span class="n">energy_errors</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">zeros_like</span><span class="p">(</span><span class="n">en_gap</span><span class="p">)</span>
<span class="n">forces_errors</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">zeros_like</span><span class="p">(</span><span class="n">forces_gap</span><span class="p">)</span>
<span class="n">GPa</span> <span class="o">=</span> <span class="mf">1.60217733e-19</span><span class="o">*</span><span class="mf">1.0e21</span>
<span class="n">stress_errors</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">zeros_like</span><span class="p">(</span><span class="n">stress_gap</span><span class="p">)</span>
<span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">nconf</span><span class="p">):</span>
    <span class="c1"># Calculate energy errors</span>
    <span class="n">energy_errors</span><span class="p">[</span><span class="n">i</span><span class="p">]</span> <span class="o">=</span> <span class="p">(</span><span class="n">atoms</span><span class="p">[</span><span class="n">i</span><span class="p">]</span><span class="o">.</span><span class="n">get_potential_energy</span><span class="p">()</span> <span class="o">-</span>\
    <span class="n">en_gap</span><span class="p">[</span><span class="n">i</span><span class="p">])</span><span class="o">/</span><span class="n">natoms</span><span class="p">[</span><span class="n">i</span><span class="p">]</span>
    <span class="c1"># Calculate errors on forces</span>
    <span class="n">forces_errors</span><span class="p">[</span><span class="n">i</span><span class="p">]</span> <span class="o">=</span> <span class="n">atoms</span><span class="p">[</span><span class="n">i</span><span class="p">]</span><span class="o">.</span><span class="n">get_forces</span><span class="p">()</span> <span class="o">-</span> <span class="n">forces_gap</span><span class="p">[</span><span class="n">i</span><span class="p">]</span>
    <span class="c1"># Calculate errors on stress</span>
    <span class="n">stress_errors</span><span class="p">[</span><span class="n">i</span><span class="p">]</span> <span class="o">=</span> <span class="n">dft_stress</span><span class="p">[</span><span class="n">i</span><span class="p">]</span> <span class="o">-</span> <span class="n">stress_gap</span><span class="p">[</span><span class="n">i</span><span class="p">]</span>

<span class="c1"># Function to plot Tersoff vs GAP results</span>
<span class="k">def</span> <span class="nf">plot_comparison</span><span class="p">(</span><span class="n">ax</span><span class="p">,</span> <span class="n">data1</span><span class="p">,</span> <span class="n">data2</span><span class="p">,</span> <span class="n">data3</span><span class="p">,</span> \
<span class="n">xlabel</span> <span class="o">=</span> <span class="s1">&#39;Original energy (eV)&#39;</span><span class="p">,</span> <span class="n">ylabel</span> <span class="o">=</span> <span class="s1">&#39;ML energy (eV)&#39;</span><span class="p">):</span>
    <span class="n">sizes</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">array</span><span class="p">(</span><span class="n">data3</span><span class="o">/</span><span class="n">np</span><span class="o">.</span><span class="n">amax</span><span class="p">(</span><span class="n">data3</span><span class="p">))</span><span class="o">*</span><span class="mf">2.0</span>
    <span class="n">ax</span><span class="o">.</span><span class="n">scatter</span><span class="p">(</span><span class="n">data1</span><span class="p">,</span> <span class="n">data2</span><span class="p">,</span> <span class="n">marker</span> <span class="o">=</span> <span class="s1">&#39;o&#39;</span><span class="p">,</span> <span class="n">s</span> <span class="o">=</span> <span class="n">sizes</span><span class="p">,</span>  <span class="n">c</span> <span class="o">=</span> <span class="s1">&#39;red&#39;</span><span class="p">)</span>
    <span class="n">lims</span> <span class="o">=</span> <span class="p">[</span><span class="n">np</span><span class="o">.</span><span class="n">min</span><span class="p">([</span><span class="n">ax</span><span class="o">.</span><span class="n">get_xlim</span><span class="p">(),</span> <span class="n">ax</span><span class="o">.</span><span class="n">get_ylim</span><span class="p">()]),</span>\
    <span class="n">np</span><span class="o">.</span><span class="n">max</span><span class="p">([</span><span class="n">ax</span><span class="o">.</span><span class="n">get_xlim</span><span class="p">(),</span> <span class="n">ax</span><span class="o">.</span><span class="n">get_ylim</span><span class="p">()])]</span>
    <span class="c1"># now plot both limits against eachother</span>
    <span class="n">ax</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="n">lims</span><span class="p">,</span> <span class="n">lims</span><span class="p">,</span> <span class="s1">&#39;k-&#39;</span><span class="p">,</span> <span class="n">alpha</span><span class="o">=</span><span class="mf">0.75</span><span class="p">,</span> <span class="n">zorder</span><span class="o">=</span><span class="mi">0</span><span class="p">)</span>
    <span class="n">ax</span><span class="o">.</span><span class="n">set_xlabel</span><span class="p">(</span><span class="n">xlabel</span><span class="p">)</span>
    <span class="n">ax</span><span class="o">.</span><span class="n">set_ylabel</span><span class="p">(</span><span class="n">ylabel</span><span class="p">)</span>

<span class="c1"># Function to plot histogram of errors</span>
<span class="k">def</span> <span class="nf">plot_error_histogram</span><span class="p">(</span><span class="n">ax</span><span class="p">,</span> <span class="n">x</span><span class="p">,</span> <span class="n">nbins</span><span class="p">,</span> <span class="n">xlabel</span><span class="p">):</span>
    <span class="kn">import</span> <span class="nn">scipy.stats</span> <span class="k">as</span> <span class="nn">st</span>
    <span class="kn">from</span> <span class="nn">scipy.stats</span> <span class="kn">import</span> <span class="n">norm</span>

    <span class="n">ax</span><span class="o">.</span><span class="n">hist</span><span class="p">(</span><span class="n">x</span><span class="p">,</span> <span class="n">density</span><span class="o">=</span><span class="kc">True</span><span class="p">,</span> <span class="n">bins</span><span class="o">=</span><span class="n">nbins</span><span class="p">)</span>
    <span class="n">mu</span><span class="p">,</span> <span class="n">std</span> <span class="o">=</span> <span class="n">norm</span><span class="o">.</span><span class="n">fit</span><span class="p">(</span><span class="n">x</span><span class="p">)</span>
    <span class="n">xmin</span><span class="p">,</span> <span class="n">xmax</span> <span class="o">=</span> <span class="n">ax</span><span class="o">.</span><span class="n">get_xlim</span><span class="p">()</span>
    <span class="n">x1</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">linspace</span><span class="p">(</span><span class="n">xmin</span><span class="p">,</span> <span class="n">xmax</span><span class="p">,</span> <span class="mi">100</span><span class="p">)</span>
    <span class="n">p</span> <span class="o">=</span> <span class="n">norm</span><span class="o">.</span><span class="n">pdf</span><span class="p">(</span><span class="n">x1</span><span class="p">,</span> <span class="n">mu</span><span class="p">,</span> <span class="n">std</span><span class="p">)</span>
    <span class="n">ax</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="n">x1</span><span class="p">,</span> <span class="n">p</span><span class="p">,</span> <span class="s1">&#39;k&#39;</span><span class="p">,</span> <span class="n">linewidth</span><span class="o">=</span><span class="mi">2</span><span class="p">)</span>
    <span class="n">ax</span><span class="o">.</span><span class="n">set_ylabel</span><span class="p">(</span><span class="s2">&quot;Probability&quot;</span><span class="p">)</span>
    <span class="n">ax</span><span class="o">.</span><span class="n">set_xlabel</span><span class="p">(</span><span class="n">xlabel</span><span class="p">)</span>

<span class="c1"># Sometimes forces arrays can be ragged list, this will flatten them</span>
<span class="k">def</span> <span class="nf">flatten</span><span class="p">(</span><span class="n">xs</span><span class="p">):</span>
    <span class="n">res</span> <span class="o">=</span> <span class="p">[]</span>
    <span class="k">def</span> <span class="nf">loop</span><span class="p">(</span><span class="n">ys</span><span class="p">):</span>
        <span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="n">ys</span><span class="p">:</span>
            <span class="k">if</span> <span class="nb">isinstance</span><span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="nb">list</span><span class="p">):</span>
                <span class="n">loop</span><span class="p">(</span><span class="n">i</span><span class="p">)</span>
            <span class="k">elif</span><span class="p">(</span><span class="nb">isinstance</span><span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="n">np</span><span class="o">.</span><span class="n">ndarray</span><span class="p">)):</span>
                <span class="n">loop</span><span class="p">(</span><span class="n">i</span><span class="o">.</span><span class="n">tolist</span><span class="p">())</span>
            <span class="k">else</span><span class="p">:</span>
                <span class="n">res</span><span class="o">.</span><span class="n">append</span><span class="p">(</span><span class="n">i</span><span class="p">)</span>
    <span class="n">loop</span><span class="p">(</span><span class="n">xs</span><span class="p">)</span>
    <span class="k">return</span> <span class="n">res</span>


<span class="c1"># Plot stuff</span>
<span class="n">plt</span><span class="o">.</span><span class="n">rcParams</span><span class="p">[</span><span class="s2">&quot;font.family&quot;</span><span class="p">]</span> <span class="o">=</span> <span class="s2">&quot;Times New Roman&quot;</span>
<span class="n">plt</span><span class="o">.</span><span class="n">rcParams</span><span class="p">[</span><span class="s1">&#39;mathtext.fontset&#39;</span><span class="p">]</span> <span class="o">=</span> <span class="s2">&quot;stix&quot;</span>
<span class="n">plt</span><span class="o">.</span><span class="n">rcParams</span><span class="o">.</span><span class="n">update</span><span class="p">({</span><span class="s1">&#39;font.size&#39;</span><span class="p">:</span> <span class="mi">16</span><span class="p">})</span>

<span class="n">fig</span> <span class="o">=</span> <span class="n">plt</span><span class="o">.</span><span class="n">figure</span><span class="p">(</span><span class="n">figsize</span><span class="o">=</span><span class="p">(</span><span class="mf">6.4</span><span class="o">*</span><span class="mf">3.0</span><span class="p">,</span> <span class="mf">4.8</span><span class="o">*</span><span class="mf">2.0</span><span class="p">))</span>
<span class="n">gs1</span> <span class="o">=</span> <span class="n">GridSpec</span><span class="p">(</span><span class="mi">2</span><span class="p">,</span> <span class="mi">3</span><span class="p">)</span>

<span class="n">ax00</span> <span class="o">=</span> <span class="n">fig</span><span class="o">.</span><span class="n">add_subplot</span><span class="p">(</span><span class="n">gs1</span><span class="p">[</span><span class="mi">0</span><span class="p">,</span> <span class="mi">0</span><span class="p">])</span>
<span class="n">plot_comparison</span><span class="p">(</span><span class="n">ax00</span><span class="p">,</span> <span class="n">dft_energies</span><span class="p">,</span> <span class="n">en_gap</span><span class="p">,</span> <span class="n">energy_errors</span><span class="p">,</span> \
<span class="n">xlabel</span> <span class="o">=</span> <span class="s1">&#39;Original energy (eV)&#39;</span><span class="p">,</span> <span class="n">ylabel</span> <span class="o">=</span> <span class="s1">&#39;ML energy (eV)&#39;</span><span class="p">)</span>
<span class="n">ax01</span> <span class="o">=</span> <span class="n">fig</span><span class="o">.</span><span class="n">add_subplot</span><span class="p">(</span><span class="n">gs1</span><span class="p">[</span><span class="mi">0</span><span class="p">,</span> <span class="mi">1</span><span class="p">])</span>
<span class="n">plot_comparison</span><span class="p">(</span><span class="n">ax01</span><span class="p">,</span> <span class="n">dft_forces</span><span class="p">,</span> <span class="n">forces_gap</span><span class="p">,</span> <span class="n">forces_errors</span><span class="p">,</span> \
<span class="n">xlabel</span> <span class="o">=</span> <span class="sa">r</span><span class="s1">&#39;Original force (eV/$\AA$)&#39;</span><span class="p">,</span> <span class="n">ylabel</span> <span class="o">=</span> <span class="sa">r</span><span class="s1">&#39;ML force (eV/$\AA$)&#39;</span><span class="p">)</span>
<span class="n">ax02</span> <span class="o">=</span> <span class="n">fig</span><span class="o">.</span><span class="n">add_subplot</span><span class="p">(</span><span class="n">gs1</span><span class="p">[</span><span class="mi">0</span><span class="p">,</span> <span class="mi">2</span><span class="p">])</span>
<span class="n">plot_comparison</span><span class="p">(</span><span class="n">ax02</span><span class="p">,</span> <span class="n">np</span><span class="o">.</span><span class="n">array</span><span class="p">(</span><span class="n">flatten</span><span class="p">(</span><span class="n">dft_stress</span><span class="p">))</span><span class="o">*</span><span class="n">GPa</span><span class="p">,</span> \
<span class="n">np</span><span class="o">.</span><span class="n">array</span><span class="p">(</span><span class="n">flatten</span><span class="p">(</span><span class="n">stress_gap</span><span class="p">))</span><span class="o">*</span><span class="n">GPa</span><span class="p">,</span> <span class="n">np</span><span class="o">.</span><span class="n">array</span><span class="p">(</span><span class="n">flatten</span><span class="p">(</span><span class="n">stress_errors</span><span class="p">))</span><span class="o">*</span><span class="n">GPa</span><span class="p">,</span> \
<span class="n">xlabel</span> <span class="o">=</span> <span class="s1">&#39;Original stress (GPa)&#39;</span><span class="p">,</span> <span class="n">ylabel</span> <span class="o">=</span> <span class="s1">&#39;ML stress (GPa)&#39;</span><span class="p">)</span>

<span class="n">ax10</span> <span class="o">=</span> <span class="n">fig</span><span class="o">.</span><span class="n">add_subplot</span><span class="p">(</span><span class="n">gs1</span><span class="p">[</span><span class="mi">1</span><span class="p">,</span> <span class="mi">0</span><span class="p">])</span>
<span class="n">plot_error_histogram</span><span class="p">(</span><span class="n">ax10</span><span class="p">,</span> <span class="n">energy_errors</span><span class="p">,</span> <span class="mi">100</span><span class="p">,</span> <span class="s1">&#39;Energy error (eV/atom)&#39;</span><span class="p">)</span>
<span class="n">ax11</span> <span class="o">=</span> <span class="n">fig</span><span class="o">.</span><span class="n">add_subplot</span><span class="p">(</span><span class="n">gs1</span><span class="p">[</span><span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">])</span>
<span class="n">flattened_forces</span> <span class="o">=</span> <span class="n">flatten</span><span class="p">(</span><span class="n">forces_errors</span><span class="p">)</span>
<span class="n">plot_error_histogram</span><span class="p">(</span><span class="n">ax11</span><span class="p">,</span> <span class="n">flattened_forces</span><span class="p">,</span> <span class="mi">100</span><span class="p">,</span> <span class="s1">&#39;Force error (eV/$\AA$)&#39;</span><span class="p">)</span>
<span class="n">ax12</span> <span class="o">=</span> <span class="n">fig</span><span class="o">.</span><span class="n">add_subplot</span><span class="p">(</span><span class="n">gs1</span><span class="p">[</span><span class="mi">1</span><span class="p">,</span> <span class="mi">2</span><span class="p">])</span>
<span class="n">plot_error_histogram</span><span class="p">(</span><span class="n">ax12</span><span class="p">,</span> <span class="n">np</span><span class="o">.</span><span class="n">array</span><span class="p">([</span><span class="n">item</span> <span class="k">for</span> <span class="n">sublist</span> <span class="ow">in</span> <span class="n">stress_errors</span> <span class="k">for</span> <span class="n">item</span> <span class="ow">in</span> <span class="n">sublist</span><span class="p">])</span><span class="o">*</span><span class="n">GPa</span><span class="p">,</span>\
<span class="mi">100</span><span class="p">,</span> <span class="s1">&#39;Stress error (GPa)&#39;</span><span class="p">)</span>

<span class="n">fig</span><span class="o">.</span><span class="n">savefig</span><span class="p">(</span><span class="s1">&#39;test.pdf&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">show</span><span class="p">()</span>
</pre></div>
</div>
<p>In the upper panel figures, ideally, we would like points to lie on the diagonal. When fitting interatomic potential we aim for normal distribution of errors centered at 0 (without bias) with as small as possible standard deviation. We should have very nice results for energies and forces. Now that we are happy with the potential let us use it to relax SSCHA at 2000 K:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">from</span> <span class="nn">quippy.potential</span> <span class="kn">import</span> <span class="n">Potential</span>
<span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span>
<span class="kn">import</span> <span class="nn">cellconstructor.Phonons</span>

<span class="c1"># Import the SSCHA engine (we will use it later)</span>
<span class="kn">import</span> <span class="nn">sscha</span><span class="o">,</span> <span class="nn">sscha.Ensemble</span><span class="o">,</span> <span class="nn">sscha.SchaMinimizer</span><span class="o">,</span> <span class="nn">sscha.Relax</span>

<span class="c1"># Declare SSCHA variables</span>
<span class="n">temperature</span> <span class="o">=</span> <span class="mf">2000.0</span>
<span class="n">nconf</span> <span class="o">=</span> <span class="mi">1000</span>
<span class="n">max_pop</span> <span class="o">=</span> <span class="mi">10000</span>

<span class="c1"># Load in the GAP potential</span>
<span class="n">gap_file</span> <span class="o">=</span> <span class="s1">&#39;./GAP.xml&#39;</span>
<span class="n">pot</span> <span class="o">=</span> <span class="n">Potential</span><span class="p">(</span><span class="s1">&#39;IP GAP&#39;</span><span class="p">,</span> <span class="n">param_filename</span><span class="o">=</span><span class="n">gap_file</span><span class="p">)</span>

<span class="c1"># Load in the SSCHA dynamical matrices</span>
<span class="n">dyn_prefix</span> <span class="o">=</span> <span class="s1">&#39;./06_the_SSCHA_with_machine_learning_potentials/start_dyn&#39;</span>
<span class="n">nqirr</span> <span class="o">=</span> <span class="mi">3</span>
<span class="n">dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="n">dyn_prefix</span><span class="p">,</span> <span class="n">nqirr</span><span class="p">)</span>

<span class="c1"># Relax the structure at 2000 K</span>
<span class="n">ensemble</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Ensemble</span><span class="o">.</span><span class="n">Ensemble</span><span class="p">(</span><span class="n">dyn</span><span class="p">,</span> <span class="n">T0</span><span class="o">=</span><span class="n">temperature</span><span class="p">,</span>
<span class="n">supercell</span> <span class="o">=</span> <span class="n">dyn</span><span class="o">.</span><span class="n">GetSupercell</span><span class="p">())</span>
<span class="n">ensemble</span><span class="o">.</span><span class="n">generate</span><span class="p">(</span><span class="n">N</span> <span class="o">=</span> <span class="n">nconf</span><span class="p">)</span>
<span class="n">minimizer</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">SchaMinimizer</span><span class="o">.</span><span class="n">SSCHA_Minimizer</span><span class="p">(</span><span class="n">ensemble</span><span class="p">)</span>
<span class="n">minimizer</span><span class="o">.</span><span class="n">min_step_dyn</span> <span class="o">=</span> <span class="mf">0.1</span>
<span class="n">minimizer</span><span class="o">.</span><span class="n">kong_liu_ratio</span> <span class="o">=</span> <span class="mf">0.5</span>
<span class="n">minimizer</span><span class="o">.</span><span class="n">meaningful_factor</span> <span class="o">=</span> <span class="mf">0.001</span>
<span class="n">minimizer</span><span class="o">.</span><span class="n">max_ka</span> <span class="o">=</span> <span class="mi">100000</span>
<span class="n">relax</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Relax</span><span class="o">.</span><span class="n">SSCHA</span><span class="p">(</span><span class="n">minimizer</span><span class="p">,</span> <span class="n">ase_calculator</span> <span class="o">=</span> <span class="n">pot</span><span class="p">,</span>
<span class="n">N_configs</span> <span class="o">=</span> <span class="n">nconf</span><span class="p">,</span> <span class="n">max_pop</span> <span class="o">=</span> <span class="n">max_pop</span><span class="p">,</span> <span class="n">save_ensemble</span> <span class="o">=</span> <span class="kc">True</span><span class="p">)</span>
<span class="n">relax</span><span class="o">.</span><span class="n">vc_relax</span><span class="p">(</span><span class="n">ensemble_loc</span><span class="o">=</span><span class="s1">&#39;Ensemble_location&#39;</span><span class="p">)</span>
<span class="n">relax</span><span class="o">.</span><span class="n">minim</span><span class="o">.</span><span class="n">dyn</span><span class="o">.</span><span class="n">save_qe</span><span class="p">(</span><span class="s1">&#39;final_dyn&#39;</span><span class="p">)</span>
<span class="c1"># We can check minimization procedure</span>
<span class="n">relax</span><span class="o">.</span><span class="n">minim</span><span class="o">.</span><span class="n">plot_results</span><span class="p">(</span><span class="n">save_filename</span> <span class="o">=</span> <span class="s1">&#39;sscha&#39;</span><span class="p">,</span> <span class="n">plot</span> <span class="o">=</span> <span class="kc">False</span><span class="p">)</span>
</pre></div>
</div>
<p>We have relaxed SSCHA at 2000 K. We can check that everything went well in “sscha” file. However, we do not know whether this is correct. We need to check how our ML potential performs at 2000 K.</p>
<div class="topic">
<p class="topic-title">Exercise:</p>
<p>Let’s create a dataset of SSCHA-generated configuration at 2000 K using GAP relaxed dynamical matrices and compute it using Tersoff potential. Next, check the performance of the GAP ML potential against this new dataset.</p>
</div>
<div class="topic">
<p class="topic-title">Excercise:</p>
<p>We should see GAP performing quite worse compared to the test.xyz case. How can we improve GAP potential? Let’s do it.</p>
</div>
<div class="topic">
<p class="topic-title">Excercise:</p>
<p>How do the Tersoff phonons compare to GAP phonons?</p>
</div>
<div class="topic">
<p class="topic-title">Excercise:</p>
<p>Does this translate to larger supercells?</p>
</div>
</section>
