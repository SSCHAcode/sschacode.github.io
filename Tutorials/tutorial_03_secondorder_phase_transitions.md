---
layout: page
title: Calculations of second-order phase transitions with the SSCHA
---

In this hands-on, we learn how to calculate second-order phase
transitions within the SSCHA.

This tutorial was prepared for the [2023 SSCHA School](http://sscha.eu/Schools/2023/home/) by Diego Martinez Gutierrez. You can see here the video os the hands-on session:

<iframe width="560" height="315" src="https://www.youtube.com/embed/JWQyO-EACVw" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>

The material needed for this tutorial can be downloaded [here](https://ehubox.ehu.eus/s/NBHTeeiM5Dmg3XF).

<h2>Structural instability: calculation of the Hessian<a class="headerlink" href="#structural-instability-calculation-of-the-hessian" title="Permalink to this headline"> </a></h2>
<p>According to Landau’s theory, a second-order phase transition occurs
when the free energy curvature around the high-symmetry structure on the direction of the order parameter becomes negative:</p>
<figure class="align-default" id="id3">
<span id="fig-second-order"></span><a class="reference internal image-reference" href="../figures_03/second_order.png"><img alt="Second order." src="../figures_03/second_order.png" style="width: 400px;" /></a>
<figcaption>
<p><span class="caption-number">Fig. 7 </span><span class="caption-text">Landau’s theory of second-order phase transitions.</span><a class="headerlink" href="#id3" title="Permalink to this image"> </a></p>
</figcaption>
</figure>
<p>For structural <em>displacive</em> phase transitions, the order parameter is associated to phonon atomic displacements:</p>
<div class="math notranslate nohighlight">
\[\frac{\partial^2 F}{\partial R_a \partial R_b}\]</div>
<p>Thus, the Free energy Hessian is the central quantity to study second-order phase transitions. The SSCHA provides an analytical equation for the free energy Hessian, derived by Raffaello Bianco in the work <a class="reference external" href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.014111">Bianco et. al. Phys. Rev. B 96, 014111</a>.
The free energy curvature can be written as:</p>
<div class="math notranslate nohighlight">
\[\frac{\partial^2 F}{\partial {R_a}\partial {R_b}} = \Phi_{ab} + \sum_{cdef} \stackrel{(3)}{\Phi}_{acd}[1 - \Lambda\stackrel{(4)}{\Phi}]^{-1}_{cdef} \stackrel{(3)}{\Phi}_{efb}\]</div>
<p>Fortunately, this complex equation can be evaluated from the ensemble with a simple function call:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="n">ensemble</span><span class="o">.</span><span class="n">get_free_energy_hessian</span><span class="p">()</span>
</pre></div>
</div>
<p>Lets see a practical example, first we calculate the SSCHA dynamical matrix for the SnTe:</p>
<p>To speedup the calculations, we will use a force-field that can mimic the physics of ferroelectric transitions in FCC lattices.</p>
<p>We begin importing some libraries:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="ch">#!/usr/bin/env python</span>
<span class="c1"># -*- coding: utf-8 -*-</span>
<span class="c1">#</span>
<span class="c1">#  SSCHA_exercise_Calculus.py</span>
<span class="c1">#</span>
<span class="c1"># Import the cellconstructor stuff</span>
<span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span>
<span class="kn">import</span> <span class="nn">cellconstructor.Phonons</span>
<span class="kn">import</span> <span class="nn">cellconstructor.ForceTensor</span>
<span class="kn">import</span> <span class="nn">cellconstructor.Structure</span>
<span class="kn">import</span> <span class="nn">cellconstructor.Spectral</span>

<span class="c1"># Import the modules of the force field</span>
<span class="kn">import</span> <span class="nn">fforces</span> <span class="k">as</span> <span class="nn">ff</span>
<span class="kn">import</span> <span class="nn">fforces.Calculator</span>

<span class="c1"># Import the modules to run the sscha</span>
<span class="kn">import</span> <span class="nn">sscha</span><span class="o">,</span> <span class="nn">sscha.Ensemble</span><span class="o">,</span> <span class="nn">sscha.SchaMinimizer</span>
<span class="kn">import</span> <span class="nn">sscha.Relax</span><span class="o">,</span> <span class="nn">sscha.Utilities</span>

<span class="kn">import</span> <span class="nn">spglib</span>
<span class="kn">from</span> <span class="nn">ase.visualize</span> <span class="kn">import</span> <span class="n">view</span>

<span class="c1"># Import Matplotlib to plot</span>
<span class="kn">import</span> <span class="nn">numpy</span> <span class="k">as</span> <span class="nn">np</span>
<span class="kn">import</span> <span class="nn">matplotlib.pyplot</span> <span class="k">as</span> <span class="nn">plt</span>
<span class="kn">from</span> <span class="nn">matplotlib</span> <span class="kn">import</span> <span class="n">cm</span>
<span class="kn">import</span> <span class="nn">timeit</span>
</pre></div>
</div>
<p>Next we set some variables for the calculation:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="c1">#Setting the variables:</span>
<span class="c1">#Setting the temperature in Kelvin:</span>
<span class="n">Temperature</span> <span class="o">=</span> <span class="mi">0</span>
<span class="c1">#Setting the number of configurations:</span>
<span class="n">configurations</span> <span class="o">=</span> <span class="mi">50</span>
<span class="c1">#Setting the names and location of the files:</span>
<span class="n">Files_dyn_SnTe</span> <span class="o">=</span> <span class="s2">&quot;ffield_dynq&quot;</span>
<span class="c1">#Set the number of irreducible q (reated to the supercell size):</span>
<span class="n">nqirr</span> <span class="o">=</span> <span class="mi">3</span>
<span class="c1">#Setting the frequencies output file:</span>
<span class="n">File_frequencies</span> <span class="o">=</span> <span class="s2">&quot;frequencies.dat&quot;</span>
<span class="c1">#Setting the dynamical matrix output filename:</span>
<span class="n">File_final_dyn</span> <span class="o">=</span> <span class="s2">&quot;final_sscha_T</span><span class="si">{}</span><span class="s2">_&quot;</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="nb">int</span><span class="p">(</span><span class="n">Temperature</span><span class="p">))</span>
</pre></div>
</div>
<p>Now we need to calculate the SSCHA dynamical matrix. For that we follow some steps:</p>
<ol class="arabic simple">
<li><p>First we prepare the Toy model force field that substitutes the usual <em>ab-initio</em> for this tutorial.
This force field needs the harmonic dynamical matrix to be initialized, and the higher order parameters.
Finally, the dynamical matrix for the minimization is loaded and readied. Since we are studying a system that has a spontaneous symmetry breaking at low temperature, the harmonic dynamical matrices will have imaginary phonons. We must enforce phonons to be positive definite to start a SSCHA minimization.</p></li>
</ol>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="c1"># Load the dynamical matrix for the force field</span>
<span class="n">ff_dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="s2">&quot;ffield_dynq&quot;</span><span class="p">,</span> <span class="mi">3</span><span class="p">)</span>

<span class="c1"># Setup the forcefield with the correct parameters</span>
<span class="n">ff_calculator</span> <span class="o">=</span> <span class="n">ff</span><span class="o">.</span><span class="n">Calculator</span><span class="o">.</span><span class="n">ToyModelCalculator</span><span class="p">(</span><span class="n">ff_dyn</span><span class="p">)</span>
<span class="n">ff_calculator</span><span class="o">.</span><span class="n">type_cal</span> <span class="o">=</span> <span class="s2">&quot;pbtex&quot;</span>
<span class="n">ff_calculator</span><span class="o">.</span><span class="n">p3</span> <span class="o">=</span> <span class="mf">0.036475</span>
<span class="n">ff_calculator</span><span class="o">.</span><span class="n">p4</span> <span class="o">=</span> <span class="o">-</span><span class="mf">0.022</span>
<span class="n">ff_calculator</span><span class="o">.</span><span class="n">p4x</span> <span class="o">=</span> <span class="o">-</span><span class="mf">0.014</span>

<span class="c1"># Initialization of the SSCHA matrix</span>
<span class="n">dyn_sscha</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="n">Files_dyn_SnTe</span><span class="p">,</span> <span class="n">nqirr</span><span class="p">)</span>
<span class="c1"># Flip the imaginary frequencies into real ones</span>
<span class="n">dyn_sscha</span><span class="o">.</span><span class="n">ForcePositiveDefinite</span><span class="p">()</span>
<span class="c1"># Apply the ASR and the symmetry group</span>
<span class="n">dyn_sscha</span><span class="o">.</span><span class="n">Symmetrize</span><span class="p">()</span>
</pre></div>
</div>
<ol class="arabic simple" start="2">
<li><p>The next step is to create the ensembles for the specified temperature. As an extra, we also look for the space group of the structure.</p></li>
</ol>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="n">ensemble</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Ensemble</span><span class="o">.</span><span class="n">Ensemble</span><span class="p">(</span><span class="n">dyn_sscha</span><span class="p">,</span>
        <span class="n">T0</span> <span class="o">=</span> <span class="n">Temperature</span><span class="p">,</span> <span class="n">supercell</span> <span class="o">=</span> <span class="n">dyn_sscha</span><span class="o">.</span><span class="n">GetSupercell</span><span class="p">())</span>
<span class="c1"># Detect space group</span>
<span class="n">symm</span><span class="o">=</span><span class="n">spglib</span><span class="o">.</span><span class="n">get_spacegroup</span><span class="p">(</span><span class="n">dyn_sscha</span><span class="o">.</span><span class="n">structure</span><span class="o">.</span><span class="n">get_ase_atoms</span><span class="p">(),</span>
        <span class="mf">0.005</span><span class="p">)</span>
<span class="nb">print</span><span class="p">(</span><span class="s1">&#39;Initial SG = &#39;</span><span class="p">,</span> <span class="n">symm</span><span class="p">)</span>
</pre></div>
</div>
<ol class="arabic" start="3">
<li><p>Next comes the minimization step. Here we can set the fourth root minimization, in which, instead of optimizing the auxiliary dynamical matrices themselves, we will optimize their fourth root.</p>
<div class="math notranslate nohighlight">
\[\Phi = \left({\sqrt[4]{\Phi}}\right)^4\]</div>
<p>This constrains the dynamical matrix to be positive definite during the minimization.
Next the automatic relaxation is set with the option here to use the Sobol sequence for the ensemble generation.</p>
<p>We also set a custom function to save the frequencies at each iteration, to see how they evolves. This is very useful to understand if the algorithm is converged or not.</p>
<p>Then the dynamical matrix of the converged minimization is saved in a file, and finally we take a look at the space group and the structure.</p>
</li>
</ol>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="n">minim</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">SchaMinimizer</span><span class="o">.</span><span class="n">SSCHA_Minimizer</span><span class="p">(</span><span class="n">ensemble</span><span class="p">)</span>

<span class="c1"># Now we setup the minimization parameters</span>
<span class="c1"># Since we are quite far from the correct solution,</span>
<span class="c1"># we will use a small optimization step</span>
<span class="n">minim</span><span class="o">.</span><span class="n">set_minimization_step</span><span class="p">(</span><span class="mf">0.25</span><span class="p">)</span>

<span class="c1"># Reduce the threshold for the gradient convergence</span>
<span class="n">minim</span><span class="o">.</span><span class="n">meaningful_factor</span> <span class="o">=</span> <span class="mf">0.01</span>

<span class="c1"># If the minimization ends with few steps (less than 10),</span>
<span class="c1"># decrease it, if it takes too much, increase it</span>

<span class="c1"># We decrease the Kong-Liu effective sample size below</span>
<span class="c1"># which the population is stopped</span>
<span class="n">minim</span><span class="o">.</span><span class="n">kong_liu_ratio</span> <span class="o">=</span> <span class="mf">0.5</span> <span class="c1"># Default 0.5</span>
<span class="c1"># We relax the structure</span>
<span class="n">relax</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Relax</span><span class="o">.</span><span class="n">SSCHA</span><span class="p">(</span><span class="n">minim</span><span class="p">,</span>
                  <span class="n">ase_calculator</span> <span class="o">=</span> <span class="n">ff_calculator</span><span class="p">,</span>
                  <span class="n">N_configs</span> <span class="o">=</span> <span class="n">configurations</span><span class="p">,</span>
                  <span class="n">max_pop</span> <span class="o">=</span> <span class="mi">50</span><span class="p">)</span>

<span class="c1"># Setup the custom function to print the frequencies</span>
<span class="c1"># at each step of the minimization</span>
<span class="n">io_func</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Utilities</span><span class="o">.</span><span class="n">IOInfo</span><span class="p">()</span>
<span class="n">io_func</span><span class="o">.</span><span class="n">SetupSaving</span><span class="p">(</span><span class="n">File_frequencies</span><span class="p">)</span>
<span class="c1"># The file that will contain the frequencies is frequencies.dat</span>

<span class="c1"># Now tell relax to call the function to save the frequencies</span>
<span class="c1"># after each iteration</span>
<span class="c1"># CFP stands for Custom Function Post (Post = after the minimization step)</span>
<span class="c1">#relax.setup_custom_functions(custom_function_post = io_func.CFP_SaveFrequencies)</span>
<span class="n">relax</span><span class="o">.</span><span class="n">setup_custom_functions</span><span class="p">(</span><span class="n">custom_function_post</span> <span class="o">=</span> <span class="n">io_func</span><span class="o">.</span><span class="n">CFP_SaveAll</span><span class="p">)</span>
<span class="c1"># Finally we do all the free energy calculations.</span>
<span class="n">relax</span><span class="o">.</span><span class="n">relax</span><span class="p">()</span>
<span class="c1">#relax.vc_relax(static_bulk_modulus=40, fix_volume = False)</span>

<span class="c1"># Save the final dynamical matrix</span>
<span class="n">relax</span><span class="o">.</span><span class="n">minim</span><span class="o">.</span><span class="n">dyn</span><span class="o">.</span><span class="n">save_qe</span><span class="p">(</span><span class="n">File_final_dyn</span><span class="p">)</span>
<span class="c1"># Detect space group</span>
<span class="n">symm</span><span class="o">=</span><span class="n">spglib</span><span class="o">.</span><span class="n">get_spacegroup</span><span class="p">(</span><span class="n">relax</span><span class="o">.</span><span class="n">minim</span><span class="o">.</span><span class="n">dyn</span><span class="o">.</span><span class="n">structure</span><span class="o">.</span><span class="n">get_ase_atoms</span><span class="p">(),</span>
            <span class="mf">0.005</span><span class="p">)</span>
<span class="nb">print</span><span class="p">(</span><span class="s1">&#39;New SG = &#39;</span><span class="p">,</span> <span class="n">symm</span><span class="p">)</span>
<span class="n">view</span><span class="p">(</span><span class="n">relax</span><span class="o">.</span><span class="n">minim</span><span class="o">.</span><span class="n">dyn</span><span class="o">.</span><span class="n">structure</span><span class="o">.</span><span class="n">get_ase_atoms</span><span class="p">())</span>
</pre></div>
</div>
<p>This code will calculate the SSCHA minimization with the <em>ff_calculator</em>. We cat use <strong>sscha-plot-data.py</strong> to take a look at the minimization.</p>
<div class="highlight-bash notranslate"><div class="highlight"><pre><span></span>python sscha-plot-data.py frequencies.dat
</pre></div>
</div>
<figure class="align-default" id="fig-plot-data">
<a class="reference internal image-reference" href="../figures_03/Figure_1_N.png"><img alt="Minimization figures" src="../figures_03/Figure_1_N.png" style="width: 400px;" /></a>
</figure>
<figure class="align-default">
<a class="reference internal image-reference" href="../figures_03/Figure_2_N.png"><img alt="Frequencies minimization" src="../figures_03/Figure_2_N.png" style="width: 400px;" /></a>
</figure>
<p>Note: this force field model is not able to compute stress, as it is defined only at fixed volume, so we cannot use it for a variable cell relaxation.</p>
<p><strong>Now we can search for instabilities.</strong></p>
<p>If we have a very small mode in the SSCHA frequencies, it means that associated to that mode we have huge fluctuations. This can indicate an instability. However, to test this we need to compute the free energy curvature along this mode. This can be obtained in one shot thanks to the theory developed in <a class="reference external" href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.014111">Bianco et. al. Phys. Rev. B 96, 014111.</a></p>
<p>For that we create another program to do the job.</p>
<p>As before, we begin importing some libraries and setting variables:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="ch">#!/usr/bin/env python</span>
<span class="c1"># -*- coding: utf-8 -*-</span>
<span class="c1">#</span>
<span class="c1">#  SSCHA_exercise_Unstable.py</span>
<span class="c1">#</span>
<span class="c1"># Import the cellconstructor stuff</span>
<span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span>
<span class="kn">import</span> <span class="nn">cellconstructor.Phonons</span>
<span class="kn">import</span> <span class="nn">cellconstructor.ForceTensor</span>
<span class="kn">import</span> <span class="nn">cellconstructor.Structure</span>
<span class="kn">import</span> <span class="nn">cellconstructor.Spectral</span>

<span class="c1"># Import the modules of the force field</span>
<span class="kn">import</span> <span class="nn">fforces</span> <span class="k">as</span> <span class="nn">ff</span>
<span class="kn">import</span> <span class="nn">fforces.Calculator</span>

<span class="c1"># Import the modules to run the sscha</span>
<span class="kn">import</span> <span class="nn">sscha</span><span class="o">,</span> <span class="nn">sscha.Ensemble</span><span class="o">,</span> <span class="nn">sscha.SchaMinimizer</span>
<span class="kn">import</span> <span class="nn">sscha.Relax</span><span class="o">,</span> <span class="nn">sscha.Utilities</span>

<span class="kn">import</span> <span class="nn">spglib</span>
<span class="kn">from</span> <span class="nn">ase.visualize</span> <span class="kn">import</span> <span class="n">view</span>

<span class="c1"># Import Matplotlib to plot</span>
<span class="kn">import</span> <span class="nn">numpy</span> <span class="k">as</span> <span class="nn">np</span>
<span class="kn">import</span> <span class="nn">matplotlib.pyplot</span> <span class="k">as</span> <span class="nn">plt</span>
<span class="kn">from</span> <span class="nn">matplotlib</span> <span class="kn">import</span> <span class="n">cm</span>
<span class="kn">import</span> <span class="nn">timeit</span>

<span class="c1">#Setting the variables:</span>
<span class="c1">#Setting the temperature in Kelvin:</span>
<span class="n">Temperature</span> <span class="o">=</span> <span class="mi">0</span>
<span class="c1">#Setting the number of configurations:</span>
<span class="n">configurations</span> <span class="o">=</span> <span class="mi">50</span>
<span class="c1">#Setting the names and location of the files:</span>
<span class="n">Files_dyn_SnTe</span> <span class="o">=</span> <span class="s2">&quot;ffield_dynq&quot;</span>
<span class="c1">#Set the number of irreducible q (reated to the supercell size):</span>
<span class="n">nqirr</span> <span class="o">=</span> <span class="mi">3</span>
<span class="c1">#Setting the frequencies output file:</span>
<span class="n">File_frequencies</span> <span class="o">=</span> <span class="s2">&quot;frequencies.dat&quot;</span>
<span class="c1">#Setting the dynamical matrix output filename:</span>
<span class="n">File_final_dyn</span> <span class="o">=</span> <span class="s2">&quot;final_sscha_T</span><span class="si">{}</span><span class="s2">_&quot;</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="nb">int</span><span class="p">(</span><span class="n">Temperature</span><span class="p">))</span>
</pre></div>
</div>
<p>Now we look for that instability:</p>
<ol class="arabic simple">
<li><p>The <em>ff_calculator</em> toy potential is defined as we have seen in the previous program.</p></li>
</ol>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="c1"># Load the dynamical matrix for the force field</span>
<span class="n">ff_dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="s2">&quot;ffield_dynq&quot;</span><span class="p">,</span> <span class="mi">3</span><span class="p">)</span>

<span class="c1"># Setup the forcefield with the correct parameters</span>
<span class="n">ff_calculator</span> <span class="o">=</span> <span class="n">ff</span><span class="o">.</span><span class="n">Calculator</span><span class="o">.</span><span class="n">ToyModelCalculator</span><span class="p">(</span><span class="n">ff_dyn</span><span class="p">)</span>
<span class="n">ff_calculator</span><span class="o">.</span><span class="n">type_cal</span> <span class="o">=</span> <span class="s2">&quot;pbtex&quot;</span>
<span class="n">ff_calculator</span><span class="o">.</span><span class="n">p3</span> <span class="o">=</span> <span class="mf">0.036475</span>
<span class="n">ff_calculator</span><span class="o">.</span><span class="n">p4</span> <span class="o">=</span> <span class="o">-</span><span class="mf">0.022</span>
<span class="n">ff_calculator</span><span class="o">.</span><span class="n">p4x</span> <span class="o">=</span> <span class="o">-</span><span class="mf">0.014</span>

<span class="c1"># Initialization of the SSCHA matrix</span>
<span class="n">dyn_sscha</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="n">Files_dyn_SnTe</span><span class="p">,</span> <span class="n">nqirr</span><span class="p">)</span>
<span class="n">dyn_sscha</span><span class="o">.</span><span class="n">ForcePositiveDefinite</span><span class="p">()</span>

<span class="c1"># Apply also the ASR and the symmetry group</span>
<span class="n">dyn_sscha</span><span class="o">.</span><span class="n">Symmetrize</span><span class="p">()</span>
</pre></div>
</div>
<ol class="arabic simple" start="2">
<li><p>Next, we will load the dynamical matrix calculated previously with the <em>ff_calculator</em> toy potential, so there is no need to calculate it again.</p></li>
</ol>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="c1"># The SSCHA dynamical matrix is needed (the one after convergence)</span>
<span class="c1"># We reload the final result (no need to rerun the sscha minimization)</span>
<span class="n">dyn_sscha_final</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="n">File_final_dyn</span><span class="p">,</span> <span class="n">nqirr</span><span class="p">)</span>
</pre></div>
</div>
<ol class="arabic simple" start="3">
<li><p>Then, as the Hessian calculation is more sensible, we generate a new ensemble with more configurations.
To compute the hessian we will use an ensemble of 10000 configurations.
Note here that we can use less if we use Sobol sequence or we can load a previously generated ensemble.</p></li>
</ol>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="c1"># We reset the ensemble</span>
<span class="n">ensemble</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Ensemble</span><span class="o">.</span><span class="n">Ensemble</span><span class="p">(</span><span class="n">dyn_sscha_final</span><span class="p">,</span> <span class="n">T0</span> <span class="o">=</span> <span class="n">Temperature</span><span class="p">,</span>
                    <span class="n">supercell</span> <span class="o">=</span> <span class="n">dyn_sscha_final</span><span class="o">.</span><span class="n">GetSupercell</span><span class="p">())</span>

<span class="c1"># We need a bigger ensemble to properly compute the hessian</span>
<span class="c1"># Here we will use 10000 configurations</span>
<span class="n">ensemble</span><span class="o">.</span><span class="n">generate</span><span class="p">(</span><span class="mi">5000</span><span class="p">,</span> <span class="n">sobol</span> <span class="o">=</span> <span class="kc">True</span><span class="p">)</span>
<span class="c1">#ensemble.generate(10000, sobol = False)</span>
<span class="c1">#We could also load the ensemble with</span>
<span class="c1"># ensemble.load(&quot;data_ensemble_final&quot;, N = 100, population = 5)</span>
</pre></div>
</div>
<ol class="arabic simple" start="4">
<li><p>We now compute forces and energies using the force field calculator.</p></li>
</ol>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="c1"># We now compute forces and energies using the force field calculator</span>
<span class="n">ensemble</span><span class="o">.</span><span class="n">get_energy_forces</span><span class="p">(</span><span class="n">ff_calculator</span><span class="p">,</span> <span class="n">compute_stress</span> <span class="o">=</span> <span class="kc">False</span><span class="p">)</span>
</pre></div>
</div>
<ol class="arabic" start="5">
<li><p>Finally the free energy hessian is calculated in the <em>hessian</em> function.
We can choose if we neglect or not in the calculation the four phonon scattering process. Four phonon scattering processes require a huge memory allocation for big systems, that scales as (3N)^4 with N the number of atoms in the supercell. Moreover, it may require also more configurations to converge.</p>
<p>In almost all the systems we studied up to now, we found this four phonon scattering at high order to be negligible. We remark, that the SSCHA minimization already includes four phonon scattering at the lowest order perturbation theory, thus neglecting this term only affects combinations of one or more four phonon scattering with two three phonon scatterings (high order diagrams). For more details, see <a class="reference external" href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.014111">Bianco et. al. Phys. Rev. B 96, 014111.</a></p>
<p>We can then print the frequencies of the hessian. If an imaginary frequency is present, then the system wants to spontaneously break the high symmetry phase.</p>
</li>
</ol>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="nb">print</span><span class="p">(</span><span class="s2">&quot;Updating the importance sampling...&quot;</span><span class="p">)</span>
<span class="c1"># If the sscha matrix was not the one used to compute the ensemble</span>
<span class="c1"># We must update the ensemble weights</span>
<span class="c1"># We can also use this function to simulate a different temperature.</span>
<span class="n">ensemble</span><span class="o">.</span><span class="n">update_weights</span><span class="p">(</span><span class="n">dyn_sscha_final</span><span class="p">,</span> <span class="n">Temperature</span><span class="p">)</span>
<span class="c1"># ----------- COMPUTE THE FREE ENERGY HESSIAN -----------</span>
<span class="nb">print</span><span class="p">(</span><span class="s2">&quot;Computing the free energy hessian...&quot;</span><span class="p">)</span>
<span class="n">dyn_hessian</span> <span class="o">=</span> <span class="n">ensemble</span><span class="o">.</span><span class="n">get_free_energy_hessian</span><span class="p">(</span><span class="n">include_v4</span> <span class="o">=</span> <span class="kc">False</span><span class="p">)</span>
<span class="c1"># We neglect high-order four phonon scattering</span>
<span class="c1">#dyn_hessian = ensemble.get_free_energy_hessian(include_v4 = True,</span>
<span class="c1">#              get_full_hessian = True,verbose = True) # Full calculus</span>
<span class="c1"># We can save the free energy hessian as a dynamical matrix</span>
<span class="c1"># in quantum espresso format</span>
<span class="n">dyn_hessian</span><span class="o">.</span><span class="n">save_qe</span><span class="p">(</span><span class="s2">&quot;hessian&quot;</span><span class="p">)</span>
<span class="c1"># -------------------------------------------------------</span>
<span class="c1"># We calculate the frequencies of the hessian:</span>
<span class="n">w_hessian</span><span class="p">,</span> <span class="n">pols_hessian</span> <span class="o">=</span> <span class="n">dyn_hessian</span><span class="o">.</span><span class="n">DiagonalizeSupercell</span><span class="p">()</span>

<span class="c1"># Print all the frequency converting them into cm-1 (They are in Ry)</span>
<span class="nb">print</span><span class="p">(</span><span class="s2">&quot;</span><span class="se">\n</span><span class="s2">&quot;</span><span class="o">.</span><span class="n">join</span><span class="p">([</span><span class="s2">&quot;</span><span class="si">{:16.4f}</span><span class="s2"> cm-1&quot;</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="n">w</span> <span class="o">*</span> <span class="n">CC</span><span class="o">.</span><span class="n">Units</span><span class="o">.</span><span class="n">RY_TO_CM</span><span class="p">)</span> <span class="k">for</span> <span class="n">w</span> <span class="ow">in</span> <span class="n">w_hessian</span><span class="p">]))</span>
</pre></div>
</div>
<p>The frequencies in the free energy hessian are temperature dependent.</p>
<p>We can look at the eigenmodes of the free energy hessian to check if we have imaginary phonons. If there are negative frequencies then we found an instability. You can check what happens if you include the fourth order.</p>
<div class="topic">
<p class="topic-title">Exercise</p>
<p>The Sobol sequences reduces the number of configurations by doing a better mapping of the gaussian than a random distribution. By uniformity spreading the samplings with a low discrepancy sequence like Sobol it is possible to reduce the number of configurations needed. Low discrepancy sequences tend to sample space “more uniformly” than random numbers. Algorithms that use such sequences may have superior convergence.
You can test this in the calculation of the hessian by changing the number of configurations and the mapping scheme in the <em>ensemble.generate()</em> function.</p>
</div>

<section id="second-order-phase-transition">
<h2>Second order phase transition<a class="headerlink" href="#second-order-phase-transition" title="Permalink to this headline"> </a></h2>
<p>Up to now we studied the system at T=0K and we found that there is an instability. However, we can repeat the minimization at many temperatures, and track the phonon frequency to see which is the temperature at which the system becomes stable.</p>
<p>Again we load and set the variables. Now the we have several temperatures so we store them in an array:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="ch">#!/usr/bin/env python</span>
<span class="c1"># -*- coding: utf-8 -*-</span>
<span class="c1">#</span>
<span class="c1">#  SSCHA_exercise_Unstable.py</span>
<span class="c1">#</span>
<span class="c1"># Import the cellconstructor stuff</span>
<span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span>
<span class="kn">import</span> <span class="nn">cellconstructor.Phonons</span>
<span class="kn">import</span> <span class="nn">cellconstructor.ForceTensor</span>
<span class="kn">import</span> <span class="nn">cellconstructor.Structure</span>
<span class="kn">import</span> <span class="nn">cellconstructor.Spectral</span>

<span class="c1"># Import the modules of the force field</span>
<span class="kn">import</span> <span class="nn">fforces</span> <span class="k">as</span> <span class="nn">ff</span>
<span class="kn">import</span> <span class="nn">fforces.Calculator</span>

<span class="c1"># Import the modules to run the sscha</span>
<span class="kn">import</span> <span class="nn">sscha</span><span class="o">,</span> <span class="nn">sscha.Ensemble</span><span class="o">,</span> <span class="nn">sscha.SchaMinimizer</span>
<span class="kn">import</span> <span class="nn">sscha.Relax</span><span class="o">,</span> <span class="nn">sscha.Utilities</span>

<span class="kn">import</span> <span class="nn">spglib</span>
<span class="kn">from</span> <span class="nn">ase.visualize</span> <span class="kn">import</span> <span class="n">view</span>

<span class="c1"># Import Matplotlib to plot</span>
<span class="kn">import</span> <span class="nn">numpy</span> <span class="k">as</span> <span class="nn">np</span>
<span class="kn">import</span> <span class="nn">matplotlib.pyplot</span> <span class="k">as</span> <span class="nn">plt</span>
<span class="kn">from</span> <span class="nn">matplotlib</span> <span class="kn">import</span> <span class="n">cm</span>
<span class="kn">import</span> <span class="nn">timeit</span>

<span class="c1">#Setting the variables:</span>
<span class="c1">#Setting the temperature in Kelvin:</span>
<span class="n">Temperature</span> <span class="o">=</span> <span class="mi">0</span>
<span class="c1">#Setting the number of configurations:</span>
<span class="n">configurations</span> <span class="o">=</span> <span class="mi">50</span>
<span class="c1">#Setting the names and location of the files:</span>
<span class="n">Files_dyn_SnTe</span> <span class="o">=</span> <span class="s2">&quot;ffield_dynq&quot;</span>
<span class="c1">#Set the number of irreducible q (reated to the supercell size):</span>
<span class="n">nqirr</span> <span class="o">=</span> <span class="mi">3</span>
<span class="c1">#Setting the frequencies output file:</span>
<span class="n">File_frequencies</span> <span class="o">=</span> <span class="s2">&quot;frequencies.dat&quot;</span>
<span class="c1">#Setting the dynamical matrix output filename:</span>
<span class="n">File_final_dyn</span> <span class="o">=</span> <span class="s2">&quot;final_sscha_T</span><span class="si">{}</span><span class="s2">_&quot;</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="nb">int</span><span class="p">(</span><span class="n">Temperature</span><span class="p">))</span>
<span class="n">sobol</span> <span class="o">=</span> <span class="kc">False</span>
<span class="n">sobol_scatter</span> <span class="o">=</span> <span class="kc">False</span>
</pre></div>
</div>
<ol class="arabic simple">
<li><p>Like in the previous program, first we prepare the Toy model force field</p></li>
</ol>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="c1"># Load the dynamical matrix for the force field</span>
<span class="n">ff_dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="s2">&quot;ffield_dynq&quot;</span><span class="p">,</span> <span class="mi">3</span><span class="p">)</span>

<span class="c1"># Setup the forcefield with the correct parameters</span>
<span class="n">ff_calculator</span> <span class="o">=</span> <span class="n">ff</span><span class="o">.</span><span class="n">Calculator</span><span class="o">.</span><span class="n">ToyModelCalculator</span><span class="p">(</span><span class="n">ff_dyn</span><span class="p">)</span>
<span class="n">ff_calculator</span><span class="o">.</span><span class="n">type_cal</span> <span class="o">=</span> <span class="s2">&quot;pbtex&quot;</span>
<span class="n">ff_calculator</span><span class="o">.</span><span class="n">p3</span> <span class="o">=</span> <span class="mf">0.036475</span>
<span class="n">ff_calculator</span><span class="o">.</span><span class="n">p4</span> <span class="o">=</span> <span class="o">-</span><span class="mf">0.022</span>
<span class="n">ff_calculator</span><span class="o">.</span><span class="n">p4x</span> <span class="o">=</span> <span class="o">-</span><span class="mf">0.014</span>
</pre></div>
</div>
<ol class="arabic simple" start="2">
<li><p>We are going to need a range of temperatures for this calculation:</p></li>
</ol>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="c1"># Define the temperatures, from 50 to 300 K, 6 temperatures</span>
<span class="n">temperatures</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">linspace</span><span class="p">(</span><span class="mi">50</span><span class="p">,</span> <span class="mi">300</span><span class="p">,</span> <span class="mi">6</span><span class="p">)</span>

<span class="n">lowest_hessian_mode</span> <span class="o">=</span> <span class="p">[]</span>
<span class="n">lowest_sscha_mode</span> <span class="o">=</span> <span class="p">[]</span>

<span class="c1"># Perform a simulation at each temperature</span>
<span class="n">t_old</span> <span class="o">=</span> <span class="n">Temperature</span>
</pre></div>
</div>
<ol class="arabic simple" start="3">
<li><p>In the next part we condense the calculation of the hessians in a loop for different temperatures. In the end, it searches for the lowest non acoustic frequency to save with the correspondent auxiliar sscha frequency.</p></li>
</ol>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="k">for</span> <span class="n">Temperature</span> <span class="ow">in</span> <span class="n">temperatures</span><span class="p">:</span>
    <span class="c1"># Load the starting dynamical matrix</span>
    <span class="n">dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="n">File_final_dyn</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="nb">int</span><span class="p">(</span><span class="n">t_old</span><span class="p">)),</span> <span class="n">nqirr</span><span class="p">)</span>

    <span class="c1"># Prepare the ensemble</span>
    <span class="n">ensemble</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Ensemble</span><span class="o">.</span><span class="n">Ensemble</span><span class="p">(</span><span class="n">dyn</span><span class="p">,</span> <span class="n">T0</span> <span class="o">=</span> <span class="n">Temperature</span><span class="p">,</span>
                                  <span class="n">supercell</span> <span class="o">=</span> <span class="n">dyn</span><span class="o">.</span><span class="n">GetSupercell</span><span class="p">())</span>

    <span class="c1"># Prepare the minimizer</span>
    <span class="n">minim</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">SchaMinimizer</span><span class="o">.</span><span class="n">SSCHA_Minimizer</span><span class="p">(</span><span class="n">ensemble</span><span class="p">)</span>
    <span class="n">minim</span><span class="o">.</span><span class="n">min_step_struc</span> <span class="o">=</span> <span class="mf">0.05</span>
    <span class="n">minim</span><span class="o">.</span><span class="n">min_step_dyn</span> <span class="o">=</span> <span class="mf">0.002</span>
    <span class="n">minim</span><span class="o">.</span><span class="n">kong_liu_ratio</span> <span class="o">=</span> <span class="mf">0.5</span>
    <span class="n">minim</span><span class="o">.</span><span class="n">meaningful_factor</span> <span class="o">=</span> <span class="mf">0.000001</span>
    <span class="c1">#minim.root_representation = &quot;root4&quot;</span>
    <span class="c1">#minim.precond_dyn = False</span>
    <span class="c1">#minim.minim_struct = True</span>
    <span class="c1">#minim.neglect_symmetries = True</span>
    <span class="n">minim</span><span class="o">.</span><span class="n">enforce_sum_rule</span> <span class="o">=</span> <span class="kc">True</span>  <span class="c1"># Lorenzo&#39;s solution to the error</span>

    <span class="c1"># Prepare the relaxer (through many population)</span>
    <span class="n">relax</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Relax</span><span class="o">.</span><span class="n">SSCHA</span><span class="p">(</span><span class="n">minim</span><span class="p">,</span> <span class="n">ase_calculator</span> <span class="o">=</span> <span class="n">ff_calculator</span><span class="p">,</span>
                  <span class="n">N_configs</span><span class="o">=</span><span class="n">configurations</span><span class="p">,</span> <span class="n">max_pop</span><span class="o">=</span><span class="mi">20</span><span class="p">)</span>

    <span class="c1"># Relax</span>
    <span class="n">relax</span><span class="o">.</span><span class="n">relax</span><span class="p">(</span><span class="n">sobol</span> <span class="o">=</span> <span class="n">sobol</span><span class="p">,</span> <span class="n">sobol_scramble</span> <span class="o">=</span> <span class="n">sobol_scatter</span><span class="p">)</span>
    <span class="c1">#relax.relax()</span>

    <span class="c1"># Save the dynamical matrix</span>
    <span class="n">relax</span><span class="o">.</span><span class="n">minim</span><span class="o">.</span><span class="n">dyn</span><span class="o">.</span><span class="n">save_qe</span><span class="p">(</span><span class="n">File_final_dyn</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="nb">int</span><span class="p">(</span><span class="n">Temperature</span><span class="p">)))</span>

    <span class="c1"># Detect space group</span>
    <span class="n">symm</span><span class="o">=</span><span class="n">spglib</span><span class="o">.</span><span class="n">get_spacegroup</span><span class="p">(</span><span class="n">relax</span><span class="o">.</span><span class="n">minim</span><span class="o">.</span><span class="n">dyn</span><span class="o">.</span><span class="n">structure</span><span class="o">.</span><span class="n">get_ase_atoms</span><span class="p">(),</span>
                                      <span class="mf">0.005</span><span class="p">)</span>
    <span class="nb">print</span><span class="p">(</span><span class="s1">&#39;Current SG = &#39;</span><span class="p">,</span> <span class="n">symm</span><span class="p">,</span><span class="s1">&#39; at T=&#39;</span><span class="p">,</span><span class="nb">int</span><span class="p">(</span><span class="n">Temperature</span><span class="p">))</span>

    <span class="c1"># Recompute the ensemble for the hessian calculation</span>
    <span class="n">ensemble</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Ensemble</span><span class="o">.</span><span class="n">Ensemble</span><span class="p">(</span><span class="n">relax</span><span class="o">.</span><span class="n">minim</span><span class="o">.</span><span class="n">dyn</span><span class="p">,</span> <span class="n">T0</span> <span class="o">=</span> <span class="n">Temperature</span><span class="p">,</span>
                            <span class="n">supercell</span> <span class="o">=</span> <span class="n">dyn</span><span class="o">.</span><span class="n">GetSupercell</span><span class="p">())</span>
    <span class="n">ensemble</span><span class="o">.</span><span class="n">generate</span><span class="p">(</span><span class="n">configurations</span><span class="p">,</span> <span class="n">sobol</span> <span class="o">=</span> <span class="n">sobol</span><span class="p">,</span>
                                      <span class="n">sobol_scramble</span> <span class="o">=</span> <span class="n">sobol_scatter</span><span class="p">)</span>
    <span class="n">ensemble</span><span class="o">.</span><span class="n">get_energy_forces</span><span class="p">(</span><span class="n">ff_calculator</span><span class="p">,</span> <span class="n">compute_stress</span> <span class="o">=</span> <span class="kc">False</span><span class="p">)</span>
    <span class="c1">#gets the energies and forces from ff_calculator</span>

    <span class="c1">#update weights!!!</span>
    <span class="n">ensemble</span><span class="o">.</span><span class="n">update_weights</span><span class="p">(</span><span class="n">relax</span><span class="o">.</span><span class="n">minim</span><span class="o">.</span><span class="n">dyn</span><span class="p">,</span> <span class="n">Temperature</span><span class="p">)</span>
    <span class="c1"># Get the free energy hessian</span>
    <span class="n">dyn_hessian</span> <span class="o">=</span> <span class="n">ensemble</span><span class="o">.</span><span class="n">get_free_energy_hessian</span><span class="p">(</span><span class="n">include_v4</span> <span class="o">=</span> <span class="kc">False</span><span class="p">)</span>
    <span class="c1">#free energy hessian as in Bianco paper 2017</span>
    <span class="n">dyn_hessian</span><span class="o">.</span><span class="n">save_qe</span><span class="p">(</span><span class="s2">&quot;hessian_T</span><span class="si">{}</span><span class="s2">_&quot;</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="nb">int</span><span class="p">(</span><span class="n">Temperature</span><span class="p">)))</span>

    <span class="c1"># Get the lowest frequencies for the sscha and the free energy hessian</span>
    <span class="n">w_sscha</span><span class="p">,</span> <span class="n">pols_sscha</span> <span class="o">=</span> <span class="n">relax</span><span class="o">.</span><span class="n">minim</span><span class="o">.</span><span class="n">dyn</span><span class="o">.</span><span class="n">DiagonalizeSupercell</span><span class="p">()</span> <span class="c1">#dynamical matrix</span>
    <span class="c1"># Get the structure in the supercell</span>
    <span class="n">superstructure</span> <span class="o">=</span> <span class="n">relax</span><span class="o">.</span><span class="n">minim</span><span class="o">.</span><span class="n">dyn</span><span class="o">.</span><span class="n">structure</span><span class="o">.</span><span class="n">generate_supercell</span><span class="p">(</span><span class="n">relax</span><span class="o">.</span><span class="n">minim</span><span class="o">.</span><span class="n">dyn</span><span class="o">.</span><span class="n">GetSupercell</span><span class="p">())</span>

    <span class="c1"># Discard the acoustic modes</span>
    <span class="n">acoustic_modes</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Methods</span><span class="o">.</span><span class="n">get_translations</span><span class="p">(</span><span class="n">pols_sscha</span><span class="p">,</span>
                                  <span class="n">superstructure</span><span class="o">.</span><span class="n">get_masses_array</span><span class="p">())</span>
    <span class="n">w_sscha</span> <span class="o">=</span> <span class="n">w_sscha</span><span class="p">[</span><span class="o">~</span><span class="n">acoustic_modes</span><span class="p">]</span>

    <span class="n">lowest_sscha_mode</span><span class="o">.</span><span class="n">append</span><span class="p">(</span><span class="n">np</span><span class="o">.</span><span class="n">min</span><span class="p">(</span><span class="n">w_sscha</span><span class="p">)</span> <span class="o">*</span> <span class="n">CC</span><span class="o">.</span><span class="n">Units</span><span class="o">.</span><span class="n">RY_TO_CM</span><span class="p">)</span> <span class="c1"># Convert from Ry to cm-1</span>

    <span class="n">w_hessian</span><span class="p">,</span> <span class="n">pols_hessian</span> <span class="o">=</span> <span class="n">dyn_hessian</span><span class="o">.</span><span class="n">DiagonalizeSupercell</span><span class="p">()</span> <span class="c1">#recomputed dyn for hessian</span>
    <span class="c1"># Discard the acoustic modes</span>
    <span class="n">acoustic_modes</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Methods</span><span class="o">.</span><span class="n">get_translations</span><span class="p">(</span><span class="n">pols_hessian</span><span class="p">,</span>
                                    <span class="n">superstructure</span><span class="o">.</span><span class="n">get_masses_array</span><span class="p">())</span>
    <span class="n">w_hessian</span> <span class="o">=</span> <span class="n">w_hessian</span><span class="p">[</span><span class="o">~</span><span class="n">acoustic_modes</span><span class="p">]</span>
    <span class="n">lowest_hessian_mode</span><span class="o">.</span><span class="n">append</span><span class="p">(</span><span class="n">np</span><span class="o">.</span><span class="n">min</span><span class="p">(</span><span class="n">w_hessian</span><span class="p">)</span> <span class="o">*</span> <span class="n">CC</span><span class="o">.</span><span class="n">Units</span><span class="o">.</span><span class="n">RY_TO_CM</span><span class="p">)</span> <span class="c1"># Convert from Ry to cm-1</span>
    <span class="c1">#print (&quot;\n&quot;.join([&quot;{:.4f} cm-1&quot;.format(w * CC.Units.RY_TO_CM) for w in pols_hessian]))</span>
    <span class="c1">#exit()</span>

    <span class="n">t_old</span> <span class="o">=</span> <span class="n">Temperature</span>
<span class="c1"># We prepare now the file to save the results</span>
<span class="n">freq_data</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">zeros</span><span class="p">(</span> <span class="p">(</span><span class="nb">len</span><span class="p">(</span><span class="n">temperatures</span><span class="p">),</span> <span class="mi">3</span><span class="p">))</span>
<span class="n">freq_data</span><span class="p">[:,</span> <span class="mi">0</span><span class="p">]</span> <span class="o">=</span> <span class="n">temperatures</span>
<span class="n">freq_data</span><span class="p">[:,</span> <span class="mi">1</span><span class="p">]</span> <span class="o">=</span> <span class="n">lowest_sscha_mode</span>
<span class="n">freq_data</span><span class="p">[:,</span> <span class="mi">2</span><span class="p">]</span> <span class="o">=</span> <span class="n">lowest_hessian_mode</span>

<span class="c1"># Save results on file</span>
<span class="n">np</span><span class="o">.</span><span class="n">savetxt</span><span class="p">(</span><span class="s2">&quot;</span><span class="si">{}</span><span class="s2">_hessian_vs_temperature.dat&quot;</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="n">configurations</span><span class="p">),</span>
            <span class="n">freq_data</span><span class="p">,</span> <span class="n">header</span> <span class="o">=</span> <span class="s2">&quot;T [K]; SSCHA mode [cm-1]; Free energy hessian [cm-1]&quot;</span><span class="p">)</span>
</pre></div>
</div>
<ol class="arabic simple" start="4">
<li><p>Finally we make a graphic output of the data.</p></li>
</ol>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="n">hessian_data</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">loadtxt</span><span class="p">(</span><span class="s2">&quot;</span><span class="si">{}</span><span class="s2">_hessian_vs_temperature.dat&quot;</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="n">configurations</span><span class="p">))</span>

<span class="n">plt</span><span class="o">.</span><span class="n">figure</span><span class="p">(</span><span class="n">dpi</span> <span class="o">=</span> <span class="mi">120</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="n">hessian_data</span><span class="p">[:,</span><span class="mi">0</span><span class="p">],</span> <span class="n">hessian_data</span><span class="p">[:,</span><span class="mi">1</span><span class="p">],</span>
                          <span class="n">label</span> <span class="o">=</span> <span class="s2">&quot;Min SCHA freq&quot;</span><span class="p">,</span> <span class="n">marker</span> <span class="o">=</span> <span class="s2">&quot;&gt;&quot;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="n">hessian_data</span><span class="p">[:,</span><span class="mi">0</span><span class="p">],</span> <span class="n">hessian_data</span><span class="p">[:,</span><span class="mi">2</span><span class="p">],</span>
                          <span class="n">label</span> <span class="o">=</span> <span class="s2">&quot;Free energy curvature&quot;</span><span class="p">,</span> <span class="n">marker</span> <span class="o">=</span> <span class="s2">&quot;o&quot;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">axhline</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="mi">0</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="n">color</span> <span class="o">=</span> <span class="s2">&quot;k&quot;</span><span class="p">,</span> <span class="n">ls</span> <span class="o">=</span> <span class="s2">&quot;dotted&quot;</span><span class="p">)</span> <span class="c1"># Draw the zero</span>
<span class="n">plt</span><span class="o">.</span><span class="n">xlabel</span><span class="p">(</span><span class="s2">&quot;Temperature [K]&quot;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">ylabel</span><span class="p">(</span><span class="s2">&quot;Frequency [cm-1]&quot;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">legend</span><span class="p">()</span>
<span class="n">plt</span><span class="o">.</span><span class="n">tight_layout</span><span class="p">()</span>
<span class="n">plt</span><span class="o">.</span><span class="n">savefig</span><span class="p">(</span><span class="s1">&#39;</span><span class="si">{}</span><span class="s1">_Temp_Freq.png&#39;</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="n">configurations</span><span class="p">))</span>
<span class="c1">#plt.show()</span>

<span class="n">plt</span><span class="o">.</span><span class="n">figure</span><span class="p">(</span><span class="n">dpi</span> <span class="o">=</span> <span class="mi">120</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="n">hessian_data</span><span class="p">[:,</span><span class="mi">0</span><span class="p">],</span> <span class="n">np</span><span class="o">.</span><span class="n">sign</span><span class="p">(</span><span class="n">hessian_data</span><span class="p">[:,</span><span class="mi">2</span><span class="p">])</span> <span class="o">*</span> <span class="n">hessian_data</span><span class="p">[:,</span><span class="mi">2</span><span class="p">]</span><span class="o">**</span><span class="mi">2</span><span class="p">,</span>
                    <span class="n">label</span> <span class="o">=</span> <span class="s2">&quot;Free energy curvature&quot;</span><span class="p">,</span> <span class="n">marker</span> <span class="o">=</span> <span class="s2">&quot;o&quot;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">axhline</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="mi">0</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="n">color</span> <span class="o">=</span> <span class="s2">&quot;k&quot;</span><span class="p">,</span> <span class="n">ls</span> <span class="o">=</span> <span class="s2">&quot;dotted&quot;</span><span class="p">)</span> <span class="c1"># Draw the zero</span>
<span class="n">plt</span><span class="o">.</span><span class="n">xlabel</span><span class="p">(</span><span class="s2">&quot;Temperature [K]&quot;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">ylabel</span><span class="p">(</span><span class="s2">&quot;$\omega^2$ [cm-2]&quot;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">legend</span><span class="p">()</span>
<span class="n">plt</span><span class="o">.</span><span class="n">tight_layout</span><span class="p">()</span>
<span class="n">plt</span><span class="o">.</span><span class="n">savefig</span><span class="p">(</span><span class="s1">&#39;</span><span class="si">{}</span><span class="s1">_Temp_Omeg.png&#39;</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="n">configurations</span><span class="p">))</span>
<span class="c1">#plt.show()</span>
</pre></div>
</div>
<p>We will simulate the temperatures up to room temperature (300 K) with steps of 50 K. Note, this will perform all the steps above 6 times, so it may take some minutes, depending on the PC (on a i3 from 2015, with one core, it took 2 hours).
If it takes too long you can reduce the number of steps by changing the temperature array in <em>Temperature_i = np.linspace(50, 300, 6)</em>.</p>
<figure class="align-default" id="id4">
<span id="fig-results1"></span><a class="reference internal image-reference" href="../figures_03/5000_Temp_Freq.png"><img alt="Freq. vs. Temps.." src="../figures_03/5000_Temp_Freq.png" style="width: 400px;" /></a>
<figcaption>
<p><span class="caption-number">Fig. 8 </span><span class="caption-text">Frequencies versus Temperatures</span><a class="headerlink" href="#id4" title="Permalink to this image"> </a></p>
</figcaption>
</figure>
<p>In <a class="reference internal" href="#fig-results1"><span class="std std-ref">Frequencies versus Temperatures</span></a> we can see that the phase transition is between 100K and 150K. We see that the data points do not drawn a linear figure. We can increase the number of Temperature points to locate the exact transition temperature, but there is another better way to find it.</p>
<figure class="align-default" id="id5">
<span id="fig-results2"></span><a class="reference internal image-reference" href="../figures_03/5000_Temp_Omeg.png"><img alt="Freqs^2 vs Temps.." src="../figures_03/5000_Temp_Omeg.png" style="width: 400px;" /></a>
<figcaption>
<p><span class="caption-number">Fig. 9 </span><span class="caption-text">squared Frequencies versus Temperatures.</span><a class="headerlink" href="#id5" title="Permalink to this image"> </a></p>
</figcaption>
</figure>
<p>For the Landau theory of phase transition, since the SSCHA is a mean-field approach, we expect that around the transition the critical exponent of the temperature goes as</p>
<div class="math notranslate nohighlight">
\[\omega \sim \sqrt{\Phi}\]</div>
<p>For this reason is better to plot the temperature versus the square of the frequency as in <a class="reference internal" href="#fig-results2"><span class="std std-ref">squared Frequencies versus Temperatures.</span></a>
This makes the graph lineal and so we can easily estimate the critic temperature by linear interpolation.</p>
<p>We are using only 50 configurations in the ensemble. Note that this makes a fast calculation but is a low number for this calculations because the free energy calculations are more noisy than the SSCHA frequencies. This is due to the fact that the computation of the free energy requires the third order force constant tensor, and that requires more configurations to converge.</p>
<div class="topic">
<p class="topic-title">Exercise</p>
<p>How the calculation of the free energy changes with the number of configurations?</p>
</div>
<figure class="align-default" id="id6">
<span id="fig-conf-freq"></span><a class="reference internal image-reference" href="../figures_03/Conf_Freq.png"><img alt="Freq. vs. Confs.." src="../figures_03/Conf_Freq.png" style="width: 400px;" /></a>
<figcaption>
<p><span class="caption-number">Fig. 10 </span><span class="caption-text">Evolution of the lowest <em>soft</em> frequency in relation to the number of configurations in the ensemble with a stable configuration. The line is the media and the shade is the standard deviation.</span><a class="headerlink" href="#id6" title="Permalink to this image"> </a></p>
</figcaption>
</figure>
<div class="topic">
<p class="topic-title">Exercise</p>
<p>Plot the Hessian phonon dispersion</p>
<figure class="align-default" id="fig-dispersion-hessian">
<a class="reference internal image-reference" href="../figures_03/dispersion.png"><img alt="Hessian phonon dispersion" src="../figures_03/dispersion.png" style="width: 400px;" /></a>
</figure>
</div>
<figure class="align-default" id="id7">
<span id="diagram"></span><a class="reference internal image-reference" href="../figures_03/Hands-on3.png"><img alt="Diagram." src="../figures_03/Hands-on3.png" style="width: 400px;" /></a>
<figcaption>
<p><span class="caption-number">Fig. 11 </span><span class="caption-text">Workflow of the SSCHA objects for: A) Free energy minimization; B) Structural instabilities search; C) Temperature loop for second order phase transition code. Dotted lines are functions within objects. The dotted and dashed lines indicate the relationship of the dynamic matrix to the objects.</span><a class="headerlink" href="#id7" title="Permalink to this image"> </a></p>
</figcaption>
</figure>
</section>

