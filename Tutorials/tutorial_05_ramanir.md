---
layout: page
title: Raman and Infrared spectra with the Time-Dependent Self Consistent Harmonic Approximation
---

This tutorial was prepared for the [2023 SSCHA School](http://sscha.eu/Schools/2023/home/) by Lorenzo Monacelli. You can see here the video os the hands-on session:

<iframe width="560" height="315" src="https://www.youtube.com/embed/_6pJ4uI0UKc" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>

<p>In the previous tutorial, you learned how to compute the spectral function by integrating the bubble in the Fourier space, with the dynamical ansatz formulated by <a class="reference external" href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.014111">Bianco et al, Physical Review B, 96 , 014111, 2017</a>.
Instead, we will employ the Lanczos algorithm within the Time-Dependent Self-Consistent Harmonic Approximation (TD-SCHA) <a class="reference external" href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.103.104305">Monacelli, Mauri, Physical Review B 103, 104305, 2021</a>.</p>
<p>For this reason, we need the package <code class="docutils literal notranslate"><span class="pre">tdscha</span></code> (it is suggested to configure it with the Julia speedup to run faster, see the installation guide).</p>
<section id="computing-the-ir-signal-in-ice">
<h2>Computing the IR signal in ICE<a class="headerlink" href="#computing-the-ir-signal-in-ice" title="Permalink to this headline"> </a></h2>
<p>We use an ensemble already computed of the phase XI of ice (low-temperature ice ad ambient pressure and prototype of standard cubic ice) to get the IR spectrum.</p>
<p>Inside the directory data, we find an already calculated ensemble of ice XI at 0K with the corresponding original dynamical matrix <em>start_dyn_ice1</em> employed to generate the ensemble and the dynamical matrix <em>final_dyn_ice1</em> after the SSCHA minimization.</p>
<section id="an-introduction">
<h3>An introduction<a class="headerlink" href="#an-introduction" title="Permalink to this headline"> </a></h3>
<p>The infrared spectrum is related to the dipole-dipole response function:</p>
<div class="math notranslate nohighlight">
\[\chi_{MM}(\omega) = \int_{-\infty}^{\infty}dt e^{-i\omega t}\left&lt;M(t) M(0)\right&gt;\]</div>
<p>where the average <span class="math notranslate nohighlight">\(\left&lt;M(t)M(0)\right&gt;\)</span> is the quantum average at finite temperature.</p>
<p>Exploiting the TD-SCHA formalism introduced in the previous lecture, this response function can be written as:</p>
<div class="math notranslate nohighlight" id="equation-eqchi">
<span class="eqno">(1)<a class="headerlink" href="#equation-eqchi" title="Permalink to this equation"> </a></span>\[\chi_{MM}(\omega) = \boldsymbol{r}(M) \boldsymbol{G}(\omega) \boldsymbol{q}(M)\]</div>
<p>where <span class="math notranslate nohighlight">\(\boldsymbol{G}(\omega)\)</span> is the TD-SCHA green function, while the <span class="math notranslate nohighlight">\(\boldsymbol{r}\)</span> and <span class="math notranslate nohighlight">\(\boldsymbol{q}\)</span> are vectors that quantify the perturbation and response, respectively.</p>
<p>In particular, if we neglect two-phonon effects (nonlinear coupling with light), we get that</p>
<div class="math notranslate nohighlight">
\[\chi_{MM}(\omega) = \sum_{ab}\frac{\mathcal Z_{\alpha a} {\mathcal Z}_{\alpha b}}{\sqrt{m_am_b}} G_{ab}(\omega)\]</div>
<p>where <span class="math notranslate nohighlight">\({\mathcal Z}_{\alpha a}\)</span> is the Born effective charge of atom <span class="math notranslate nohighlight">\(a\)</span>, with polarization <span class="math notranslate nohighlight">\(\alpha\)</span>, and <span class="math notranslate nohighlight">\(G_{ab}(\omega)\)</span> is the one-phonon green function, (its imaginary part is precisely the spectral function).</p>
<p>Indeed, we need to compute the effective charges. This can be done directly by quantum espresso using linear response theory (ph.x).</p>
<div class="topic">
<p class="topic-title">Exercise</p>
<p>Use the knowledge of cellconstructor to extract a structure file from the final dynamical matrix to submit the calculation of the dielectric tensor, Effective charges, and Raman tensor in quantum espresso.</p>
<p>Hint. The structure is the attribute <em>structure</em> of the Phonons object.
The structure in the SCF file can be saved with the <em>save_scf</em> method of the Structure object.</p>
<p>You can then attach the structure to the header of the espresso
<em>ir_raman_header.pwi</em>.</p>
<p>Notice that we are using norm-conserving pseudo-potentials and LDA exchange-correlation functional,
as the Raman Tensor in quantum espresso is implemented only with them.
However, it is usually an excellent approximation.
<em>ir_raman_header.pwi</em>.</p>
<p>You must run the pw.x code and the ph.x code (<em>ir_raman_complete.phi</em>), which performs the phonon calculation.</p>
<p>We provide the final output file in <em>ir_raman_complete.pho</em></p>
</div>
</section>
<section id="prepare-the-infrared-response">
<h3>Prepare the infrared response<a class="headerlink" href="#prepare-the-infrared-response" title="Permalink to this headline"> </a></h3>
<p>We need to attach the Raman Tensor and effective charges computed inside <em>ir_raman_complete.pho</em> to the final dynamical matrix,
we will use this to initialize the response function calculation, as in <a class="reference internal" href="#equation-eqchi">Eq.1</a>.</p>
<p>To attach the content of an espresso ph calculation (only Dielectric tensor, Raman Tensor, and Born effective charges)
to a specific dynamical matrix, use</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="n">dyn</span><span class="o">.</span><span class="n">ReadInfoFromESPRESSO</span><span class="p">(</span><span class="s2">&quot;ir_raman_complete.pho&quot;</span><span class="p">)</span>
</pre></div>
</div>
<p>If you save the dynamical matrix in quantum espresso format, before the frequencies and the diagonalization, there will be
the Dielectric tensor</p>
<div class="highlight-text notranslate"><div class="highlight"><pre><span></span>Dielectric Tensor:

     1.890128098000           0.000000000000           0.000000000000
     0.000000000000           1.912811137000           0.000000000000
     0.000000000000           0.000000000000           1.916728724000
</pre></div>
</div>
<p>Followed by the effective charges and the Raman tensor.</p>
</section>
<section id="submitting-the-ir-calculation">
<h3>Submitting the IR calculation<a class="headerlink" href="#submitting-the-ir-calculation" title="Permalink to this headline"> </a></h3>
<p>With the following script, we submit a TD-SCHA calculation for the IR.</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">numpy</span> <span class="k">as</span> <span class="nn">np</span>
<span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span><span class="o">,</span> <span class="nn">cellconstructor.Phonons</span>
<span class="kn">import</span> <span class="nn">sscha</span><span class="o">,</span> <span class="nn">sscha.Ensemble</span>
<span class="kn">import</span> <span class="nn">tdscha</span><span class="o">,</span> <span class="nn">tdscha.DynamicalLanczos</span> <span class="k">as</span> <span class="nn">DL</span>

<span class="c1"># Load the starting dynamical matrix</span>
<span class="n">dyn_start</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="s2">&quot;start_dyn_ice&quot;</span><span class="p">)</span>

<span class="c1"># Load the ensemble</span>
<span class="n">temperature</span> <span class="o">=</span> <span class="mi">0</span> <span class="c1"># K</span>
<span class="n">population</span> <span class="o">=</span> <span class="mi">2</span>
<span class="n">n_configs</span> <span class="o">=</span> <span class="mi">10000</span>

<span class="n">ensemble</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Ensemble</span><span class="o">.</span><span class="n">Ensemble</span><span class="p">(</span><span class="n">dyn_start</span><span class="p">,</span> <span class="n">temperature</span><span class="p">)</span>
<span class="n">ensemble</span><span class="o">.</span><span class="n">load</span><span class="p">(</span><span class="s2">&quot;data&quot;</span><span class="p">,</span> <span class="n">population</span><span class="p">,</span> <span class="n">n_configs</span><span class="p">)</span>

<span class="c1"># Load the final dynamical matrix</span>
<span class="n">final_dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="s2">&quot;final_dyn_ice&quot;</span><span class="p">)</span>
<span class="n">final_dyn</span><span class="o">.</span><span class="n">ReadInfoFromESPRESSO</span><span class="p">(</span><span class="s2">&quot;ir_raman_complete.pho&quot;</span><span class="p">)</span>

<span class="c1"># Update the ensemble weights for the final dynamical matrix</span>
<span class="n">ensemble</span><span class="o">.</span><span class="n">update_weights</span><span class="p">(</span><span class="n">final_dyn</span><span class="p">,</span> <span class="n">temperature</span><span class="p">)</span>

<span class="c1"># Setup the TD-SCHA calculation with the Lanczos algorithm</span>
<span class="n">lanczos</span> <span class="o">=</span> <span class="n">DL</span><span class="o">.</span><span class="n">Lanczos</span><span class="p">(</span><span class="n">ensemble</span><span class="p">)</span>
<span class="n">lanczos</span><span class="o">.</span><span class="n">ignore_v3</span> <span class="o">=</span> <span class="kc">True</span>
<span class="n">lanczos</span><span class="o">.</span><span class="n">ignore_v4</span> <span class="o">=</span> <span class="kc">True</span>

<span class="c1"># If you have julia-enabled tdscha installed uncomment</span>
<span class="c1"># lanczos.mode = DL.MODE_FAST_JULIA</span>
<span class="c1"># for a x10-x15 speedup.</span>

<span class="n">lanczos</span><span class="o">.</span><span class="n">init</span><span class="p">()</span>


<span class="c1"># Setup the IR response</span>
<span class="n">polarization</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">array</span><span class="p">([</span><span class="mi">1</span><span class="p">,</span><span class="mi">0</span><span class="p">,</span><span class="mi">0</span><span class="p">])</span>  <span class="c1"># Polarization of light</span>
<span class="n">lanczos</span><span class="o">.</span><span class="n">prepare_ir</span><span class="p">(</span><span class="n">pol_vec</span> <span class="o">=</span> <span class="n">polarization</span><span class="p">)</span>


<span class="c1"># Run the algorithm</span>
<span class="n">n_iterations</span> <span class="o">=</span> <span class="mi">1000</span>
<span class="n">lanczos</span><span class="o">.</span><span class="n">run_FT</span><span class="p">(</span><span class="n">n_iterations</span><span class="p">)</span>
<span class="n">lanczos</span><span class="o">.</span><span class="n">save_status</span><span class="p">(</span><span class="s2">&quot;ir_xpol&quot;</span><span class="p">)</span>
</pre></div>
</div>
<p><strong>Congratulations!</strong> You ran your first TD-SCHA calculation.
You can plot the results by using:</p>
<div class="highlight-console notranslate"><div class="highlight"><pre><span></span><span class="go">tdscha-plot.py ir_xpol.npz</span>
</pre></div>
</div>
<p>The script <code class="docutils literal notranslate"><span class="pre">tdscha-plot.py</span></code> is automatically installed with the tdscha package.</p>
<figure class="align-default" id="ir-scha">
<a class="reference internal image-reference" href="../figures_05/IR_v2.png"><img alt="../figures_05/IR_v2.png" src="../figures_05/IR_v2.png" style="width: 50%;" /></a>
<figcaption>
<p><span class="caption-number">Fig. 12 </span><span class="caption-text">IR spectrum with both <em>include_v3</em> and <em>include_v4</em> set to False.</span><a class="headerlink" href="#ir-scha" title="Permalink to this image"> </a></p>
</figcaption>
</figure>
<p>Additionally, <code class="docutils literal notranslate"><span class="pre">tdscha-plot.py</span></code> takes three more parameters: the range of the frequencies to be displayed and the smearing.</p>
</section>
<section id="deep-dive-into-the-calculation">
<h3>Deep dive into the calculation<a class="headerlink" href="#deep-dive-into-the-calculation" title="Permalink to this headline"> </a></h3>
<p>Let us dive a bit into the calculation. The beginning of the script should
be almost self-explanatory, as we are just loading dynamical matrices, dielectric tensors, and effective charges.</p>
<p>The line</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="n">ensemble</span><span class="o">.</span><span class="n">update_weights</span><span class="p">(</span><span class="n">final_dyn</span><span class="p">,</span> <span class="n">temperature</span><span class="p">)</span>
</pre></div>
</div>
<p>deserves special attention. Here, we are changing the weights of the configurations inside the ensemble to simulate the specified dynamical matrix and temperature, even if they differ from those used to generate the ensemble.
This is useful to compute the spectrum at several temperatures without extracting and calculating a new ensemble each time.</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="c1"># Setup the TD-SCHA calculation with the Lanczos algorithm</span>
<span class="n">lanczos</span> <span class="o">=</span> <span class="n">DL</span><span class="o">.</span><span class="n">Lanczos</span><span class="p">(</span><span class="n">ensemble</span><span class="p">)</span>
<span class="n">lanczos</span><span class="o">.</span><span class="n">ignore_v3</span> <span class="o">=</span> <span class="kc">True</span>
<span class="n">lanczos</span><span class="o">.</span><span class="n">ignore_v4</span> <span class="o">=</span> <span class="kc">True</span>
<span class="n">lanczos</span><span class="o">.</span><span class="n">init</span><span class="p">()</span>
</pre></div>
</div>
<p>Then we initialize the Lanczos algorithm for the tdscha, passing the ensemble.</p>
<p>The ignore_v3 and ignore_v4 are flags that, if set to True, the 3-phonon and 4-phonon scattering will be ignored during the calculation.</p>
<p>As you can see from the output, our IR signal had very sharp peaks because we ignored any phonon-phonon scattering process that may give rise to a finite lifetime.</p>
<p>By setting only ignore_v4 to True, we reproduce the behavior of the bubble approximation.
Notably, while the four-phonon scattering is exceptionally computationally and memory demanding in free energy hessian calculations, within the Lanczos algorithm, accounting for the four-phonon scattering is only a factor two more expensive than using just the third order, without requiring any additional memory.</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="c1"># Setup the IR response</span>
<span class="n">polarization</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">array</span><span class="p">([</span><span class="mi">1</span><span class="p">,</span><span class="mi">0</span><span class="p">,</span><span class="mi">0</span><span class="p">])</span>  <span class="c1"># Polarization of light</span>
<span class="n">lanczos</span><span class="o">.</span><span class="n">prepare_ir</span><span class="p">(</span><span class="n">pol_vec</span> <span class="o">=</span> <span class="n">polarization</span><span class="p">)</span>
</pre></div>
</div>
<p>Here we tell the Lanczos which kind of calculation we want to do. In other words, we set the <span class="math notranslate nohighlight">\(\boldsymbol{r}\)</span> and <span class="math notranslate nohighlight">\(\boldsymbol{q}\)</span> vectors in <a class="reference internal" href="#equation-eqchi">Eq.1</a> for the Lanczos calculation.
- prepare_ir
- prepare_raman
- prepare_mode
- prepare_perturbation</p>
<p>The names are intuitive; besides the Raman and IR, prepare_mode allows you to study the response function of a specific phonon mode, and prepare_perturbation enables defining a custom perturbation function.</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="c1"># Run the algorithm</span>
<span class="n">n_iterations</span> <span class="o">=</span> <span class="mi">1000</span>
<span class="n">lanczos</span><span class="o">.</span><span class="n">run_FT</span><span class="p">(</span><span class="n">n_iterations</span><span class="p">)</span>
<span class="n">lanczos</span><span class="o">.</span><span class="n">save_status</span><span class="p">(</span><span class="s2">&quot;ir_xpol&quot;</span><span class="p">)</span>
</pre></div>
</div>
<p>Here we start the calculation of the response function.
The number of iterations indicates how many Lanczos steps are required. Each step adds a new pole to the green function. Therefore, many steps are necessary to converge broad spectrum features, while much less if the peaks are sharp.
We save the status in such a way that we can get it back later.</p>
<p>Last, the commented line</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="n">lanczos</span><span class="o">.</span><span class="n">mode</span> <span class="o">=</span> <span class="n">DL</span><span class="o">.</span><span class="n">MODE_FAST_JULIA</span>
</pre></div>
</div>
<p>This line only works if Julia and PyCall are correctly set up in the PC; in that case, run the script with <em>python-jl</em> instead of python. It will exploit a massive speedup of a factor between 10x and 15x.
The calculation can also be run in parallel using <em>mpirun</em> before calling the Python executable (or python-jl). In this case, to work correctly, you should have mpi4py installed and working.</p>
<div class="topic">
<p class="topic-title">Exercise</p>
<p>Compute the Lanczos with the bubble approximation and without any approximation, and check the differences.</p>
</div>
<figure class="align-default" id="irv3">
<a class="reference internal image-reference" href="../figures_05/IR_v3.png"><img alt="../figures_05/IR_v3.png" src="../figures_05/IR_v3.png" style="width: 50%;" /></a>
<figcaption>
<p><span class="caption-number">Fig. 13 </span><span class="caption-text">IR signal accounting for the three-phonon scattering</span><a class="headerlink" href="#irv3" title="Permalink to this image"> </a></p>
</figcaption>
</figure>
<figure class="align-default" id="irv4">
<a class="reference internal image-reference" href="../figures_05/IR_v4.png"><img alt="../figures_05/IR_v4.png" src="../figures_05/IR_v4.png" style="width: 50%;" /></a>
<figcaption>
<p><span class="caption-number">Fig. 14 </span><span class="caption-text">IR signal accounting for all anharmonic scattering.
The peaks that appear slightly below 2500 cm-1 is a
combination mode known to be present in ice. See <a class="reference external" href="https://pubs.aip.org/aip/jcp/article-abstract/155/18/184502/199619/The-microscopic-origin-of-the-anomalous-isotopic?redirectedFrom=fulltext">Cherubini et al, J Chem Phys 155, 184502, 2021</a></span><a class="headerlink" href="#irv4" title="Permalink to this image"> </a></p>
</figcaption>
</figure>
<div class="topic" id="exercize-polarization-ir">
<p class="topic-title">Exercise</p>
<p>Try to see how different polarization of the light affects the result.</p>
</div>
</section>
<section id="analyze-the-output">
<h3>Analyze the output<a class="headerlink" href="#analyze-the-output" title="Permalink to this headline"> </a></h3>
<p>In the last part, we employed the script <code class="docutils literal notranslate"><span class="pre">tdscha-plot.py</span></code> to display the simulation result. This is a quick way to show the results of a calculation.</p>
<p>Here, we will dive deeper into the calculation output file to extract the response function and get the results.</p>
<p>The Lanczos algorithm provides a set of coefficients <span class="math notranslate nohighlight">\(a_i\)</span>, <span class="math notranslate nohighlight">\(b_i\)</span>, and <span class="math notranslate nohighlight">\(c_i\)</span> through which the green function is evaluated thanks to a continued fraction:</p>
<div class="math notranslate nohighlight">
\[G(\omega) = \frac{1}{a_1 - (\omega + i\eta)^2 +  \frac{b_1c_1}{a_2 - (\omega+i\eta)^2 + \frac{c_2b_2}{a_3 - \cdots}}}\]</div>
<p>Each iteration of the algorithm adds a new set of coefficients written in the standard output.
Thanks to this expression, we only need the series of coefficients to compute the dynamical Green function at any frequency and with any smearing.
The Green function can be computed with:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="n">green_function</span> <span class="o">=</span> <span class="n">lanczos</span><span class="o">.</span><span class="n">get_green_function_continued_fraction</span><span class="p">(</span><span class="n">frequencies</span><span class="p">,</span> <span class="n">smearing</span><span class="o">=</span><span class="n">smearing</span><span class="p">)</span>
</pre></div>
</div>
<p>Here <code class="docutils literal notranslate"><span class="pre">frequencies</span></code> is an array in Rydberg.
The response function is the opposite of the imaginary part of the green function; thus, to reproduce the plot, we have:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">tdscha</span><span class="o">,</span> <span class="nn">tdscha.DynamicalLanczos</span>
<span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span><span class="o">,</span> <span class="nn">cellconstructor.Units</span>
<span class="kn">import</span> <span class="nn">numpy</span> <span class="k">as</span> <span class="nn">np</span>
<span class="kn">import</span> <span class="nn">matplotlib.pyplot</span> <span class="k">as</span> <span class="nn">plt</span>

<span class="c1"># Load the result of the previous calculation</span>
<span class="n">lanczos</span> <span class="o">=</span> <span class="n">tdscha</span><span class="o">.</span><span class="n">DynamicalLanczos</span><span class="o">.</span><span class="n">Lanczos</span><span class="p">()</span>
<span class="n">lanczos</span><span class="o">.</span><span class="n">load_status</span><span class="p">(</span><span class="s2">&quot;ir_xpol_v4&quot;</span><span class="p">)</span>

<span class="c1"># Get the green function</span>
<span class="n">W_START</span> <span class="o">=</span> <span class="mi">0</span>
<span class="n">W_END</span> <span class="o">=</span> <span class="mi">3700</span>
<span class="n">N_W</span> <span class="o">=</span> <span class="mi">10000</span>
<span class="n">SMEARING</span> <span class="o">=</span> <span class="mi">10</span>

<span class="n">frequencies</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">linspace</span><span class="p">(</span><span class="n">W_START</span><span class="p">,</span> <span class="n">W_END</span><span class="p">,</span> <span class="n">N_W</span><span class="p">)</span>

<span class="c1"># Convert in RY</span>
<span class="n">frequencies_ry</span> <span class="o">=</span> <span class="n">frequencies</span> <span class="o">/</span> <span class="n">CC</span><span class="o">.</span><span class="n">Units</span><span class="o">.</span><span class="n">RY_TO_CM</span>
<span class="n">smearing_ry</span> <span class="o">=</span> <span class="n">SMEARING</span> <span class="o">/</span> <span class="n">CC</span><span class="o">.</span><span class="n">Units</span><span class="o">.</span><span class="n">RY_TO_CM</span>

<span class="c1"># Compute the green function</span>
<span class="n">green_function</span> <span class="o">=</span> <span class="n">lanczos</span><span class="o">.</span><span class="n">get_green_function_continued_fraction</span><span class="p">(</span><span class="n">frequencies_ry</span><span class="p">,</span>
        <span class="n">smearing</span><span class="o">=</span><span class="n">smearing_ry</span><span class="p">)</span>

<span class="c1"># Get the response function</span>
<span class="n">ir_response_function</span> <span class="o">=</span> <span class="o">-</span> <span class="n">np</span><span class="o">.</span><span class="n">imag</span><span class="p">(</span><span class="n">green_function</span><span class="p">)</span>

<span class="c1"># Plot the data</span>
<span class="n">plt</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="n">frequencies</span><span class="p">,</span> <span class="n">ir_response_function</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">show</span><span class="p">()</span>
</pre></div>
</div>
<p>The previous script plots the data, precisely like <em>plot-tdscha.py</em>;
however, now you have full access to the response function, both its imaginary and real parts.</p>
<div class="topic">
<p class="topic-title">Exercise</p>
<p>Plot the IR data at various smearing and as a function of the number of steps (50, 100, 200, 300, and 1000). How does the signal change with smearing and the number of steps? When is it converged?</p>
</div>
</section>
</section>
<section id="raman-response">
<h2>Raman response<a class="headerlink" href="#raman-response" title="Permalink to this headline"> </a></h2>
<p>The Raman response is very similar to the IR.
Raman probes the fluctuations of the polarizability instead of those of the polarization,
and it occurs when the samples interact with two light sources: the incoming electromagnetic
radiation and the outcoming one.
The outcoming radiation has a frequency that is shifted with respect to the incoming one
by the energy of the scattering phonons.
The signal on the red side of the pump is called Stokes, while the signal on the blue side is the Antistokes.
Since the outcoming radiation has higher energy than the incoming one in the Antistokes, it is generated only by existing (thermally excited phonons) inside the sample. Therefore it has a lower intensity than the Stokes.</p>
<p>On the Stokes side, the intensity of the scattered light with a frequency redshift of <span class="math notranslate nohighlight">\(\omega\)</span> is</p>
<div class="math notranslate nohighlight">
\[I(\omega) \propto \left&lt;\alpha_{xy}(\omega)\alpha_{xy}(0)\right&gt; (n(\omega) + 1)\]</div>
<p>where <span class="math notranslate nohighlight">\(\alpha\)</span> is the polarizability along the <span class="math notranslate nohighlight">\(xy\)</span> axis.
We can do a linear expansion around the equilibrium position of the polarizability, and we get:</p>
<div class="math notranslate nohighlight">
\[ \begin{align}\begin{aligned}\alpha_{xy}(\omega) = \sum_{a = 1}^{3N}\frac{\partial \alpha_{xy}}{\partial R_a(\omega)} (R_a(\omega) - \mathcal R_a)\\\alpha_{xy}(\omega) = \sum_{a = 1}^{3N}\Xi_{xya} (R_a(\omega) - \mathcal R_a)\end{aligned}\end{align} \]</div>
<p>If we insert it in the expression of the intensity, the average between the positions is the atomic green function divided by the square root of the masses, and we get</p>
<div class="math notranslate nohighlight">
\[I(\omega) \propto \sum_{ab} \frac{\Xi_{xy a} \Xi_{xy b}}{\sqrt{m_a m_b}} G_{ab}(\omega)(n(\omega) + 1)\]</div>
<p>where <span class="math notranslate nohighlight">\(G_{ab}(\omega)\)</span> is the atomic green function on atoms <span class="math notranslate nohighlight">\(a\)</span> and <span class="math notranslate nohighlight">\(b\)</span>,
while <span class="math notranslate nohighlight">\(\Xi_{xy a}\)</span> is the Raman tensor along the electric fields directed in <span class="math notranslate nohighlight">\(x\)</span> and <span class="math notranslate nohighlight">\(y\)</span> and on atom <span class="math notranslate nohighlight">\(a\)</span>.</p>
<p>The multiplication factor <span class="math notranslate nohighlight">\(n(\omega) + 1\)</span> comes from the observation of the Stokes nonresonant Raman (it would be just <span class="math notranslate nohighlight">\(n(\omega)\)</span> for the antistokes).</p>
<p>As we did for the IR signal, we can prepare the calculation of the Raman raman scattering by computing the polarizability-polarizability.</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="c1"># Setup the polarized Raman response</span>
<span class="n">polarization_in</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">array</span><span class="p">([</span><span class="mi">1</span><span class="p">,</span><span class="mi">0</span><span class="p">,</span><span class="mi">0</span><span class="p">])</span>
<span class="n">polarization_out</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">array</span><span class="p">([</span><span class="mi">1</span><span class="p">,</span><span class="mi">0</span><span class="p">,</span><span class="mi">0</span><span class="p">])</span>
<span class="n">lanczos</span><span class="o">.</span><span class="n">prepare_raman</span><span class="p">(</span><span class="n">pol_vec_in</span><span class="o">=</span><span class="n">polarization_in</span><span class="p">,</span>
        <span class="n">pol_vec_out</span><span class="o">=</span><span class="n">polarization_out</span><span class="p">)</span>
</pre></div>
</div>
<p>Note that here we have to specify two polarization of the light, the incoming radiation, and the outcoming radiation.</p>
<div class="topic">
<p class="topic-title">Exercise</p>
<p>Compute and plot the Intensity of the Raman in the Stokes and antistokes configurations.
Try with different polarization and even orthogonal polarization; what does it change?</p>
<p>The Bose-Einstein factor <span class="math notranslate nohighlight">\(n(\omega)\)</span> can be computed with the following function:</p>
</div>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="c1"># n(w) Bose-Einstein occupation number:</span>
<span class="c1"># w is in Ry, T is in K</span>
<span class="n">n_w</span> <span class="o">=</span> <span class="n">tdscha</span><span class="o">.</span><span class="n">DynamicalLanczos</span><span class="o">.</span><span class="n">bose_occupation</span><span class="p">(</span><span class="n">w</span><span class="p">,</span> <span class="n">T</span><span class="p">)</span>
</pre></div>
</div>
<section id="unpolarize-raman-and-ir">
<h3>Unpolarize Raman and IR<a class="headerlink" href="#unpolarize-raman-and-ir" title="Permalink to this headline"> </a></h3>
<p>In the previous section, we saw how to compute Raman and IR with specific polarization of the incoming and
outcoming radiation, and on oriented crystals (single crystals).
However, the most common situation is a powder sample probed with unpolarized light.</p>
<p>In this case, we need to look at the Raman and IR response for unpolarized samples.
While this is just the average of the IR signal’s x, y, and z, the Raman is more complex. In particular, unpolarized Raman signal can be computed from the so-called <em>invariants</em>, where the perturbations in the polarizations are the following:</p>
<div class="math notranslate nohighlight">
\[ \begin{align}\begin{aligned}I_A = \frac{1}{3}( {xx} + {yy} + {zz})^2/9\\I_{B_1} = ({xx} - {yy})^2 / 2\\I_{B_2} = ({xx} - {zz})^2/2\\I_{B_3} = ({yy} - {zz})^2/2\\I_{B_4} = 3({xy})^2\\I_{B_5} = 3({yz})^2\\I_{B_6} = 3({xz})^2\end{aligned}\end{align} \]</div>
<p>The total Intensity of unpolarized Raman is:</p>
<div class="math notranslate nohighlight">
\[I_{unpol}(\omega) = 45 \cdot I_a(\omega) + 7 \cdot \sum_{i=1}^6 I_{B_i}(\omega)\]</div>
<p>The tdscha code implements a way to compute each perturbation separately. For example, the Raman response related to <span class="math notranslate nohighlight">\(I_A\)</span> is calculated with</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="n">lanczos</span><span class="o">.</span><span class="n">prepare_raman</span><span class="p">(</span><span class="n">unpolarized</span><span class="o">=</span><span class="mi">0</span><span class="p">)</span>
</pre></div>
</div>
<p>While the <span class="math notranslate nohighlight">\(I_{B_i}\)</span> is computed using index <span class="math notranslate nohighlight">\(i\)</span>. For example, to compute <span class="math notranslate nohighlight">\(I_{B_5}\)</span>:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="c1"># To compute I_B5 we do</span>
<span class="n">lanczos</span><span class="o">.</span><span class="n">prepare_raman</span><span class="p">(</span><span class="n">unpolarized</span><span class="o">=</span><span class="mi">5</span><span class="p">)</span>
</pre></div>
</div>
<p>To get the total spectrum, you need to add the scattering factor <span class="math notranslate nohighlight">\(n(\omega) + 1\)</span> and sum all these perturbation with the correct prefactor (45 for <span class="math notranslate nohighlight">\(I_A\)</span> and 7 for the sum of all <span class="math notranslate nohighlight">\(I_B\)</span>).</p>
<p>To reset a calculation and start a new one, you can use</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="n">lanczos</span><span class="o">.</span><span class="n">reset</span><span class="p">()</span>
</pre></div>
</div>
<p>which may be called before preparing the perturbation.</p>
<div class="topic">
<p class="topic-title">Exercise</p>
<p>Compute the unpolarized Raman spectrum of ice and plot the results.</p>
</div>
<figure class="align-default">
<a class="reference internal image-reference" href="../figures_05/raman_unpolarized.png"><img alt="../figures_05/raman_unpolarized.png" src="../figures_05/raman_unpolarized.png" style="width: 50%;" /></a>
</figure>
<p>You should employ a supercell size sufficiently big to converge the simulation properly. In this case, the 1x1x1 supercell is too tiny to converge the calculation and get meaningful results.</p>
</section>
</section>
