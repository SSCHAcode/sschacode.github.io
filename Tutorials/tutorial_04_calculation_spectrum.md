---
layout: page
title: Calculation of spectral properties with the Self Consistent Harmonic Approximation
---

This tutorial was prepared for the [2023 SSCHA School](http://sscha.eu/Schools/2023/home/) by Raffaello Bianco. You can see here the video os the hands-on session:

<iframe width="560" height="315" src="https://www.youtube.com/embed/e-F-EKHA_3M" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>

<section id="theoretical-introduction">
<h2>Theoretical introduction<a class="headerlink" href="#theoretical-introduction" title="Permalink to this headline"> </a></h2>
<p>The SCHA phonons are non-interacting quasiparticles that already include anharmonic effects (i.e. interaction between standard harmonic phonons) at some level. However, anharmoncity causes interaction between  the SCHA phonons too. The interactions between phonons causes a change of their energy spectrum: from the overlap of simple Dirac-delta functions  centered around the SSCHA phonon frequencies, to the overlap of Lorentzians with finite width (i.e. the quasiparticles have finite lifetime) and centered around shifted energies, or structures even more complex  (when the anharmoncity is so strong that the quasiparticle picture has to be abandoned).
For each <span class="math notranslate nohighlight">\(\boldsymbol{q}\)</span> of the Brillouin zone, the SCHA phonons energy spectrum <span class="math notranslate nohighlight">\(\sigma(\boldsymbol{q},\Omega)\)</span> is given by Eq.(70)</p>
<figure class="align-default">
<a class="reference internal image-reference" href="../figures_04/Fig1.png"><img alt="../figures_04/Fig1.png" src="../figures_04/Fig1.png" style="width: 100%;" /></a>
</figure>
<p>where <span class="math notranslate nohighlight">\(\boldsymbol{G}\)</span> is the SCHA phonons Green function, given by
<span class="math notranslate nohighlight">\(\boldsymbol{G}=\boldsymbol{G}^{(0)}+\boldsymbol{G}^{(0)}\boldsymbol{\Pi}\,\,\boldsymbol{G}\)</span>, with
<span class="math notranslate nohighlight">\(\boldsymbol{G}^{(0)}\)</span> the Green function of the noninteracting SCHA phonons, and <span class="math notranslate nohighlight">\(\boldsymbol{\Pi}\)</span> their selfenergy
taking into account the interaction (in order to make easier the comparison with the literature, here and in the subsequent equations,
the equation numbers refer to the paper here <a class="reference external" href="https://arxiv.org/abs/2103.03973">https://arxiv.org/abs/2103.03973</a>). The SSCHA code allows to compute these quantities, even if in the current implementation the selfenergy can be computed only in the bubble approximation (Eq. (75))</p>
<figure class="align-default">
<a class="reference internal image-reference" href="../figures_04/Fig2.png"><img alt="../figures_04/Fig2.png" src="../figures_04/Fig2.png" style="width: 100%;" /></a>
</figure>
<p>i.e. the self-energy terms including the 4th order FCs are discarded. In this equation <span class="math notranslate nohighlight">\(\delta_{\scriptscriptstyle{\text{se}}}\)</span> is an infinitely small positive number (a smearing parameter), that in actual calculations has to be chosen, toegher with the integration <span class="math notranslate nohighlight">\(\boldsymbol{k}\)</span>-grid, in order to find converged results.</p>
<p>The code, in addition to providing the ability to compute the spectral function through the full formula Eq. (70), allows for various approximations to be used in order to both compute the spectral function with reduced computational cost and conduct an analysis of the different contributions to the spectral function provided by each mode. The main approximation is to negeclet the off-diagonal terms in the self-energy written in the SCHA-modes basis set. In other words, we can negelct the possibility that the interaction mixes different SCHA phonons. In that case, the total spectrum  is given by the the sum of the spectrum of each mode, <span class="math notranslate nohighlight">\(\sigma_{\mu}(\boldsymbol{q},\Omega)\)</span>, as shown in Eq.(78),</p>
<figure class="align-default">
<a class="reference internal image-reference" href="../figures_04/Fig3.png"><img alt="../figures_04/Fig3.png" src="../figures_04/Fig3.png" style="width: 100%;" /></a>
</figure>
<p>with  <span class="math notranslate nohighlight">\(\sigma_{\mu}(\boldsymbol{q},\Omega)\)</span> having a generalized Lorentzian-like expression shown in Eq.(79)</p>
<figure class="align-default">
<a class="reference internal image-reference" href="../figures_04/Fig4.png"><img alt="../figures_04/Fig4.png" src="../figures_04/Fig4.png" style="width: 100%;" /></a>
</figure>
<p>with <span class="math notranslate nohighlight">\(\mathcal{Z}_{\mu}(\boldsymbol{q},\Omega)\)</span> defined in Eq.(80)</p>
<figure class="align-default">
<a class="reference internal image-reference" href="../figures_04/Fig5.png"><img alt="../figures_04/Fig5.png" src="../figures_04/Fig5.png" style="width: 100%;" /></a>
</figure>
<p>where <span class="math notranslate nohighlight">\(\omega_{\mu}(\boldsymbol{q})\)</span> is the frequency (energy) of the SCHA phonon <span class="math notranslate nohighlight">\((\boldsymbol{q},\mu)\)</span>, and
<span class="math notranslate nohighlight">\(\boldsymbol{\Pi}_{\mu\mu}(\boldsymbol{q},\Omega)\)</span> is the corresponding diagonal element of the self-energy. This is the spectrum in the
so called no mode-mixing approximation. At this level, the single-mode spectral functions resemble Lorentzian functions, but they are not true Lorentzians as they
have kind of frequency-dependent center and width. As a matter of fact, in general, the spectrum of a mode can be very different from a true Lorentzian function,
meaning that the quasiparticle picture for that mode is not appropriate. However, there are cases where the interaction  between the SCHA phonons does not affect the
quasiparticle picture but causes only a shift in the quasiparticle energy and the appearance of a finite linediwth (i.e. finite lifetime) with respect to the non-interacting case. In that case, we can write the spectral function of the mode as a true Lorentzian, Eq.(81),</p>
<figure class="align-default">
<a class="reference internal image-reference" href="../figures_04/Fig6.png"><img alt="../figures_04/Fig6.png" src="../figures_04/Fig6.png" style="width: 100%;" /></a>
</figure>
<p>i.e., the SCHA phonon <span class="math notranslate nohighlight">\((\boldsymbol{q},\mu)\)</span> is a quasiparticle with definite energy <span class="math notranslate nohighlight">\(\Omega_{\mu}(\boldsymbol{q})\)</span> (
<span class="math notranslate nohighlight">\(\Delta_{\mu}(\boldsymbol{q})=\Omega_{\mu}(\boldsymbol{q})-\omega_{\mu}(\boldsymbol{q})\)</span> is called the energy shift)
and lifetime <span class="math notranslate nohighlight">\(\tau_{\mu}(\boldsymbol{q})=1/2\Gamma_{\mu}(\boldsymbol{q})\)</span>, where <span class="math notranslate nohighlight">\(\Gamma_{\mu}(\boldsymbol{q})\)</span>
is the Lorentzian half width at half maximum (HWHM). The quantities  <span class="math notranslate nohighlight">\(\Omega_{\mu}(\boldsymbol{q})\)</span> and
<span class="math notranslate nohighlight">\(\Gamma_{\mu}(\boldsymbol{q})\)</span> satisfy the relations given in Eqs.(82),(83)</p>
<figure class="align-default">
<a class="reference internal image-reference" href="../figures_04/Fig7.png"><img alt="../figures_04/Fig7.png" src="../figures_04/Fig7.png" style="width: 100%;" /></a>
</figure>
<p>Notice that the first one is a self-consistent equation. Instead of solving the self-consistent equation to evaluate <span class="math notranslate nohighlight">\(\Omega_{\mu}(\boldsymbol{q})\)</span> two approximated  approaches can be adpoted, both implemented in the SSCHA. One, that we call “one-shot”, evaluates the r.h.s of Eq.(82) at the SCHA frequency (Eqs.(84))</p>
<figure class="align-default">
<a class="reference internal image-reference" href="../figures_04/Fig8.png"><img alt="../figures_04/Fig8.png" src="../figures_04/Fig8.png" style="width: 100%;" /></a>
</figure>
<p>This approximation is  reasonable as long as the energy shift <span class="math notranslate nohighlight">\(\Delta_{\mu}(\boldsymbol{q})=\Omega_{\mu}(\boldsymbol{q})-\omega_{\mu}(\boldsymbol{q})\)</span> is small.
In particular, this is true if the SCHA self-energy is a (small) perturbation of the SCHA free propagator (not meaning that we are in a perturbative regime with respect to the harmonic approximation). In thas case, perturbation theory can be employed to evaluate the spectral function. If in Eq.(80) we keep only the first-order term in the self-energy, we get Eqs.(86),(87):</p>
<figure class="align-default">
<a class="reference internal image-reference" href="../figures_04/Fig9.png"><img alt="../figures_04/Fig9.png" src="../figures_04/Fig9.png" style="width: 100%;" /></a>
</figure>
<p>This concludes the overview on the quantities that we are going to compute for PbTe.</p>
</section>

<section id="calculations-on-pbte">
<h2>Calculations on PbTe<a class="headerlink" href="#calculations-on-pbte" title="Permalink to this headline"> </a></h2>
<p>We will perform calculations on PbTe in the rock-salt structure and, in order to speed-up the calculation, we will employ the force-field model already used in previous tutorials (the force-field model can be downloaded and installed from here <a class="reference external" href="https://github.com/SSCHAcode/F3ToyModel">https://github.com/SSCHAcode/F3ToyModel</a>). The calculations that we are going to perform are heavily underconverged and have to be consiered just as a guide to use of the SSCHA code. We will use a 2x2x2 supercell for PbTe. In order to define the force-field model
on this superell we need three FCs, <em>PbTe.ff.2x2x2.dyn1</em>, <em>PbTe.ff.2x2x2.dyn2</em>, and <em>PbTe.ff.2x2x2.dyn3</em></p>
<p>First, we need to do the SSCHA minimization. Create a directory <em>minim</em>, and go into it. We take the force-field FCs as starting point of the SSCHA mnimization too, with this input file <em>min.py</em></p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span>
<span class="kn">import</span> <span class="nn">cellconstructor.Phonons</span>
<span class="kn">import</span> <span class="nn">fforces</span> <span class="k">as</span> <span class="nn">ff</span>
<span class="kn">import</span> <span class="nn">fforces.Calculator</span>
<span class="kn">import</span> <span class="nn">sscha</span><span class="o">,</span> <span class="nn">sscha.Ensemble</span><span class="o">,</span> <span class="nn">sscha.SchaMinimizer</span>
<span class="kn">import</span> <span class="nn">sscha.Relax</span><span class="o">,</span> <span class="nn">sscha.Utilities</span>

<span class="c1"># ========================= TOY MODEL DEFINITION ===========================</span>
<span class="c1"># Dynamical matrices that set up the harmonic part of the force-field</span>
<span class="n">ff_dyn_name</span><span class="o">=</span><span class="s2">&quot;04_spectral_calculations/toy_matrices_2x2x2/PbTe.ff.2x2x2.dyn&quot;</span>
<span class="c1"># Paramters that set up the anharmonic part of the force-field</span>
<span class="n">p3</span> <span class="o">=</span> <span class="o">-</span><span class="mf">0.01408</span>
<span class="n">p4</span> <span class="o">=</span> <span class="o">-</span><span class="mf">0.01090</span>
<span class="n">p4x</span> <span class="o">=</span> <span class="mf">0.00254</span>
<span class="c1"># ====================================================================</span>

<span class="c1"># ==========================================================</span>
<span class="c1"># dynamical matrices to be used as starting guess</span>
<span class="n">dyn_sscha_name</span><span class="o">=</span><span class="s2">&quot;04_spectral_calculations/toy_matrices_2x2x2/PbTe.ff.2x2x2.dyn&quot;</span>
<span class="c1"># temperature</span>
<span class="n">T</span><span class="o">=</span><span class="mi">300</span>
<span class="c1"># minimization parameters</span>
<span class="n">N_CONFIGS</span> <span class="o">=</span> <span class="mi">50</span>
<span class="n">MAX_ITERATIONS</span> <span class="o">=</span> <span class="mi">10</span>
<span class="c1"># ====================================================================</span>
<span class="c1"># Setup the harmonic part of the force-field</span>
<span class="n">ff_dynmat</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="n">ff_dyn_name</span><span class="p">,</span> <span class="mi">3</span><span class="p">)</span>
<span class="n">ff_calculator</span> <span class="o">=</span> <span class="n">ff</span><span class="o">.</span><span class="n">Calculator</span><span class="o">.</span><span class="n">ToyModelCalculator</span><span class="p">(</span><span class="n">ff_dynmat</span><span class="p">)</span>
<span class="c1"># Setup the anharmonic part of the force-field</span>
<span class="n">ff_calculator</span><span class="o">.</span><span class="n">type_cal</span> <span class="o">=</span> <span class="s2">&quot;pbtex&quot;</span>
<span class="n">ff_calculator</span><span class="o">.</span><span class="n">p3</span> <span class="o">=</span> <span class="n">p3</span>
<span class="n">ff_calculator</span><span class="o">.</span><span class="n">p4</span> <span class="o">=</span> <span class="n">p4</span>
<span class="n">ff_calculator</span><span class="o">.</span><span class="n">p4x</span> <span class="o">=</span> <span class="n">p4x</span>
<span class="c1"># Load matrices</span>
<span class="n">dyn_sscha</span><span class="o">=</span><span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span> <span class="n">dyn_sscha_name</span><span class="p">,</span><span class="mi">3</span><span class="p">)</span>
<span class="n">dyn_sscha</span><span class="o">.</span><span class="n">Symmetrize</span><span class="p">()</span>
<span class="c1"># Generate the ensemble</span>
<span class="n">supercell</span><span class="o">=</span><span class="n">dyn_sscha</span><span class="o">.</span><span class="n">GetSupercell</span><span class="p">()</span>
<span class="n">ens</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Ensemble</span><span class="o">.</span><span class="n">Ensemble</span><span class="p">(</span><span class="n">dyn_sscha</span><span class="p">,</span> <span class="n">T</span><span class="p">,</span> <span class="n">supercell</span><span class="p">)</span>
<span class="n">ens</span><span class="o">.</span><span class="n">generate</span><span class="p">(</span><span class="n">N_CONFIGS</span><span class="p">)</span>
<span class="c1"># Compute energy and forces for the ensemble elements</span>
<span class="n">ens</span><span class="o">.</span><span class="n">get_energy_forces</span><span class="p">(</span><span class="n">ff_calculator</span> <span class="p">,</span> <span class="n">compute_stress</span> <span class="o">=</span> <span class="kc">False</span><span class="p">)</span>
<span class="c1"># Set up minimizer</span>
<span class="n">minimizer</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">SchaMinimizer</span><span class="o">.</span><span class="n">SSCHA_Minimizer</span><span class="p">(</span><span class="n">ens</span><span class="p">)</span>
<span class="c1"># Ignore the structure minimization (is fixed by symmetry)</span>
<span class="n">minimizer</span><span class="o">.</span><span class="n">minim_struct</span> <span class="o">=</span> <span class="kc">False</span>
<span class="c1"># max number steps (negative infinite)</span>
<span class="n">minimizer</span><span class="o">.</span><span class="n">max_ka</span><span class="o">=-</span><span class="mi">1</span>
<span class="c1"># Setup the minimization parameter for the covariance matrix</span>
<span class="n">minimizer</span><span class="o">.</span><span class="n">set_minimization_step</span><span class="p">(</span><span class="mf">1.0</span><span class="p">)</span>
<span class="c1"># Setup the threshold for the ensemble wasting</span>
<span class="n">minimizer</span><span class="o">.</span><span class="n">kong_liu_ratio</span> <span class="o">=</span><span class="mf">0.8</span>     <span class="c1"># Usually 0.5 is a good value</span>
<span class="n">minimizer</span><span class="o">.</span><span class="n">meaningful_factor</span><span class="o">=</span><span class="mf">1e-5</span>  <span class="c1"># meaningul factor</span>
<span class="c1"># Initialize the simulation</span>
<span class="n">relax</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Relax</span><span class="o">.</span><span class="n">SSCHA</span><span class="p">(</span><span class="n">minimizer</span><span class="p">,</span>
                          <span class="n">ff_calculator</span><span class="p">,</span>
                          <span class="n">N_configs</span> <span class="o">=</span> <span class="n">N_CONFIGS</span><span class="p">,</span>
                          <span class="n">max_pop</span> <span class="o">=</span> <span class="n">MAX_ITERATIONS</span><span class="p">)</span>
<span class="c1"># Define the I/O operations</span>
<span class="c1"># To save info about the free energy minimization after each step</span>
<span class="n">ioinfo</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Utilities</span><span class="o">.</span><span class="n">IOInfo</span><span class="p">()</span>
<span class="n">ioinfo</span><span class="o">.</span><span class="n">SetupSaving</span><span class="p">(</span><span class="s2">&quot;minim_info&quot;</span><span class="p">)</span>
<span class="n">relax</span><span class="o">.</span><span class="n">setup_custom_functions</span><span class="p">(</span><span class="n">custom_function_post</span> <span class="o">=</span> <span class="n">ioinfo</span><span class="o">.</span><span class="n">CFP_SaveAll</span><span class="p">)</span>
<span class="c1"># Start the minimization</span>
<span class="n">relax</span><span class="o">.</span><span class="n">relax</span><span class="p">()</span>
<span class="c1"># Print in stdout the info about the minimization</span>
<span class="c1"># and save the final dynamical matrix</span>
<span class="n">relax</span><span class="o">.</span><span class="n">minim</span><span class="o">.</span><span class="n">finalize</span><span class="p">()</span>
<span class="n">relax</span><span class="o">.</span><span class="n">minim</span><span class="o">.</span><span class="n">ensemble</span><span class="o">.</span><span class="n">save_bin</span><span class="p">(</span><span class="s2">&quot;./data_pop&quot;</span><span class="p">,</span> <span class="n">population_id</span><span class="o">=</span><span class="mi">1</span><span class="p">)</span>
<span class="n">relax</span><span class="o">.</span><span class="n">minim</span><span class="o">.</span><span class="n">dyn</span><span class="o">.</span><span class="n">save_qe</span><span class="p">(</span><span class="s2">&quot;SSCHA.T</span><span class="si">{}</span><span class="s2">.dyn&quot;</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="n">T</span><span class="p">))</span>
</pre></div>
</div>
<p>launching</p>
<div class="highlight-bash notranslate"><div class="highlight"><pre><span></span>$ python min.py &gt; min.out
</pre></div>
</div>
<p>after a few second the minimization has concluded. The three SSCHA FCs matrices have been saved as <em>SSCHA.T300.dyn#q</em>. Now we need to compute the third order FCs (FC3s) (and the Hessian FCs).
At the end of the SSCHA minimization, we saved the last population and the dynamical matrices that generated it too.
Using them, we compute the Hessian matrices and the FC3s. Exit from <em>minim</em>, create a directory <em>hessian</em>, go into it, and use this input file <cite>hessian.py</cite></p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span>
<span class="kn">import</span> <span class="nn">sscha.Ensemble</span>
<span class="c1">#</span>
<span class="n">NQIRR</span> <span class="o">=</span> <span class="mi">3</span>
<span class="n">Tg</span> <span class="o">=</span> <span class="mi">300</span>
<span class="n">T</span> <span class="o">=</span> <span class="mi">300</span>
<span class="c1">#</span>
<span class="n">pop_dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="s1">&#39;minim/data_pop/dyn_gen_pop1_&#39;</span><span class="p">,</span> <span class="mi">3</span><span class="p">)</span>
<span class="n">sscha_dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="s1">&#39;minim/SSCHA/SSCHA.T300.dyn&#39;</span><span class="p">,</span> <span class="mi">3</span><span class="p">)</span>
<span class="c1">#</span>
<span class="n">ens</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Ensemble</span><span class="o">.</span><span class="n">Ensemble</span><span class="p">(</span><span class="n">pop_dyn</span><span class="p">,</span> <span class="n">Tg</span><span class="p">)</span>
<span class="n">ens</span><span class="o">.</span><span class="n">load_bin</span><span class="p">(</span><span class="s1">&#39;minim/data_pop&#39;</span><span class="p">,</span> <span class="n">population_id</span> <span class="o">=</span> <span class="mi">1</span><span class="p">)</span>
<span class="n">ens</span><span class="o">.</span><span class="n">update_weights</span><span class="p">(</span><span class="n">sscha_dyn</span><span class="p">,</span> <span class="n">T</span><span class="p">)</span>
<span class="c1">#</span>
<span class="n">hessian_dyn</span><span class="p">,</span> <span class="n">d3</span> <span class="o">=</span> <span class="n">ens</span><span class="o">.</span><span class="n">get_free_energy_hessian</span><span class="p">(</span><span class="n">include_v4</span> <span class="o">=</span> <span class="kc">False</span><span class="p">,</span>
                                              <span class="n">return_d3</span> <span class="o">=</span> <span class="kc">True</span><span class="p">)</span>

<span class="n">hessian_dyn</span><span class="o">.</span><span class="n">save_qe</span><span class="p">(</span><span class="s1">&#39;Hessian.dyn&#39;</span><span class="p">)</span>

<span class="c1">############################### FC3 part #############################################</span>

<span class="n">tensor3</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">ForceTensor</span><span class="o">.</span><span class="n">Tensor3</span><span class="p">(</span><span class="n">dyn</span><span class="o">=</span><span class="n">sscha_dyn</span><span class="p">)</span>          <span class="c1"># initialize 3rd order tensor</span>
<span class="n">tensor3</span><span class="o">.</span><span class="n">SetupFromTensor</span><span class="p">(</span><span class="n">d3</span><span class="p">)</span>                              <span class="c1"># assign values</span>
<span class="n">tensor3</span><span class="o">.</span><span class="n">Center</span><span class="p">()</span>                                         <span class="c1"># center it</span>
<span class="n">tensor3</span><span class="o">.</span><span class="n">Apply_ASR</span><span class="p">()</span>                                      <span class="c1"># apply ASR</span>
<span class="n">tensor3</span><span class="o">.</span><span class="n">WriteOnFile</span><span class="p">(</span><span class="n">fname</span><span class="o">=</span><span class="s2">&quot;FC3&quot;</span><span class="p">,</span><span class="n">file_format</span><span class="o">=</span><span class="s1">&#39;D3Q&#39;</span><span class="p">)</span>       <span class="c1"># write on file</span>
</pre></div>
</div>
<p>giving</p>
<div class="highlight-bash notranslate"><div class="highlight"><pre><span></span>$ python hessian.py &gt; hessian.out
</pre></div>
</div>
<p>We already did an Hessian calculation in a previous tutorial. The new part is the creation of the 3rd order FCs (FC3s), which we wrote in the file <em>FC3</em>. Some comments about the input file. In <code class="docutils literal notranslate"><span class="pre">get_free_energy_hessian</span></code> we set <code class="docutils literal notranslate"><span class="pre">return_d3</span> <span class="pre">=</span> <span class="pre">True</span></code> because we need these informations to set up the FC3s. Moreover, since we are not going to use forth order FCs (we will work within the “bubble approximation”), we set <code class="docutils literal notranslate"><span class="pre">include_v4</span> <span class="pre">=</span> <span class="pre">False</span></code>, which in general saves a lot of computation time. Before writing the FC3s on file, we center it (a step necessary to perform Fourier interpolation), and apply the acoustic sum rule (ASR), since the centering spoils it.s</p>
<p>The format chosen here to write the FC3 file is the same used in the d3q.x code <a class="reference external" href="https://anharmonic.github.io/d3q/">https://anharmonic.github.io/d3q/</a>.
To be precise: exploiting the lattice translation symmetry, the third order FCs
can be written as <span class="math notranslate nohighlight">\(\Phi^{\alpha_1,\alpha_2,\alpha_3}_{a_1 a_2 a_3}(0,\boldsymbol{R},\boldsymbol{S})\)</span>, where
<span class="math notranslate nohighlight">\(\alpha_1,\alpha_2,\alpha_3\)</span> are cartesian indices, <span class="math notranslate nohighlight">\(a_1, a_2, a_3\)</span> atomic indices in the unit cell, and <span class="math notranslate nohighlight">\(\boldsymbol{R},\boldsymbol{S}\)</span> lattice vectors.
The FC3 file in D3Q format is</p>
<div class="highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">alpha_1</span> <span class="n">alpha_2</span> <span class="n">alpha_3</span> <span class="n">at_1</span> <span class="n">at_2</span> <span class="n">at_3</span>
<span class="n">N_RS</span>
<span class="n">R_x</span> <span class="n">R_y</span> <span class="n">R_z</span> <span class="n">S_x</span> <span class="n">S_y</span> <span class="n">S_z</span> <span class="n">phi</span><span class="p">(</span><span class="n">alpha_1</span><span class="p">,</span><span class="n">at_1</span><span class="p">,</span><span class="n">alpha_2</span><span class="p">,</span><span class="n">at_2</span><span class="p">,</span><span class="n">alpha_3</span><span class="p">,</span><span class="n">at_3</span><span class="p">)</span>

<span class="o">...</span>

<span class="n">alpha_1</span> <span class="n">alpha_2</span> <span class="n">alpha_3</span> <span class="n">at_1</span> <span class="n">at_2</span> <span class="n">at_3</span>
<span class="n">N_RS</span>
<span class="n">R_x</span> <span class="n">R_y</span> <span class="n">R_z</span> <span class="n">S_x</span> <span class="n">S_y</span> <span class="n">S_z</span> <span class="n">phi</span><span class="p">(</span><span class="n">alpha_1</span><span class="p">,</span><span class="n">at_1</span><span class="p">,</span><span class="n">alpha_2</span><span class="p">,</span><span class="n">at_2</span><span class="p">,</span><span class="n">alpha_3</span><span class="p">,</span><span class="n">at_3</span><span class="p">)</span>

<span class="o">...</span>
</pre></div>
</div>
<p>For each <code class="docutils literal notranslate"><span class="pre">alpha_1</span> <span class="pre">alpha_2</span> <span class="pre">alpha_3</span> <span class="pre">at_1</span> <span class="pre">at_2</span> <span class="pre">at_3</span></code> we have a block where: the first line is <code class="docutils literal notranslate"><span class="pre">N_RS</span></code>, which is
the number of <span class="math notranslate nohighlight">\(\boldsymbol{R},\boldsymbol{S}\)</span> considered. Each subsequent line refers to a couple <span class="math notranslate nohighlight">\(\boldsymbol{R},\boldsymbol{S}\)</span>, with <code class="docutils literal notranslate"><span class="pre">R_x</span> <span class="pre">R_y</span> <span class="pre">R_z</span></code> and <code class="docutils literal notranslate"><span class="pre">S_x</span> <span class="pre">S_y</span> <span class="pre">S_z</span></code> the crystal coordinates of <span class="math notranslate nohighlight">\(\boldsymbol{R}\)</span> and <span class="math notranslate nohighlight">\(\boldsymbol{S}\)</span>, respectively, and <code class="docutils literal notranslate"><span class="pre">phi(alpha_1,at_1,alpha_2,at_2,alpha_3,at_3)</span></code> the corresponding FCs value  <span class="math notranslate nohighlight">\(\Phi^{\alpha_1\alpha_2\alpha_3}_{a_1 a_2 a_3}(0,\boldsymbol{R},\boldsymbol{S})\)</span>.</p>
<p>Equipped with the third order SSCHA FCs ,written in real space in the <cite>FC3</cite> file, and the second-order SSCHA FCs, written in reciprocal space in the <em>SSCHA.T300.dyn#q</em> files, we have all the ingredient to compute the spectral functions. As first calculation, we compute the spectral function using Eq.(70) but within the “static approximation”, this meaning that we keep the selfenergy blocked with <span class="math notranslate nohighlight">\(\Omega=0\)</span>, as shown in Eq.(66) (within the bubble approximation)</p>
<figure class="align-default">
<a class="reference internal image-reference" href="../figures_04/Fig10.png"><img alt="../figures_04/Fig10.png" src="../figures_04/Fig10.png" style="width: 100%;" /></a>
</figure>
<p>In order to do that, exit from the current directory <em>hessian</em>, create a directory <em>spectral_static</em>, enter into it, and use this input file <em>spectral_static.py</em> to compute the spectral function in the static approximation, for the special point <span class="math notranslate nohighlight">\(X\)</span></p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span>
<span class="kn">import</span> <span class="nn">cellconstructor.ForceTensor</span>

<span class="n">dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="s2">&quot;minim/SSCHA.T300.dyn&quot;</span><span class="p">,</span><span class="mi">3</span><span class="p">)</span>
<span class="n">FC3</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">ForceTensor</span><span class="o">.</span><span class="n">Tensor3</span><span class="p">(</span><span class="n">dyn</span><span class="o">=</span><span class="n">dyn</span><span class="p">)</span>
<span class="n">FC3</span><span class="o">.</span><span class="n">SetupFromFile</span><span class="p">(</span><span class="n">fname</span><span class="o">=</span><span class="s2">&quot;hessian/FC3&quot;</span><span class="p">,</span><span class="n">file_format</span><span class="o">=</span><span class="s1">&#39;D3Q&#39;</span><span class="p">)</span>


<span class="c1"># integration grid</span>
<span class="n">k_grid</span><span class="o">=</span><span class="p">[</span><span class="mi">2</span><span class="p">,</span><span class="mi">2</span><span class="p">,</span><span class="mi">2</span><span class="p">]</span>

<span class="c1"># X in 2pi/Angstrom</span>
<span class="n">points</span><span class="o">=</span><span class="p">[</span><span class="mf">0.0</span><span class="p">,</span><span class="o">-</span><span class="mf">0.1547054</span><span class="p">,</span> <span class="mf">0.0</span><span class="p">]</span>

<span class="n">CC</span><span class="o">.</span><span class="n">Spectral</span><span class="o">.</span><span class="n">get_full_dynamic_correction_along_path</span><span class="p">(</span><span class="n">dyn</span><span class="o">=</span><span class="n">dyn</span><span class="p">,</span>
                                                <span class="n">tensor3</span><span class="o">=</span><span class="n">FC3</span><span class="p">,</span>
                                                <span class="n">k_grid</span><span class="o">=</span><span class="n">k_grid</span><span class="p">,</span>
                                                <span class="n">e1</span><span class="o">=</span><span class="mi">100</span><span class="p">,</span> <span class="n">de</span><span class="o">=</span><span class="mf">0.1</span><span class="p">,</span> <span class="n">e0</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span>     <span class="c1"># energy grid</span>
                                                <span class="n">T</span><span class="o">=</span><span class="mi">300</span><span class="p">,</span>
                                                <span class="n">q_path</span><span class="o">=</span><span class="n">points</span><span class="p">,</span>
                                                <span class="n">static_limit</span> <span class="o">=</span> <span class="kc">True</span><span class="p">,</span>
                                                <span class="n">filename_sp</span><span class="o">=</span><span class="s1">&#39;full_spectral_func_X&#39;</span><span class="p">)</span>
</pre></div>
</div>
<p>We used as integration <span class="math notranslate nohighlight">\(\boldsymbol{k}\)</span>-grid (i.e. the <span class="math notranslate nohighlight">\(\boldsymbol{k}\)</span>-grid of the summation in Eqs.(66), (75)) the grid commensurate with the supercell, i.e. a 2x2x2 grid. In that case, using the centering is irrelevant, as there is no Fourier interpolation.
With</p>
<div class="highlight-bash notranslate"><div class="highlight"><pre><span></span>$ mpirun -np <span class="m">4</span> python spectral_static.py &gt; spectral_static.out
</pre></div>
</div>
<p>we run the code with MPI with 4 parallel processes. In output we have the file  <em>full_spectral_func_X_static.dat</em> that contains the spectral function</p>
<div class="highlight-default notranslate"><div class="highlight"><pre><span></span><span class="c1"># -------------------------------------------------------------</span>
<span class="c1"># len (2pi/Angstrom), energy (cm-1), spectral function (1/cm-1)</span>
<span class="c1"># -------------------------------------------------------------</span>
 <span class="mf">0.000000</span>        <span class="mf">0.0000000</span>       <span class="mf">0.0000000</span>
 <span class="mf">0.000000</span>        <span class="mf">0.1000000</span>       <span class="mf">0.0000000</span>
 <span class="mf">0.000000</span>        <span class="mf">0.2000000</span>       <span class="mf">0.0000000</span>
 <span class="mf">0.000000</span>        <span class="mf">0.3000000</span>       <span class="mf">0.0000001</span>
 <span class="mf">0.000000</span>        <span class="mf">0.4000000</span>       <span class="mf">0.0000002</span>
 <span class="mf">0.000000</span>        <span class="mf">0.5000000</span>       <span class="mf">0.0000003</span>

                   <span class="o">...</span>
</pre></div>
</div>
<p>where the first line indicates the lenght of the path. This would be relevant in case we had not just a single <span class="math notranslate nohighlight">\(\boldsymbol{q}\)</span> point, but a path of <span class="math notranslate nohighlight">\(\boldsymbol{q}\)</span>-points. In that case, we would have several blocks, one for each <span class="math notranslate nohighlight">\(\boldsymbol{q}\)</span> point, and the first column of each block would indicate the lenght of the <span class="math notranslate nohighlight">\(\boldsymbol{q}\)</span>-path. In this case, since we have just one point, we have just one block with the first column equal to zero. Plotting the 3rd vs 2nd column we obtain this result:</p>
<figure class="align-default">
<a class="reference internal image-reference" href="../figures_04/static_X.png"><img alt="../figures_04/static_X.png" src="../figures_04/static_X.png" style="width: 100%;" /></a>
</figure>
<p>We have Dirac deltas ( to be precise, extremely narrow Lorentzians whose width is given only by the choice of the finite size of the used energy grid) around values that coincides with the Hessian frequency values (plotted here with vertical lines), that you can find in the <code class="docutils literal notranslate"><span class="pre">Hessian.dyn3</span></code> file obtained in the previous run. Indeed, the Hessian calculation corresponds exactly to a calculation done with the static self-energy. Two observations. The height of the spikes is proportional to the degeneracy of the modes. The yellow line indicates the integral function <span class="math notranslate nohighlight">\(\int_0^{\Omega}\sigma(\Omega',\boldsymbol{q})\,d\Omega'\)</span>, which at the end returns the value: [number of modes]/2 (therefore 3 in this case). This is a general sum rule fulfilled by the spectral function (not only in the static approximation).</p>
<p>Now we do a full calculation (no static approximation anymore). In this case, we need to specify the smearing parameter
<span class="math notranslate nohighlight">\(\delta_{\scriptscriptstyle{\text{se}}}\)</span> to compute the dynamic selfenergy from Eq.(75). In order to be tidy, let us do this calculation in another directory <em>spectral</em> (and let us do the same for all the subsequent calculations, new calculations in new directories). Using this <em>spectral.py</em>
input file</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span>
<span class="kn">import</span> <span class="nn">cellconstructor.ForceTensor</span>

<span class="n">dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="s2">&quot;minim/SSCHA.T300.dyn&quot;</span><span class="p">,</span><span class="mi">3</span><span class="p">)</span>
<span class="n">FC3</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">ForceTensor</span><span class="o">.</span><span class="n">Tensor3</span><span class="p">(</span><span class="n">dyn</span><span class="o">=</span><span class="n">dyn</span><span class="p">)</span>
<span class="n">FC3</span><span class="o">.</span><span class="n">SetupFromFile</span><span class="p">(</span><span class="n">fname</span><span class="o">=</span><span class="s2">&quot;hessian/FC3&quot;</span><span class="p">,</span><span class="n">file_format</span><span class="o">=</span><span class="s1">&#39;D3Q&#39;</span><span class="p">)</span>


<span class="c1"># integration grid</span>
<span class="n">k_grid</span><span class="o">=</span><span class="p">[</span><span class="mi">20</span><span class="p">,</span><span class="mi">20</span><span class="p">,</span><span class="mi">20</span><span class="p">]</span>

<span class="c1"># X in 2pi/Angstrom</span>
<span class="n">points</span><span class="o">=</span><span class="p">[</span><span class="mf">0.0</span><span class="p">,</span><span class="o">-</span><span class="mf">0.1547054</span><span class="p">,</span> <span class="mf">0.0</span><span class="p">]</span>

<span class="n">CC</span><span class="o">.</span><span class="n">Spectral</span><span class="o">.</span><span class="n">get_full_dynamic_correction_along_path</span><span class="p">(</span><span class="n">dyn</span><span class="o">=</span><span class="n">dyn</span><span class="p">,</span>
                                           <span class="n">tensor3</span><span class="o">=</span><span class="n">FC3</span><span class="p">,</span>
                                           <span class="n">k_grid</span><span class="o">=</span><span class="n">k_grid</span><span class="p">,</span>
                                           <span class="n">e1</span><span class="o">=</span><span class="mi">100</span><span class="p">,</span> <span class="n">de</span><span class="o">=</span><span class="mf">0.1</span><span class="p">,</span> <span class="n">e0</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span>     <span class="c1"># energy grid</span>
                                           <span class="n">sm1</span><span class="o">=</span><span class="mf">10.0</span><span class="p">,</span><span class="n">sm0</span><span class="o">=</span><span class="mf">1.0</span><span class="p">,</span><span class="n">nsm</span><span class="o">=</span><span class="mi">3</span><span class="p">,</span>    <span class="c1"># smearing values</span>
                                           <span class="n">T</span><span class="o">=</span><span class="mi">300</span><span class="p">,</span>
                                           <span class="n">q_path</span><span class="o">=</span><span class="n">points</span><span class="p">,</span>
                                           <span class="n">static_limit</span> <span class="o">=</span> <span class="kc">True</span><span class="p">,</span>
                                           <span class="n">filename_sp</span><span class="o">=</span><span class="s1">&#39;full_spectral_func_X&#39;</span><span class="p">)</span>
</pre></div>
</div>
<p>where now we have specified that we want to do the calculation with 3 smearing values (equally spaced) between 1.0 and 10.0
cm-1 (thus we will have 1.0, 5.5, and 10.0 cm-1) With</p>
<div class="highlight-bash notranslate"><div class="highlight"><pre><span></span>$ mpirun -np <span class="m">4</span> python spectral.py &gt; spectral.out
</pre></div>
</div>
<p>in output we have three files with the spectral functions, one for each smearing value. In general, convergence must be studied with respect to the integration <span class="math notranslate nohighlight">\(\boldsymbol{k}\)</span>-grid  and smearing used. Plotting the static result and the dynamic result for sm=1.0 cm-1, both computed with 20x20x20 <span class="math notranslate nohighlight">\(\boldsymbol{k}\)</span>-grid,  we see this result</p>
<figure class="align-default">
<a class="reference internal image-reference" href="../figures_04/dynamic_X.png"><img alt="../figures_04/dynamic_X.png" src="../figures_04/dynamic_X.png" style="width: 100%" /></a>
</figure>
<p>Therefore, we can conclude that in <span class="math notranslate nohighlight">\(X\)</span> the SSCHA phonons are barely affected by the interaction. However, the situation is different if analogous calculation is done in <span class="math notranslate nohighlight">\(\Gamma\)</span>.</p>
<div class="topic">
<p class="topic-title">Exercise</p>
<p>In <span class="math notranslate nohighlight">\(\Gamma\)</span>, do the same calculations previously done in <span class="math notranslate nohighlight">\(X\)</span>, and plot the results</p>
</div>
<p>This is the result you should obtain</p>
<figure class="align-default">
<a class="reference internal image-reference" href="../figures_04/G_spectrum.png"><img alt="../figures_04/G_spectrum.png" src="../figures_04/G_spectrum.png" style="width: 100%;" /></a>
</figure>
<p>Notice that here the triple-degenerate optical mode of the Hessian dynamical matrix is splitted into due different peaks of the static spectral function (LO and double degenerate TO). This is due to the LO-TO splitting occurring in PbTe. The frequencies in the Hessian dynamical matrix in <span class="math notranslate nohighlight">\(\Gamma\)</span> refer only to the short-range part of the FCs. However, the long-range dipole-dipole contribution coming from the Effective Charges (nonanalytic contribution), which is at the origin of the LO-TO splitting, is taken into account when the spectral function is computed. Moreover, notice that when dynamic spectral function is considered, the double-degenerate TO mode gets smeared, showing a strong non-Lorentzian character. When 4x4x4 supercell calculations are performed, it clearly appears a satellite peak.</p>
<p>Before continuing the spectral analysis, let us spend some time to investigate the static correction. As said, the static spectral
function is nothing but a collection of Dirac-deltas centered around the Hessian eigenvalues. Therefore, in the static case the only information
are the eigenvalues, there is not a complex spectrum to be analyzed. Indeed, there is a routine that we can use to compute the static correction for any q point, and for any integration grid. With this input file <em>static.py</em></p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">cellconstructor.ForceTensor</span>
<span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span>
<span class="kn">import</span> <span class="nn">numpy</span> <span class="k">as</span> <span class="nn">np</span> <span class="c1"># will be used just to create a path</span>

<span class="n">dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="s2">&quot;minim/SSCHA.T300.dyn&quot;</span><span class="p">,</span><span class="mi">3</span><span class="p">)</span> <span class="c1"># SSCHA matrices</span>
<span class="n">FC3</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">ForceTensor</span><span class="o">.</span><span class="n">Tensor3</span><span class="p">(</span><span class="n">dyn</span><span class="o">=</span><span class="n">dyn</span><span class="p">)</span>
<span class="n">FC3</span><span class="o">.</span><span class="n">SetupFromFile</span><span class="p">(</span><span class="n">fname</span><span class="o">=</span><span class="s2">&quot;hessian/FC3&quot;</span><span class="p">,</span><span class="n">file_format</span><span class="o">=</span><span class="s1">&#39;D3Q&#39;</span><span class="p">)</span>


<span class="c1"># integration grid</span>
<span class="n">k_grid</span><span class="o">=</span><span class="p">[</span><span class="mi">4</span><span class="p">,</span><span class="mi">4</span><span class="p">,</span><span class="mi">4</span><span class="p">]</span>

<span class="n">Xcoord</span><span class="o">=</span><span class="mf">0.1547054</span>
<span class="n">points</span><span class="o">=</span><span class="p">[[</span><span class="mf">0.0</span><span class="p">,</span><span class="n">z</span><span class="p">,</span><span class="mf">0.0</span><span class="p">]</span> <span class="k">for</span> <span class="n">z</span> <span class="ow">in</span> <span class="n">np</span><span class="o">.</span><span class="n">linspace</span><span class="p">(</span><span class="o">-</span><span class="n">Xcoord</span><span class="p">,</span><span class="n">Xcoord</span><span class="p">,</span><span class="mi">100</span><span class="p">)]</span> <span class="c1"># create the path</span>
                                   <span class="c1"># you can also download the path from a file</span>

<span class="n">CC</span><span class="o">.</span><span class="n">Spectral</span><span class="o">.</span><span class="n">get_static_correction_along_path</span><span class="p">(</span><span class="n">dyn</span><span class="o">=</span><span class="n">dyn</span><span class="p">,</span>
                                            <span class="n">tensor3</span><span class="o">=</span><span class="n">FC3</span><span class="p">,</span>
                                            <span class="n">k_grid</span><span class="o">=</span><span class="n">k_grid</span><span class="p">,</span>
                                            <span class="n">T</span><span class="o">=</span><span class="mi">300</span><span class="p">,</span>
                                            <span class="n">q_path</span><span class="o">=</span><span class="n">points</span><span class="p">)</span>
</pre></div>
</div>
<p>with</p>
<div class="highlight-bash notranslate"><div class="highlight"><pre><span></span>$ mpirun -np <span class="m">4</span> python static.py &gt; static.out
</pre></div>
</div>
<p>we obtain the file <em>v2+d3static_freq.dat</em>, done like this</p>
<div class="highlight-default notranslate"><div class="highlight"><pre><span></span><span class="c1"># ------------------------------------------------------------------------</span>
<span class="c1"># len (2pi/Angstrom), sscha freq (cm-1), sscha + static bubble freq (cm-1)</span>
<span class="c1"># ------------------------------------------------------------------------</span>
 <span class="mf">0.000000</span>    <span class="mf">22.6699858</span>  <span class="mf">22.6699858</span> <span class="o">...</span>        <span class="mf">22.1602770</span>  <span class="mf">22.1717260</span> <span class="o">...</span>
 <span class="mf">0.003125</span>    <span class="mf">22.7847829</span>  <span class="mf">22.7847829</span> <span class="o">...</span>        <span class="mf">29.9763683</span>  <span class="mf">80.8414572</span> <span class="o">...</span>

                              <span class="o">...</span>
</pre></div>
</div>
<p>where the first column is the lenght of the path in 2:math:<cite>pi</cite>/Angstrom, the next (6, in this case) columns are the
SSCHA frequencies, and the next (6, in this case) columns are the SSCHA+static bubble self-energy-corrected frequencies. This is the plot obtained with this result</p>
<figure class="align-default">
<a class="reference internal image-reference" href="../figures_04/static_path.png"><img alt="../figures_04/static_path.png" src="../figures_04/static_path.png" style="width: 100%;" /></a>
</figure>
<p>Therefore, as long as one is interested only in the static correction, e.g. because one wants to study the structural instability,
the routine <code class="docutils literal notranslate"><span class="pre">get_static_correction_along_path</span></code> is the one that has to be employed. Indeed, notice that this is the proper way to
detect instabilities (imaginary frequencies) in points of the Brillouin zone that do not belong to the grid used to compute the Hessian.
One should not Fourier interpolate the Hessian matrices computed on a grid, in order to obtain the frequency dispersion along a path, but rather should interpolate the correction
and add it to the SSCHA frequency, point by point (which is what we are doing here). Moreover, in this way we can increase the integration <span class="math notranslate nohighlight">\(k\)</span>-grid to reach the convergence. Notice that the LO-TO splitting has been properly taken into account.</p>
<p>We can now go back to the spectral calculations. In general, calculations done with Eq.(70) can be heavy, but often the off-diagonal terms of the phonon self-energy
in the mode basis set can be neglected and use Eqs. (78), (79), (80). This is the case of PbTe. With this input file <em>nomm_spectral.py</em></p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span>
<span class="kn">import</span> <span class="nn">cellconstructor.ForceTensor</span>

<span class="n">dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="s2">&quot;minim/SSCHA.T300.dyn&quot;</span><span class="p">,</span><span class="mi">3</span><span class="p">)</span>
<span class="n">FC3</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">ForceTensor</span><span class="o">.</span><span class="n">Tensor3</span><span class="p">(</span><span class="n">dyn</span><span class="o">=</span><span class="n">dyn</span><span class="p">)</span>
<span class="n">FC3</span><span class="o">.</span><span class="n">SetupFromFile</span><span class="p">(</span><span class="n">fname</span><span class="o">=</span><span class="s2">&quot;hessian/FC3&quot;</span><span class="p">,</span><span class="n">file_format</span><span class="o">=</span><span class="s1">&#39;D3Q&#39;</span><span class="p">)</span>


<span class="c1"># integration grid</span>
<span class="n">k_grid</span><span class="o">=</span><span class="p">[</span><span class="mi">20</span><span class="p">,</span><span class="mi">20</span><span class="p">,</span><span class="mi">20</span><span class="p">]</span>

<span class="n">G</span><span class="o">=</span><span class="p">[</span><span class="mf">0.0</span><span class="p">,</span><span class="mf">0.0</span><span class="p">,</span><span class="mf">0.0</span><span class="p">]</span>
<span class="n">points</span><span class="o">=</span><span class="n">G</span>

<span class="n">CC</span><span class="o">.</span><span class="n">Spectral</span><span class="o">.</span><span class="n">get_diag_dynamic_correction_along_path</span><span class="p">(</span><span class="n">dyn</span><span class="p">,</span>
                                           <span class="n">tensor3</span><span class="o">=</span><span class="n">FC3</span><span class="p">,</span>
                                           <span class="n">k_grid</span><span class="o">=</span><span class="n">k_grid</span><span class="p">,</span>
                                           <span class="n">e1</span><span class="o">=</span><span class="mi">150</span><span class="p">,</span> <span class="n">de</span><span class="o">=</span><span class="mf">0.1</span><span class="p">,</span> <span class="n">e0</span><span class="o">=</span><span class="mf">0.0</span><span class="p">,</span>
                                           <span class="n">sm1</span><span class="o">=</span><span class="mf">1.0</span><span class="p">,</span> <span class="n">sm0</span><span class="o">=</span><span class="mf">1.0</span><span class="p">,</span>
                                           <span class="n">nsm</span><span class="o">=</span><span class="mi">1</span><span class="p">,</span>
                                           <span class="n">q_path</span><span class="o">=</span><span class="n">points</span><span class="p">,</span>
                                           <span class="n">T</span><span class="o">=</span><span class="mf">300.0</span><span class="p">)</span>
</pre></div>
</div>
<p>and</p>
<div class="highlight-bash notranslate"><div class="highlight"><pre><span></span>$ mpirun -np <span class="m">4</span> python nomm_spectral.py &gt; nomm_spectral.out
</pre></div>
</div>
<p>we obtain several files:</p>
<ul>
<li><p><em>spectral_func_1.00.dat</em></p>
<blockquote>
<div><p>with the structure</p>
<div class="highlight-default notranslate"><div class="highlight"><pre><span></span><span class="c1"># ---------------------------------------------------------------------------------</span>
<span class="c1"># len (2pi/Ang), ene. (cm-1), spec. func. (1/cm-1), spec. fun. mode comp. (1/cm-1)</span>
<span class="c1"># ---------------------------------------------------------------------------------</span>
  <span class="mf">0.000</span>         <span class="mf">0.000</span>         <span class="mf">0.000</span>           <span class="o">-</span><span class="mf">0.000</span>          <span class="mf">0.000</span>   <span class="o">...</span>
  <span class="mf">0.000</span>         <span class="mf">0.100</span>         <span class="mf">1.527</span>            <span class="mf">0.509</span>          <span class="mf">0.509</span>   <span class="o">...</span>
  <span class="mf">0.000</span>         <span class="mf">0.200</span>         <span class="mf">2.387</span>            <span class="mf">0.795</span>          <span class="mf">0.795</span>   <span class="o">...</span>
                           <span class="o">...</span>
</pre></div>
</div>
<p>The file is made by several blocks, one for each <span class="math notranslate nohighlight">\(\boldsymbol{q}\)</span> point of the path considered
(now we have just one     point, thus one block). The first column of the block
gives the length of the path (at that point). The second and third column give the energy and total
spectrum,    respectively. The subsequent columns give the contribution to the spectrum given by each mode
(Cfr. Eqs.(78),(79)). Plotting the 3rd vs 2nd column you can verify that in this case the spectrum is essentially
identical to the one already obtained with the full formula (i.e. to the spectrum obtained considering the
off-diagonal terms of the self-energy too). Notice that in this file we have the spectrum given by the
three acoustic modes in <span class="math notranslate nohighlight">\(\Gamma\)</span>. We did not consider the translational modes (because trivial)
in the calculation done with <em>get_full_dynamic_correction_along_path</em>, where the flag <code class="docutils literal notranslate"><span class="pre">notransl</span></code>
by default was set equal to <code class="docutils literal notranslate"><span class="pre">True</span></code>.</p>
</div></blockquote>
</li>
<li><p><em>spectral_func_lorentz_one_shot_1.00.dat</em>, <em>spectral_func_lorentz_perturb_1.00.dat</em></p>
<blockquote>
<div><p>These files have the same structure of <em>spectral_func_1.00.dat</em>. However, now the spectral functions are computed in the Lorenztian approximation, Eq.(81),
using the one-shot, Eqs.(84),(85) and the perturbative, Eqs.(86),(87),  values of the energy and HWHM. The codes offers also the possibility to use a still
very tentative approach to solve the self-consistent relation Eqs.(82),(83), and produce the relative Lorentzian spectral functions.
Plotting the spectral functions of the TO mode from <em>spectral_func_1.00.dat</em>  and <em>spectral_func_lorentz_one_shot_1.00.dat</em> we obtain this</p>
<figure class="align-default">
<a class="reference internal image-reference" href="../figures_04/NoMM_Spectral.png"><img alt="../figures_04/NoMM_Spectral.png" src="../figures_04/NoMM_Spectral.png" style="width: 100%;" /></a>
</figure>
<p>This confirms the strong non-Lorentzian character of this mode</p>
</div></blockquote>
</li>
<li><p><em>v2_freq_shift_hwhm_one_shot_1.00.dat, v2_freq_shift_hwhm_perturb_1.00.dat</em></p>
<blockquote>
<div><p>They have the structure</p>
<div class="highlight-default notranslate"><div class="highlight"><pre><span></span><span class="c1"># -----------------------------------------------------------------</span>
<span class="c1"># len (2pi/Angstrom), SSCHA freq (cm-1), shift (cm-1) , HWHM (cm-1)</span>
<span class="c1"># -----------------------------------------------------------------</span>

                      <span class="o">....</span>
</pre></div>
</div>
<p>For each <span class="math notranslate nohighlight">\(\boldsymbol{q}\)</span> point there is a line, with the lenght along the path, the SSCHA frequencies for that point, <span class="math notranslate nohighlight">\(\omega_{\mu}(\boldsymbol{q})\)</span>, the energy shift <span class="math notranslate nohighlight">\(\Delta_\mu(\boldsymbol{q})=\Omega_\mu(\boldsymbol{q})-\omega_\mu(\boldsymbol{q})\)</span> and the HWHM <span class="math notranslate nohighlight">\(\Gamma_\mu(\boldsymbol{q})\)</span>.</p>
</div></blockquote>
</li>
<li><p><em>freq_dynamic_one_shot_1.00.dat, freq_dynamic_perturb_1.00.dat</em></p>
<blockquote>
<div><p>They have the structure</p>
<div class="highlight-default notranslate"><div class="highlight"><pre><span></span><span class="c1"># ------------------------------------------------------------</span>
<span class="c1"># len (2pi/Angstrom), SSCHA+shift (sorted) (cm-1), HWHM (cm-1)</span>
<span class="c1"># ------------------------------------------------------------</span>

                      <span class="o">....</span>
</pre></div>
</div>
<p>For each <span class="math notranslate nohighlight">\(\boldsymbol{q}\)</span> point there is a line, with the lenght along the path, the shifted frequencies <span class="math notranslate nohighlight">\(\Omega_\mu(\boldsymbol{q})=\omega_\mu(\boldsymbol{q})+\Delta_\mu(\boldsymbol{q})\)</span> (sorted in incresing value, in principle different from the SSCHA-frequency increasing order) and the corresponding HWHMs <span class="math notranslate nohighlight">\(\Gamma_\mu(\boldsymbol{q})\)</span>.  This is the file that has to be used to plot the correct phonon dispersion, together with the linewidth.</p>
</div></blockquote>
</li>
</ul>
<p>There is a dedicated routine to compute, with less computational time and more accuracy, only the energy shift (i.e. the corrected frequency) and the linewidth of the modes in the one-shot and perturbative no-mode-mixing Lorentzian approach. It is the <code class="docutils literal notranslate"><span class="pre">get_os_perturb_dynamic_correction_along_path</span></code> routine, using the
<em>os_perturb_correction.py</em> input file</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span>
<span class="kn">import</span> <span class="nn">cellconstructor.ForceTensor</span>

<span class="n">dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="s2">&quot;minim/SSCHA.T300.dyn&quot;</span><span class="p">,</span><span class="mi">3</span><span class="p">)</span>
<span class="n">FC3</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">ForceTensor</span><span class="o">.</span><span class="n">Tensor3</span><span class="p">(</span><span class="n">dyn</span><span class="o">=</span><span class="n">dyn</span><span class="p">)</span>
<span class="n">FC3</span><span class="o">.</span><span class="n">SetupFromFile</span><span class="p">(</span><span class="n">fname</span><span class="o">=</span><span class="s2">&quot;hessian/FC3&quot;</span><span class="p">,</span><span class="n">file_format</span><span class="o">=</span><span class="s1">&#39;D3Q&#39;</span><span class="p">)</span>


<span class="c1"># integration grid</span>
<span class="n">k_grid</span><span class="o">=</span><span class="p">[</span><span class="mi">10</span><span class="p">,</span><span class="mi">10</span><span class="p">,</span><span class="mi">10</span><span class="p">]</span>

<span class="n">points</span><span class="o">=</span><span class="p">[</span><span class="mf">0.0</span><span class="p">,</span><span class="mf">0.0</span><span class="p">,</span><span class="mf">0.0</span><span class="p">]</span>

<span class="n">CC</span><span class="o">.</span><span class="n">Spectral</span><span class="o">.</span><span class="n">get_os_perturb_dynamic_correction_along_path</span><span class="p">(</span><span class="n">dyn</span><span class="p">,</span>
                                           <span class="n">tensor3</span><span class="o">=</span><span class="n">FC3</span><span class="p">,</span>
                                           <span class="n">k_grid</span><span class="o">=</span><span class="n">k_grid</span><span class="p">,</span>
                                           <span class="n">sm1</span><span class="o">=</span><span class="mf">1.0</span><span class="p">,</span> <span class="n">sm0</span><span class="o">=</span><span class="mf">1.0</span><span class="p">,</span>
                                           <span class="n">nsm</span><span class="o">=</span><span class="mi">1</span><span class="p">,</span>
                                           <span class="n">q_path</span><span class="o">=</span><span class="n">points</span><span class="p">,</span>
                                           <span class="n">T</span><span class="o">=</span><span class="mf">300.0</span><span class="p">)</span>
</pre></div>
</div>
<div class="topic">
<p class="topic-title">Exercise</p>
<p>Do this calculation to obtain the dispersion along the path <span class="math notranslate nohighlight">\(X-\Gamma-X\)</span>.</p>
</div>
<p>After that, using this script to extract the results</p>
<div class="highlight-bash notranslate"><div class="highlight"><pre><span></span><span class="k">for</span> i <span class="k">in</span> <span class="sb">`</span>seq <span class="m">1</span> <span class="m">6</span><span class="sb">`</span>
<span class="k">do</span>
 awk -v <span class="nv">i</span><span class="o">=</span><span class="s2">&quot;</span><span class="si">${</span><span class="nv">i</span><span class="si">}</span><span class="s2">&quot;</span> <span class="s1">&#39;{ if ( NR  &gt; 3 ) printf&quot;%22.11f%22.11f%22.11f%22.11f\n&quot;, \</span>
<span class="s1"> $1,$(i+1),$(i+1)+$(i+7),$(i+1)-$(i+7)}&#39;</span> freq_dynamic_1.0.os.dat &gt; freq_<span class="si">${</span><span class="nv">i</span><span class="si">}</span>
<span class="k">done</span>
</pre></div>
</div>
<p>and then using this gnuplot script <em>plot.gp</em></p>
<div class="highlight-gnuplot notranslate"><div class="highlight"><pre><span></span><span class="k">set</span> <span class="nb">terminal</span> <span class="n">pngcairo</span> <span class="n">size</span> <span class="mi">2048</span><span class="o">,</span><span class="mi">1536</span> <span class="n">enhanced</span>
<span class="k">set</span> <span class="nb">output</span> <span class="s">&#39;freq.png&#39;</span>

<span class="nv">xmin</span><span class="o">=</span><span class="mi">0</span>
<span class="nv">xmax</span><span class="o">=</span><span class="mf">0.30941100000</span>
<span class="k">set</span> <span class="nb">xrange</span><span class="p">[</span><span class="n">xmin</span><span class="o">:</span><span class="n">xmax</span><span class="p">]</span>

<span class="k">set</span> <span class="nb">multiplot</span>
<span class="k">set</span> <span class="nb">object</span> <span class="n">rectangle</span> <span class="n">from</span> <span class="n">graph</span> <span class="mi">0</span><span class="o">,</span><span class="mi">0</span> <span class="n">to</span> <span class="n">graph</span> <span class="mi">1</span><span class="o">,</span><span class="mi">1</span> \
<span class="n">behind</span> <span class="n">fillcolor</span> <span class="n">rgb</span> <span class="s">&#39;black&#39;</span> <span class="n">fillstyle</span> <span class="n">transparent</span> <span class="n">solid</span> <span class="mf">0.100</span> <span class="n">noborder</span>

<span class="k">plot</span> <span class="n">for</span> <span class="p">[</span><span class="nb">i</span><span class="o">=</span><span class="mi">1</span><span class="o">:</span><span class="mi">6</span><span class="p">]</span> <span class="s">&#39;freq_&#39;</span><span class="o">.</span><span class="nb">i</span> <span class="nb">using</span> <span class="mi">1</span><span class="o">:</span><span class="mi">3</span><span class="o">:</span><span class="mi">4</span> <span class="nb">with</span> <span class="n">filledcurves</span> \
    <span class="n">fs</span> <span class="n">transparent</span> <span class="n">solid</span> <span class="mf">0.125</span> <span class="n">noborder</span> <span class="n">lc</span> <span class="n">rgb</span> <span class="s">&quot;blue&quot;</span> <span class="nb">notitle</span>
<span class="k">plot</span> <span class="n">for</span> <span class="p">[</span><span class="nb">i</span><span class="o">=</span><span class="mi">1</span><span class="o">:</span><span class="mi">6</span><span class="p">]</span> <span class="s">&#39;freq_&#39;</span><span class="o">.</span><span class="nb">i</span> <span class="nb">using</span> <span class="mi">1</span><span class="o">:</span><span class="mi">2</span>   <span class="nb">with</span> <span class="n">lines</span> \
    <span class="n">lc</span> <span class="n">rgb</span> <span class="s">&quot;black&quot;</span> <span class="n">lw</span> <span class="mf">3.5</span> <span class="nb">notitle</span>
</pre></div>
</div>
<div class="highlight-bash notranslate"><div class="highlight"><pre><span></span>$ gnuplot plot.gp
</pre></div>
</div>
<p>you obtain this figure</p>
<blockquote>
<div><figure class="align-default">
<a class="reference internal image-reference" href="../figures_04/freq.png"><img alt="../figures_04/freq.png" src="../figures_04/freq.png" style="width: 100%;" /></a>
</figure>
</div></blockquote>
<p>Here you have the plot of the shifted SSCHA phonon frequencies with the linewidth. However, it must be remembered that this picture is appropriate as long
as the Lorentzian picture is valid. We already know that at least in <span class="math notranslate nohighlight">\(\Gamma\)</span> this is not really the case (as said, this is even more evident if the calculation
is done with a 4x4x4 supercell). In that case, the best thing to do is a spectral calculation (full, or in the no-mode-mixing approximation), and use the
three columns (lenght of the path &amp; energy &amp; spectral value) to do a colorplot. For example, this is the kind of result that you obtain for PbTe with the 4x4x4 supercell</p>
<figure class="align-default">
<a class="reference internal image-reference" href="../figures_04/spec_path.png"><img alt="../figures_04/spec_path.png" src="../figures_04/spec_path.png" style="width: 100%;" /></a>
</figure>
<div class="topic">
<p class="topic-title">Exercise</p>
<p>Do the no-mode-mixing calculation along the path <span class="math notranslate nohighlight">\(X-\Gamma-X\)</span> with smearing <span class="math notranslate nohighlight">\(10.0\,\text{cm}^{-1}\)</span></p>
</div>
<p>We conclde this tutorial stressing that a convergence analysis in terms of integration grid and smearing has to be done in order to obtain reliable results. The <code class="docutils literal notranslate"><span class="pre">get_os_perturb_dynamic_correction_along_path</span></code> routine is the best tool to do that. With this <em>input.py</em>  input file you can do
the calculation for several integration grids</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span>
<span class="kn">import</span> <span class="nn">cellconstructor.ForceTensor</span>


<span class="n">dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="s2">&quot;minim/SSCHA.T300.dyn&quot;</span><span class="p">,</span><span class="mi">3</span><span class="p">)</span>
<span class="n">FC3</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">ForceTensor</span><span class="o">.</span><span class="n">Tensor3</span><span class="p">(</span><span class="n">dyn</span><span class="o">=</span><span class="n">dyn</span><span class="p">)</span>
<span class="n">FC3</span><span class="o">.</span><span class="n">SetupFromFile</span><span class="p">(</span><span class="n">fname</span><span class="o">=</span><span class="s2">&quot;hessian/FC3&quot;</span><span class="p">,</span><span class="n">file_format</span><span class="o">=</span><span class="s1">&#39;D3Q&#39;</span><span class="p">)</span>


<span class="k">for</span> <span class="n">kval</span> <span class="ow">in</span> <span class="p">[</span><span class="mi">4</span><span class="p">,</span><span class="mi">8</span><span class="p">,</span><span class="mi">16</span><span class="p">,</span><span class="mi">32</span><span class="p">]:</span>

   <span class="nb">print</span><span class="p">(</span><span class="s2">&quot;COMPUTING </span><span class="si">{}</span><span class="s2">&quot;</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="n">kval</span><span class="p">))</span>

   <span class="n">CC</span><span class="o">.</span><span class="n">Spectral</span><span class="o">.</span><span class="n">get_os_perturb_dynamic_correction_along_path</span><span class="p">(</span><span class="n">dyn</span><span class="p">,</span>
                                  <span class="n">tensor3</span><span class="o">=</span><span class="n">FC3</span><span class="p">,</span>
                                  <span class="n">k_grid</span><span class="o">=</span><span class="p">[</span><span class="n">kval</span><span class="p">,</span><span class="n">kval</span><span class="p">,</span><span class="n">kval</span><span class="p">],</span>
                                  <span class="n">sm1</span><span class="o">=</span><span class="mf">20.0</span><span class="p">,</span> <span class="n">sm0</span><span class="o">=</span><span class="mf">0.1</span><span class="p">,</span>
                                  <span class="n">nsm</span><span class="o">=</span><span class="mi">80</span><span class="p">,</span>
                                  <span class="n">q_path</span><span class="o">=</span><span class="p">[</span><span class="mf">0.0</span><span class="p">,</span><span class="mf">0.0</span><span class="p">,</span><span class="mf">0.0</span><span class="p">],</span>
                                  <span class="n">T</span><span class="o">=</span><span class="mf">300.0</span><span class="p">,</span>
                                  <span class="n">filename_shift_lw</span>  <span class="o">=</span> <span class="s1">&#39;shift_hwhm_</span><span class="si">{}</span><span class="s1">&#39;</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="n">kval</span><span class="p">),</span>
                                  <span class="n">filename_freq_dyn</span> <span class="o">=</span> <span class="s1">&#39;freq_</span><span class="si">{}</span><span class="s1">&#39;</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="n">kval</span><span class="p">))</span>
</pre></div>
</div>
<p>giving</p>
<div class="highlight-bash notranslate"><div class="highlight"><pre><span></span>$ mpirun -np <span class="m">4</span> python input.py &gt; output
</pre></div>
</div>
<p>From <em>output</em> you can take the list of smearing values and write them in a file <em>sm.dat</em>.
After that, you can collect the results with this <em>extract.sh</em> script</p>
<div class="highlight-bash notranslate"><div class="highlight"><pre><span></span><span class="k">for</span> grid <span class="k">in</span> <span class="m">4</span> <span class="m">8</span> <span class="m">16</span> <span class="m">32</span>
<span class="k">do</span>
 &gt; <span class="si">${</span><span class="nv">grid</span><span class="si">}</span>x<span class="si">${</span><span class="nv">grid</span><span class="si">}</span>x<span class="si">${</span><span class="nv">grid</span><span class="si">}</span>.dat
 <span class="k">while</span> <span class="nb">read</span> sm
 <span class="k">do</span>
  tail -1  shift_hwhm_<span class="si">${</span><span class="nv">grid</span><span class="si">}</span>_<span class="si">${</span><span class="nv">sm</span><span class="si">}</span>.os.dat <span class="se">\</span>
  <span class="p">|</span> awk -v <span class="nv">sm</span><span class="o">=</span><span class="s2">&quot;</span><span class="si">${</span><span class="nv">sm</span><span class="si">}</span><span class="s2">&quot;</span> <span class="s1">&#39;{printf&quot;%11.7f\t%22.11f\n&quot;,sm,$(NF-1)}&#39;</span> <span class="se">\</span>
  &gt;&gt; <span class="si">${</span><span class="nv">grid</span><span class="si">}</span>x<span class="si">${</span><span class="nv">grid</span><span class="si">}</span>x<span class="si">${</span><span class="nv">grid</span><span class="si">}</span>.dat
 <span class="k">done</span> &lt; sm.dat
<span class="k">done</span>
</pre></div>
</div>
<p>giving</p>
<div class="highlight-bash notranslate"><div class="highlight"><pre><span></span>$ bash extract.sh
</pre></div>
</div>
<p>so as to obtain the files <em>4x4x4.dat</em>, <em>8x8x8.dat</em>, <em>16x16x16.dat</em>, and <em>32x32x32.dat</em>. Plotting them you obtain</p>
<figure class="align-default">
<a class="reference internal image-reference" href="../figures_04/convergence.png"><img alt="../figures_04/convergence.png" src="../figures_04/convergence.png" style="width: 100%;" /></a>
</figure>
<p>From the plot, you can see that you need at least a 8x8x8 grid to obtain a converged value, using smearing equal to 1.0 cm-1.</p>
</section>

