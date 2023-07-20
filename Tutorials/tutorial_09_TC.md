---
layout: page
title: Thermal conductivity calculations with the SSCHA
---

This tutorial was prepared for the [2023 SSCHA School](http://sscha.eu/Schools/2023/home/) by Đorđe Dangić. You can see here the video os the hands-on session:

<iframe width="560" height="315" src="https://www.youtube.com/embed/FjhvCO_dWdY" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>

<p>In previous lessons we saw how to calculate vibrational properties of material using SSCHA. Now we will use this acquired knowledge to calculate lattice thermal conductivity of materials. We will need dynamical matrices (<strong>auxiliary ones, not hessians</strong>) and the third order force constants (we already calculated them when we checked the dynamical stability of the system). With these we can calculate materials’ harmonic (phonon frequencies and phonon group velocities) and anharmonic properties (phonon lifetimes and spectral functions) which is all we need to calculate lattice thermal conductivity.</p>
<section id="lattice-thermal-conductivity-of-silicon">
<h2>Lattice thermal conductivity of silicon<a class="headerlink" href="#lattice-thermal-conductivity-of-silicon" title="Permalink to this headline"> </a></h2>
<p>As a first exercise let’s calculate lattice thermal conductivity of silicon. Silicon is very harmonic material which means it’s lattice thermal conductivity is very high. This also makes it a good test case to check the equivalence of Green-Kubo and Boltzmann transport equation approaches in the limit of vanishing anharmonicity. To speed up the calculation we will use Tersoff potential to obtain the second and third order force constants. We will do this with this simple script:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">numpy</span> <span class="k">as</span> <span class="nn">np</span>
<span class="kn">from</span> <span class="nn">quippy.potential</span> <span class="kn">import</span> <span class="n">Potential</span>
<span class="kn">from</span> <span class="nn">ase</span> <span class="kn">import</span> <span class="n">Atoms</span>
<span class="kn">import</span> <span class="nn">ase.io</span>
<span class="kn">from</span> <span class="nn">ase.eos</span> <span class="kn">import</span> <span class="n">calculate_eos</span>
<span class="kn">from</span> <span class="nn">ase.units</span> <span class="kn">import</span> <span class="n">kJ</span>
<span class="kn">from</span> <span class="nn">ase.phonons</span> <span class="kn">import</span> <span class="n">Phonons</span> <span class="k">as</span> <span class="n">AsePhonons</span>
<span class="kn">from</span> <span class="nn">ase.constraints</span> <span class="kn">import</span> <span class="n">ExpCellFilter</span>
<span class="kn">from</span> <span class="nn">ase.optimize</span> <span class="kn">import</span> <span class="n">BFGS</span><span class="p">,</span> <span class="n">QuasiNewton</span>
<span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span>
<span class="kn">import</span> <span class="nn">cellconstructor.Phonons</span>
<span class="kn">import</span> <span class="nn">cellconstructor.Structure</span>
<span class="kn">import</span> <span class="nn">sscha</span><span class="o">,</span> <span class="nn">sscha.Ensemble</span><span class="o">,</span> <span class="nn">sscha.SchaMinimizer</span><span class="o">,</span> <span class="nn">sscha.Relax</span>


<span class="c1"># This function will use ASE to give us a starting</span>
<span class="c1"># guess for dynamical matrices</span>
<span class="k">def</span> <span class="nf">get_starting_dynamical_matrices</span><span class="p">(</span><span class="n">structure_filename</span><span class="p">,</span>
        <span class="n">potential</span><span class="p">,</span> <span class="n">supercell</span><span class="p">):</span>
    <span class="n">atoms</span> <span class="o">=</span> <span class="n">ase</span><span class="o">.</span><span class="n">io</span><span class="o">.</span><span class="n">read</span><span class="p">(</span><span class="n">structure_filename</span><span class="p">)</span>
    <span class="n">atoms</span><span class="o">.</span><span class="n">set_calculator</span><span class="p">(</span><span class="n">potential</span><span class="p">)</span>
    <span class="n">ecf</span> <span class="o">=</span> <span class="n">ExpCellFilter</span><span class="p">(</span><span class="n">atoms</span><span class="p">)</span>
    <span class="n">qn</span> <span class="o">=</span> <span class="n">QuasiNewton</span><span class="p">(</span><span class="n">ecf</span><span class="p">)</span>
    <span class="n">qn</span><span class="o">.</span><span class="n">run</span><span class="p">(</span><span class="n">fmax</span><span class="o">=</span><span class="mf">0.0005</span><span class="p">)</span>

    <span class="n">structure</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Structure</span><span class="o">.</span><span class="n">Structure</span><span class="p">()</span>
    <span class="n">structure</span><span class="o">.</span><span class="n">generate_from_ase_atoms</span><span class="p">(</span><span class="n">atoms</span><span class="p">,</span> <span class="n">get_masses</span> <span class="o">=</span> <span class="kc">True</span><span class="p">)</span>
    <span class="n">dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">compute_phonons_finite_displacements</span><span class="p">(</span><span class="n">structure</span><span class="p">,</span>
    <span class="n">potential</span><span class="p">,</span> <span class="n">supercell</span> <span class="o">=</span> <span class="n">supercell</span><span class="p">)</span>
    <span class="n">dyn</span><span class="o">.</span><span class="n">Symmetrize</span><span class="p">()</span>
    <span class="n">dyn</span><span class="o">.</span><span class="n">ForcePositiveDefinite</span><span class="p">()</span>

    <span class="n">eos</span> <span class="o">=</span> <span class="n">calculate_eos</span><span class="p">(</span><span class="n">atoms</span><span class="p">)</span>
    <span class="n">v0</span><span class="p">,</span> <span class="n">e0</span><span class="p">,</span> <span class="n">B</span> <span class="o">=</span> <span class="n">eos</span><span class="o">.</span><span class="n">fit</span><span class="p">()</span>
    <span class="n">bulk</span> <span class="o">=</span> <span class="n">B</span> <span class="o">/</span> <span class="n">kJ</span> <span class="o">*</span> <span class="mf">1.0e24</span>

    <span class="k">return</span> <span class="n">dyn</span><span class="p">,</span> <span class="n">bulk</span>


<span class="c1"># Our input variables</span>
<span class="n">temperature</span> <span class="o">=</span> <span class="mf">100.0</span>
<span class="n">nconf</span> <span class="o">=</span> <span class="mi">1000</span>
<span class="n">max_pop</span> <span class="o">=</span> <span class="mi">1000</span>

<span class="c1"># Load in Tersoff potential</span>
<span class="n">pot</span> <span class="o">=</span> <span class="n">Potential</span><span class="p">(</span><span class="s1">&#39;IP Tersoff&#39;</span><span class="p">,</span> <span class="n">param_filename</span><span class="o">=</span>
<span class="s1">&#39;../06_the_SSCHA_with_machine_learning_potentials/ip.parms.Tersoff.xml&#39;</span><span class="p">)</span>
<span class="n">supercell</span> <span class="o">=</span> <span class="nb">tuple</span><span class="p">((</span><span class="mi">4</span><span class="o">*</span><span class="n">np</span><span class="o">.</span><span class="n">ones</span><span class="p">(</span><span class="mi">3</span><span class="p">,</span> <span class="n">dtype</span><span class="o">=</span><span class="nb">int</span><span class="p">))</span><span class="o">.</span><span class="n">tolist</span><span class="p">())</span>
<span class="n">dyn</span><span class="p">,</span> <span class="n">bulk</span> <span class="o">=</span> <span class="n">get_starting_dynamical_matrices</span><span class="p">(</span>
<span class="s1">&#39;../06_the_SSCHA_with_machine_learning_potentials/POSCAR&#39;</span><span class="p">,</span> <span class="n">pot</span><span class="p">,</span> <span class="n">supercell</span><span class="p">)</span>


<span class="c1"># Generate the ensemble and the minimizer objects</span>
<span class="n">ensemble</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Ensemble</span><span class="o">.</span><span class="n">Ensemble</span><span class="p">(</span><span class="n">dyn</span><span class="p">,</span> <span class="n">T0</span><span class="o">=</span><span class="n">temperature</span><span class="p">,</span>
<span class="n">supercell</span> <span class="o">=</span> <span class="n">dyn</span><span class="o">.</span><span class="n">GetSupercell</span><span class="p">())</span>
<span class="n">ensemble</span><span class="o">.</span><span class="n">generate</span><span class="p">(</span><span class="n">N</span> <span class="o">=</span> <span class="n">nconf</span><span class="p">)</span>
<span class="n">minimizer</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">SchaMinimizer</span><span class="o">.</span><span class="n">SSCHA_Minimizer</span><span class="p">(</span><span class="n">ensemble</span><span class="p">)</span>
<span class="n">minimizer</span><span class="o">.</span><span class="n">min_step_dyn</span> <span class="o">=</span> <span class="mf">0.1</span>
<span class="n">minimizer</span><span class="o">.</span><span class="n">kong_liu_ratio</span> <span class="o">=</span> <span class="mf">0.5</span>
<span class="n">minimizer</span><span class="o">.</span><span class="n">meaningful_factor</span> <span class="o">=</span> <span class="mf">0.001</span>
<span class="n">minimizer</span><span class="o">.</span><span class="n">max_ka</span> <span class="o">=</span> <span class="mi">1000</span>

<span class="c1"># Relax structure</span>
<span class="n">relax</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Relax</span><span class="o">.</span><span class="n">SSCHA</span><span class="p">(</span><span class="n">minimizer</span><span class="p">,</span> <span class="n">ase_calculator</span> <span class="o">=</span> <span class="n">pot</span><span class="p">,</span>
<span class="n">N_configs</span> <span class="o">=</span> <span class="n">nconf</span><span class="p">,</span> <span class="n">max_pop</span> <span class="o">=</span> <span class="n">max_pop</span><span class="p">,</span> <span class="n">save_ensemble</span> <span class="o">=</span> <span class="kc">True</span><span class="p">)</span>
<span class="n">relax</span><span class="o">.</span><span class="n">vc_relax</span><span class="p">(</span><span class="n">static_bulk_modulus</span> <span class="o">=</span> <span class="n">bulk</span><span class="p">,</span>
<span class="n">ensemble_loc</span> <span class="o">=</span> <span class="s2">&quot;directory_of_the_ensemble&quot;</span><span class="p">)</span>

<span class="c1"># Generate ensemble for third-order FC with the relaxed dynamical matrices</span>
<span class="n">new_ensemble</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">Ensemble</span><span class="o">.</span><span class="n">Ensemble</span><span class="p">(</span><span class="n">relax</span><span class="o">.</span><span class="n">minim</span><span class="o">.</span><span class="n">dyn</span><span class="p">,</span> <span class="n">T0</span><span class="o">=</span><span class="n">temperature</span><span class="p">,</span>
<span class="n">supercell</span> <span class="o">=</span> <span class="n">relax</span><span class="o">.</span><span class="n">minim</span><span class="o">.</span><span class="n">dyn</span><span class="o">.</span><span class="n">GetSupercell</span><span class="p">())</span>
<span class="n">new_ensemble</span><span class="o">.</span><span class="n">generate</span><span class="p">(</span><span class="n">N</span> <span class="o">=</span> <span class="n">nconf</span><span class="o">*</span><span class="mi">5</span><span class="p">)</span>
<span class="n">new_ensemble</span><span class="o">.</span><span class="n">compute_ensemble</span><span class="p">(</span><span class="n">pot</span><span class="p">,</span> <span class="n">compute_stress</span> <span class="o">=</span> <span class="kc">True</span><span class="p">,</span>
<span class="n">stress_numerical</span> <span class="o">=</span> <span class="kc">False</span><span class="p">,</span> <span class="n">cluster</span> <span class="o">=</span> <span class="kc">None</span><span class="p">,</span> <span class="n">verbose</span> <span class="o">=</span> <span class="kc">True</span><span class="p">)</span>
<span class="c1"># We minimize the free energy with this new ensemble</span>
<span class="n">new_minimizer</span> <span class="o">=</span> <span class="n">sscha</span><span class="o">.</span><span class="n">SchaMinimizer</span><span class="o">.</span><span class="n">SSCHA_Minimizer</span><span class="p">(</span><span class="n">new_ensemble</span><span class="p">)</span>
<span class="n">new_minimizer</span><span class="o">.</span><span class="n">minim_struct</span> <span class="o">=</span> <span class="kc">False</span>
<span class="n">new_minimizer</span><span class="o">.</span><span class="n">set_minimization_step</span><span class="p">(</span><span class="mf">0.1</span><span class="p">)</span>
<span class="n">new_minimizer</span><span class="o">.</span><span class="n">meaningful_factor</span> <span class="o">=</span> <span class="mf">0.001</span>
<span class="n">new_minimizer</span><span class="o">.</span><span class="n">max_ka</span> <span class="o">=</span> <span class="mi">10000</span>
<span class="n">new_minimizer</span><span class="o">.</span><span class="n">init</span><span class="p">()</span>
<span class="n">new_minimizer</span><span class="o">.</span><span class="n">run</span><span class="p">()</span>
<span class="n">new_minimizer</span><span class="o">.</span><span class="n">dyn</span><span class="o">.</span><span class="n">save_qe</span><span class="p">(</span><span class="s1">&#39;final_dyn&#39;</span><span class="p">)</span>
<span class="c1"># Update weights with a new dynamical matrice</span>
<span class="n">new_ensemble</span><span class="o">.</span><span class="n">update_weights</span><span class="p">(</span><span class="n">new_minimizer</span><span class="o">.</span><span class="n">dyn</span><span class="p">,</span> <span class="n">temperature</span><span class="p">)</span>

<span class="c1"># Calculate Hessian and the third order tensor (return_d3 = True)</span>
<span class="n">dyn_hessian</span><span class="p">,</span> <span class="n">d3_tensor</span> <span class="o">=</span> <span class="n">new_ensemble</span><span class="o">.</span><span class="n">get_free_energy_hessian</span><span class="p">(</span><span class="n">include_v4</span> <span class="o">=</span> <span class="kc">False</span><span class="p">,</span>
    <span class="n">get_full_hessian</span> <span class="o">=</span> <span class="kc">True</span><span class="p">,</span> <span class="n">return_d3</span> <span class="o">=</span> <span class="kc">True</span><span class="p">)</span>
<span class="n">np</span><span class="o">.</span><span class="n">save</span><span class="p">(</span><span class="s2">&quot;d3.npy&quot;</span><span class="p">,</span> <span class="n">d3_tensor</span><span class="p">)</span>
<span class="n">dyn_hessian</span><span class="o">.</span><span class="n">save_qe</span><span class="p">(</span><span class="s1">&#39;hessian_dyn&#39;</span><span class="p">)</span>
</pre></div>
</div>
<p>Here we used 4x4x4 supercell. <strong>You will need to converge results with respect to the size of the supercell</strong>. A good check for the convergence could be the decay of the second and third order force constants with the distance. Now that we have second and third order force constants, we can calculate lattice thermal conductivity. For this we provide following script:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">from</span> <span class="nn">__future__</span> <span class="kn">import</span> <span class="n">print_function</span>
<span class="kn">from</span> <span class="nn">__future__</span> <span class="kn">import</span> <span class="n">division</span>

<span class="kn">import</span> <span class="nn">numpy</span> <span class="k">as</span> <span class="nn">np</span>
<span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span>
<span class="kn">import</span> <span class="nn">cellconstructor.Phonons</span>
<span class="kn">import</span> <span class="nn">cellconstructor.ForceTensor</span>
<span class="kn">import</span> <span class="nn">cellconstructor.ThermalConductivity</span>
<span class="kn">import</span> <span class="nn">time</span>

<span class="n">dyn_prefix</span> <span class="o">=</span> <span class="s1">&#39;final_dyn&#39;</span>
<span class="n">nqirr</span> <span class="o">=</span> <span class="mi">8</span>

<span class="n">SSCHA_TO_MS</span> <span class="o">=</span> <span class="n">cellconstructor</span><span class="o">.</span><span class="n">ThermalConductivity</span><span class="o">.</span><span class="n">SSCHA_TO_MS</span>
<span class="n">RY_TO_THZ</span> <span class="o">=</span> <span class="n">cellconstructor</span><span class="o">.</span><span class="n">ThermalConductivity</span><span class="o">.</span><span class="n">SSCHA_TO_THZ</span>
<span class="n">dyn</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">Phonons</span><span class="o">.</span><span class="n">Phonons</span><span class="p">(</span><span class="n">dyn_prefix</span><span class="p">,</span> <span class="n">nqirr</span><span class="p">)</span>

<span class="n">supercell</span> <span class="o">=</span> <span class="n">dyn</span><span class="o">.</span><span class="n">GetSupercell</span><span class="p">()</span>

<span class="n">fc3</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">ForceTensor</span><span class="o">.</span><span class="n">Tensor3</span><span class="p">(</span><span class="n">dyn</span><span class="o">.</span><span class="n">structure</span><span class="p">,</span>
<span class="n">dyn</span><span class="o">.</span><span class="n">structure</span><span class="o">.</span><span class="n">generate_supercell</span><span class="p">(</span><span class="n">supercell</span><span class="p">),</span> <span class="n">supercell</span><span class="p">)</span>

<span class="n">d3</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">load</span><span class="p">(</span><span class="s2">&quot;d3.npy&quot;</span><span class="p">)</span>
<span class="n">fc3</span><span class="o">.</span><span class="n">SetupFromTensor</span><span class="p">(</span><span class="n">d3</span><span class="p">)</span>
<span class="n">fc3</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">ThermalConductivity</span><span class="o">.</span><span class="n">centering_fc3</span><span class="p">(</span><span class="n">fc3</span><span class="p">)</span>

<span class="n">mesh</span> <span class="o">=</span> <span class="p">[</span><span class="mi">10</span><span class="p">,</span><span class="mi">10</span><span class="p">,</span><span class="mi">10</span><span class="p">]</span>
<span class="n">smear</span> <span class="o">=</span> <span class="mf">0.03</span><span class="o">/</span><span class="n">RY_TO_THZ</span>

<span class="n">tc</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">ThermalConductivity</span><span class="o">.</span><span class="n">ThermalConductivity</span><span class="p">(</span><span class="n">dyn</span><span class="p">,</span> <span class="n">fc3</span><span class="p">,</span>
<span class="n">kpoint_grid</span> <span class="o">=</span> <span class="n">mesh</span><span class="p">,</span> <span class="n">scattering_grid</span> <span class="o">=</span> <span class="n">mesh</span><span class="p">,</span> <span class="n">smearing_scale</span> <span class="o">=</span> <span class="kc">None</span><span class="p">,</span>
<span class="n">smearing_type</span> <span class="o">=</span> <span class="s1">&#39;constant&#39;</span><span class="p">,</span> <span class="n">cp_mode</span> <span class="o">=</span> <span class="s1">&#39;quantum&#39;</span><span class="p">,</span> <span class="n">off_diag</span> <span class="o">=</span> <span class="kc">False</span><span class="p">)</span>

<span class="n">temperatures</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">linspace</span><span class="p">(</span><span class="mi">200</span><span class="p">,</span><span class="mi">1200</span><span class="p">,</span><span class="mi">10</span><span class="p">,</span><span class="n">dtype</span><span class="o">=</span><span class="nb">float</span><span class="p">)</span>
<span class="n">start_time</span> <span class="o">=</span> <span class="n">time</span><span class="o">.</span><span class="n">time</span><span class="p">()</span>
<span class="n">tc</span><span class="o">.</span><span class="n">setup_harmonic_properties</span><span class="p">(</span><span class="n">smear</span><span class="p">)</span>
<span class="n">tc</span><span class="o">.</span><span class="n">write_harmonic_properties_to_file</span><span class="p">()</span>

<span class="n">tc</span><span class="o">.</span><span class="n">calculate_kappa</span><span class="p">(</span><span class="n">mode</span> <span class="o">=</span> <span class="s1">&#39;SRTA&#39;</span><span class="p">,</span> <span class="n">temperatures</span> <span class="o">=</span> <span class="n">temperatures</span><span class="p">,</span>
<span class="n">write_lifetimes</span> <span class="o">=</span> <span class="kc">True</span><span class="p">,</span> <span class="n">gauss_smearing</span> <span class="o">=</span> <span class="kc">True</span><span class="p">,</span> <span class="n">offdiag_mode</span> <span class="o">=</span> <span class="s1">&#39;wigner&#39;</span><span class="p">,</span>
<span class="n">kappa_filename</span> <span class="o">=</span> <span class="s1">&#39;Thermal_conductivity_SRTA&#39;</span><span class="p">,</span> <span class="n">lf_method</span> <span class="o">=</span> <span class="s1">&#39;fortran-LA&#39;</span><span class="p">)</span>

<span class="nb">print</span><span class="p">(</span><span class="s1">&#39;Calculated SSCHA kappa in: &#39;</span><span class="p">,</span> <span class="n">time</span><span class="o">.</span><span class="n">time</span><span class="p">()</span> <span class="o">-</span> <span class="n">start_time</span><span class="p">)</span>

<span class="n">tc</span><span class="o">.</span><span class="n">calculate_kappa</span><span class="p">(</span><span class="n">mode</span> <span class="o">=</span> <span class="s1">&#39;GK&#39;</span><span class="p">,</span> <span class="n">write_lineshapes</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span>
<span class="n">ne</span> <span class="o">=</span> <span class="mi">1000</span><span class="p">,</span> <span class="n">temperatures</span> <span class="o">=</span> <span class="n">temperatures</span><span class="p">,</span>
<span class="n">kappa_filename</span> <span class="o">=</span> <span class="s1">&#39;Thermal_conductivity_GK&#39;</span><span class="p">)</span>

<span class="nb">print</span><span class="p">(</span><span class="s1">&#39;Calculated SSCHA kappa in: &#39;</span><span class="p">,</span> <span class="n">time</span><span class="o">.</span><span class="n">time</span><span class="p">()</span> <span class="o">-</span> <span class="n">start_time</span><span class="p">)</span>
<span class="c1"># Save ThermalConductivity object for postprocessing.</span>
<span class="n">tc</span><span class="o">.</span><span class="n">save_pickle</span><span class="p">()</span>
</pre></div>
</div>
<p>Important parts of the script are:</p>
<blockquote>
<div><ul>
<li><p>We define mesh on which we calculate phonon properties to be the same as the mesh we are calculating scattering processes (variable <code class="code docutils literal notranslate"><span class="pre">mesh</span></code>). This does not have to be true. In most cases <code class="code docutils literal notranslate"><span class="pre">scattering_grid</span></code> can be much smaller than <code class="code docutils literal notranslate"><span class="pre">kpoint_grid</span></code>. <strong>Converge your results with respect to both grids.</strong></p></li>
<li><p>We use smearing approach to satisfy energy conservation laws. There are two ways: constant and adaptive. In the case of <code class="code docutils literal notranslate"><span class="pre">smearing_type</span> <span class="pre">=</span> <span class="pre">'constant'</span></code> we have to provide smearing value in <strong>Ry</strong> as the argument to <code class="code docutils literal notranslate"><span class="pre">setup_harmonic_properties</span></code> function. In case we choose adaptive smearing, the smearing constant will be different for different phonon modes. We still can define global variable <code class="code docutils literal notranslate"><span class="pre">smearing_scale</span></code> with which we multiply precomputed smearing constants. <code class="code docutils literal notranslate"><span class="pre">smearing_scale</span> <span class="pre">=</span> <span class="pre">1.0</span></code> works pretty well in most cases. <strong>Converge your results with respect to smearing variables.</strong></p></li>
<li><p><code class="code docutils literal notranslate"><span class="pre">off_diag</span></code> variable defines whether we are doing calculation with what was termed as <em>coherent transport</em>. This will be important for highly anharmonic materials with large bunching of phonon modes.</p></li>
<li><p>Function <code class="code docutils literal notranslate"><span class="pre">calculate_kappa</span></code> does most of the work. Here we will describe main options:</p>
<blockquote>
<div><ul class="simple">
<li><p><code class="code docutils literal notranslate"><span class="pre">mode</span></code> defines which method to use to calculate lattice thermal conductivity. Options are <strong>SRTA</strong> which is Boltzmann transport equation solution in single relaxation time approximation and <strong>GK</strong> (<a class="reference external" href="https://www.nature.com/articles/s41524-021-00523-7">Dangic et al.</a>) which is Green-Kubo method that uses phonon spectral functions instead of phonon lifetimes. These two modes should give similar results in low anharmonicity materials, but different in strongly anharmonic ones.</p></li>
<li><p><code class="code docutils literal notranslate"><span class="pre">gauss_smearing</span></code> defines how we treat energy conservation in the calculation of self energy. If <strong>True</strong> it will use Gaussian functions, if <strong>False</strong> it will use Lorentzian functions. In case of Gaussian smearing real part of the self energy is calculated using Kramers-Kronig transformation.</p></li>
<li><p><code class="code docutils literal notranslate"><span class="pre">offdiag_mode</span></code> defines how we calculate coherent transport if <code class="code docutils literal notranslate"><span class="pre">mode</span> <span class="pre">=</span> <span class="pre">'SRTA'</span></code>. Two options: <strong>wigner</strong> (<a class="reference external" href="https://www.nature.com/articles/s41567-019-0520-x">Simoncelli et al.</a>) and <strong>gk</strong> (<a class="reference external" href="https://www.nature.com/articles/s41467-019-11572-4">Isaeva et al.</a>). If <code class="code docutils literal notranslate"><span class="pre">mode</span></code> is <strong>GK</strong>, coherent transport is included naturally.</p></li>
<li><p><code class="code docutils literal notranslate"><span class="pre">lf_method</span></code> defines how lifetimes are calculated in case <code class="code docutils literal notranslate"><span class="pre">mode</span> <span class="pre">=</span> <span class="pre">'SRTA'</span></code>. In short you want to keep <em>fortan-</em>, and then add <em>LA</em> or <em>P</em>. These should give more or less same results. Additional option is <em>SC</em> where we solve phonon lifetimes self-consistently, meaning we account for the phonon lineshifts.</p></li>
<li><p><code class="code docutils literal notranslate"><span class="pre">ne</span></code> defines the number of frequency steps if we are calculating phonon lineshapes. Also important in case of <code class="code docutils literal notranslate"><span class="pre">lf_method</span> <span class="pre">=</span> <span class="pre">'SC'</span></code> because we solve self-consistent equation on a grid of frequency values linearly interpolating real and imaginary part. Larger is better. <strong>Converge your results with respect to ne.</strong></p></li>
</ul>
</div></blockquote>
</li>
</ul>
</div></blockquote>
<p>This calculation should take a few minutes. The results are save in the <strong>kappa_filename</strong>.</p>
<div class="topic">
<p class="topic-title">Question:</p>
<p>If we check results we see that <em>SRTA</em> and <em>GK</em> results are different. Why? How can we improve this calculation?</p>
</div>
</section>
<section id="lattice-thermal-conductivity-of-gete">
<h2>Lattice thermal conductivity of GeTe<a class="headerlink" href="#lattice-thermal-conductivity-of-gete" title="Permalink to this headline"> </a></h2>
<p>As a second example we will calculate lattice thermal conductivity of GeTe. GeTe is a highly anharmonic material with a phase transition from rhombohedral to cubic phase at around 700 K. This means its lattice thermal conductivity is very low. Additionally, it should show difference between <em>SRTA</em> and <em>GK</em> methods.</p>
<p>For SSCHA minimization we can calculate atomic properties using <a class="reference external" href="https://archive.materialscloud.org/record/2021.42">Gaussian Approximation Potential</a> developed for this material. However, in the interest of time we provided the dynamical matrices calculated at 0 K and the third order force constantsin the folder 09_Thermal_conductivity_calculations_with_the_SSCHA.</p>
<div class="topic">
<p class="topic-title">Exercise:</p>
<p>Calculate lattice thermal conductivity of GeTe up to 1200 K (sample temperature from 300 K every 200 K). Is there a difference between <em>GK</em> and <em>SRTA</em> methods?</p>
</div>
<div class="topic">
<p class="topic-title">Exercise:</p>
<p>Check if coherent transport has an influence on thermal conductivity in this material system.</p>
</div>
<p>Finally, in case we want to do some postprocessing we can load in the previously saved <code class="code docutils literal notranslate"><span class="pre">ThermalConductivity</span></code> object and access all previously calculated data. For example, we can calculate phonon density of states calculated with auxiliary force constants and the one calculated with phonon lineshapes:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="kn">from</span> <span class="nn">__future__</span> <span class="kn">import</span> <span class="n">print_function</span>
<span class="kn">from</span> <span class="nn">__future__</span> <span class="kn">import</span> <span class="n">division</span>

<span class="c1"># Import the modules to read the dynamical matrix</span>
<span class="kn">import</span> <span class="nn">numpy</span> <span class="k">as</span> <span class="nn">np</span>
<span class="kn">import</span> <span class="nn">cellconstructor</span> <span class="k">as</span> <span class="nn">CC</span>
<span class="kn">import</span> <span class="nn">cellconstructor.Phonons</span>
<span class="kn">import</span> <span class="nn">cellconstructor.ForceTensor</span>
<span class="kn">import</span> <span class="nn">cellconstructor.ThermalConductivity</span>
<span class="kn">import</span> <span class="nn">time</span>
<span class="kn">import</span> <span class="nn">matplotlib.pyplot</span> <span class="k">as</span> <span class="nn">plt</span>
<span class="kn">import</span> <span class="nn">matplotlib.gridspec</span> <span class="k">as</span> <span class="nn">gridspec</span>

<span class="n">tc</span> <span class="o">=</span> <span class="n">CC</span><span class="o">.</span><span class="n">ThermalConductivity</span><span class="o">.</span><span class="n">load_thermal_conductivity</span><span class="p">()</span>

<span class="c1"># See at which temperatures we calculated stuff</span>
<span class="n">tc</span><span class="o">.</span><span class="n">what_temperatures</span><span class="p">()</span>

<span class="n">key</span> <span class="o">=</span> <span class="nb">list</span><span class="p">(</span><span class="n">tc</span><span class="o">.</span><span class="n">lineshapes</span><span class="o">.</span><span class="n">keys</span><span class="p">())</span> <span class="c1"># Get Ts for lineshapes</span>

<span class="c1"># DOS calculated from auxiliary force constants</span>
<span class="n">harm_dos</span> <span class="o">=</span> <span class="n">tc</span><span class="o">.</span><span class="n">get_dos</span><span class="p">()</span>
<span class="c1"># Temperature dependent DOS calculated from lineshapes</span>
<span class="c1"># first two arrays are raw data</span>
<span class="c1"># second two is gaussian smoothed results \</span>
<span class="c1">#for the distance between energy points de</span>
<span class="n">anharm_dos</span> <span class="o">=</span> <span class="n">tc</span><span class="o">.</span><span class="n">get_dos_from_lineshapes</span><span class="p">(</span><span class="nb">float</span><span class="p">(</span><span class="n">key</span><span class="p">[</span><span class="o">-</span><span class="mi">1</span><span class="p">]),</span> <span class="n">de</span> <span class="o">=</span> <span class="mf">0.1</span><span class="p">)</span>


<span class="c1"># Plot results</span>
<span class="n">fig</span> <span class="o">=</span> <span class="n">plt</span><span class="o">.</span><span class="n">figure</span><span class="p">(</span><span class="n">figsize</span><span class="o">=</span><span class="p">(</span><span class="mf">6.4</span><span class="p">,</span> <span class="mf">4.8</span><span class="p">))</span>
<span class="n">gs1</span> <span class="o">=</span> <span class="n">gridspec</span><span class="o">.</span><span class="n">GridSpec</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">)</span>
<span class="n">ax</span> <span class="o">=</span> <span class="n">fig</span><span class="o">.</span><span class="n">add_subplot</span><span class="p">(</span><span class="n">gs1</span><span class="p">[</span><span class="mi">0</span><span class="p">,</span> <span class="mi">0</span><span class="p">])</span>
<span class="n">ax</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="n">harm_dos</span><span class="p">[</span><span class="mi">0</span><span class="p">],</span> <span class="n">harm_dos</span><span class="p">[</span><span class="mi">1</span><span class="p">],</span> <span class="s1">&#39;k-&#39;</span><span class="p">,</span>
<span class="n">zorder</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span> <span class="n">label</span> <span class="o">=</span> <span class="s1">&#39;Harmonic&#39;</span><span class="p">)</span>
<span class="n">ax</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="n">anharm_dos</span><span class="p">[</span><span class="mi">0</span><span class="p">],</span> <span class="n">anharm_dos</span><span class="p">[</span><span class="mi">1</span><span class="p">],</span> <span class="s1">&#39;r-&#39;</span><span class="p">,</span>
<span class="n">zorder</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span> <span class="n">label</span> <span class="o">=</span> <span class="s1">&#39;Anharmonic raw @ &#39;</span> <span class="o">+</span> <span class="n">key</span><span class="p">[</span><span class="o">-</span><span class="mi">1</span><span class="p">]</span> <span class="o">+</span> <span class="s1">&#39; K&#39;</span><span class="p">)</span>
<span class="n">ax</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="n">anharm_dos</span><span class="p">[</span><span class="mi">2</span><span class="p">],</span> <span class="n">anharm_dos</span><span class="p">[</span><span class="mi">3</span><span class="p">],</span> <span class="s1">&#39;b-&#39;</span><span class="p">,</span>
<span class="n">zorder</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span> <span class="n">label</span> <span class="o">=</span> <span class="s1">&#39;Anharmonic smooth @ &#39;</span> <span class="o">+</span> <span class="n">key</span><span class="p">[</span><span class="o">-</span><span class="mi">1</span><span class="p">]</span> <span class="o">+</span> <span class="s1">&#39; K&#39;</span><span class="p">)</span>
<span class="n">ax</span><span class="o">.</span><span class="n">set_xlabel</span><span class="p">(</span><span class="s1">&#39;Frequency (THz)&#39;</span><span class="p">)</span>
<span class="n">ax</span><span class="o">.</span><span class="n">set_ylabel</span><span class="p">(</span><span class="s1">&#39;Density of states&#39;</span><span class="p">)</span>
<span class="n">ax</span><span class="o">.</span><span class="n">legend</span><span class="p">(</span><span class="n">loc</span> <span class="o">=</span> <span class="s1">&#39;upper right&#39;</span><span class="p">)</span>
<span class="n">ax</span><span class="o">.</span><span class="n">set_ylim</span><span class="p">(</span><span class="n">bottom</span> <span class="o">=</span> <span class="mf">0.0</span><span class="p">)</span>
<span class="n">fig</span><span class="o">.</span><span class="n">savefig</span><span class="p">(</span><span class="s1">&#39;test.pdf&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">show</span><span class="p">()</span>
</pre></div>
</div>
<p>Additionally, if we want to check a specific phonon lineshape (for example at <span class="math notranslate nohighlight">\(\Gamma\)</span> point) we can do it with a bit of hacking:</p>
<div class="highlight-python notranslate"><div class="highlight"><pre><span></span><span class="k">for</span> <span class="n">iqpt</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">tc</span><span class="o">.</span><span class="n">nkpt</span><span class="p">):</span>
    <span class="k">if</span><span class="p">(</span><span class="n">np</span><span class="o">.</span><span class="n">linalg</span><span class="o">.</span><span class="n">norm</span><span class="p">(</span><span class="n">tc</span><span class="o">.</span><span class="n">k_points</span><span class="p">[</span><span class="n">iqpt</span><span class="p">])</span> <span class="o">==</span> <span class="mf">0.0</span><span class="p">):</span>
        <span class="k">break</span>

<span class="n">energies</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">arange</span><span class="p">(</span><span class="nb">len</span><span class="p">(</span><span class="n">tc</span><span class="o">.</span><span class="n">lineshapes</span><span class="p">[</span><span class="n">key</span><span class="p">[</span><span class="o">-</span><span class="mi">1</span><span class="p">]][</span><span class="mi">0</span><span class="p">,</span><span class="mi">0</span><span class="p">]),</span>
<span class="n">dtype</span><span class="o">=</span><span class="nb">float</span><span class="p">)</span><span class="o">*</span><span class="n">tc</span><span class="o">.</span><span class="n">delta_omega</span> <span class="o">+</span> <span class="n">tc</span><span class="o">.</span><span class="n">delta_omega</span>
<span class="n">tc</span><span class="o">.</span><span class="n">write_lineshape</span><span class="p">(</span><span class="s1">&#39;Lineshape_at_Gamma&#39;</span><span class="p">,</span>
<span class="n">tc</span><span class="o">.</span><span class="n">lineshapes</span><span class="p">[</span><span class="n">key</span><span class="p">[</span><span class="o">-</span><span class="mi">1</span><span class="p">]][</span><span class="n">iqpt</span><span class="p">],</span> <span class="n">iqpt</span><span class="p">,</span> <span class="n">energies</span><span class="p">,</span> <span class="s1">&#39;no&#39;</span><span class="p">)</span>
</pre></div>
</div>
</section>
</section>
