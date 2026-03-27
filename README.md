# Metastability of the active Potts model (hydrodynamics)

Codes on the model studied in the scientific publication: S. Chatterjee, M. Karmakar, M. Mangeat, H. Rieger, and R. Paul, <i>Stability of discrete-symmetry flocks: Sandwich state, traveling domains and motility-induced pinning</i>, <a href='https://journals.aps.org/pre/abstract/10.1103/1r19-ryx9'>Phys. Rev. E <b>112</b>, 064115 (2025)</a>. A preprint is available on <a href='https://arxiv.org/abs/2507.08187'>arXiv</a>. The equations solved here are slightly different from Eq. (5) of this reference.</br></br>
A C++ code to compute the numerical solutions of the hydrodynamic equations is available in this repository, as well as a Python code to generate movies of the system dynamics, published on <a href='https://www.youtube.com/playlist?list=PLid9FhdFLQVn_GB4a_3MrzeTDo9EsABKP'>youtube</a>.</br></br>
<b>Exportations:</b> magnetization snapshots.</br>
<b>Compile:</b> g++ APM_meta_hydro_omp.cpp -fopenmp -lgsl -lgslcblas -lm -O3 -s -o APM_meta_hydro_omp.out.</br>
<b>Run:</b> ./APM_meta_hydro_omp.out -parameter=value.</br>
<b>Parameters:</b> beta, epsilon, rho0, Rd, rhod, LX, LY, dt, tmax, dx, init, threads (details as comments in the code).</br>

<b>Generate the movie:</b> python figure_APM_meta_dynamics2d.py -parameter=value.</br>
<b>Parameters:</b> beta, epsilon, rho0, Rd, rhod, LX, LY, tmax, init, NCPU, movie (details as comments in the code).</br>
