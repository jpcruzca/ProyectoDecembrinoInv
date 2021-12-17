Readme file: readme_0.txt                  Last modified:        June 10, 2014

This is the readme file for code "He4_0.java" that contains all the routines
for the variational and diffusion quantum Monte Carlo simulations of helium-4
clusters on a graphite surface.  The Aziz (HFD-B3-FCI1) potential [PRL 74, 1586
(1995)] is adopted for the interaction between any two helium atoms.  The JLB
potential [Surface Sci. 264, 419 (1992)] is used for the interaction between
any helium atom and the graphite surface and A linear combination of local
orbitals around the virtual lattice sites, modulated by a surface state, is
used as the single-particle part of the variational/guide wavefunction that
also contains two-body Jastrow factors for the correlation effects.  The code
He4_0.java creates the class of the entire simulation.  File "fort.15" stores
the additional input data.

I. Methods/functions: 

   initiate: Initialize the simulation.
   cqmc: Gather and analyze the autocorrelation of the data points.
   vqmc: Gather and analyze the variational Monte Carlo data.
   dqmc: Gather and analyze the diffusion Monte Carlo data.
   configure: Initiate a configuration.
   wave: Set up the wavefunction and corresponding weight.
   energy: Evaluate the local energy of a given configuration. 
   drift: Create the drift velocity.
   ensemble: Initiate an ensemble of configurations.
   metropolis: Move a configuration one Metropolis step.
   move: Move the ensemble one time step.
   pair: Create histogram for the pair-distribution.
   auto: Evaluate the autocorrelation function.
   save: Save the data for continuation runs.
   ranf: Create a uniform random number.
   rang: Generate two Gaussian random numbers.
   xe: exponential function exp(-t).
   u: function u in the Jastrow factor. 
   u1: du/dr of the Jastrow factor. 
   u2: du^2/dr^2 of the Jastrow factor. 
   up: Laplace operation of u.
   v: The Aziz potential.
   ue: The JLB surface potential.


Each routine here is written as a method in Java, for example, configure()
for configure method.  There are common blocks stored under class He4_0 before
any method is initiated and these variable are common to all the methods.
Functions needed for constructing the guiding wavefunction and the external
and interaction potential are created at the end of the file.

A combination of constants hbar^2/(m k)=12.119232 K angstroms^2 is used, where
m is the mass of a helium-4 atom, hbar is the Planck constant, and k is the
Boltzmann constant. Under such a choice, the energy is measured in K, length
in angstroms, imaginary time in 1/K.

If pair-distribution function is calculated (ig=2 or 3), a data file "gr.dat"
is created to store both g(r) and r.


II. Global Constants, Parameters, and Variables:

(1) Constants:

    hm: Combination hbar^2/(m k).
    rm: The distance of the minimum Aziz potential r_m.
    e0: Interaction energy scale epsilon/k.
    u0: Coefficient for the exponent term in the surface potential.
    g0: Coefficient gamma in the surface potential.
    a3: Coefficient for the 1/y^3 term in the surface potential.
    a4: Coefficient for the 1/y^4 term in the surface potential.
    a0: Coefficient alpha in the Aziz potential.
    b0: Coefficient beta in the Aziz potential.
    v0 Coefficient v_0 in the Aziz potential.
    c6 Coefficient a_6 in the Aziz potential.
    c8 Coefficient a_8 in the Aziz potential.
    c10 Coefficient a_{10} in the Aziz potential.
    d: Coefficient d in the Aziz potential.

    ax: x-increment for setting up the initial configuration.
    axh: ax/2.
    ay: y-increment (sqrt{3}ax/2) for setting up the initial configuration.
    xh: Half of the initial length of the cluster in x-direction.
    yh: Half of the initial length of the cluster in y-direction.
    dx: Twice of the Metropolis step size in x-direction.
    dy: Twice of the Metropolis step size in y-direction.
    dz: Twice of the Metropolis step size in z-direction.
    dr: Spacing size in the pair-distribution function.
    dt: Time step in the diffusion move.
    st: chi*sqrt(2)=sqrt(hm*dt), diffusion step size.


(2) Integer Parameters and Variables:

    lx: Number of particles in a row of the initial unit block.
    ly: Number of rows of the initial unit block.
    N: Number of particles.
    L: Number of virtual lattice sites.
    nv: Total number of coordinates of the particles. 
    nc: Number of data points for autocorrelation calculation.
    ne: Maximum number of configurations allowed in the ensemble. 
    bmax: The total number of points in the pair-distribution function.

    seed: Random number generator seed.
    iq: Control parameter for variational (1) or diffusion (2) simulation.
    is: Control parameter for saving (1 or 2) or reading in (2) data.
    ic: Control parameter for auto correlation if not equal to 0;
    ig: Control parameter for pair-distribution if not equal to 0;
    ia: The total number of the accepted Metropolis steps.

    np: The starting number of configurations in the ensemble. 
    dn: The change allowed in the number of configurations in the ensemble. 
    ng: Number of groups used is gathering data.
    nf: Number of Monte Carlo steps skipped between two adjacent data points.
    ns: Number of data points in each group of data.
    nr: Number of Monte Carlo steps in the equilibration runs.
    step: Total number of time steps taken in each individual run.


(3) Real Parameters and Variables:

    z0: Parameter z_0 in the single-particle wavefunction.
    ze: parameter z_e in the single-particle wavefunction.
    r0: Parameter r_0 in the single-particle wavefunction.
    a: Parameter a in the Jastrow factor.
    b: Parameter b in the Jastrow factor.
    c: Parameter c in the Jastrow factor.
    rho: Density of bulk helium-4 for reference in pair-distribution function.
    rho2D: Density of sqrt(3)*sqrt(3) commensurate lattice density.
    norm: 4*pi*rho/3 or pi*rho2D.
    hideal: Relative g(r)-1 for a corresponding ideal gas.
    ka: Parameter kappa used in branching.

    e: Total (or ensemble average) energy of the system.
    ek: Total (or ensemble average) kinetic energy of the system.
    ep: Total (or ensemble average) potential energy of the system.
    er: Adjusted reference energy of the system.
    w: Exponent of the guide/trial wavefunction squared, -2[J(R)-X(R)].
    ws: Sum of the birth/death rate.

(4) Variable Arrays:
    phi[N]: Single-particle wavefunction.
    x[N]: x coordinates of the particles in the system.
    y[N]: y coordinates of the particles in the system.
    z[N]: z coordinates of the particles in the system.
    zb[N]: z[i]-z_e.
    za[N]: zb[i]*zb[i]/(z0*z0).
    xx[N][N]: x_{ij} for all the pairs of particles in the system.
    yy[N][N]: y_{ij} for all the pairs of particles in the system.
    zz[N][N]: z_{ij} for all the pairs of particles in the system.
    zz[N][N]: z_{ij} for all the pairs of particles in the system.
    rr[N][N]: x_{ij}^2+y_{ij}^2+z_{ij}^2.
    sx[L]: x coordinates of the virtual lattice sites.
    sy[L]: y coordinates of the virtual lattice sites.
    sz[L]: z coordinates of the virtual lattice sites.
    xs[N][L]: x_i-sx_j.
    ys[N][L]: y_i-sy_j.
    zs[N][L]: z_i-sz_j.
    rs[N][L]: xs_{ij}^2+ys_{ij}^2+zs_{ij}^2.
    conf[nv]: Current configuration.
    ensm[nv][ne]: Current ensemble.
    vix[nv]: x component of the single-particle part of the drift velocity.
    viy[nv]: y component of the single-particle part of the drift velocity.
    viz[nv]: z component of the single-particle part of the drift velocity.
    vd[nv]: Drift velocity for the current configuration.
    eold[ne]: Energies for the previous ensemble.
    ecrl[nc]: Autocorrelation function for the energies.
    gr[bmax]: pair-distribution.
    hist[bmax]: Histogram for pair-distribution.


III. Local Variables Used in Methods/Functions

(1) initiate(): 
    t1 through t6: Current time (second to year).
    Other integers used as counters in the loops.

(2) vqmc() and dqmc():
    sum0: The acceptance rate of the Metropolis steps.
    sum1, sum2, sum3: Various summations for total energy.
    kum1, kum2, kum3: Various summations for kinetic energy.
    sig1: Variance on total energy calculated from all the data points.
    sig2: Variance on total energy by treating group average as a data point.
    kig1: Variance on kinetic energy calculated from all the data points.
    kig2: Variance on kinetic energy by treating group average as a data point.
    rl and ru: r and r+dr.
    Other integers used as counters in the loops.

(3) wave():
    psi: Product of single-particle wavefunctions.
    chi: X(R).
    ut: J(R).
    Other integers used as counters in the loops.

(4) energy():
    ei[N] and ez[N]: Arrays used to store parts of kinetic energy.
    ed, ek11, ek12, ek1 and ek2: Parts of kinetic energy.
    ed: |F_i|^2.
    et: T_i.

(5) drift():
    vx[N]: x component of a part of drift velocity.
    vy[N]: y component of a part of drift velocity.
    vz[N]: z component of a part of drift velocity.
    Other integers used as counters in the loops.

(6) metropolis():
    cfsv[nv]: Temporary storage of the current configuration.
    wold: The previous weight.
    Other integers used as counters in the loops.

(7) move():
    xg[2]: Two Gaussian random numbers.
    gn[nv]: Array of Gaussian random numbers.
    cfsv[nv]: Temporary storage of the current configuration.
    enew[ne]: Energies for the current ensemble.
    weight[ne]: Weights for the current ensemble.
    ensv[nv][ne]: Temporary ensemble.
    wbar: Variable for summing up weights.
    ebar: Variable for calculating average total energy.
    ekav: Variable for calculating average kinetic energy.
    Other integers used as counters in the loops.

(8) pair():
    rij: r_{ij}^2.
    Other integers used as counters in the loops.

(9) save():
    data.new: Output data file
    data.old: Input data file
    Other integers used as counters in the loops.

(10) auto():
     kmax: Number of points to be calculated for the autocorrelation function.
     corr: Autocorrelation function.
     sum11, sum12, sum21, sum22, sum30: Summations carried out.
     Other integers used as counters in the loops.

(11) ranf():
     a, c, q, and r: Temporary integer variables used.
     cd, Double value of c. 
     Other integers used in creating a random number.

(12) rang():
     x[2]; Array for two Gaussian random numbers.
     Other variables used in creating two Gaussian random numbers.


