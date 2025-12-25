#include "LGTDescriptor.hpp"

LGTDescriptor::LGTDescriptor(int argc, char* argv[]) {
    // Constructor implementation
    lgt_params.add_options()
        (                    "LS",       po::value<int>(    &(LS)                   )->default_value(         3), "Spatial lattice size")
        (                    "LT",       po::value<int>(    &(LT)                   )->default_value(        20), "Temporal lattice size")
        (                     "g",       po::value<double>( &(g)                    )->default_value(       1.0), "Coupling constant (inverse of the 'g' used in the paper)")
        (                 "kappa",       po::value<double>( &(kappa)                )->default_value(       1.0), "Hopping parameter")
        (                  "beta",       po::value<double>( &(beta)                 )->default_value(       1.0), "Inverse temperature")
        (                   "spc",       po::value<double>( &(spatial_plaq_coeff)   )->default_value(       0.0), "Spatial plaquette coefficient")
        (                    "NW",       po::value<int>(    &(NW)                   )->default_value(       200), "Number of bins for spectral function calculation")
        (                  "wmax",       po::value<double>( &(wmax)                 )->default_value(       2.0), "Frequency cutoff for spectral function calculation")
        (              "fermions",       po::bool_switch(   &(fermions)             )->default_value(     false), "Include fermions in the model")
        (    "spectral_functions",       po::bool_switch(   &spectral_functions     )->default_value(     false), "Enable spectral functions")
        (  "static_approximation",       po::bool_switch(   &static_approximation   )->default_value(     false), "Enable static projection")
        (        "jj_correlators",       po::bool_switch(   &jj_correlators         )->default_value(     false), "Enable JJ correlators")
        (        "dd_correlators",       po::bool_switch(   &dd_correlators         )->default_value(     false), "Enable DD correlators")
        (       "internal_checks",       po::bool_switch(   &internal_checks        )->default_value(     false), "Perform internal checks (e.g. eigensystem errors etc.)")
        ( "real_time_correlators",       po::bool_switch(   &real_time_correlators  )->default_value(     false), "Calculate real-time correlators")
        (                   "dt",        po::value<double>( &(dt)                   )->default_value(        dt), "Time step size")
        (                 "tmax",        po::value<double>( &(tmax)                 )->default_value(      tmax), "Maximum time for real-time correlators");

	//Reading parameters from the command line
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(lgt_params).allow_unregistered().run(), vm);
	po::notify(vm);

    suffix0 = new char[256];
    sprintf(suffix0, "LS%d_LT%d_g%.2f_kappa%.2f_beta%.2f_spc%.2f%s", LS, LT, g, kappa, beta, spatial_plaq_coeff, (fermions ? "_fermions" : "_pure_lgt"));

	// Check if help is needed
	if(vm.count( "help" )){std::cout<<lgt_params<<std::endl;};

	init_parameters();
	print_parameters();
}

void LGTDescriptor::init_parameters() 
{
    // Initialize any additional parameters if needed
    dtau  = beta / LT; // Euclidean time step size
    LT2 = LT / 2;
    dw  = wmax/(double)NW;
    nsteps = (int)std::ceil(tmax/dt);
}

void LGTDescriptor::print_parameters() 
{
    std::cout << ansi::cyan << "Parameters of 2D U(1) Lattice Gauge Theory" << (fermions ? " (with fermions)" : " (without fermions)") << ansi::reset << std::endl;

    std::cout << ansi::green << "\tSpatial lattice size (LS):                      " << ansi::yellow << LS                 << std::endl;
    std::cout << ansi::green << "\tTemporal lattice size (LT):                     " << ansi::yellow << LT                 << std::endl;
    std::cout << ansi::green << "\tCoupling constant (g):                          " << ansi::yellow << g                  << std::endl;
    if(fermions)
    std::cout << ansi::green << "\tHopping parameter (kappa):                      " << ansi::yellow << kappa              << std::endl;
    std::cout << ansi::green << "\tInverse temperature (beta):                     " << ansi::yellow << beta               << std::endl;
    std::cout << ansi::green << "\tEuclidean time step size (dtau):                " << ansi::yellow << dtau                 << std::endl;
    std::cout << ansi::green << "\tSpatial plaquette coefficient (spc):            " << ansi::yellow << spatial_plaq_coeff << std::endl;
    std::cout << ansi::green << "\tSuffix for output files:                        " << ansi::yellow << suffix0            << std::endl;

    if(spectral_functions)
    {
    std::cout << ansi::green << "\tFrequency cutoff for the spectral function:     " << ansi::yellow << wmax               << std::endl;
    std::cout << ansi::green << "\tNumber of bins for spectral function histogram: " << ansi::yellow << NW                 << std::endl;
    std::cout << ansi::green << "\tBin size for spectral function histogram: "       << ansi::yellow << dw                 << std::endl;
    };

    if(real_time_correlators)
    {
    std::cout << ansi::green << "\tMaximum time for real-time correlators (tmax):  " << ansi::yellow << tmax               << std::endl;
    std::cout << ansi::green << "\tMinkowski time step size (dt):                  " << ansi::yellow << dt                 << std::endl;
    std::cout << ansi::green << "\tNumber of time steps for real-time correlators: " << ansi::yellow << nsteps             << std::endl;
    };

    std::cout << ansi::green << "\tCalculating the spectral functions?             " << ansi::yellow << (spectral_functions?    "YES" : "NO") << std::endl;
    std::cout << ansi::green << "\tCalculating real-time correlators?              " << ansi::yellow << (real_time_correlators? "YES" : "NO") << std::endl;
    std::cout << ansi::green << "\tPerforming the static projection?               " << ansi::yellow << (static_approximation?  "YES" : "NO") << std::endl;
    std::cout << ansi::green << "\tCalculating the current-current correlators?    " << ansi::yellow << (jj_correlators?        "YES" : "NO") << std::endl;
    std::cout << ansi::green << "\tCalculating the D-D correlators?                " << ansi::yellow << (dd_correlators?        "YES" : "NO") << std::endl;
    std::cout << ansi::green << "\tPerforming internal checks?                     " << ansi::yellow << (internal_checks?       "YES" : "NO") << std::endl;

    std::cout << ansi::reset << std::endl;
}

LGTDescriptor::~LGTDescriptor() {
    // Destructor implementation
    delete [] suffix0;
}
// Add method implementations here

double LGTDescriptor::Z_pure_lgt()
{
    double Z, nZ = 1.0;
    int n = 1;
    do
    {
       Z = nZ; 
       double W = exp(-beta*LS*n*n/(2.0*g*g));
       nZ = Z + 2.0*W;
       n++;
    } while(nZ > Z && n < 1000);
    return Z;
}

double LGTDescriptor::E2_pure_lgt()
{
    double Z, E2, nZ = 1.0, nE2 = 0.0;
    int n = 1;
    do
    {
       Z   = nZ; E2 = nE2; 
       double W   = exp(-beta*LS*n*n/(2.0*g*g));
       nZ  = Z  + 2.0*W;
       nE2 = E2 + 2.0*W*(n*n);
       n++;
    } 
    while((nZ > Z || nE2 > E2) && n < 1000); // Continue until the terms become negligible
    return E2/Z;
}