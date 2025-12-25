#ifndef LGT_DESCRIPTOR_HPP
#define LGT_DESCRIPTOR_HPP

#include <iostream>
#include <boost/program_options.hpp>
#include "ansi_colors.hpp"

namespace po = boost::program_options;

class LGTDescriptor 
{
public:
    LGTDescriptor(int argc, char* argv[]);
    ~LGTDescriptor();
    po::options_description lgt_params{"Lattice Gauge Theory Parameters"};
    // Add public methods and members here
    int LS = 3; // Spatial lattice size
    int LT = 20; // Temporal lattice size
    int LT2 = 10; // Half of the temporal lattice size
    double g = 1.0; // Coupling constant. For historical/technical reasons, this is the inverse of the "g" used in the paper
    double kappa = 1.0; // Hopping parameter
    double beta = 1.0; // Inverse temperature
    double dtau = 0.05; // Euclidean time step size
    double dt = 0.1;     // Minkowski time step size for real-time correlators
    double tmax = 100.0; // Maximum time for real-time correlators
    int    nsteps = 1000; // Number of time steps for real-time correlators, calculated as tmax/dt
    double spatial_plaq_coeff = 0.0; // Spatial plaquette coefficient
    bool fermions = false; // Whether fermions are included
    bool spectral_functions = false;
    bool static_approximation = false;
    bool jj_correlators = false;
    bool dd_correlators = false;
    bool real_time_correlators = false;
    bool internal_checks = false;
    //Spectral function params
    double wmax = 0.0;
    int    NW   = 0;
    double dw   = 0.0;
    void init_parameters();
	void print_parameters();
    char* suffix0 = NULL; //A suffix for output files based on the parameters
    //Some quick basic calculations
    double Z_pure_lgt();
    double E2_pure_lgt();
};

#endif // LGT_DESCRIPTOR_HPP
