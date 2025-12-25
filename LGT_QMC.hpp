#ifndef LGT_QMC_HPP
#define LGT_QMC_HPP

#include <fftw3.h>
#include <numbers>
#include <random>
#include <algorithm>
#include "LGTDescriptor.hpp"
#include "linalg.hpp"

class LGT_QMC : public LGTDescriptor 
{
public:
    LGT_QMC(int argc, char* argv[]);
    ~LGT_QMC();
    po::options_description qmc_params{"QMC Simulation Parameters"};
    //Random number generators
    std::random_device rd;
    std::ranlux48 gen;
    std::uniform_real_distribution<double> angle_dist{-std::numbers::pi, std::numbers::pi};
    std::uniform_real_distribution<double> unit_dist{0.0, 1.0};
    std::normal_distribution<double>       unit_normal_dist{0.0, 1.0};
    // Add more methods or members as needed
    int    NMC      = 100000; // Number of Monte Carlo steps
    int    dmeas    = 10;   // Measurement step for the observables
    void print_parameters();
    //QMC variables
    double*      phi = NULL;
    double*    theta = NULL;
    double        fs = 0.0; //Winding number in the Villain model
    double        pa = 0.0; //Plaquette action
    double      fdet = 1.0; //Fermionic determinant
    double sigma_nw2 = 0.0; // Variance of the winding number distribution
    //Probability distribution of the winding number
    int nwmax; // Maximum winding number 
    double* nwp = NULL; // Probability distribution of the winding number
    //Exponentiated fields for the fermionic determinant calculation
    t_complex* ephi = NULL;
    t_complex* etheta = NULL;
    //Products of transfer matrices
    t_complex* T0Tt    = NULL;
    t_complex* TtT0    = NULL;
    void fermionic_hamiltonian(int it, t_complex* res);
    void transfer_matrix(int it, t_complex* res);
    void current_matrix(int it, t_complex* res);
    void d_matrix(t_complex* res);
    double fermion_det();
    void   fermion_det(double* out){ out[0] = fermion_det(); };
    void update_exps();
    //Fermionic observables
    void JJ_func(double* out);
    void DD_func(double* out);
    void real_time_DD_func(double* out);
    void diagonalize_single_particle_hamiltonian();
    void OO_spf(t_complex* O, double* hist); //Update the spectral function
    t_complex* F0 = NULL; // Temporary storage for the fermionic fermi factors
    t_complex* Ft = NULL; // Temporary storage for the fermionic fermi factors
    t_complex* J0 = NULL; // Temporary storage for the fermionic current operator
    t_complex* Jt = NULL; // Temporary storage for the fermionic current operator
    //Temporary storage for transfer matrix eigensystem calculation
    /*t_complex* revecs = NULL;
    t_complex* levecs = NULL;
    t_complex* evals  = NULL;*/
    t_complex* evecs = NULL;
    double*    evals = NULL;
    //Backup of old fields in case of update rejection + temporary storage
    double*     phi_backup    = NULL;
    double*     theta_backup  = NULL;
    t_complex*  ephi_backup   = NULL;
    t_complex*  etheta_backup = NULL;
    t_complex*  T0Tt_backup   = NULL;
    t_complex*  TtT0_backup   = NULL;
    t_complex*  tmp_matrix    = NULL;
    t_complex*  tmp_matrix1   = NULL;
    t_complex*  tmp_matrix2   = NULL;
    t_complex*  tmp_matrix3   = NULL;
    void backup_fields();
    void restore_fields();
    //FFTW vars
    fftw_complex* fx = NULL;
    fftw_plan fftw_spatial, fftw_temporal;
    //QMC generator functions
    void random_scalar_field(double* out, int N, fftw_plan& plan, int stride=1);
    void gauge_field_update();
    void static_gauge();
    void static_projection();
    //Diagnostic variables
    int accepted = 0; // Number of accepted Metropolis updates
    //Various bits in the action function
    double plaquette_action();
    //Observables functions
    void electric_field_correlator(double* out);
    //Test functions
    void test_random_scalar_field(int ntrials, int N, fftw_plan& plan, int reporting_interval=0);
    void rand_spatial_matrix(t_complex* out);
};

#endif // LGT_QMC_HPP
