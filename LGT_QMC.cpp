#include "LGT_QMC.hpp"
#include <algorithm>

LGT_QMC::LGT_QMC(int argc, char* argv[]) : LGTDescriptor(argc, argv), gen(rd())
{
    // Initialize QMC parameters
    qmc_params.add_options()
        ("NMC",   po::value<int>(&NMC  )->default_value(200000), "Number of Monte Carlo steps")
        ("dmeas", po::value<int>(&dmeas)->default_value(5     ), "Measurement step for the observables");
    
    po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(qmc_params).allow_unregistered().run(), vm);
	po::notify(vm);
	
	// Check if help is needed
	if(vm.count( "help" )){std::cout<<qmc_params<<std::endl;};

    // Allocate memory for QMC variables
    phi           = new double[LT*LS];           std::fill(phi,     phi    + LT*LS,        0.0);
    theta         = new double[LT*LS];           std::fill(theta,   theta  + LT*LS,        0.0);
    ephi          = new t_complex[LT*LS];
    etheta        = new t_complex[LT*LS];
    T0Tt          = new t_complex[(LT+1)*LS*LS]; std::fill(T0Tt,    T0Tt   + (LT+1)*LS*LS,  t_complex(0.0, 0.0));
    TtT0          = new t_complex[(LT+1)*LS*LS]; std::fill(TtT0,    TtT0   + (LT+1)*LS*LS,  t_complex(0.0, 0.0));
    //Allocate memory for the eigensystem calculation
    /*revecs        = new t_complex[LS*LS]; 
    levecs        = new t_complex[LS*LS];
    evals         = new t_complex[LS];*/ //TODO: deprecated, should eventually be removed
    evecs = new t_complex[LS*LS];
    evals = new double[LS];

    //Allocate memory for the backups
    phi_backup    = new double[LT*LS];
    theta_backup  = new double[LT*LS];
    ephi_backup   = new t_complex[LT*LS];
    etheta_backup = new t_complex[LT*LS];
    T0Tt_backup   = new t_complex[(LT+1)*LS*LS];
    TtT0_backup   = new t_complex[(LT+1)*LS*LS];

    tmp_matrix    = new t_complex[LS*LS];
    tmp_matrix1   = new t_complex[LS*LS];
    tmp_matrix2   = new t_complex[LS*LS];
    tmp_matrix3   = new t_complex[LS*LS];
    F0           = new t_complex[LS*LS];
    Ft           = new t_complex[LS*LS];
    J0           = new t_complex[LS*LS];
    Jt           = new t_complex[LS*LS];

    fs = 0.0; //Initial winding number is zero
    accepted = 0;

    fx  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * std::max(LS, LT));
    fftw_spatial  = fftw_plan_dft_1d(LS, fx, fx, FFTW_FORWARD,  FFTW_MEASURE);
    fftw_temporal = fftw_plan_dft_1d(LT, fx, fx, FFTW_FORWARD,  FFTW_MEASURE);

    //Initializing the probability distribution for the winding number
    //Probability distribution of the winding number is exp(-g^2*T(2*pi*k)^2/(2*LS))
    sigma_nw2 = (beta*(double)LS)/SQR(2.0*pi*g); // Variance of the winding number distribution 
    nwmax = 0;
    while(exp(-(nwmax*nwmax)/(2.0*sigma_nw2)) > 0.0) nwmax++;
    nwp = new double[nwmax];
    nwp[0] = 1.0; // Probability of zero winding number
    double normalization = 1.0; // Normalization factor
    for(int nw = 1; nw < nwmax; nw++) 
    {
        nwp[nw] = 2.0*exp(-(nw*nw)/(2.0*sigma_nw2));
        normalization += nwp[nw];
    };
    for(int nw = 0; nw<nwmax; nw++) 
        nwp[nw] /= normalization; // Normalize the probabilities

    print_parameters();
}

void LGT_QMC::gauge_field_update()
{
    //Back up all the core variables in case of update rejection
    double opa   = pa;
    double ofdet = fdet;
    double ofs   = fs; 
    backup_fields();

    //Randomly choosing the winding number, the probability is stored in probs
    int    nw = 0;
    double s  = nwp[0];
    double q  = unit_dist(gen);
    while((q>s)&&(nw<nwmax)) {nw++; s += nwp[nw];};
    if(unit_dist(gen)<0.5) nw = -nw; // Randomly choosing the sign of the winding number
    fs = 2.0*pi*(double)nw/(double)LS;

    std::fill(phi, phi + LS*LT, 0.0);
    random_scalar_field(phi + LS*(LT - 1), LS, fftw_spatial);

    double ff = sqrt(beta)/g;
    double phi_global = 8.0*angle_dist(gen); // Randomly choosing the global phase 
    for(int is = 0; is<LS; is++) 
        phi[LS*(LT-1) + is] = ff*phi[LS*(LT-1) + is] + fs*is + phi_global;

    //Now generate theta's for each is = 0 to LS-1
    ff = sqrt(dtau)/g;
    phi_global = 8.0*angle_dist(gen); // Randomly choosing the global phase
    for(int is = 0; is<LS; is++) 
    {
        random_scalar_field(theta + is, LT, fftw_temporal, LS);
        double tfs = phi[LS*(LT-1) + is] - phi[LS*(LT-1) + (is + 1)%LS] - (is==LS-1? fs*(double)LS : 0.0); //Theta field from the spatial plaquette
        tfs /= (double)LT; 
        for(int it = 0; it < LT; it++) 
            theta[LS*it + is] = ff*theta[LS*it + is] + tfs*it + phi_global;
    };
    update_exps(); //Update the complex exponentials

    pa = plaquette_action();
    double acceptance = exp(-pa + opa); //Acceptance ratio for the Metropolis update
    if(fermions)
    {
        fdet = fermion_det();
        acceptance *= std::min(1.0, fdet/ofdet); //Acceptance ratio for the fermionic determinant
    };
    if(unit_dist(gen) > acceptance) //Reject the proposal
    {
        restore_fields();
        fs = ofs; //Restore the old winding number
        pa = opa; //Restore the old plaquette action
        fdet = ofdet; //Restore the old fermionic determinant        //update_exps(); //Revert the complex exponentials to old values
    }
    else 
        accepted++; //Count the accepted Metropolis updates
}

void LGT_QMC::backup_fields()
{
    std::copy(phi,     phi    + LS*LT,        phi_backup);
    std::copy(theta,   theta  + LS*LT,        theta_backup);
    std::copy(ephi,    ephi   + LS*LT,        ephi_backup);
    std::copy(etheta,  etheta + LS*LT,        etheta_backup);
    std::copy(T0Tt,    T0Tt   + LS*LS*(LT+1), T0Tt_backup);
    std::copy(TtT0,    TtT0   + LS*LS*(LT+1), TtT0_backup);
}

void LGT_QMC::restore_fields()
{
    std::copy(phi_backup,    phi_backup    + LS*LT,        phi);
    std::copy(theta_backup,  theta_backup  + LS*LT,        theta);
    std::copy(ephi_backup,   ephi_backup   + LS*LT,        ephi);
    std::copy(etheta_backup, etheta_backup + LS*LT,        etheta);
    std::copy(T0Tt_backup,   T0Tt_backup   + LS*LS*(LT+1), T0Tt);
    std::copy(TtT0_backup,   TtT0_backup   + LS*LS*(LT+1), TtT0);
}

LGT_QMC::~LGT_QMC()
{
    // Free allocated memory
    delete[] phi;
    delete[] theta;
    delete[] ephi;
    delete[] etheta;
    delete[] T0Tt; 
    delete[] TtT0;

    /*delete[] revecs;
    delete[] levecs;*/ //TODO: deprecated, should eventually be removed
    delete[] evecs;
    delete[] evals;

    delete[] phi_backup;
    delete[] theta_backup;
    delete[] ephi_backup;
    delete[] etheta_backup;
    delete[] T0Tt_backup;
    delete[] TtT0_backup;

    delete[] tmp_matrix;
    delete[] tmp_matrix1;
    delete[] tmp_matrix2;
    delete[] tmp_matrix3;

    delete[] F0;
    delete[] Ft;
    delete[] J0;
    delete[] Jt;

    delete[] nwp;

    fftw_free(fx);
    fftw_destroy_plan(fftw_spatial);
    fftw_destroy_plan(fftw_temporal);
}

void LGT_QMC::print_parameters()
{
    std::cout << ansi::cyan << "Parameters of QMC Simulation: " << ansi::reset << std::endl;

    std::cout << ansi::green << "\tNumber of Monte Carlo steps:                    " << ansi::yellow << NMC       << std::endl;
    std::cout << ansi::green << "\tMeasurement step for the observables:           " << ansi::yellow << dmeas     << std::endl;
    std::cout << ansi::green << "\tVariance of the winding number distribution:    " << ansi::yellow << sigma_nw2 << std::endl;
    std::cout << ansi::reset << std::endl;
}

void LGT_QMC::random_scalar_field(double* out, int N, fftw_plan& plan, int stride)
{
    fx[0][0] = 0.0; fx[0][1] = 0.0;
    int nm = (N%2 == 0) ? N/2 : N/2 + 1;
    for(int n=1; n<nm; n++) 
    {
        double scale = 1.0/(2.0*sqrt(2.0*(double)N)*sin(pi*n/(double)N));
        fx[n][0]   = scale*unit_normal_dist(gen); fx[n][1]   = scale*unit_normal_dist(gen);
        fx[N-n][0] =  fx[n][0];                   fx[N-n][1] = -fx[n][1];
    };
    if(N%2 == 0)
    {
        fx[N/2][0] = 1.0/(2.0*sqrt((double)N))*unit_normal_dist(gen); fx[N/2][1] = 0.0;
    };

    fftw_execute(plan);
    
    for(int n = 0; n<N; n++) 
        out[n*stride] = fx[n][0];
}

void LGT_QMC::electric_field_correlator(double* out)
{
    double prefactor = SQR((g*g)/dtau);
    double E0 = theta[LS*1 + 0] - theta[LS*0 + 0] - phi[LS*0 + 1] + phi[LS*0 + 0]; //Electric field at the origin
    out[0] += g*g/dtau;
    for(int it = 0; it<LT; it++)
        for(int is = 0; is<LS; is++)
        {
            double Et = theta[LS*((it+1)%LT) + is] - theta[LS*it + is] - phi[LS*it + (is + 1)%LS] + phi[LS*it + is];
            if(is==LS-1 && it==LT-1) 
                Et -= fs*(double)LS; //Electric field at time it
            out[LS*it + is] -= prefactor*E0*Et;
        };
}

void LGT_QMC::test_random_scalar_field(int ntrials, int N, fftw_plan& plan, int reporting_interval)
{
    double* out = new double[N];
    double* res = new double[N*N]; std::fill(res, res + N*N, 0.0);

    if(reporting_interval == 0) reporting_interval = ntrials/50;

    for(int itrial = 0; itrial < ntrials; itrial++)
    {
        random_scalar_field(out, N, plan);
        for(int i=0; i<N; i++)
            for(int j=0; j<N; j++)
                res[i*N + j] += out[i]*out[j];

        if(itrial%reporting_interval == 0 && itrial > 0)
        {
            double err = 0.0;
            for(int i=0; i<N; i++)
                for(int j=0; j<N; j++)
                {
                    int ifw = (i+1)%N, ibw = (i-1+N)%N;
                    double loc_err = (2.0*res[i*N + j] - res[ifw*N + j] - res[ibw*N + j])/(double)itrial + 1.0/(double)N;
                    if(i==j) loc_err -= 1.0;
                    err += loc_err*loc_err;
                };
            std::cout << "Trial " << itrial << " Error: " << sqrt(err) << std::endl;
        };
    };

    delete[] out;
    delete[] res;
}

double LGT_QMC::plaquette_action()
{
    double S = 0.0;
    for(int it = 0; it<LT; it++)
    {
        double theta_total = 0.0;
        for(int is = 0; is<LS; is++)
            theta_total += theta[LS*it + is];
        S += (2.0 - 2.0*cos(theta_total));
    };
    return 0.5*g*g*dtau*spatial_plaq_coeff*S;
}

/* Here goes the fermionic stuff */

void LGT_QMC::fermionic_hamiltonian(int it, t_complex* res) 
{
    std::fill(res, res + LS*LS, t_complex(0.0, 0.0));
    for(int is = 0; is < LS; is++)
    {
        int isf = (is + 1)%LS; //Forward spatial index
        int isb = (is - 1 + LS)%LS; //Backward spatial index
        res[LS*is + isf] = kappa*etheta[LS*it + is]; //Forward hopping
        res[LS*is + isb] = kappa*std::conj(etheta[LS*it + isb]); //Backward hopping
    };
}

void LGT_QMC::current_matrix(int it, t_complex* res) 
{
    std::fill(res, res + LS*LS, t_complex(0.0, 0.0));
    for(int is = 0; is < LS; is++)
    {
        int isf = (is + 1)%LS; //Forward spatial index
        int isb = (is - 1 + LS)%LS; //Backward spatial index
        res[LS*is + isf] = t_complex(0.0,  1.0)*kappa*etheta[LS*it + is]; //Forward hopping
        res[LS*is + isb] = t_complex(0.0, -1.0)*kappa*std::conj(etheta[LS*it + isb]); //Backward hopping
    };
}

void LGT_QMC::transfer_matrix(int it, t_complex* res) 
{
    identity_matrix(res, LS);
    fermionic_hamiltonian(it, tmp_matrix1);
    A_pluseq_bB(res, t_complex(-dtau, 0.0), tmp_matrix1, LS*LS);
    A_eq_B_mult_C(tmp_matrix2, tmp_matrix1, tmp_matrix1, LS); //tmp_matrix2 = (H_F)^2
    A_pluseq_bB(res, t_complex(0.5*dtau*dtau, 0.0), tmp_matrix2, LS*LS);
    for(int is=0; is<LS; is++)
        for(int js=0; js<LS; js++)
            res[LS*is + js] *= ephi[LS*it + js];
}

double LGT_QMC::fermion_det()
{
    std::copy(T0Tt + LS*LS*LT, T0Tt + LS*LS*(LT+1), tmp_matrix);
    for(int i=0; i<LS; i++) tmp_matrix[i*LS + i] += 1.0; // Add the identity matrix
    double r = abs_det(tmp_matrix, LS);
    return r*r;
}

void LGT_QMC::update_exps()
{
    for(int it = 0; it < LT; it++)
        for(int is = 0; is < LS; is++)
        { 
            etheta[LS*it + is] = std::exp(t_complex(0.0, theta[LS*it + is]));
            ephi[LS*it + is]   = std::exp(t_complex(0.0, phi[LS*it + is]));
        };
    std::fill(T0Tt, T0Tt + LS*LS, t_complex(0.0, 0.0));
    for(int is=0; is<LS; is++) T0Tt[is*LS + is] = t_complex(1.0, 0.0); //Initialize T0Tt[0] to the identity matrix
    for(int it=1; it<=LT; it++)
    {
        transfer_matrix(it-1, tmp_matrix);
        A_eq_B_mult_C(T0Tt + LS*LS*it, T0Tt + LS*LS*(it-1), tmp_matrix, LS);
    };

    std::fill(TtT0 + LS*LS*LT, TtT0 + LS*LS*(LT+1), t_complex(0.0, 0.0));
    for(int is=0; is<LS; is++) TtT0[LS*LS*LT + is*LS + is] = t_complex(1.0, 0.0); //Initialize TtT0[LT] to the identity matrix
    for(int it=LT-1; it>=0; it--)
    {
        transfer_matrix(it, tmp_matrix);
        A_eq_B_mult_C(TtT0 + LS*LS*it, tmp_matrix, TtT0 + LS*LS*(it+1), LS);
    };
}

void LGT_QMC::rand_spatial_matrix(t_complex* out)
{
    for(int i=0; i<LS*LS; i++)
        out[i] = t_complex(unit_normal_dist(gen), unit_normal_dist(gen));
}

void LGT_QMC::JJ_func(double* out)
{
    std::copy(T0Tt + LS*LS*LT, T0Tt + LS*LS*(LT+1), tmp_matrix1); //F0 = T0Tt[NT-1];
    for(int is=0; is<LS; is++) tmp_matrix1[is*LS + is] += t_complex(1.0, 0.0); // Add the identity matrix
    inverse_matrix(tmp_matrix1, F0, LS); //F0 = (1 + T0Tt[NT-1])^-1

    current_matrix(0, J0);
    A_eq_B_mult_C(tmp_matrix1, F0, J0, LS); //tmp1 = F0*J0
    A_eq_B_mult_C(J0, tmp_matrix1, F0, LS); //J0 -> F0*J0*F0 // Current operator at time slice 0

    for(int it=0; it<LT; it++) 
    {
        current_matrix(it, Jt);
        std::copy(J0, J0 + LS*LS, tmp_matrix1); //tmp1 = J0
        A_eq_B_mult_C(tmp_matrix2, tmp_matrix1, T0Tt + LS*LS*it, LS); //tmp2 = J0*T0Tt[it]
        A_eq_B_mult_C(tmp_matrix1, tmp_matrix2, Jt, LS); //tmp1 = J0*T0Tt[it]*Jt
        A_eq_B_mult_C(tmp_matrix2, tmp_matrix1, TtT0 + LS*LS*it, LS); //tmp2 = J0*T0Tt[it]*Jt*TtT0[it]
       
        for(int is=0; is<LS; is++)
            out[it] += 2.0*tmp_matrix2[is*LS + is].real(); //Trace of the product of the matrices
    };

    current_matrix(0, J0);
    A_eq_B_mult_C(tmp_matrix1, J0, F0, LS);

    double tr1 = 0.0;
    for(int is=0; is<LS; is++)
        tr1 += 2.0*tmp_matrix1[is*LS + is].real(); //Trace of the product of the matrices

    for(int it=0; it<LT; it++) 
    {
        current_matrix(it, Jt);
        A_eq_B_mult_C(tmp_matrix1, TtT0 + LS*LS*it, T0Tt + LS*LS*it, LS);
        for(int is=0; is<LS; is++) tmp_matrix1[is*LS + is] += t_complex(1.0, 0.0); // Add the identity matrix
        inverse_matrix(tmp_matrix1, Ft, LS); //Ft = (1 + T0Tt[NT-1])^-1
        A_eq_B_mult_C(tmp_matrix1, Jt, Ft, LS);

        double tr2 = 0.0;
        for(int is=0; is<LS; is++)
            tr2 += 2.0*tmp_matrix1[is*LS + is].real(); //Trace of the product of the matrices
        
        out[it] += tr1*tr2; //Disconnected correlator
    };
}

void LGT_QMC::d_matrix(t_complex* res)
{
    std::fill(res, res + LS*LS, t_complex(0.0, 0.0));
    res[0*LS + 0] = t_complex(+1.0, 0.0);
    res[1*LS + 1] = t_complex(-1.0, 0.0);
}

void LGT_QMC::DD_func(double* out)
{
    //std::cout << "Entering DD_func" << std::endl << std::flush;
    std::copy(T0Tt + LS*LS*LT, T0Tt + LS*LS*(LT+1), tmp_matrix1); //F0 = T0Tt[NT-1];
    for(int is=0; is<LS; is++) tmp_matrix1[is*LS + is] += t_complex(1.0, 0.0); // Add the identity matrix
    inverse_matrix(tmp_matrix1, F0, LS); //F0 = (1 + T0Tt[NT-1])^-1
    //std::cout << "Inverted F0 matrix for DD_func" << std::endl << std::flush;

    d_matrix(J0);
    A_eq_B_mult_C(tmp_matrix1, F0, J0, LS); //tmp1 = F0*J0
    A_eq_B_mult_C(J0, tmp_matrix1, F0, LS); //J0 -> F0*J0*F0 // Current operator at time slice 0

    for(int it=0; it<LT; it++) 
    {
        d_matrix(Jt);
        std::copy(J0, J0 + LS*LS, tmp_matrix1); //tmp1 = J0
        A_eq_B_mult_C(tmp_matrix2, tmp_matrix1, T0Tt + LS*LS*it, LS); //tmp2 = J0*T0Tt[it]
        A_eq_B_mult_C(tmp_matrix1, tmp_matrix2, Jt, LS); //tmp1 = J0*T0Tt[it]*Jt
        A_eq_B_mult_C(tmp_matrix2, tmp_matrix1, TtT0 + LS*LS*it, LS); //tmp2 = J0*T0Tt[it]*Jt*TtT0[it]
       
        for(int is=0; is<LS; is++)
            out[it] += 2.0*tmp_matrix2[is*LS + is].real(); //Trace of the product of the matrices
    };

    d_matrix(J0);
    A_eq_B_mult_C(tmp_matrix1, J0, F0, LS);

    double tr1 = 0.0;
    for(int is=0; is<LS; is++)
        tr1 += 2.0*tmp_matrix1[is*LS + is].imag(); //Trace of the product of the matrices

    for(int it=0; it<LT; it++) 
    {
        d_matrix(Jt);
        A_eq_B_mult_C(tmp_matrix1, TtT0 + LS*LS*it, T0Tt + LS*LS*it, LS);
        for(int is=0; is<LS; is++) tmp_matrix1[is*LS + is] += t_complex(1.0, 0.0); // Add the identity matrix
        //std::cout << "Inverting matrix for time slice " << it << " in DD_func" << std::flush;
        inverse_matrix(tmp_matrix1, Ft, LS); //Ft = (1 + T0Tt[NT-1])^-1
        //std::cout << "... done!" << std::endl << std::flush;
        A_eq_B_mult_C(tmp_matrix1, Jt, Ft, LS);

        double tr2 = 0.0;
        for(int is=0; is<LS; is++)
            tr2 += 2.0*tmp_matrix1[is*LS + is].imag(); //Trace of the product of the matrices
        
        out[it] -= tr1*tr2; //Disconnected correlator
    };
}

void LGT_QMC::static_gauge()
{
    double* phi_avg = new double[LS]; std::fill(phi_avg, phi_avg + LS, 0.0); //Averages of phi1 and phi2
    double* omega   = new double[LS]; std::fill(  omega,   omega + LS, 0.0); //Gauge transformation function

    for(int is=0; is<LS; is++)
    {
        for(int it=0; it<LT; it++) phi_avg[is] += phi[LS*it + is];
        phi_avg[is] /= (double)LT;
    };

    for(int it=0; it<LT; it++) 
    {
        for(int is=0; is<LS; is++) theta[LS*it + is] += omega[(is+1)%LS] - omega[is];
        for(int is=0; is<LS; is++) omega[is]         += phi_avg[is] - phi[LS*it + is];
        for(int is=0; is<LS; is++) phi[  LS*it + is]  = phi_avg[is];
    };

    delete [] phi_avg; delete [] omega;
}

void LGT_QMC::diagonalize_single_particle_hamiltonian()
{
    fermionic_hamiltonian(0, evecs);
    int info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'U', LS, evecs, LS, evals);
    if(info != 0) { std::cerr << ansi::red << "LAPACKE_zheev failed, info = " << ansi::yellow << info << ansi::reset << std::endl; return; };
    transpose(evecs, LS);
    sort_eigensystem(evals, evecs, LS, LS);
    //TODO: check ortogonality and normalization of evecs here
}

void LGT_QMC::OO_spf(t_complex* O, double* hist)
{
    for(int m=0; m<LS; m++)
    {
        psi_eq_A_mult_chi(tmp_matrix1, O, evecs + LS*m, LS); //tmp1 = O |m>, here the first entries of tmp_matrix1 are used as a vector
        for(int n=m+1; n<LS; n++)
        {
            double w = evals[n] - evals[m];
            int ibin = (int)std::floor(w/dw);
            if(ibin<0 || ibin>=NW) continue;
            double Onm2 = std::norm(scalar_prod(evecs + LS*n, tmp_matrix1, LS)); //|<n|O|m>|^2

            double phi_avg = 0.0;
            for(int ist=0; ist<LS*LT; ist++) phi_avg += phi[ist];
            phi_avg /= (double)(LS);
            t_complex exp_iphi = std::exp(t_complex(0.0, phi_avg)); 
            
            t_complex spf_contrib = 2.0*Onm2*exp_iphi*(exp(-beta*evals[m]) + exp(-beta*evals[n]));
            spf_contrib /= (1.0 + exp(-beta*evals[m])*exp_iphi);
            spf_contrib /= (1.0 + exp(-beta*evals[n])*exp_iphi);
            hist[ibin] += spf_contrib.real(); //Only the real part contributes
        };
    };
}

void LGT_QMC::static_projection()
{
    static_gauge();

    for(int is=0; is<LS; is++)
    {
        double theta_avg = 0.0;
        for(int it=0; it<LT; it++) theta_avg += theta[LS*it + is];
        theta_avg /= (double)LT;
        for(int it=0; it<LT; it++) theta[LS*it + is] = theta_avg;
    };

    for(int it=0; it<LT; it++)
    {
        double phi_avg = 0.0;
        for(int is=0; is<LS; is++) phi_avg += phi[LS*it + is];
        phi_avg /= (double)LS;
        for(int is=0; is<LS; is++) phi[LS*it + is] = phi_avg;
    };
}

void LGT_QMC::real_time_DD_func(double* out)
{
    d_matrix(J0);

    double phi_avg = 0.0;
    for(int ist=0; ist<LS*LT; ist++) phi_avg += phi[ist];
    phi_avg /= (double)(LS);
    t_complex exp_iphi = std::exp(t_complex(0.0, phi_avg)); 

    for(int m=0; m<LS; m++)
    {
        psi_eq_A_mult_chi(tmp_matrix1, J0, evecs + LS*m, LS); //tmp1 = O |m>, here the first entries of tmp_matrix1 are used as a vector
        for(int n=0; n<LS; n++)
        {
            double Onm2 = std::norm(scalar_prod(evecs + LS*n, tmp_matrix1, LS)); //|<n|O|m>|^2
            
            t_complex gr_contrib_c = exp_iphi/((1.0 + exp(-beta*evals[m])*exp_iphi)*(1.0 + exp(-beta*evals[n])*exp_iphi));
            double gr_contrib = 2.0*Onm2*exp(-beta*evals[m])*gr_contrib_c.real();

            for(int it=0; it<nsteps; it++)
            {
                double t = dt*(double)it;
                out[2*it + 0] += gr_contrib*cos((evals[n] - evals[m])*t); //Real part
                out[2*it + 1] += gr_contrib*sin((evals[n] - evals[m])*t); //Imaginary part
            };

            //Disconnected part
            if(m==n)
            {
                t_complex disc_contrib_c = 1.0/(1.0 + exp(-beta*evals[m])*exp_iphi);
                double disc_contrib = 4.0*Onm2*SQR(disc_contrib_c.imag());
                for(int it=0; it<nsteps; it++)
                    out[2*it + 0] -= disc_contrib;
            };

        };
    };
}