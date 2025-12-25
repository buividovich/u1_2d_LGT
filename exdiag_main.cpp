#include <iostream>
#include <fstream>
#include <iomanip>
#include "LGT_ExDiag.hpp"

int main(int argc, char* argv[])
{
    LGT_ExDiag exdiag(argc, argv);
    if(exdiag.internal_checks)
    {
        exdiag.test_fermionic_indices();
        exdiag.test_hopping();
    };

    exdiag.diagonalize_hamiltonian();

    if(exdiag.internal_checks)
    {
        test_hermiticity(exdiag.HamiltonianMatrix, exdiag.NS);
        test_eigensystem(exdiag.HamiltonianMatrix, exdiag.E, exdiag.psi, exdiag.NS);
        exdiag.test_current_hamiltonian_algebra();
    };

    double* E0Matrix = exdiag.construct_efield_matrix(0);
    test_hermiticity(E0Matrix, exdiag.NS);
    double* EE = exdiag.EuclideanCorrelator(E0Matrix);

    //Placeholders for observables
    double* JJ = NULL; double* DD = NULL; double* JMatrix = NULL; double* DMatrix = NULL;

    if(exdiag.jj_correlators)
    {
        JMatrix = exdiag.construct_current_matrix();
        if(exdiag.internal_checks) test_hermiticity(JMatrix, exdiag.NS, -1);
        JJ = exdiag.EuclideanCorrelator(JMatrix);
    }
    else { JJ = new double[exdiag.LT]; std::fill(JJ, JJ + exdiag.LT, 0.0); };

    if(exdiag.dd_correlators)
    {
        DMatrix = exdiag.construct_d_matrix();
        if(exdiag.internal_checks) test_hermiticity(DMatrix, exdiag.NS);
        DD = exdiag.EuclideanCorrelator(DMatrix);
    } else { DD = new double[exdiag.LT]; std::fill(DD, DD + exdiag.LT, 0.0); };

    char fname_OO[256];
    snprintf(fname_OO, 256, "./data/OO_ExDiag_%s_NM%i.dat", exdiag.suffix0, exdiag.NM);
    
    std::cout << "Output file for Euclidean-time correlators: " << fname_OO << std::endl;
    std::ofstream outfile_OO(fname_OO);
    if (!outfile_OO.is_open()) {std::cerr << ansi::red << "Failed to open " << ansi::cyan << fname_OO << ansi::red << " for writing." << std::endl;  return 1;};

    std::cout  << "Euclidean-time correlators:" << std::endl;
    std::cout  << std::fixed << std::setprecision(4);
    outfile_OO << std::fixed << std::setprecision(4);
    for (int it = 0; it < exdiag.LT; it++)
    {
        double tau = it*exdiag.beta/(double)exdiag.LT;
        std::cout  << tau << "\t" << JJ[it] << "\t" << EE[it] << "\t" << DD[it] << std::endl;
        outfile_OO << tau << "\t" << JJ[it] << "\t" << EE[it] << "\t" << DD[it] << std::endl;
    }
    outfile_OO.close();

    if(exdiag.spectral_functions)
    {
        SpectralFunction* SPF_JJ = NULL;
        SpectralFunction* SPF_DD = NULL;
        if(exdiag.jj_correlators) SPF_JJ = exdiag.spectral_function(JMatrix);
        if(exdiag.dd_correlators) SPF_DD = exdiag.spectral_function(DMatrix);
        char fname_SPF[256];
        snprintf(fname_SPF, 256, "./data/SPF_ExDiag_%s_NM%i.dat", exdiag.suffix0, exdiag.NM);
        std::cout << "Output file for spectral function: " << fname_SPF << std::endl;
        std::ofstream outfile_SPF(fname_SPF);
        if (!outfile_SPF.is_open()) {std::cerr << ansi::red << "Failed to open " << ansi::cyan << fname_SPF << ansi::red << " for writing." << std::endl;  return 1;};

        for (int iw = 0; iw < exdiag.NW; iw++)
        {
            outfile_SPF << std::fixed << std::setprecision(6) << (iw + 0.5)*exdiag.dw << "\t";
            if(SPF_JJ!=NULL) outfile_SPF  << SPF_JJ->rho[iw]; else outfile_SPF  << 0.0;
            outfile_SPF << "\t";
            if(SPF_DD!=NULL) outfile_SPF  << SPF_DD->rho[iw]; else outfile_SPF  << 0.0;
            outfile_SPF << std::endl;
        };
        outfile_SPF.close();

        if(SPF_JJ!=NULL) delete SPF_JJ;
        if(SPF_DD!=NULL) delete SPF_DD;
    };

    if(exdiag.real_time_correlators)
    {
        t_complex* GR_JJ = NULL;
        t_complex* GR_DD = NULL;
        if(exdiag.jj_correlators) GR_JJ = exdiag.MinkowskiCorrelator(JMatrix);
        if(exdiag.dd_correlators) GR_DD = exdiag.MinkowskiCorrelator(DMatrix);
        char fname_GR[256];
        snprintf(fname_GR, 256, "./data/GR_ExDiag_%s_NM%i.dat", exdiag.suffix0, exdiag.NM);
        std::cout << "Output file for Minkowski-time correlators: " << fname_GR << std::endl;
        std::ofstream outfile_GR(fname_GR);
        if (!outfile_GR.is_open()) {std::cerr << ansi::red << "Failed to open " << ansi::cyan << fname_GR << ansi::red << " for writing." << std::endl;  return 1;};

        for (int it = 0; it < exdiag.nsteps; it++)
        {
            double t = it*exdiag.dt;
            outfile_GR << std::fixed << std::setprecision(4) << t << "\t";
            if(GR_JJ!=NULL) outfile_GR  << GR_JJ[it].real() << "\t" << GR_JJ[it].imag(); else outfile_GR  << 0.0 << "\t" << 0.0;
            outfile_GR << "\t";
            if(GR_DD!=NULL) outfile_GR  << GR_DD[it].real() << "\t" << GR_DD[it].imag(); else outfile_GR  << 0.0 << "\t" << 0.0;
            outfile_GR << std::endl;
        };
        outfile_GR.close();
        if(GR_JJ!=NULL) delete [] GR_JJ;
        if(GR_DD!=NULL) delete [] GR_DD;     
    }

    if(JMatrix!=NULL) delete [] JMatrix; delete [] JJ;
    if(DMatrix!=NULL) delete [] DMatrix; delete [] DD;
    delete [] E0Matrix; delete [] EE;

    return 0;
}