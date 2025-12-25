#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>

#include "LGT_QMC.hpp"
#include "observables.hpp"

int main(int argc, char* argv[]) 
{
    LGT_QMC QMC(argc, argv);    

    Observable eeq(QMC.LT*QMC.LS, QMC.dmeas, &QMC, &LGT_QMC::electric_field_correlator);
    Observable ees(QMC.LT*QMC.LS, QMC.dmeas, &QMC, &LGT_QMC::electric_field_correlator);
    Observable jjq(QMC.LT,        QMC.dmeas, &QMC, &LGT_QMC::JJ_func);
    Observable jjs(QMC.LT,        QMC.dmeas, &QMC, &LGT_QMC::JJ_func);
    Observable ddq(QMC.LT,        QMC.dmeas, &QMC, &LGT_QMC::DD_func, true); //true here is a switch to calculate also the covariance matrix
    Observable dds(QMC.LT,        QMC.dmeas, &QMC, &LGT_QMC::DD_func);
    Observable  fd(1,             QMC.dmeas, &QMC, &LGT_QMC::fermion_det);
    Observable ddr(2*QMC.nsteps,  QMC.dmeas, &QMC, &LGT_QMC::real_time_DD_func);

    //Variables for spectral func histogramming
    t_complex* Jmatrix = new t_complex[QMC.LS*QMC.LS]; double* Jhist = new double[QMC.NW]; std::fill(Jhist, Jhist + QMC.NW, 0.0);
    t_complex* Dmatrix = new t_complex[QMC.LS*QMC.LS]; double* Dhist = new double[QMC.NW]; std::fill(Dhist, Dhist + QMC.NW, 0.0);
    int nhist = 0;

    double max_evec_err = 0.0, max_ortho_err = 0.0;
    t_complex* hsp = NULL;
    if(QMC.internal_checks) hsp = new t_complex[QMC.LS*QMC.LS];

    std::cout << std::fixed << std::setprecision(3);
    for(int itrial=0; itrial<QMC.NMC; itrial++) 
    {
        QMC.gauge_field_update();
        if(itrial < 10*QMC.dmeas) continue; // Skip the first 10*dmeas steps to allow the system to equilibrate
        if(itrial%10000 == 0) std::cout << ansi::save_cursor << " " << 100.0*(double)itrial/(double)QMC.NMC << "% " << ansi::restore_cursor << std::flush; // Print a dot every 100*dmeas steps
        eeq.measure();
        if(QMC.dd_correlators) ddq.measure();
        if(QMC.jj_correlators) jjq.measure();
        fd.measure();

        if(QMC.static_approximation)
        {
            QMC.backup_fields();
            QMC.static_projection();
            QMC.update_exps();

            QMC.current_matrix(0, Jmatrix);
            QMC.d_matrix(Dmatrix);

            ees.measure();
            if(QMC.jj_correlators) jjs.measure();
            if(QMC.dd_correlators) dds.measure();

            if(QMC.spectral_functions || QMC.real_time_correlators)
            {
                QMC.diagonalize_single_particle_hamiltonian();
                if(QMC.internal_checks) 
                {
                    double evec_err, ortho_err;
                    QMC.fermionic_hamiltonian(0, hsp);
                    test_eigensystem(hsp, QMC.evals, QMC.evecs, QMC.LS, &evec_err, &ortho_err);
                    max_evec_err  = std::max(max_evec_err, evec_err);
                    max_ortho_err = std::max(max_ortho_err, ortho_err);
                };
            };

            if(QMC.spectral_functions)
            {
                if(QMC.jj_correlators) QMC.OO_spf(Jmatrix, Jhist);
                if(QMC.dd_correlators) QMC.OO_spf(Dmatrix, Dhist);
                nhist++;
            };

            if(QMC.real_time_correlators && QMC.dd_correlators) ddr.measure();

            QMC.restore_fields();
        };
    };
    std::cout << ansi::reset << std::endl;
    std::cout << "Acceptance ratio: " << (double)QMC.accepted/(double)QMC.NMC << std::endl;
    std::cout << "Max. eigenvector error: "   << std::scientific << max_evec_err << std::endl;
    std::cout << "Max. orthogonality error: " << std::scientific << max_ortho_err << std::endl;

    double *afd, *efd;
    fd.statistical_averages(afd, efd); // Just to print the fermionic determinant at the end of the simulation
    std::cout << "Average fermionic determinant: " << afd[0] << " +/- " << efd[0] << std::endl;

    double *aeeq = NULL; double *eeeq = NULL;
    eeq.statistical_averages(aeeq, eeeq);

    double *aees = NULL; double *eees = NULL;
    ees.statistical_averages(aees, eees);

    double *ajjq = NULL; double *ejjq = NULL;
    jjq.statistical_averages(ajjq, ejjq);

    double* ajjs = NULL; double *ejjs = NULL;
    jjs.statistical_averages(ajjs, ejjs);

    double *addq = NULL; double *eddq = NULL; double* cddq = NULL;
    ddq.statistical_averages(addq, eddq, cddq);

    double* adds = NULL; double *edds = NULL;
    dds.statistical_averages(adds, edds);

    double* addr = NULL; double* eddr = NULL;
    ddr.statistical_averages(addr, eddr);

    std::cout << "QMC.suffix0 = " << QMC.suffix0 << std::endl;

    char fname_OO[512];
    std::sprintf(fname_OO, "./data/OO_QMC_%s.dat", QMC.suffix0);
    std::cout << "Output file for operator-operator correlators: " << fname_OO << std::endl;
    std::ofstream outfile_OO(fname_OO);
    if (!outfile_OO.is_open()) {std::cerr << "Failed to open " << fname_OO << " for writing." << std::endl; return 1;};

    std::cout << "Electric field correlator <E(t)E(0)> and current-current correlator <J(t) J(0)> (E2 should be "<< QMC.E2_pure_lgt() << " in pure LGT):" << std::endl;
    for(int it=0; it<QMC.LT; it++)
    {
        double tau = it*QMC.beta/(double)QMC.LT;
        std::cout  << std::fixed << std::setprecision(2) << tau << "\t";
        outfile_OO << std::fixed << std::setprecision(8) << tau << "\t";

        outfile_OO << ajjq[it] << "\t" << ejjq[it] << "\t";
        outfile_OO << ajjs[it] << "\t" << ejjs[it] << "\t";
        std::cout << ajjq[it] << " +/- " << ejjq[it] << "\t|";
        std::cout << ajjs[it] << " +/- " << ejjs[it] << "\t|";
        for(int is=0; is<QMC.LS; is++)
            outfile_OO << aeeq[QMC.LS*it + is] << "\t" << eeeq[QMC.LS*it + is] << "\t";
        for(int is=0; is<QMC.LS; is++)
            outfile_OO << aees[QMC.LS*it + is] << "\t" << eees[QMC.LS*it + is] << "\t";
        outfile_OO << addq[it] << "\t" << eddq[it] << "\t";
        outfile_OO << adds[it] << "\t" << edds[it] << "\n";

        int is0 = 0;
        std::cout << "(" << aeeq[QMC.LS*it + is0] << " +/- " << eeeq[QMC.LS*it + is0] << ")\t";
        std::cout << "(" << aees[QMC.LS*it + is0] << " +/- " << eees[QMC.LS*it + is0] << ")\t";

        std::cout << std::endl;
    };
    std::cout << std::endl;
    outfile_OO.close();

    char fname_cov[512];
    std::sprintf(fname_cov, "./data/DD_QMC_%s.cov", QMC.suffix0);
    std::cout << "Output file for DD operator covariance matrix: " << fname_cov << std::endl;
    std::ofstream outfile_cov(fname_cov, std::ios::out | std::ios::binary);
    if (!outfile_cov.is_open()) {std::cerr << "Failed to open " << fname_cov << " for writing." << std::endl; return 1;};
    outfile_cov.write((char*)(cddq), QMC.LT*QMC.LT*sizeof(double));
    outfile_cov.close();

    if(QMC.spectral_functions)
    {
        // Output spectral functions to a binary file
        char fname_spf[512];
        std::sprintf(fname_spf, "./data/SPF_QMC_%s.bin", QMC.suffix0);
        std::ofstream outfile_spf(fname_spf, std::ios::out | std::ios::binary);
        if (!outfile_spf.is_open()) { std::cerr << ansi::red << "Failed to open " << fname_spf << " for writing." << ansi::reset << std::endl; return 1; }

        double norm_factor = (double)(nhist)*QMC.dw;
        rescale(Dhist, 1.0/norm_factor, QMC.NW);
        rescale(Jhist, 1.0/norm_factor, QMC.NW);

        outfile_spf.write((char*)(Dhist), QMC.NW*sizeof(double));
        outfile_spf.write((char*)(Jhist), QMC.NW*sizeof(double));

        outfile_spf.close();
        std::cout << "Output file for spf_hist (binary): " << fname_spf << ", dw = " << QMC.dw << std::endl;
    };

    if(QMC.real_time_correlators)
    {
        char fname_rt[512];
        std::sprintf(fname_rt, "./data/GR_QMC_%s.dat", QMC.suffix0);
        std::cout << "Output file for real-time DD correlator: " << fname_rt << std::endl;
        std::ofstream outfile_rt(fname_rt);
        if (!outfile_rt.is_open()) {std::cerr << "Failed to open " << fname_rt << " for writing." << std::endl; return 1;};

        for(int it=0; it<QMC.nsteps; it++)
        {
            double t = it*QMC.dt;
            outfile_rt << std::fixed << std::setprecision(8) << t << "\t";
            outfile_rt << addr[2*it+0] << "\t" << eddr[2*it+0] << "\t";
            outfile_rt << addr[2*it+1] << "\t" << eddr[2*it+1] << "\n";
            std::cout << std::fixed << std::setprecision(4) << t << "\t";
            std::cout << addr[2*it+0] << " +/- " << eddr[2*it+0] << "\t";
            std::cout << addr[2*it+1] << " +/- " << eddr[2*it+1] << std::endl;
        };
        outfile_rt.close();
    };

    std::cout << std::endl;
    return 0;
}