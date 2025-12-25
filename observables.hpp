#ifndef OBSERVABLES_HPP
#define OBSERVABLES_HPP

#include <cstddef>
#include <algorithm>
#include <cmath>
#include "LGT_QMC.hpp"

class Observable
{
    // Make Observable a friend of LGT_QMC
    friend class LGT_QMC;

    public:
        Observable(int NO, int dmeas, LGT_QMC* lgt_qmc, void (LGT_QMC::*func)(double*), bool calculate_covariance=false) : NO_(NO), dmeas_(dmeas), lgt_qmc_(lgt_qmc), func_(func) 
        {
            O_loc = new double[NO_];
            aO = new double[NO_];
            aO2 = new double[NO_];
            if(calculate_covariance) cO = new double[NO_*NO_];
            std::fill(O_loc, O_loc + NO_, 0.0);
            std::fill(aO, aO + NO_, 0.0);
            std::fill(aO2, aO2 + NO_, 0.0);
            if(cO != NULL) std::fill(cO, cO + NO_*NO_, 0.0);
            nmeas_loc = 0; // Initialize the local measurement counter
            nmeas = 0; // Initialize the global measurement counter
        }

        void measure()
        {
            (lgt_qmc_->*func_)(O_loc);
            nmeas_loc ++;
            if(nmeas_loc==dmeas_)
            { 
                for(int i=0; i<NO_; i++) 
                    O_loc[i] /= (double)nmeas_loc; // Normalize the local observabl
                for(int i=0; i<NO_; i++) 
                {
                    aO[i]  += O_loc[i];
                    aO2[i] += O_loc[i]*O_loc[i];
                    if(cO != NULL)
                        for(int j=0; j<NO_; j++)
                            cO[i*NO_ + j] += O_loc[i]*O_loc[j];
                };

                std::fill(O_loc, O_loc + NO_, 0.0); // Reset the local observable

                //Reset the local measurement counter and increase the measurement counter
                nmeas_loc = 0;
                nmeas++;
            };
        };

        void statistical_averages(double*& mean, double*& err)
        {
            // aO contains the sum of observables, aO2 contains the sum of squares
            // Compute averages and statistical uncertainty (standard error of the mean)
            mean = aO; err = aO2;
            if (nmeas <1 ){std::fill(aO, aO + NO_, 0.0); std::fill(aO2, aO2 + NO_, 0.0); return;};
            for(int i = 0; i < NO_; i++) 
            {
                aO[i]  /= (double)nmeas;
                aO2[i] /= (double)nmeas;
                if(nmeas> 1)
                    aO2[i] = std::sqrt((aO2[i] - aO[i]*aO[i])/(double)(nmeas - 1)); // Second moment of the observable
            };
            if(nmeas == 1) err = NULL;
        };

        void statistical_averages(double*& mean, double*& err, double*& covariance)
        {
            // aO contains the sum of observables, aO2 contains the sum of squares
            // Compute averages and statistical uncertainty (standard error of the mean)
            mean = aO; err = aO2; covariance = cO;
            if (nmeas <1 ){std::fill(aO, aO + NO_, 0.0); std::fill(aO2, aO2 + NO_, 0.0); std::fill(cO, cO + NO_*NO_, 0.0); return;};
            for(int i = 0; i < NO_; i++) 
            {
                aO[i]  /= (double)nmeas;
                aO2[i] /= (double)nmeas;
                if(nmeas> 1)
                    aO2[i] = std::sqrt((aO2[i] - aO[i]*aO[i])/(double)(nmeas - 1)); // Second moment of the observable
            };
            if(nmeas == 1) { err = NULL; covariance = NULL; return; };
            if(nmeas>1 && cO != NULL)
            {
                //Compute covariance matrix
                for(int i=0; i<NO_; i++)
                    for(int j=0; j<NO_; j++)
                        cO[i*NO_ + j] = (cO[i*NO_ + j]/(double)nmeas - aO[i]*aO[j])/(double)(nmeas - 1);
            };
        };

        ~Observable()
        {
            delete[] O_loc;
            delete[] aO;
            delete[] aO2;
            if(cO != NULL) delete[] cO;
        };

        int nmeas = 0;
        int nmeas_loc = 0; // Number of measurements
            
        int NO_;
        int dmeas_;
        LGT_QMC* lgt_qmc_;
        void (LGT_QMC::*func_)(double*);
        double *O_loc = NULL;
        double *aO  = NULL; //Expectation values
        double *aO2 = NULL; //Expectation values of the squares
        double* cO  = NULL; //Covariance matrix storage
};

#endif // OBSERVABLES_HPP