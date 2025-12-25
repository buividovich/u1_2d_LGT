#ifndef LGT_EXDIAG_HPP
#define LGT_EXDIAG_HPP

#include <boost/program_options.hpp>
#include <iostream>
#include <algorithm>
#include <complex>
#include "linalg.hpp"
#include "LGTDescriptor.hpp"
#include "timing.hpp"

namespace po = boost::program_options;

class FermionicState
{
    public:
    FermionicState(int LS) : LS(LS)
    {
        na = new int[LS];
        nb = new int[LS];
    };
    ~FermionicState()
    {
        delete [] na;
        delete [] nb;
    };
    int state_index() const
    {
        int f = 1; int res = 0;
        for(int is=0; is<LS; is++){ res += f*na[is]; f *= 2; };
        for(int is=0; is<LS; is++){ res += f*nb[is]; f *= 2; };
        return res;
    };
    int total_charge() const
    {
        int res = 0;
        for(int is=0; is<LS; is++) res += (na[is] - nb[is]);
        return res;
    };
    int LS;
    int* na;
    int* nb;
};

class LGTState : public FermionicState
{
    public:
    int* nf; //Flux through each of the links
    LGTState(int LS) : FermionicState(LS) { nf = new int[LS]; std::fill(nf, nf+LS, 0); };
    ~LGTState() { delete [] nf; };
    void flux_dress();
    bool is_valid_state();
    LGTState& operator=(const LGTState& other);
    LGTState(const LGTState& other);
    friend std::ostream& operator<<(std::ostream& Stream, const LGTState& LS);
};

class SpectralFunction
{
    public:
    double* rho = nullptr; //Spectral function values at the frequency points
    int NW = 0; //Number of frequency points
    double dw = 0.0; //Frequency resolution
    double wmax = 0.0; //Maximum frequency
    ~SpectralFunction(){ if(rho!=nullptr) delete [] rho; };
};

void IndexToFermionicState(int index, FermionicState* S);

class LGT_ExDiag : public LGTDescriptor 
{
public:
    LGT_ExDiag(int argc, char* argv[]);
    ~LGT_ExDiag();
    int NM           = 0; // Number of matrix elements (example default)
    int NFS_ungauged = 0; //Number of possible fermion states without gauging
    int NFS          = 0; //Number of charge-zero fermion states
    int NS           = 0; //Total number of gauge-invariant states
    po::options_description exdiag_params{"Exact Diagonalization Parameters"};
    //"Dictionaries" for fermionic states
    int* ParentStateDirectory = NULL;
    LGTState** ParentStates   = NULL;
    // Fermion state enumeration and indexing
    void IndexToLGTState(int index, LGTState* S);
    int LGTStateToIndex(const LGTState* S);
    //Hopping term in the Hamiltonian
    double* HamiltonianMatrix = NULL; //Hamiltonian matrix now real-valued
    double* E = NULL; //Energy eigenvalues of the Hamiltonian
    double* psi = NULL; //Eigenvectors of the Hamiltonian
    int hopping(const LGTState* in, int is, int dir, int charge, LGTState* out) const; //Moves a fermion from site number "is" to "(is + dir)%LS" in the state "in", puts the result in "out". Returns 0 if not possible, +1 or -1 depending on fermionic sign otherwise.
    void construct_hamiltonian_matrix(); //Builds the Hamiltonian matrix in the basis of LGT states
    void diagonalize_hamiltonian(); //Diagonalizes the Hamiltonian matrix and prints the eigenvalues
    void print_eigenvalues() const;
    //Operators of physical observables
    double* construct_current_matrix();
    double* construct_d_matrix();
    double* construct_efield_matrix(int is=0); //Constructs the electric field operator matrix for link "is"
    double* EuclideanCorrelator(double* O); //Computes the Euclidean correlator of the operator O (given as a matrix in the LGT basis)
    t_complex* MinkowskiCorrelator(double* O); //Computes the Minkowski correlator of the operator O (given as a matrix in the LGT basis)
    SpectralFunction* spectral_function(double* O); //Computes the spectral function of the operator O (given as a matrix in the LGT basis). On output, NW is set to the number of frequency points, dw is the frequency resolution, wmax is the maximum frequency.  
    //Thermodynamics
    double partition_function(); //Partition function
    //Test routines
    void test_fermionic_indices();
    void print_parent_states();
    void print_parameters();
    void test_hopping();
    void test_current_hamiltonian_algebra();
};

std::ostream& operator<<(std::ostream& Stream, const FermionicState& FS); 

#endif // LGT_EXDIAG_HPP
