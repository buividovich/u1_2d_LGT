#include "LGT_ExDiag.hpp"

void LGTState::flux_dress()
{
    // Simple example: set flux to be the difference in occupation numbers
    std::fill(nf, nf+LS, 0);
    for(int is=1; is<LS; is++)
        nf[is] = nf[is-1] + (na[is] - nb[is]);
}

bool LGTState::is_valid_state()
{
    bool res = true;
    for(int is=0; is<LS; is++)
    {
        int isb = (LS + is - 1) % LS;
        res = res && (nf[is] - nf[isb] == na[is] - nb[is]);
    };
    return res;
}

LGTState& LGTState::operator=(const LGTState& other)
{
    if(this != &other)
    { 
        this->LS = other.LS;
        this->na = new int[other.LS];
        this->nb = new int[other.LS];
        this->nf = new int[other.LS];
        std::copy(other.na, other.na+other.LS, this->na);
        std::copy(other.nb, other.nb+other.LS, this->nb);
        std::copy(other.nf, other.nf+other.LS, this->nf);
    };
    return *this;
};

LGTState::LGTState(const LGTState& other) : FermionicState(other.LS) 
{
    nf = new int[other.LS];
    std::copy(other.na, other.na + other.LS, this->na);
    std::copy(other.nb, other.nb + other.LS, this->nb);
    std::copy(other.nf, other.nf + other.LS, this->nf);
}

LGT_ExDiag::LGT_ExDiag(int argc, char* argv[]) : LGTDescriptor(argc, argv)
{
    // Initialize Exact Diagonalization parameters
    exdiag_params.add_options()
        ("NM", po::value<int>(&NM)->default_value(10), "Number of matrix elements");

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(exdiag_params).allow_unregistered().run(), vm);
    po::notify(vm);

    // Check if help is needed
    if(vm.count("help")) { std::cout << exdiag_params << std::endl; }

    NFS_ungauged = 1;
    for(int is=0; is<LS; is++) NFS_ungauged *= 4;

    ParentStateDirectory = new int[NFS_ungauged];
    ParentStates         = new LGTState*[NFS_ungauged];
    std::fill(ParentStateDirectory, ParentStateDirectory+NFS_ungauged, -1);
    std::fill(ParentStates, ParentStates+NFS_ungauged, nullptr);

    NFS = 0;
    LGTState tmp(LS);
    for(int idx=0; idx<NFS_ungauged; idx++)
    {
        IndexToFermionicState(idx, &tmp);
        if(tmp.total_charge()==0)
        {
            tmp.flux_dress();
            ParentStateDirectory[idx] = NFS;
            ParentStates[NFS]         = new LGTState(tmp);
            NFS++;
        };
    };

    NS = (2*NM + 1)*NFS;

    //print_parent_states();
    print_parameters();
}

void LGT_ExDiag::print_parameters()
{
    std::cout << ansi::cyan << "Exact Diagonalization Parameters:" << ansi::reset << std::endl;
    std::cout << ansi::green << "\tNM           = " << ansi::yellow << NM           << ansi::reset << std::endl;
    std::cout << ansi::green << "\tNFS_ungauged = " << ansi::yellow << NFS_ungauged << ansi::reset << std::endl;
    std::cout << ansi::green << "\tNFS          = " << ansi::yellow << NFS          << ansi::reset << std::endl;
    std::cout << ansi::green << "\tNS           = " << ansi::yellow << NS           << ansi::reset << std::endl;
    std::cout << std::endl;
}

void LGT_ExDiag::print_parent_states()
{
    std::cout << ansi::cyan << "Sparse array of parent states:" << ansi::reset << std::endl;
    for(int idx=0; idx<NFS_ungauged; idx++)
        if(ParentStateDirectory[idx]!=-1)
        {
            std::cout << ansi::white << "\t" << idx << " = " << ansi::magenta << *(ParentStates[ParentStateDirectory[idx]]);
            if(ParentStates[ParentStateDirectory[idx]]->is_valid_state())
                std::cout << ansi::green << " V " << ansi::reset << std::endl;
            else
                std::cout << ansi::red   << " X " << ansi::reset << std::endl;
        }
        else
            std::cout << ansi::blue << "\t" << idx << " = NULL" << std::endl;
    std::cout << ansi::cyan << "Compact array of parent states:" << ansi::reset << std::endl;
    for(int idx=0; idx<NFS; idx++)
        std::cout << ansi::white << "\t" << idx << " = " << ansi::magenta << *(ParentStates[idx]) << std::endl;
    std::cout << std::endl;    
}

std::ostream& operator<<(std::ostream& Stream, const LGTState& LS)
{
    Stream << static_cast<const FermionicState&>(LS);
    Stream << " nf = (";
    for(int is=0; is<LS.LS-1; is++) Stream << LS.nf[is] << ", ";
    Stream << LS.nf[LS.LS-1] << ") ";
    return Stream;
}

LGT_ExDiag::~LGT_ExDiag() 
{
    // Destructor implementation
    delete [] ParentStateDirectory;

    for(int idx=0; idx<NFS; idx++) delete ParentStates[idx];
    delete [] ParentStates;

    if(HamiltonianMatrix != NULL) delete [] HamiltonianMatrix;
    if(E != NULL)  delete [] E;
    if(psi != NULL) delete [] psi;
}
// Add more method implementations as needed

void LGT_ExDiag::test_fermionic_indices()
{
    LGTState* tmp = new LGTState(LS);
    bool res = true;
    for(int idx=0; idx<NFS_ungauged; idx++)
    {
        IndexToFermionicState(idx, tmp);
        bool check = (tmp->state_index()==idx);
        res = res && check;
        if(!check)
            std::cout << "Error with the state " << *tmp << std::endl;
    };
    if(res)
        std::cout << ansi::green << "Test of fermionic state indices passed!" << ansi::reset << std::endl;
    else
        std::cout << ansi::red   << "Test of fermionic state indices failed!" << ansi::reset << std::endl;  

    res = true;
    for(int idx=0; idx<NS; idx++)
    {
        IndexToLGTState(idx, tmp);
        int idx_check = LGTStateToIndex(tmp);
        bool check = (idx_check==idx);
        res = res && check;
        if(!check)
           std::cout << ansi::red << "Error with the state " << *tmp << ": got " << idx_check << " instead of " << idx << ansi::reset << std::endl;
    };

    if(res)
        std::cout << ansi::green << "Test of multi-flux LGT states indices passed!" << ansi::reset << std::endl;
    else
        std::cout << ansi::red   << "Test of multi-flux LGT states indices failed!" << ansi::reset << std::endl;  
    delete tmp;
}

void LGT_ExDiag::IndexToLGTState(int index, LGTState* S)
{
    *S = *ParentStates[(index % NFS)];
    for(int is=0; is<LS; is++) S->nf[is] += ((index / NFS) - NM);
}

int LGT_ExDiag::LGTStateToIndex(const LGTState* S)
{
    int idx = ParentStateDirectory[S->state_index()];
    if(idx == -1)
    {
        std::cerr << ansi::red << "Error: Attempting to index a non-parent state!" << ansi::reset << std::endl;
        return -1;
    };
    return (S->nf[0] + NM)*NFS + idx;
}

void IndexToFermionicState(int index, FermionicState* S)
{
    int f = 1;
    for(int is=0; is<S->LS; is++){S->na[is] = (index/f)%2; f*=2;};
    for(int is=0; is<S->LS; is++){S->nb[is] = (index/f)%2; f*=2;};
}

std::ostream& operator<<(std::ostream& Stream, const FermionicState& FS) 
{
  Stream << " na = (";
  for(int is=0; is<FS.LS-1; is++) Stream << FS.na[is] << ", ";
  Stream << FS.na[FS.LS-1] << "), nb = (";
  for(int is=0; is<FS.LS-1; is++) Stream << FS.nb[is] << ", ";
  Stream << FS.nb[FS.LS-1] << ") ";
  return Stream;
}

int LGT_ExDiag::hopping(const LGTState* in, int is, int dir, int charge, LGTState* out) const
{
    int is1 = (LS + dir + is)%LS;
    *out = *in; //Copy the content of "in" to "out"
    int* nq = (charge==1 ? out->na : out->nb);
    if(nq[is]==1 && nq[is1]==0)
    {
        nq[is]  = 0; nq[is1] = 1;
        int ia = std::min(is, is1);  int ib = std::max(is, is1); int res = 0;
        for(int i=ia+1; i<ib; i++) res += nq[i];
        int ifs = (dir>0 ? is : is1);
        out->nf[ifs] -= charge*dir;
        if(std::abs(out->nf[0]) > NM) return 0; //Flux too large, hopping not allowed
        return (res%2==0 ? 1 : -1); //No fermionic sign change for 1D nearest-neighbor hopping
    }
    else
        return 0;
}

void LGT_ExDiag::test_hopping()
{
    LGTState* from = new LGTState(LS);
    LGTState* to   = new LGTState(LS);
    bool res = true;
    for(int idx=0; idx<NS; idx++)
    {
        IndexToLGTState(idx, from);
        if(!from->is_valid_state()){ std::cout << ansi::red << "Error: State " << *from << " is not a valid state!" << ansi::reset << std::endl; };
        for(int is=0; is<LS; is++)
            for(int dir=-1; dir<=1; dir+=2)
                for(int charge=-1; charge<=1; charge+=2)
                {
                    int hop_sign = hopping(from, is, dir, charge, to);
                    if(hop_sign == 0) continue; //Hopping not possible

                    if(!to->is_valid_state())
                    {   
                        res = false; 
                        std::cout << ansi::red << "Error - an invalid state is produced! : " << ansi::reset << std::endl;
                        std::cout << ansi::cyan << "From state: " << ansi::yellow << *from << ansi::reset << std::endl;
                        std::cout << ansi::cyan << "(is, dir, charge) = (" << ansi::yellow << is << ", " << dir << ", " << charge << ansi::reset << ")" << std::endl;
                        std::cout << ansi::cyan << "To state:   " << ansi::yellow << *to   << ansi::reset << std::endl;
                        std::cout << ansi::cyan << "Hopping sign: " << ansi::yellow << hop_sign << ansi::reset << std::endl;
                        continue;
                    };

                    int idx1 = LGTStateToIndex(to);
                    if(idx1<0 || idx1>=NS)
                    {
                        res = false; 
                        std::cout << ansi::red << "Error - index of a \"to\" state is invalid: " << ansi::reset << std::endl;
                        std::cout << ansi::cyan << "From state: " << ansi::yellow << *from << ansi::reset << std::endl;
                        std::cout << ansi::cyan << "(is, dir, charge) = (" << ansi::yellow << is << ", " << dir << ", " << charge << ansi::reset << ")" << std::endl;
                        std::cout << ansi::cyan << "To state:   " << ansi::yellow << *to   << ansi::reset << std::endl;
                        std::cout << ansi::cyan << "Hopping sign: " << ansi::yellow << hop_sign << ansi::reset << std::endl;
                        std::cout << ansi::cyan << "Index of to state: " << ansi::yellow << idx1   << ansi::reset << std::endl;
                    };
                }; //End of loop over is, dir, charge
    }; //End of loop over all basis states
    if(res)
        std::cout << ansi::green << "Test of hopping term passed!" << ansi::reset << std::endl;
    else
        std::cout << ansi::red   << "Test of hopping term failed!" << ansi::reset << std::endl;  
    delete from;
    delete to;
}

void LGT_ExDiag::construct_hamiltonian_matrix()
{
    HamiltonianMatrix = new double[NS*NS];
    std::fill(HamiltonianMatrix, HamiltonianMatrix+NS*NS, 0.0);
    LGTState* from = new LGTState(LS);
    LGTState* to   = new LGTState(LS);
    for(int jdx=0; jdx<NS; jdx++)
    {
        IndexToLGTState(jdx, from);
        //Kinetic terms in the Hamiltonian - diagonal elements
        double EK = 0.0;
        for(int is=0; is<LS; is++) EK += (double)(from->nf[is]*from->nf[is]);
        HamiltonianMatrix[jdx*NS + jdx] += (1.0/(2*g*g))*EK;
        //Hopping terms in the Hamiltonian - off-diagonal elements
        for(int is=0; is<LS; is++)
            for(int dir=-1; dir<=1; dir+=2)
                for(int charge=-1; charge<=1; charge+=2)
                {
                    int hop_sign = hopping(from, is, dir, charge, to);
                    if(hop_sign == 0) continue; //Hopping not possible
                    int idx = LGTStateToIndex(to);
                    HamiltonianMatrix[idx*NS + jdx] -= kappa*hop_sign;
                }; //End of loop over is, dir, charge
    }
    delete from;
    delete to;
}

void LGT_ExDiag::diagonalize_hamiltonian()
{
    TIMING_INIT;
    if (HamiltonianMatrix == NULL) construct_hamiltonian_matrix();
    E   = new double[NS];
	psi = new double[NS*NS]; //Eigenvectors stored column-wise
    std::copy(HamiltonianMatrix, HamiltonianMatrix+NS*NS, psi); //Copy the Hamiltonian matrix to psi, which will be overwritten by LAPACK_dsyev

	std::cout << "Running LAPACK_dsyev to find all " << ansi::yellow << NS << ansi::reset << " eigenstates of the Hamiltonian matrix " << std::endl << std::flush;
	TIMING_START;
	int res = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', NS, psi, NS, E);
	TIMING_FINISH;
	std::cout << ansi::green << " ... Done in " << ansi::magenta << a_time << ansi::reset << " sec.\n" << std::endl;
	if(res!=0)
	{
		cerr << ansi::red << "Something went wrong, LAPACKE_dsyev returned " << ansi::cyan << res << ansi::red << " !!!\n" << ansi::reset << std::endl << std::flush;
		std::exit(EXIT_FAILURE);
	};

	std::cout << "Transposing the eigensystem... " << std::endl << std::flush;
	TIMING_START;
	transpose(psi, NS);
	TIMING_FINISH;
	std::cout << ansi::green << " ... Done in " << ansi::magenta << a_time << ansi::reset << " sec.\n" << std::endl;
	
	std::cout << "Sorting the eigensystem... " << std::endl << std::flush;
	TIMING_START;
	sort_eigensystem(E, psi, NS, NS);
	TIMING_FINISH;
	std::cout << ansi::green << " ... Done in " << ansi::magenta << a_time << ansi::reset << " sec.\n" << endl;
}

void LGT_ExDiag::print_eigenvalues() const
{
    if (E == NULL) {
        std::cout << ansi::red << "Eigenvalue array E is not initialized." << ansi::reset << std::endl;
        return;
    }
    std::cout << ansi::cyan << "Eigenvalues:" << ansi::reset << std::endl;
    for (int i = 0; i < NS; ++i) {
        std::cout << E[i] << std::endl;
    }
}

double* LGT_ExDiag::construct_current_matrix()
{
    double* CurrentMatrix = new double[NS*NS];
    std::fill(CurrentMatrix, CurrentMatrix + NS*NS, 0.0);

    LGTState* from = new LGTState(LS);
    LGTState* to   = new LGTState(LS);

    for (int jdx = 0; jdx < NS; jdx++)
    {
        IndexToLGTState(jdx, from);
        // Example: Off-diagonal elements (replace with your current operator logic)
        for (int is = 0; is < LS; is++)
            for (int dir = -1; dir <= 1; dir += 2)
                for (int charge = -1; charge <= 1; charge += 2)
                {
                    int hop_sign = hopping(from, is, dir, charge, to);
                    if (hop_sign == 0) continue;
                    int idx = LGTStateToIndex(to);
                    CurrentMatrix[idx*NS + jdx] += hop_sign*charge*dir*kappa;
                }
    }
    delete from;
    delete to;
    return CurrentMatrix;
}

double* LGT_ExDiag::construct_d_matrix()
{
    double* Dmatrix = new double[NS*NS];
    std::fill(Dmatrix, Dmatrix + NS*NS, 0.0);

    LGTState* from = new LGTState(LS);
    
    for (int jdx = 0; jdx < NS; jdx++)
    {
        IndexToLGTState(jdx, from);
        int q0 = from->na[0] - from->nb[0];
        int q1 = from->na[1] - from->nb[1];
        Dmatrix[jdx*NS + jdx] = (double)(q0 - q1);
    }
    delete from;
    return Dmatrix;
}

double* LGT_ExDiag::construct_efield_matrix(int is)
{
    double* EFieldMatrix = new double[NS*NS];
    std::fill(EFieldMatrix, EFieldMatrix + NS*NS, 0.0);

    LGTState* from = new LGTState(LS);

    for (int jdx = 0; jdx < NS; jdx++)
    {
        IndexToLGTState(jdx, from);
        EFieldMatrix[jdx*NS + jdx] += from->nf[is];
    };
    delete from;
    return EFieldMatrix;
}

double* LGT_ExDiag::EuclideanCorrelator(double* O)
{
    TIMING_INIT;
    if(E == NULL || psi == NULL) { std::cout << ansi::red << "Eigenvalues or eigenvectors are not initialized." << ansi::reset << std::endl; return nullptr; };
    if(O == NULL) { std::cout << ansi::red << "Operator matrix O is not initialized." << ansi::reset << std::endl; return nullptr; };
    double* GE = new double[LT];
    std::fill(GE, GE+LT, 0.0);

    double Z = 0.0;

    std::cout << ansi::white << "Calculating the Euclidean correlators ... " << ansi::reset << std::endl << std::flush;
	TIMING_START;
    double* tmp   = new double[NS];
	for(int m=0; m<NS; m++)
	{
        psi_eq_A_mult_chi(tmp, O, psi + m*NS, NS);
        Z += exp(-beta*E[m]);
		for(int n=0; n<NS; n++) 
		{
			double Omn2 = std::norm(scalar_prod(psi + n*NS, tmp, NS)); //Squared matrix element
			//Current-current correlator
			for(uint it=0; it<LT; it++)
			{
				double tau = it*beta/(double)LT;
				GE[it] += Omn2*exp(-tau*E[n] - (beta - tau)*E[m]);
			};
		};
    };
    for(int it=0; it<LT; it++) GE[it] /= Z;
	TIMING_FINISH;
	std::cout << ansi::green << " ... Done in " << ansi::magenta << a_time << ansi::reset << " sec.\n" << endl;
    delete [] tmp;
    return GE;
}

t_complex* LGT_ExDiag::MinkowskiCorrelator(double* O) //Computes the Minkowski correlator of the operator O (given as a matrix in the LGT basis)
{
    TIMING_INIT;
    if(E == NULL || psi == NULL) { std::cout << ansi::red << "Eigenvalues or eigenvectors are not initialized." << ansi::reset << std::endl; return nullptr; };
    if(O == NULL) { std::cout << ansi::red << "Operator matrix O is not initialized." << ansi::reset << std::endl; return nullptr; };
    t_complex* GR = new t_complex[nsteps];
    std::fill(GR, GR+nsteps, t_complex(0.0, 0.0));

    double Z = 0.0;

    std::cout << ansi::white << "Calculating the Minkowski correlators for " << nsteps << " time steps ... " << ansi::reset << std::endl << std::flush;

    TIMING_START;
    double* tmp   = new double[NS];
    for(int m=0; m<NS; m++)
    {
        psi_eq_A_mult_chi(tmp, O, psi + m*NS, NS);
        Z += exp(-beta*E[m]);
        for(int n=0; n<NS; n++) 
        {
            double Omn2 = std::norm(scalar_prod(psi + n*NS, tmp, NS)); //Squared matrix element
            #pragma omp parallel for
            for(uint it=0; it<nsteps; it++)
                GR[it] += Omn2*std::exp(t_complex(-beta*E[m], (E[n] - E[m])*it*dt));
        };
    };
    for(int it=0; it<nsteps; it++) GR[it] /= Z;
    TIMING_FINISH;
    std::cout << ansi::green << " ... Done in " << ansi::magenta << a_time << ansi::reset << " sec.\n" << endl;
    delete [] tmp;
    return GR;
}

SpectralFunction* LGT_ExDiag::spectral_function(double* O)
{
    TIMING_INIT;
    if(E == NULL || psi == NULL) { std::cout << ansi::red << "Eigenvalues or eigenvectors are not initialized." << ansi::reset << std::endl; return nullptr; };
    if(O == NULL) { std::cout << ansi::red << "Operator matrix O is not initialized." << ansi::reset << std::endl; return nullptr; };
    SpectralFunction* SPF = new SpectralFunction();

    if(wmax<=0.0)
    {
        for(int i=0; i<NS; i++)
            for(int j=i+1; j<NS; j++)
            {
                double dE = std::abs(E[i] - E[j]);
                SPF->wmax = std::max(SPF->wmax, dE);
            }
    }
    else
        SPF->wmax = wmax;

    SPF->NW = NW;
    SPF->dw = SPF->wmax/(double)(NW);

    if(SPF->rho!=nullptr) delete [] SPF->rho;
    SPF->rho = new double[NW];
    std::fill(SPF->rho, SPF->rho+NW, 0.0);

    double Z = 0.0;
    
    std::cout << ansi::white << "Calculating the spectral function ... " << ansi::reset << std::endl << std::flush;
    TIMING_START;
    double* tmp   = new double[NS];
    for(int m=0; m<NS; m++)
    {
        psi_eq_A_mult_chi(tmp, O, psi + m*NS, NS);
        Z += exp(-beta*E[m]);
        for(int n=0; n<NS; n++) 
        {
            double Omn2 = std::norm(scalar_prod(psi + n*NS, tmp, NS)); //Squared matrix element
            double w = E[n] - E[m];
            int iw = (int)std::floor(w/SPF->dw);
            double SPF_contrib = Omn2*(exp(-beta*E[m]) + exp(-beta*E[n]));
			if(iw>=0 && iw < NW) SPF->rho[iw] += SPF_contrib;
        };
    };
    for(int iw=0; iw<NW; iw++) SPF->rho[iw] /= (Z*dw);
    TIMING_FINISH;
    std::cout << ansi::green << " ... Done in " << ansi::magenta << a_time << ansi::reset << " sec.\n" << endl;
    delete [] tmp;
    return SPF;
}

double LGT_ExDiag::partition_function()
{
    if(E == NULL) { std::cout << ansi::red << "Eigenvalues are not initialized." << ansi::reset << std::endl; return 0.0; };
    double Z = 0.0;
    for(int m=0; m<NS; m++) Z += exp(-beta*E[m]);
    return Z;
}

void LGT_ExDiag::test_current_hamiltonian_algebra()
{
    if(HamiltonianMatrix == NULL) construct_hamiltonian_matrix();
    double* CurrentMatrix = construct_current_matrix();

    double* EFieldMatrix = new double[NS*NS];
    std::fill(EFieldMatrix, EFieldMatrix + NS*NS, 0.0);
    for(int is=0; is<LS; is++)
    {
        double* tmp = construct_efield_matrix(is);
        A_pluseq_bB(EFieldMatrix, 1.0, tmp, NS*NS);
        delete [] tmp;
    };
        
    
    // Compute H*E
    double* HE = A_eq_B_mult_C(HamiltonianMatrix, EFieldMatrix, NS);
    double* EH = A_eq_B_mult_C(EFieldMatrix, HamiltonianMatrix, NS);
    A_pluseq_bB(EH, -1.0, HE, NS*NS); // HE now contains H*E - E*H
    
    double err = norm_diff(EH, CurrentMatrix, NS*NS);
    std::cout << std::endl << ansi::cyan << "Norm of the difference [H,E] - J: " << ansi::yellow << err << ansi::reset << std::endl;

    delete [] HE;
    delete [] EH;
    delete [] CurrentMatrix;
    delete [] EFieldMatrix;
}