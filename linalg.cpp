#include "linalg.hpp"

double test_hermiticity(const double* O, int N, int sign, bool verbose)
{
    if (O == NULL) {
        std::cout << ansi::red << "Operator Matrix is not initialized." << ansi::reset << std::endl;
        return -1.0;
    }
    double max_err = 0.0;
    for(int i=0; i<N; i++)
        for(int j=0; j<N; j++) 
        {
            double err = std::abs(O[i*N + j] - sign*O[j*N + i]);
            max_err = std::max(max_err, err);
        }
    if(verbose)
        std::cout << ansi::cyan << "Maximum error in Hermiticity test: " << ansi::yellow << max_err << ansi::reset << std::endl;
    return max_err;
}

void        rescale(double* A, double a, uint n)
{
	cblas_dscal(n, a, A, 1);
}

void		rescale(t_complex* A, t_complex a, uint n)
{
	cblas_zscal(n, &a, A, 1);
}

void  A_pluseq_bB(double* A, double b, double *B, uint n)
{
	cblas_daxpy(n,b,B,1,A,1);
}

void  A_pluseq_bB(t_complex* A, t_complex b, t_complex *B, uint n)
{
	cblas_zaxpy(n,&b,B,1,A,1);
}

void  A_eq_B_mult_C(double* A, double* B, double* C, uint n) //Matrix multiplication C = A*B
{
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, B, n, C, n, 0.0, A, n);
}

void  A_eq_B_mult_C(t_complex* A, t_complex* B, t_complex* C, uint n) //Matrix multiplication C = A*B
{
	t_complex alpha = 1.0 + 0.0i;
	t_complex beta  = 0.0 + 0.0i; 
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, &alpha, B, n, C, n, &beta, A, n);
}

void  commutator(t_complex* C, t_complex* A, t_complex* B, uint n) //C = [A, B] = A*B - B*A
{
	t_complex alpha = 1.0 + 0.0i; t_complex beta  = 0.0 + 0.0i; 
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, &alpha, A, n, B, n, &beta, C, n);
	alpha = -1.0 + 0.0i; beta = 1.0 + 0.0i;
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, &alpha, B, n, A, n, &beta, C, n);
}


void		A_pluseq_a_psi_dirprod_chi(double* A, double a, double* psi, double* chi, uint n)
{
	for(uint ij=0; ij<n*n; ij++)
	{
		uint i = ij/n; uint j = ij - i*n;
		A[ij] += a*psi[i]*chi[j];
	};
}

double  scalar_prod(double* psi1, double* psi2, uint n)
{
	return cblas_ddot(n, psi1, 1, psi2, 1);
}

t_complex   scalar_prod(t_complex* psi1, t_complex* psi2, uint n)
{
	t_complex res  = 0.0 + 0.0i;
	cblas_zdotc_sub(n, psi1, 1, psi2, 1, &res);
	return res;
}

double      norm(double*    psi, uint n)
{
	return cblas_dnrm2(n, psi, 1);
};

double   norm(t_complex* psi, uint n)
{
	return cblas_dznrm2(n, psi, 1);
};

double  norm_diff(double* psi1, double* psi2, uint n)
{
	double res = 0.0;
	for(uint i=0; i<n; i++)
	{
		double d = psi1[i] - psi2[i];
		res += d*d;
	};
	return sqrt(res);
}

double    norm_diff(t_complex* psi1, t_complex* psi2, uint n)
{
	double res = 0.0;
	for(uint i=0; i<n; i++)
	{
		t_complex d = (psi1[i] - psi2[i]);
		res += real(conj(d)*d);
	};
	return sqrt(res);
}

double      unitarity_norm(t_complex* U, uint n)
{
	double unorm = 0.0;
	for(uint i=0; i<n; i++)
		for(uint j=0; j<n; j++)
		{
			t_complex r = (i==j? -1.0 + 0.0i : 0.0 + 0.0i);
			for(uint k=0; k<n; k++)
				r += U[i*n + k]*conj(U[j*n + k]);
			unorm += real(r*conj(r));
		}
	return sqrt(unorm);
}

void test_eigensystem(t_complex* A, t_complex* evals, t_complex* revecs, uint n, uint nev, double* evec_err)
{
	t_complex* tmp = new t_complex[n];
	double max_evec_err = 0.0;
	if(evec_err==NULL) std::cout << std::endl << "Checking the eigensystem... " << std::flush;
	//Eigenvector error
	for(uint ie=0; ie<nev; ie++)
	{
		psi_eq_A_mult_chi(tmp, A, revecs + n*ie, n);
		A_pluseq_bB(tmp, -evals[ie], revecs + n*ie, n);
		max_evec_err = std::max(max_evec_err, norm(tmp, n));
	};
	delete [] tmp;
	if(evec_err!=NULL)
		*evec_err = max_evec_err;
	else
		std::cout << "\t Max. eigenstate error: " << max_evec_err << std::endl; 
}

void        double2complex(t_complex* out, double* in, uint n)
{
	for(uint i=0; i<n; i++)
		out[i] = in[i] + 0.0i;
}

t_complex*  double2complex(double* in, uint n)
{
	t_complex* res = new t_complex[n];
	double2complex(res, in, n);
	return res;
}

double      unitarity_err(t_complex* U, uint n)
{
	double max_err = 0.0;
	for(uint i=0; i<n; i++)
		for(uint j=0; j<n; j++)
		{
			t_complex UUd = 0.0 + 0.0i;
			for(uint k=0; k<n; k++)
				UUd += U[i*n + k]*conj(U[j*n + k]);
			UUd -= (i==j? 1.0 + 0.0i : 0.0i);
			max_err = std::max(max_err, std::abs(UUd));
		};
	return max_err;
}

void hermitian_conjugate(t_complex* out, t_complex* in, uint n)
{
	for(uint i=0; i<n; i++)
		for(uint j=0; j<n; j++)
		{
			out[i*n + j] = conj(in[j*n + i]);
		};
}

t_complex* hermitian_conjugate(t_complex* A, uint n)
{
	t_complex* r = new t_complex[n*n];
	hermitian_conjugate(r, A, n);
	return r;
}

double abs_det(t_complex* A, uint n)
{
	if(n==1) return std::abs(A[0]);
	if(n==2) return std::abs(A[0]*A[3] - A[1]*A[2]);
	if(n==3)
	{
		return std::abs(-A[2]*A[4]*A[6] + A[1]*A[5]*A[6] + A[2]*A[3]*A[7] - A[0]*A[5]*A[7] - A[1]*A[3]*A[8] + A[0]*A[4]*A[8]);
	}
	else
	{
		t_complex* LU = new t_complex[n*n];
		std::copy(A, A + n*n, LU);
		int* piv = new int[n];
		int info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, LU, n, piv);
		if(info!=0)
		{
			std::cerr << "Error in LU decomposition for determinant calculation, info = " << info << std::endl;
			for(int i=0; i<n; i++)
			{
				for(int j=0; j<n; j++)
					std::cerr << LU[i*n + j] << " ";
				std::cerr << std::endl;
			}

			exit(1);
		};
		double r = 1.0;
		for(uint i=0; i<n; i++)
			r *= std::abs(LU[i*n + i]);
				
		delete[] LU;
		delete[] piv;
		return r;
	};
	return 0.0;
}

void identity_matrix(t_complex* A, uint n)
{
	std::fill(A, A + n*n, t_complex(0.0, 0.0));
	for(uint i=0; i<n; i++)
		A[i*n + i] = t_complex(1.0, 0.0);
}

void inverse_matrix(t_complex* A, t_complex* A_inv, uint n)
{
	std::copy(A, A + n*n, A_inv);
	int* piv = new int[n];
	int info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, A_inv, n, piv);
	if(info!=0)
	{
		std::cerr << "Error in LU decomposition for matrix inversion, info = " << info << std::endl;
		exit(1);
	};
	info = LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, A_inv, n, piv);
	if(info!=0)
	{
		std::cerr << "Error in matrix inversion, info = " << info << std::endl;
		exit(1);
	};
	delete[] piv;
}

void transpose(double* A, int n)
{
    for(int i=0; i<n; i++)
        for(int j=i+1; j<n; j++)
        {
            double tmp = A[i*n + j];
            A[i*n + j] = A[j*n + i];
            A[j*n + i] = tmp;
        }
}

void transpose(t_complex* A, int n)
{
    for(int i=0; i<n; i++)
        for(int j=i+1; j<n; j++)
        {
            t_complex tmp = A[i*n + j];
            A[i*n + j] = A[j*n + i];
            A[j*n + i] = tmp;
        }
}

void conjugate_transpose(t_complex* A, int n)
{
	for(int i=0; i<n; i++)
	{
		A[i*n + i] = conj(A[i*n + i]);
		for(int j=i+1; j<n; j++)
		{
			t_complex tmp = conj(A[i*n + j]);
			A[i*n + j] = conj(A[j*n + i]);
			A[j*n + i] = tmp;
		}
	}
}