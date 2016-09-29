#include "mex.h"
#include <Eigen/Core>

template <class Matrix, class Matrix2, typename R>
void computeLipschitz(Matrix &L, Matrix2& idxOfMaxL, const R *vr, const R *vi, 
    const R *fzr, const R *fzi, double *Lfz)
{
    const size_t m = L.rows();
    const size_t n = L.cols();
    using complex = std::complex < R > ;

#pragma omp parallel for
    for (int i = 0; i < m; i++){
        size_t idx0 = size_t(idxOfMaxL(i, 0) - 1);
        size_t idx1 = size_t(idxOfMaxL(i, 1) - 1);
        complex c_fz( fzr[idx0] - fzr[idx1], fzi[idx0] - fzi[idx1] );
        c_fz /= complex(vr[idx1] - vr[idx0], vi[idx1] - vi[idx0]);

        complex d_fz = complex(-fzr[idx0], -fzi[idx0]) - c_fz*complex(vr[idx0], vi[idx0]);

        Lfz[i] = 0;
        for (size_t j = 0; j < n; j++){
            Lfz[i] += L(i, j)*std::abs( d_fz+complex(fzr[j], fzi[j]) + c_fz*complex(vr[j], vi[j]) );
        }
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    if (nrhs != 4) {
        mexErrMsgTxt("Incorrect inputs!\n");
        return;
    }

	const size_t m = mxGetM(prhs[0]);
	const size_t n = mxGetN(prhs[0]);

	//const mwSize *dim = mxGetDimensions(m);
    using Mat = Eigen::Map < const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> > ;

    //Mat L(mxGetData(prhs[0]), m, n);
    Mat L(mxGetPr(prhs[0]), m, n);
    Mat idxOfMaxL(mxGetPr(prhs[1]), m, 2);

    plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);

    const mxArray *v = prhs[2];
    const mxArray *fz = prhs[3];
    computeLipschitz(L, idxOfMaxL, mxGetPr(v), mxGetPi(v), mxGetPr(fz), mxGetPi(fz), mxGetPr(plhs[0]));
}
