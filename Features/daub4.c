/* daub4.c - Daubechies 4-coeff wavelet transform
 * 
 * Applies the daubechies 4-coefficient multiscale wavelet transform on  
 * columns of the input matrix.
 *
 * USAGE
 *
 *      [W] = daub4(D, [isign])
 *
 * Compile with:
 *
 *      mex daub4.c
 *
 * INPUT
 *  D - signal matrix (each column of D is transformed independently)
 *  isign - (optional) scalar >=0 for forward wavelet transform (default), 
 *           or <0 for inverse wavelet transform
 *
 * OUTPUT
 *  W - wavelet coefficients (Nsamples-by-Nsignals)
 *
 * This is a MEX file for MATLAB.
 *
 * daub4() and pdaub4() adapted from "Numerical Recipies in C"
 *
 * Brendan Hasz
 * haszx010@umn.edu
 * June 2016
 * David Redish Lab, University of Minnesota, Twin Cities
 */


#include "mex.h"
#define D0 0.4829629131445341  //Daubechies coeffs
#define D1 0.8365163037378079
#define D2 0.2241438680420134
#define D3 -0.1294095225512604


/* Daubechies 4-coeff partial (single scale) wavelet transform
 * 
 * INPUT
 *  a - the data (must be length power of 2)
 *  n - length of a
 *  isign - whether to do a forward (isign>=0) or inverse (isign<0) wavelet transform
 */
void pdaub4(double * a, unsigned int n, int isign)
{
    double * wksp;
    double * owksp;
    unsigned int nh, nh1, i, j;

    if (n < 4) return;
    owksp = malloc(n*sizeof(double));
    wksp = owksp-1; //for 1-offset
    nh = n>>1; //midpoint
    nh1 = nh+1; //midpoint + 1
    if (isign >= 0) { //forward transform
        for (i=1,j=1; j<=n-3; j+=2,i++) {
            wksp[i] = D0*a[j] + D1*a[j+1] + D2*a[j+2] + D3*a[j+3];
            wksp[i+nh] = D3*a[j] - D2*a[j+1] + D1*a[j+2] - D0*a[j+3];
        }
        wksp[i] = D0*a[n-1] + D1*a[n] + D2*a[1] + D3*a[2]; //wrap-arounds
        wksp[i+nh] = D3*a[n-1] - D2*a[n] + D1*a[1] - D0*a[2];
    } else { //inverse transform
        wksp[1] = D2*a[nh] + D1*a[n] + D0*a[1] + D3*a[nh1];
        wksp[2] = D3*a[nh] - D0*a[n] + D1*a[1] - D2*a[nh1];
        for (i=1,j=3; i<nh; i++) {
            wksp[j++] = D2*a[i] + D1*a[i+nh] + D0*a[i+1] + D3*a[i+nh1];
            wksp[j++] = D3*a[i] - D0*a[i+nh] + D1*a[i+1] - D2*a[i+nh1];
        }
    }
    for (i=1; i<=n; i++) a[i] = wksp[i];
    free(owksp);
}


/* Daubechies 4-coeff multiscale wavelet transform
 * 
 * INPUT
 *  a - the data (must be length power of 2)
 *  n - length of a
 *  isign - whether to do a forward (isign>=0) or inverse (isign<0) wavelet transform
 */
void daub4(double * a, unsigned int n, int isign)
{
    unsigned int nn;
    if (n < 4) return;
    if (isign >= 0) {
        for (nn=n; nn>=4; nn>>=1) pdaub4(a-1, nn, isign); //a-1 for 1-offset
    } else {
        for (nn=4; nn<=n; nn<<=1) pdaub4(a-1, nn, isign);
    }
}


/* MATLAB MEX FUNCTION: [WC] = daub4_mex(D) */
void mexFunction( int nlhs, mxArray *plhs[], 
                  int nrhs, const mxArray * prhs[])
{

    // Declare variables
    double * d; //signal and the wavelet-transformed signal
    unsigned int N; //number of elements in the signal
    unsigned int Ns; //number of signals
    const mwSize * dims; //dimensions of input array
    unsigned int i, j, k;
    double * outPr;
    double * inPr;
    int isign;


    // Get inputs
    inPr = (double *) mxGetPr(prhs[0]);
    dims = mxGetDimensions(prhs[0]); //dimensions of input array
    N = (unsigned int) dims[0];
    if (mxGetNumberOfDimensions(prhs[0]) == 1) { //only one signal
        Ns = 1;
    } else { //>1 signal, find out how many
        Ns = (unsigned int) dims[1]; //second dim is # signals
    }
    if (nrhs>1) { //passed more than 1 args
        isign = (int) mxGetScalar(prhs[1]); // >=0 for forward transform, <0 for inverse transform
    } else { //otherwise assume forward transform
        isign = 1;
    }


    // Copy data into local array
    d = malloc(N*Ns*sizeof(double));
    memcpy(d, inPr, N*Ns*sizeof(double));

    // Run the C function on each signal
    for (i=0; i<Ns; i++) {
        daub4(d+i*N, N, isign); //wavelet transform on each signal
    }

    // Allocate output array
    plhs[0] = mxCreateDoubleMatrix(N, Ns, mxREAL); //allocate output array in matlab
    outPr = mxGetPr(plhs[0]);
    memcpy(outPr, d, N*Ns*sizeof(double)); //copy data -> output array
    free(d);

}



