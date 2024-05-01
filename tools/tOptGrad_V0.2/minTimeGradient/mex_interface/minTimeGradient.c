/*
 * [C,time,g,s,k, phi, sta, stb] = minTimeGradient(C, RIV/RV, g0, gfin, gmax, smax,T, ds)
 * [C,time,g,s,k, phi, sta, stb] = minTimeGradient(C)
 *
 * Given a k-space trajectory C = [x y z], gradient and slew constraints,
 * This function will return a new parameterization that will meet these constraints
 * while getting from one point to the other in minimum time.
 * The constraints can be either magnitude (|g| < Gmax and |s| < Smax) or as a set of
 * individual for each direction (|gx|, |gy|, |gz| < Gmax and |sx|, |sy|, |sz| < Smax).
 *
 * Input Values :
 * C        -  The curve in k-space given in any parameterization [1/cm].
 *             Accepts 2D and 3D trajectories as 2D or 3D arrays
 *             C = [x, y] or C = [x, y, z].
 * RIV/RV   -  0 for rotationally invariant solution (magnitude constraints),
               1 for rotationally variant solution (individualgradient/slew constraints).
               Default is 0 for the rotationally invariant solution.
 * g0       -  Initial gradient amplitude.
 * gfin     -  Gradient value at the end of the trajectory.
               If given value is not possible
               the result would be the largest possible amplitude.
               (Enter -1 if you do not care to get maximum gradient)
 * gmax     -  Maximum gradient [G/cm] (4 default)
 * smax     -  Maximum slew [G/cm/ms] (15 default)
 * T        -  Sampling time intervale [ms] (4e-3 default)
 * ds       -  step size for ODE integration, enter -1 to use default value
 *
 *
 * Return Values :
 *
 * C        -  The curve reparameterized and sampled at T[ms].
 *              Returned as Nx3 array.
 * time     -  total time to get to the end
 * g        -  gradient waveform [G/cm] Nx3 array.
 * s        -  slew rate [G/cm/ms] , Nx3 array.
 * k        -  exact k-space corresponding to the gradient g
 *              (This function reparametrizes   C, then takes a derivative.
 *              Numerical errors in the derivative can lead to
 *              deviation.
 * phi       -  Geometry constrains on the amplitude vs. arclength
 * sta       -  Solution for the forward ODE
 * stb       -  Solution for the backward ODE
 *
 * Examples :
 * t = linspace(-1,1, 500)';
 * C = (5+ t.*cos(5.5*pi*t)).*exp(i*pi*t);
 * C = [real(C), imag(C)];
 * [C, time, g, s, k] = minTimeGradient(C,0);
 *
 * L = length(s);
 * figure, subplot(2,2,1), plot(C(:,1), C(:,2)); axis([-6.5,6.5,-6.5,6.5]), title('k-space')
 * subplot(2,2,2), plot(g(:,1)); axis([0,L,-4.5,4.5]); title('gradient waveforms')
 * hold on, subplot(2,2,2), plot(g(:,2), 'r');
 * legend('gx', 'gy', 'Location', 'NorthEast');
 * subplot(2,2,3), plot((g(:,1).^2 + g(:,2).^2).^0.5);  axis([0 L 0 6]); title('gradient magnitude')
 * subplot(2,2,4), plot((s(:,1).^2 + s(:,2).^2).^0.5); axis([0,L,0,20]); title('slew-rate magnitude')
 */

#include <float.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include "matrix.h"
#include "mex.h"
#include "mtg_functions.c"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ) {
    
    int i = 0;
    double rv=0;               /* 0 for rotationally invariant (default) and 1 for rotationally variant solution */
    int gfin_empty = 0;        /* Indicates whether or not a gfin value was specified*/
    int ds_empty = 0;          /*   Indicates whether or not a ds value was specified */
    
    mwSize mrows, ncols;
    
    double g0, gfin, smax, gmax, T, ds;
    double *ptr;
    
    if (nrhs < 1) {
        mexErrMsgTxt("You gotta give me a curve!");
    }
    
    /* The curve input */
    double *Ci;
    Ci = mxGetPr(prhs[0]);
    
    mrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    
    /* Dimension of the input C */
    int r = mrows;
    int c = ncols;
    
    if((ncols < 2) || (ncols > 3 )) {
        mexErrMsgTxt("Curve should be either Nx2 or Nx3 double array");
    }
    
    /* Parse input curve. Cin should be given as [x y z] Nx3 double array. */
    double *x, *y, *z;
    x = double_malloc_check(mrows * sizeof(double));
    y = double_malloc_check(mrows * sizeof(double));
    z = double_malloc_check(mrows * sizeof(double));
    for(i=0; i < mrows; i++) {
        x[i] = Ci[i];
        y[i] = Ci[i+mrows];
        if(ncols != 3) {
            z[i] = 0;
        } else {
            z[i] = Ci[i+2*mrows];
        }
    }
    if (nrhs < 2) {
        rv = 0; /* Default is rotationally invariant solution */
    }else {
        mrows = mxGetM(prhs[1]);
        ncols = mxGetN(prhs[1]);
        if((mrows != 1) || (ncols != 1)) {
            mexErrMsgTxt("second input must be either 0 or 1 for type of solution");
        } else {
            ptr = mxGetPr(prhs[1]);
            rv = ptr[0];
        }
    }
    if(nrhs < 3) {
        g0 = 0;
    } else {
        mrows = mxGetM(prhs[2]);
        ncols = mxGetN(prhs[2]);
        if((mrows != 1) || (ncols != 1)) {
            mexErrMsgTxt("g0 must be a 1x1 double value");
        } else {
            ptr = mxGetPr(prhs[2]);
            g0 = ptr[0];
        }
    }
    if (nrhs < 4) {
        gfin_empty = 1; /* No gfin is specified */
        gfin = -1;
    } else {
        mrows = mxGetM(prhs[3]);
        ncols = mxGetN(prhs[3]);
        if((mrows != 1) || (ncols != 1)) {
            mexErrMsgTxt("gfin must be a 1x1 double value");
        } else {
            
            ptr = mxGetPr(prhs[3]);
            gfin = ptr[0];
            if (gfin < 0 ) {    /* negative value means gfin input ignored */
                gfin_empty = 1;
            }
        }
    }
    if (nrhs < 5) {
        gmax = 4;
    } else {
        mrows = mxGetM(prhs[4]);
        ncols = mxGetN(prhs[4]);
        if((mrows != 1) || (ncols != 1)) {
            mexErrMsgTxt("gmax must be a 1x1 double value");
        } else {
            ptr = mxGetPr(prhs[4]);
            gmax = ptr[0];
        }
    }
    if (nrhs < 6) {
        smax = 15;
    } else {
        mrows = mxGetM(prhs[5]);
        ncols = mxGetN(prhs[5]);
        if((mrows != 1) || (ncols != 1)) {
            mexErrMsgTxt("smax must be a 1x1 double value");
        } else {
            ptr = mxGetPr(prhs[5]);
            smax = ptr[0];
        }
    }
    if (nrhs < 7) {
        T = 4e-3;
    } else {
        mrows = mxGetM(prhs[6]);
        ncols = mxGetN(prhs[6]);
        if((mrows != 1) || (ncols != 1)) {
            mexErrMsgTxt("T must be a 1x1 double value");
        } else {
            ptr = mxGetPr(prhs[6]);
            T = ptr[0];
        }
    }
    if (nrhs < 8) {
        ds_empty = 1;
    } else {
        mrows = mxGetM(prhs[7]);
        ncols = mxGetN(prhs[7]);
        if((mrows != 1) || (ncols != 1)) {
            mexErrMsgTxt("ds must be a 1x1 double value");
        } else {
            ptr = mxGetPr(prhs[7]);
            ds = ptr[0];
            if (ds < 0) {
                ds_empty = 1;   /* ignore inputed ds and use default */
            }
        }
    }
    
    /*pointers to output arrays*/
    double **Cx, **Cy, **Cz, **gx, **gy, **gz, **sx, **sy, **sz, **kx, **ky, **kz, **sdot, **sta, **stb;
    Cx = doublep_malloc_check(1);
    Cy = doublep_malloc_check(1);
    Cz = doublep_malloc_check(1);
    gx = doublep_malloc_check(1);
    gy = doublep_malloc_check(1);
    gz = doublep_malloc_check(1);
    sx = doublep_malloc_check(1);
    sy = doublep_malloc_check(1);
    sz = doublep_malloc_check(1);
    kx = doublep_malloc_check(1);
    ky = doublep_malloc_check(1);
    kz = doublep_malloc_check(1);
    sdot = doublep_malloc_check(1);
    sta = doublep_malloc_check(1);
    stb = doublep_malloc_check(1);
    
    double *timep;
    int *size_interpolatedp, *size_sdotp, *size_stp;
    double time = 0;
    timep = &time;
    int size_interpolated = 0;
    size_interpolatedp = &size_interpolated;
    int size_sdot = 0;
    size_sdotp = &size_sdot;
    int size_st = 0;
    size_stp = &size_st;
        
    if (rv ==1 ) {
        /*Compute the rotationally variant solution*/
        printf("Computing the rotationally variant solution\n");
        minTimeGradientRV(Ci, r, c, g0, gfin, gmax, smax, T, ds,
                Cx, Cy, Cz, gx, gy, gz, sx, sy, sz,
                kx, ky, kz, sdot, sta, stb, timep,
                size_interpolatedp, size_sdotp, size_stp, gfin_empty, ds_empty);
    }
    else {
        /*Compute the rotationally invariant solution*/
        printf("Computing the rotationally invariant solution\n");
        minTimeGradientRIV(Ci, r, c, g0, gfin, gmax, smax, T, ds,
                Cx, Cy, Cz, gx, gy, gz, sx, sy, sz,
                kx, ky, kz, sdot, sta, stb, timep,
                size_interpolatedp, size_sdotp, size_stp, gfin_empty, ds_empty);
    }
    
    /*Convert returned arrays into mex outputs*/
    plhs[0] = mxCreateDoubleMatrix(size_interpolated, 3, mxREAL);
    double *C_ptr;
    C_ptr = mxGetPr(plhs[0]);
    
    for(i=0; i < size_interpolated; i ++) {
        C_ptr[i] = Cx[0][i];
        C_ptr[i+size_interpolated] = Cy[0][i];
        C_ptr[i+2*size_interpolated] = Cz[0][i];
    }
    
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *time_ptr;
    time_ptr = mxGetPr(plhs[1]);
    time_ptr[0] = time;
    
    plhs[2] = mxCreateDoubleMatrix(size_interpolated, 3, mxREAL);
    double *g_ptr;
    g_ptr = mxGetPr(plhs[2]);
    for(i=0; i < size_interpolated; i ++) {
        g_ptr[i] = gx[0][i];
        g_ptr[i+size_interpolated] = gy[0][i];
        g_ptr[i+2*size_interpolated] = gz[0][i];
    }
    
    plhs[3] = mxCreateDoubleMatrix(size_interpolated, 3, mxREAL);
    double *s_ptr;
    s_ptr = mxGetPr(plhs[3]);
    
    for(i=0; i < size_interpolated; i ++) {
        s_ptr[i] = sx[0][i];
        s_ptr[i+size_interpolated] = sy[0][i];
        s_ptr[i+2*size_interpolated] = sz[0][i];
    }
    
    plhs[4] = mxCreateDoubleMatrix(size_interpolated, 3, mxREAL);
    double *k_ptr;
    k_ptr = mxGetPr(plhs[4]);
    
    for(i=0; i < size_interpolated; i ++) {
        k_ptr[i] = kx[0][i];
        k_ptr[i+size_interpolated] = ky[0][i];
        k_ptr[i+2*size_interpolated] = kz[0][i];
    }
    
    
    plhs[5] = mxCreateDoubleMatrix(size_sdot, 1, mxREAL);
    double *sdot_ptr;
    sdot_ptr = mxGetPr(plhs[5]);
    for(i=0; i < size_sdot; i++) {
        sdot_ptr[i] = sdot[0][i];
    }
    
    plhs[6] = mxCreateDoubleMatrix(size_st, 1, mxREAL);
    plhs[7] = mxCreateDoubleMatrix(size_st, 1, mxREAL);
    
    double *sta_ptr, *stb_ptr;
    sta_ptr = mxGetPr(plhs[6]);
    stb_ptr = mxGetPr(plhs[7]);
    
    for(i=0; i < size_st; i++) {
        sta_ptr[i] = sta[0][i];
        stb_ptr[i] = stb[0][i];
    }
    
    mxFree(Cx[0]); mxFree(Cy[0]); mxFree(Cz[0]);
    mxFree(gx[0]); mxFree(gy[0]); mxFree(gz[0]);
    mxFree(sx[0]); mxFree(sy[0]); mxFree(sz[0]);
    mxFree(kx[0]); mxFree(ky[0]); mxFree(kz[0]);
    mxFree(sdot[0]); mxFree(sta[0]); mxFree(stb[0]);
    
    mxFree(x);  mxFree(y);  mxFree(z);
    mxFree(Cx); mxFree(Cy); mxFree(Cz);
    mxFree(gx); mxFree(gy); mxFree(gz);
    mxFree(sx); mxFree(sy); mxFree(sz);
    mxFree(kx); mxFree(ky); mxFree(kz);
    mxFree(sdot); mxFree(sta); mxFree(stb);
}
