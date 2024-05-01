/*This file contains all the functions used to calculated the time optimal gradient waveforms.

minTimeGradientRIV   -   Computes the rotationally invariant solution
     RungeKutte_riv  -   Used to solve the ODE using RK4
     beta            -   calculates sqrt (gamma^2 * smax^2 - k^2 * st^4) in the ODE. Used in RungeKutte_riv
minTimeGradientRV    -   Computes the rotationally variant solution
     RungeKutte_rv   -   Used to solve the ODE using RK4
     sdotdot         -   calculates the maximum possible value for d^2s/dt^t, used in RugeKutte_rv */

#include <float.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "spline.c"
#include <time.h>
#include <sys/types.h>

double* double_malloc_check(int L) {
    /* Double array malloc */
    double *rtn = (double*) mxMalloc(L*sizeof(double));
    if(!rtn) {
        printf("malloc failed");
        exit(1);
    }
    return rtn;
}

double** doublep_malloc_check(int L) {
    /* Double array pointer malloc  */
    double **rtn = (double**) mxMalloc(L*sizeof(double *));
    if(!rtn) {
        printf("malloc failed");
        exit(1);
    }
    return rtn;
}

double sdotdot(double xpp, double ypp, double zpp, double xp, double yp, double zp, double st, double smax) {
    /* Computes maximum possible value for sdotdot , the second time derivative of s, the arc-length.
     xp, yp, zp are the first derivatives of the curve derivatives of the curve in the arc-length(s) parameterization.
     xpp, ypp, and zpp are the second derivatives */
    double gamma = 4.257;
    double sx, sy, sz;      /* constraints for sdotdot in x, y, and z directions */
                            /* maximum sdotdot will be limited by the smallest of these */
    sx = (-xpp*st*st + gamma * smax) / xp;
    sy = (-ypp*st*st + gamma * smax) / yp;
    sz = (-zpp*st*st + gamma * smax) / zp;
    double rtn = sx;
    if (sy < rtn) {
        rtn = sy;
    }
    if(sz < rtn) {
        rtn = sz;
    }
    return rtn;
}

double beta(double k, double st, double smax) {
    /* calculates sqrt (gamma^2 * smax^2 - k^2 * st^4) used in RK4 method for rotationally invariant ODE solver */
    double gamma = 4.257;
    return  sqrt(sqrt((gamma*gamma*smax*smax - k*k*st*st*st*st)*(gamma*gamma*smax*smax - k*k*st*st*st*st)));
}

double RungeKutte_riv(double ds, double st, double k[], double smax) {
    /*  Solves ODE for rotationally invariant solution using Runge-Kutte method*/
    double k1 = ds * (1/st) * beta(k[0], st, smax);
    double k2 = ds * 1 / (st + k1/2) * beta(k[1], st + k1/2, smax);
    double k3 = ds * 1 / (st + k2/2) * beta(k[1], st + k2/2, smax);
    double k4 = ds * 1 / (st + k3/2) * beta(k[2], st + k3/2, smax);
    double rtn = k1/6 + k2/3 + k3/3 + k4/6;
    return rtn;
}

double RungeKutte_rv(double ds, double st, double xpp[], double ypp[], double zpp[],
        double xp[], double yp[], double zp[], double smax) {
    /*  Solves ODE for rotationally variant solution using Runge-Kutte method */
    double k1, k2, k3, k4;
    k1 = ds * (1/st) * sdotdot(xpp[0], ypp[0], zpp[0], xp[0], yp[0], zp[0], st, smax);
    k2 = ds * (1/(st+k1/2)) * sdotdot(xpp[1], ypp[1], zpp[1], xp[1], yp[1], zp[1], st+k1/2, smax);
    k3 = ds * (1/(st+k2/2)) * sdotdot(xpp[1], ypp[1], zpp[1], xp[1], yp[1], zp[1], st+k2/2, smax);
    k4 = ds * (1/(st+k3/2)) * sdotdot(xpp[2], ypp[2], zpp[2], xp[2], yp[2], zp[2], st+k3/2, smax);
    double rtn = k1/6 + k2/3 + k3/3 + k4/6;
    return rtn;
}


void minTimeGradientRV(double *Ci, int Cr, int Cc, double g0, double gfin, double gmax, double smax, double T, double ds,
        double **Cx, double **Cy, double **Cz, double **gx, double **gy, double **gz,
        double **sx, double **sy, double **sz, double **kx, double **ky, double **kz, double **sdot, double **sta, double **stb, double *time,
        int *size_interpolated, int *size_sdot, int *size_st, int gfin_empty, int ds_empty) {
    
    /* Finds the time optimal gradient waveforms for the rotationally variant constraints case.
    
    Ci           -       The input curve (Nx3 double array)
    Cr           -       row dimension of Ci
    Cc           -       column dimension of Ci
    g0           -       Initial gradient amplitude.
    gfin         -       Gradient value at the end of the trajectory.
                         If given value is not possible
                         the result would be the largest possible amplitude.
    gmax         -       Maximum gradient [G/cm] (4 default)
    smax         -       Maximum slew [G/cm/ms] (15 default)
    T            -       Sampling time intervale [ms] (4e-3 default)
    gx, gy, gz   -       pointers to gradient waveforms to be returned
    kx, ky, kz   -       pointers to k-space trajectory after interpolation
    sx, sy, sz   -       pointers to slew waveforms to be returned
    sdot         -       geometry constrains on the amplitude vs. arclength
    sta          -       pointer to solution for the forward ODE to be returned
    stb          -       pointer to solution for the backward ODE to be returned
    size_interpolated -  Dimension of interpolated k-space trajectory (kx, ky, kz) needed for creating mex return arrays.
    size_sdot    -       Dimension of sdot needed for creating mex return arrays.
    size_st      -       Dimension of sta and stb, needed for creating mex return array.
    gfin_empty   -       Indicates whether or not a final gradient amplitude was specified by the user */
    
    int i = 0;
    /* iflag used in spline method to signal error */
    int *iflag;
    int iflagp;
    iflag = &iflagp;
    double dt = T;
    double gamma = 4.257;
    int Lp = Cr;            /* Length of the curve in p-parameterization */
    
    /* Ci given as Nx3 array, parse into x, y, z components: */
    double *x, *y, *z;
    x = double_malloc_check(Cr * sizeof(double));
    y = double_malloc_check(Cr * sizeof(double));
    z = double_malloc_check(Cr * sizeof(double));
    for(i=0; i < Cr; i++) {
        x[i] = Ci[i];
        y[i] = Ci[i+Cr];
        if(Cc == 2) {       /* if Nx2 curve was inputed, z = 0; */
            z[i] = 0;
        } else {
            z[i] = Ci[i+2*Cr];
        }
    }
    /* Representing the curve with parameter p */
    double p [Lp];
    for (i = 0; i < Lp; i++) {
        p[i] = i;
    }
    /* Interpolation of curve in p-parameterization for gradient accuracy, using cubic spline interpolation */
    double *c1x, *c2x, *c3x,
            *c1y, *c2y, *c3y,
            *c1z, *c2z, *c3z;
    /* arrays used by spline function for spline coeficcients: */
    c1z = double_malloc_check(Lp * sizeof(double));
    c2z = double_malloc_check(Lp * sizeof(double));
    c3z = double_malloc_check(Lp * sizeof(double));
    c1x = double_malloc_check(Lp * sizeof(double));
    c2x = double_malloc_check(Lp * sizeof(double));
    c3x = double_malloc_check(Lp * sizeof(double));
    c1y = double_malloc_check(Lp * sizeof(double));
    c2y = double_malloc_check(Lp * sizeof(double));
    c3y = double_malloc_check(Lp * sizeof(double));
    
    spline(Lp, 0, 0, 1, 1, p, x, c1x, c2x, c3x, iflag);
    spline(Lp, 0, 0, 1, 1, p, y, c1y, c2y, c3y, iflag);
    spline(Lp, 0, 0, 1, 1, p, z, c1z, c2z, c3z, iflag);
    
    double dp = 0.1;
    int num_evals = (int) floor((Lp-1) / dp)+ 1;
    double *CCx, *CCy, *CCz;
    CCx = double_malloc_check(num_evals * sizeof(double));
    CCy = double_malloc_check(num_evals * sizeof(double));
    CCz = double_malloc_check(num_evals * sizeof(double));
    
    double toeval = 0;  /* used by spline eval function, seval */
    int *last;
    int holder = 0;
    last = &holder;
    
    double *Cpx, *Cpy, *Cp_abs, *Cpz;       /* interpolated curve in p-parameterization */
    Cpx = double_malloc_check(num_evals * sizeof(double));
    Cpy =  double_malloc_check(num_evals * sizeof(double));
    Cpz = double_malloc_check(num_evals * sizeof(double));
    Cp_abs =  double_malloc_check(num_evals * sizeof(double));
    
    for (i = 0; i < num_evals; i++) {
        toeval = (double) i * dp;
        CCx[i] = seval(Lp, toeval, p, x, c1x, c2x, c3x, last);
        CCy[i] = seval(Lp, toeval, p, y, c1y, c2y, c3y, last);
        CCz[i] = seval(Lp, toeval, p, z, c1z, c2z, c3z, last);
        Cpx[i] = deriv(Lp, toeval, p, c1x, c2x, c3x, last);
        Cpy[i] = deriv(Lp, toeval, p, c1y, c2y, c3y, last);
        Cpz[i] = deriv(Lp, toeval, p, c1z, c2z, c3z, last);
        Cp_abs[i] = sqrt(Cpx[i]*Cpx[i] + Cpy[i]*Cpy[i] + Cpz[i]*Cpz[i]);
    }
    mxFree(Cpx);    mxFree(Cpy);    mxFree(Cpz);
    
    /* converting to arc-length parameterization from p, using trapezoidal integration */
    double *s_of_p;
    s_of_p = double_malloc_check(num_evals * sizeof(double));
    s_of_p[0] = 0;
    double sofar = 0;
    
    for (i=0; i < num_evals; i++) {
        sofar += (Cp_abs[i]+ Cp_abs[i-1])/2;
        s_of_p[i] =  dp * sofar;
    }
    mxFree(Cp_abs);
    
    /* length of the curve */
    double L = s_of_p[num_evals-1];
    
    /* decide ds and compute st for the first point */
    double stt0 = gamma*smax;   /* always assumes first point is max slew */
    double st0 = (stt0*dt)/2;   /* start at half the gradient for accuracy close to g=0 */
    double s0 = st0*dt;
    if(ds_empty == 1 ) {        /* if a ds value was not specified */
        ds = s0/1.5;            /* smaller step size for numerical accuracy */
    }
    
    int length_of_s =  (int) floor(L/ds);
    int half_ls = (int) floor(L/(ds/2));
    
    double *s;
    s = double_malloc_check(length_of_s * sizeof(double));
    *size_st= length_of_s;
    
    for (i = 0; i< length_of_s; i++) {
        s[i] = i*ds;
    }
    
    double *s_half =  double_malloc_check(half_ls*sizeof(double));
    for (i=0; i < half_ls; i++) {
        s_half[i] = (double)i*(ds/2);
    }
    
    double *p_of_s_half;
    p_of_s_half = double_malloc_check(half_ls * sizeof(double));
    
    /* Convert from s(p) to p(s) and interpolate for accuracy */
    double *a1x, *a2x, *a3x;
    a1x = double_malloc_check(num_evals * sizeof(double));
    a2x = double_malloc_check(num_evals * sizeof(double));
    a3x = double_malloc_check(num_evals * sizeof(double));
    
    double sop_num[num_evals];
    for (i=0; i<num_evals; i++) {
        sop_num[i] = i * dp;
    }
    spline(num_evals, 0, 0, 1, 1, s_of_p, sop_num, a1x, a2x, a3x, iflag);
    
    for (i=0; i < half_ls; i++) {
        p_of_s_half[i] = seval(num_evals, s_half[i], s_of_p, sop_num, a1x, a2x, a3x, last);
    }
    
    mxFree(a1x); mxFree(a2x); mxFree(a3x);
    mxFree(s_of_p);
    
    int size_p_of_s = half_ls/2;
    double *p_of_s;
    p_of_s = double_malloc_check(size_p_of_s * sizeof(double));
    for (i=0; i<size_p_of_s; i++) {
        p_of_s[i] = p_of_s_half[2*i];
    }
    
    printf("Compute geometry dependent constraints \n");
    
    double *k;  /* k is the curvature along the curve */
    k = double_malloc_check(half_ls*sizeof(double));
    
    double *Cspx, *Cspy, *Cspz;
    /* Csp is C(s(p)) = [Cx(p(s)) Cy(p(s)) Cz(p(s))]  */
    Cspx =  double_malloc_check(length_of_s*sizeof(double));
    Cspy =  double_malloc_check(length_of_s*sizeof(double));
    Cspz =  double_malloc_check(length_of_s*sizeof(double));
    
    for (i=0; i<length_of_s; i++) {
        Cspx[i] = seval(Lp, p_of_s[i], p, x, c1x, c2x, c3x, last);
        Cspy[i] = seval(Lp, p_of_s[i], p, y, c1y, c2y, c3y, last);
        Cspz[i] = seval(Lp, p_of_s[i], p, z, c1z, c2z, c3z, last);
    }
    
    double *Csp1x, *Csp2x, *Csp3x, *Csp1y, *Csp2y, *Csp3y,  *Csp1z, *Csp2z, *Csp3z; /* arrays used for coefficients in spline interpolation function */
    Csp1x =  double_malloc_check(length_of_s*sizeof(double));
    Csp2x =  double_malloc_check(length_of_s*sizeof(double));
    Csp3x =  double_malloc_check(length_of_s*sizeof(double));
    Csp1y =  double_malloc_check(length_of_s*sizeof(double));
    Csp2y =  double_malloc_check(length_of_s*sizeof(double));
    Csp3y =  double_malloc_check(length_of_s*sizeof(double));
    Csp1z =  double_malloc_check(length_of_s*sizeof(double));
    Csp2z =  double_malloc_check(length_of_s*sizeof(double));
    Csp3z =  double_malloc_check(length_of_s*sizeof(double));
    
    /* interpolation for C(s(p)) */
    spline(length_of_s, 0, 0, 1, 1, s, Cspx, Csp1x, Csp2x, Csp3x, iflag);
    spline(length_of_s, 0, 0, 1, 1, s, Cspy, Csp1y, Csp2y, Csp3y, iflag);
    spline(length_of_s, 0, 0, 1, 1, s, Cspz, Csp1z, Csp2z, Csp3z, iflag);
    
    double *xpp, *ypp, *zpp;        /* xpp, ypp, zpp are d^2/ds^2 - second derivative in s parameterization */
    xpp = double_malloc_check(half_ls*sizeof(double));
    ypp = double_malloc_check(half_ls*sizeof(double));
    zpp = double_malloc_check(half_ls*sizeof(double));
    
    double *xp, *yp, *zp;          /* d/ds - first derivative in s parameterization */
    xp = double_malloc_check(half_ls*sizeof(double));
    yp = double_malloc_check(half_ls*sizeof(double));
    zp = double_malloc_check(half_ls*sizeof(double));
    
    /* Computing the curvature along the curve */
    for (i=0; i<half_ls; i++) {
        double kx = (deriv2(length_of_s, s_half[i], s, Csp1x, Csp2x, Csp3x, last));
        double ky = (deriv2(length_of_s, s_half[i], s, Csp1y, Csp2y, Csp3y, last));
        double kz = (deriv2(length_of_s, s_half[i], s, Csp1z, Csp2z, Csp3z, last));
        zpp[i] = sqrt(kz*kz);
        xpp[i] = sqrt(kx*kx);
        ypp[i] = sqrt(ky*ky);
        xp[i] = sqrt(deriv(length_of_s, s_half[i], s, Csp1x, Csp2x, Csp3x, last)*deriv(length_of_s, s_half[i], s, Csp1x, Csp2x, Csp3x, last));
        yp[i] = sqrt(deriv(length_of_s, s_half[i], s, Csp1y, Csp2y, Csp3y, last)*deriv(length_of_s, s_half[i], s, Csp1y, Csp2y, Csp3y, last));
        zp[i] = sqrt(deriv(length_of_s, s_half[i], s, Csp1z, Csp2z, Csp3z, last)*deriv(length_of_s, s_half[i], s, Csp1z, Csp2z, Csp3z, last));
        k[i] =  sqrt(kx*kx + ky*ky + kz*kz);        /* the curvature, magnitude of the second derivative of the curve in arc-length parameterization */
    }
    
    mxFree(CCx);    mxFree(CCy);    mxFree(CCz);
    mxFree(p_of_s_half);
    
    /* Compute the geometry dependant constraints */
    double *sdot1, *sdot2, *sdot3;
    sdot1 =  double_malloc_check(half_ls*sizeof(double));
    sdot2 =  double_malloc_check(half_ls*sizeof(double));
    sdot3 =  double_malloc_check(half_ls*sizeof(double));
    sdot[0] =  double_malloc_check(half_ls * sizeof(double));
    
    *size_sdot = half_ls;
    for (i=0; i< half_ls; i++) {
        sdot1[i] = gamma*gmax/xp[i];
        sdot2[i] = gamma*gmax/yp[i];
        sdot3[i] = gamma*gmax/zp[i];
        sdot[0][i] = sdot1[i];
        if (sdot2[i] < sdot[0][i]) {
            sdot[0][i] = sdot2[i];
        }
        if (sdot3[i] < sdot[0][i]){
            sdot[0][i] = sdot3[i];
        }
    }
    
    mxFree(s_half);
    mxFree(Cspx); mxFree(Cspy); mxFree(Cspz);
    mxFree(Csp1x); mxFree(Csp2x); mxFree(Csp3x);
    mxFree(Csp1y); mxFree(Csp2y); mxFree(Csp3y);
    mxFree(Csp1z); mxFree(Csp2z); mxFree(Csp3z);
    mxFree(sdot1);
    mxFree(sdot2);
    mxFree(sdot3);
    
    /* Extend k, xp, yp, zp, xpp, ypp, zpp for RK4 end points */
    int size_k2 = half_ls+2;
    double *k2;
    k2 = double_malloc_check(size_k2*sizeof(double));
    
    for(i=0; i < half_ls; i++) {
        k2[i] = k[i];
    }
    double *xpp2, *ypp2, *zpp2, *xp2, *yp2, *zp2;
    xpp2 = double_malloc_check(size_k2*sizeof(double));
    ypp2 = double_malloc_check(size_k2*sizeof(double));
    zpp2 = double_malloc_check(size_k2*sizeof(double));
    xp2 = double_malloc_check(size_k2*sizeof(double));
    yp2 = double_malloc_check(size_k2*sizeof(double));
    zp2 = double_malloc_check(size_k2*sizeof(double));
    
    for (i=0; i<half_ls; i++) {
        xpp2[i] = xpp[i];
        ypp2[i] = ypp[i];
        zpp2[i] = zpp[i];
        xp2[i] = xp[i];
        yp2[i] = yp[i];
        zp2[i] = zp[i];
    }
    xpp2[size_k2-2] = xpp2[size_k2-3];
    xpp2[size_k2-1] = xpp2[size_k2-3];
    
    ypp2[size_k2-2] = ypp2[size_k2-3];
    ypp2[size_k2-1] = ypp2[size_k2-3];
    
    zpp2[size_k2-2] = zpp2[size_k2-3];
    zpp2[size_k2-1] = zpp2[size_k2-3];
    
    xp2[size_k2-2] = xp2[size_k2-3];
    xp2[size_k2-1] = xp2[size_k2-3];
    
    yp2[size_k2-2] = yp2[size_k2-3];
    yp2[size_k2-1] = yp2[size_k2-3];
    
    zp2[size_k2-2] = zp2[size_k2-3];
    zp2[size_k2-1] = zp2[size_k2-3];
    
    k2[size_k2-2] = k2[size_k2-3];
    k2[size_k2-1] = k2[size_k2-3];
    
    /* Solving the ODE */
    sta[0] =  double_malloc_check(length_of_s * sizeof(double));        /* forward solution */
    stb[0] =  double_malloc_check(length_of_s * sizeof(double));        /* backward solution */
    
    /* Forward solution  */
    /* Initial cond.  */
    double g0gamma = g0*gamma + st0;
    double gammagmax = gamma *gmax;
    if (g0gamma < gammagmax) {
        sta[0][0] = g0gamma;
    } else {
        sta[0][0] = gammagmax;
    }
    
    printf("Solving ODE Forward\n");
    for (i=1; i<length_of_s; i++) {
        double k_rk[3], xpp_rk[3], ypp_rk[3], zpp_rk[3], xp_rk[3], yp_rk[3], zp_rk[3];
        k_rk[0] = k2[2*i-2];
        k_rk[1] = k2[2*i-1];
        k_rk[2] = k2[2*i];
        
        xpp_rk[0] = xpp2[2*i-2];
        xpp_rk[1] = xpp2[2*i-1];
        xpp_rk[2] = xpp2[2*i];
        
        ypp_rk[0] = ypp2[2*i-2];
        ypp_rk[1] = ypp2[2*i-1];
        ypp_rk[2] = ypp2[2*i];
        
        zpp_rk[0] = zpp2[2*i-2];
        zpp_rk[1] = zpp2[2*i-1];
        zpp_rk[2] = zpp2[2*i];
        
        xp_rk[0] = xp2[2*i-2];
        xp_rk[1] = xp2[2*i-1];
        xp_rk[2] = xp2[2*i];
        
        yp_rk[0] = yp2[2*i-2];
        yp_rk[1] = yp2[2*i-1];
        yp_rk[2] = yp2[2*i];
        
        zp_rk[0] = zp2[2*i-2];
        zp_rk[1] = zp2[2*i-1];
        zp_rk[2] = zp2[2*i];
        
        double dstds = RungeKutte_rv(ds, sta[0][i-1], xpp_rk , ypp_rk, zpp_rk, xp_rk, yp_rk, zp_rk, smax);
        double tmpst = sta[0][i-1] + dstds;
        if (sdot[0][2*i] < tmpst) {
            sta[0][i] = sdot[0][2*i];
        } else {
            sta[0][i] = tmpst;
        }
    }
    
    mxFree(k2);
    
    printf("Solving ODE Backwards\n");
    
    if (gfin_empty == 1) {  /* if gfin is not provided */
        stb[0][length_of_s-1] = sta[0][length_of_s - 1];
    } else {
        double max;
        if (gfin * gamma > st0) {
            max = gfin*gamma;
        } else {
            max = st0;
        }
        if (gamma * gmax < max) {
            stb[0][length_of_s-1] = gamma * gmax;
        } else {
            stb[0][length_of_s-1] = max;
        }
    }
    for (i=length_of_s-2; i>-1; i--) {
        double k_rk[3], xpp_rk[3], ypp_rk[3], zpp_rk[3], xp_rk[3], yp_rk[3], zp_rk[3];
        k_rk[0] = k[2*i+2];
        k_rk[1] = k[2*i+1];
        k_rk[2] = k[2*i];
        
        xpp_rk[0] = xpp2[2*i+2];
        xpp_rk[1] = xpp2[2*i+1];
        xpp_rk[2] = xpp2[2*i];
        
        ypp_rk[0] = ypp2[2*i+2];
        ypp_rk[1] = ypp2[2*i+1];
        ypp_rk[2] = ypp2[2*i];
        
        zpp_rk[0] = zpp2[2*i+2];
        zpp_rk[1] = zpp2[2*i+1];
        zpp_rk[2] = zpp2[2*i];
        
        xp_rk[0] = xp2[2*i+2];
        xp_rk[1] = xp2[2*i+1];
        xp_rk[2] = xp2[2*i];
        
        yp_rk[0] = yp2[2*i+2];
        yp_rk[1] = yp2[2*i+1];
        yp_rk[2] = yp2[2*i];
        
        zp_rk[0] = zp2[2*i+2];
        zp_rk[1] = zp2[2*i+1];
        zp_rk[2] = zp2[2*i];
        
        double dstds = RungeKutte_rv(ds, stb[0][i+1], xpp_rk, ypp_rk, zpp_rk, xp_rk, yp_rk, zp_rk, smax);
        double tmpst = stb[0][i+1] + dstds;
        
        if (sdot[0][2*i] < tmpst) {
            stb[0][i] = sdot[0][2*i];
        } else {
            stb[0][i] = tmpst;
        }
    }
    
    mxFree(k);
    mxFree(xp); mxFree(yp); mxFree(zp);
    mxFree(xp2); mxFree(yp2); mxFree(zp2);
    mxFree(xpp);    mxFree(ypp);    mxFree(zpp);
    mxFree(xpp2);   mxFree(ypp2);   mxFree(zpp2);
    
    double *st_of_s, *st_ds_i;
    st_of_s =  double_malloc_check(length_of_s*sizeof(double));
    st_ds_i =  double_malloc_check(length_of_s*sizeof(double));
    
    for (i=0; i<length_of_s; i++) {
        if (sta[0][i] < stb[0][i]) {
            st_of_s[i] = sta[0][i];             /* use the minimum of the two solutions, sta and stb */
        }
        else {
            st_of_s[i] = stb[0][i];
        }
        st_ds_i[i] = ds*(1/st_of_s[i]);         /* ds * 1/st(s) used in below calculation of t(s) */
    }
    
    /* st is v(s)
    
    printf("Final interpolation\n");
    /* Converting to the time parameterization, t(s) using trapezoidal integration. t(s) = integral (1/st) ds  */
    double *t_of_s= double_malloc_check(length_of_s*sizeof(double));
    t_of_s[0] = 0;
    for (i=1; i < length_of_s; i++) {
        t_of_s[i] =  t_of_s[i-1] + (st_ds_i[i]+ st_ds_i[i-1])/2;
    }
    
    int l_t =  (int) floor(t_of_s[length_of_s-1]/dt);   /* size of the interpolated trajectory  */
    *size_interpolated = l_t;
    
    double t[l_t];
    for (i=0; i<l_t; i++) {
        t[i] = i*dt;            /* time array  */
    }
    
    double *t1x, *t2x, *t3x;    /* coefficient arrays for spline interpolation of t(s) to get s(t)  */
    t1x = double_malloc_check(length_of_s*sizeof(double));
    t2x = double_malloc_check(length_of_s*sizeof(double));
    t3x = double_malloc_check(length_of_s*sizeof(double));
    
    double *s_of_t;
    s_of_t = double_malloc_check(l_t * sizeof(double));
    spline(length_of_s, 0, 0, 1, 1, t_of_s, s, t1x, t2x, t3x, iflag);
    
    for (i=0; i < l_t; i++){
        s_of_t[i] = seval(length_of_s, t[i], t_of_s, s, t1x, t2x, t3x, last);
    }
    mxFree(t1x);    mxFree(t2x);    mxFree(t3x);
    mxFree(st_ds_i);
    mxFree(st_of_s);
    mxFree(t_of_s);
    
    double *p1x, *p2x, *p3x;    /* coefficient arrays for spline interpolation of p(s) with s(t) to get p(s(t)) = p(t) */
    p1x = double_malloc_check(length_of_s*sizeof(double));
    p2x = double_malloc_check(length_of_s*sizeof(double));
    p3x = double_malloc_check(length_of_s*sizeof(double));
    
    spline(length_of_s, 0, 0, 1, 1, s, p_of_s, p1x, p2x, p3x, iflag);
    
    double *p_of_t;
    p_of_t = double_malloc_check(l_t*sizeof(double));
    
    for (i=0; i < l_t; i++){
        p_of_t[i] = seval(length_of_s, s_of_t[i], s, p_of_s, p1x, p2x, p3x, last);
    }
    
    mxFree(s);
    mxFree(p_of_s);
    mxFree(p1x);    mxFree(p2x);    mxFree(p3x);
    mxFree(s_of_t);
    
    /*  Reparameterized curve C, sampled at T */
    Cx[0] =  double_malloc_check(l_t * sizeof(double));
    Cy[0] =  double_malloc_check(l_t * sizeof(double));
    Cz[0] =  double_malloc_check(l_t * sizeof(double));
    
    for (i=0; i<l_t; i++) {
        Cx[0][i] = seval(Lp, p_of_t[i], p, x, c1x, c2x, c3x, last);
        Cy[0][i] = seval(Lp, p_of_t[i], p, y, c1y, c2y, c3y, last);
        Cz[0][i] = seval(Lp, p_of_t[i], p, z, c1z, c2z, c3z, last);
    }
    mxFree(p_of_t);
    mxFree(x);  mxFree(y);  mxFree(z);
    mxFree(c1x);    mxFree(c2x);    mxFree(c3x);
    mxFree(c1y);    mxFree(c2y);    mxFree(c3y);
    mxFree(c1z);    mxFree(c2z);    mxFree(c3z);
    
    /* Final gradient waveforms to be returned */
    gx[0] =  double_malloc_check(l_t * sizeof(double));
    gy[0] =  double_malloc_check(l_t * sizeof(double));
    gz[0] =  double_malloc_check(l_t * sizeof(double));
    
    for (i=0; i< l_t -1; i++) {
        gx[0][i] = (Cx[0][i+1] - Cx[0][i]) / (gamma * dt);
        gy[0][i] = (Cy[0][i+1] - Cy[0][i]) / (gamma * dt);
        gz[0][i] = (Cz[0][i+1] - Cz[0][i]) / (gamma * dt);
    }
    
    gx[0][l_t-1] = gx[0][l_t-2] + gx[0][l_t-2] - gx[0][l_t-3];
    gy[0][l_t-1] = gy[0][l_t-2] + gy[0][l_t-2] - gy[0][l_t-3];
    gz[0][l_t-1] = gz[0][l_t-2] + gz[0][l_t-2] - gz[0][l_t-3];
    
    /* k-space trajecoty to be returned (calculated by integrating gradient waveforms by trapezoidal integration) */
    kx[0] =  double_malloc_check(l_t * sizeof(double));
    ky[0] =  double_malloc_check(l_t * sizeof(double));
    kz[0] =  double_malloc_check(l_t * sizeof(double));
    
    double sofarx = 0;
    double sofary = 0;
    double sofarz = 0;
    
    kx[0][0] = 0;
    ky[0][0] = 0;
    kz[0][0] = 0;
    
    for (i=1; i < l_t; i++) {
        sofarx += (gx[0][i] + gx[0][i-1]) / 2;
        sofary += (gy[0][i] + gy[0][i-1]) / 2;
        sofarz += (gz[0][i] + gz[0][i-1]) / 2;
        kx[0][i] = sofarx * dt * gamma;
        ky[0][i] = sofary * dt * gamma;
        kz[0][i] = sofarz * dt * gamma;
    }
    
    /* slew waveforms to be returned */
    sx[0] =  double_malloc_check(l_t * sizeof(double));
    sy[0] =  double_malloc_check(l_t * sizeof(double));
    sz[0] =  double_malloc_check(l_t * sizeof(double));
    
    for (i=0; i < l_t-1; i++) {
        sx[0][i] = (gx[0][i+1] - gx[0][i])/dt;
        sy[0][i] = (gy[0][i+1] - gy[0][i])/dt;
        sz[0][i] = (gz[0][i+1] - gz[0][i])/dt;
    }
    sx[0][l_t-1] = sx[0][l_t-2];
    sy[0][l_t-1] = sy[0][l_t-2];
    sz[0][l_t-1] = sz[0][l_t-2];
    
    /* total traversal time */
    *time = t[l_t-1];
    printf("Done\n");
}

void minTimeGradientRIV(double *Ci, int Cr, int Cc, double g0, double gfin, double gmax, double smax, double T, double ds,
        double **Cx, double **Cy, double **Cz, double **gx, double **gy, double **gz,
        double **sx, double **sy, double **sz, double **kx, double **ky, double **kz, double **sdot, double **sta, double **stb, double *time,
        int *size_interpolated, int *size_sdot, int *size_st, int gfin_empty, int ds_empty) {
    
    /*Finds the time optimal gradient waveforms for the rotationally invariant constraints case.
    
    Ci           -       The input curve (Nx3 double array)
    Cr           -       row dimension of Ci
    Cc           -       column dimension of Ci
    g0           -       Initial gradient amplitude.
    gfin         -       Gradient value at the end of the trajectory.
                         If given value is not possible
                         the result would be the largest possible amplitude.
    gmax         -       Maximum gradient [G/cm] (4 default)
    smax         -       Maximum slew [G/cm/ms] (15 default)
    T            -       Sampling time intervale [ms] (4e-3 default)
    gx, gy, gz   -       pointers to gradient waveforms to be returned
    kx, ky, kz   -       pointers to k-space trajectory after interpolation
    sx, sy, sz   -       pointers to slew waveforms to be returned
    sdot         -       geometry constrains on the amplitude vs. arclength
    sta          -       pointer to solution for the forward ODE to be returned
    stb          -       pointer to solution for the backward ODE to be returned
    size_interpolated -  Dimension of interpolated k-space trajectory (kx, ky, kz) needed for creating mex return arrays.
    size_sdot    -       Dimension of sdot needed for creating mex return arrays.
    size_st      -       Dimension of sta and stb, need for creating mex return array.
    gfin_empty   -       Indicats wheter or not the final gradient amplitude was specifed */
    
    int i = 0;
    /* iflag used in spline method to signal error */
    int *iflag;
    int iflagp;
    iflag = &iflagp;
    
    double dt = T;
    double gamma = 4.257;
    
    /* Length of the curve in p-parameterization */
    int Lp = Cr;
    
    double *x, *y, *z;
    /* Ci given as Nx3 array, parse into x, y, z components */
    x = double_malloc_check(Cr * sizeof(double));
    y = double_malloc_check(Cr * sizeof(double));
    z = double_malloc_check(Cr * sizeof(double));
    
    for(i=0; i < Cr; i++) {
        x[i] = Ci[i];
        y[i] = Ci[i+Cr];
        if(Cc == 2) {
            z[i] = 0;       /* if inputed curve is Nx2, z = 0 */
        } else {
            z[i] = Ci[i+2*Cr];
        }
    }
    
    double p [Lp];
    
    /* Representing the curve with parameter p */
    
    for (i = 0; i < Lp; i++) {
        p[i] = i;
    }
    
    /* Interpolation of curve for gradient accuracy, using cubic spline interpolation */
    double *c1x, *c2x, *c3x,
            *c1y, *c2y, *c3y,
            *c1z, *c2z, *c3z;       /* arrays used by spline function to store coefficients. */
    
    c1z = double_malloc_check(Lp * sizeof(double));
    c2z = double_malloc_check(Lp * sizeof(double));
    c3z = double_malloc_check(Lp * sizeof(double));
    
    c1x = double_malloc_check(Lp * sizeof(double));
    c2x = double_malloc_check(Lp * sizeof(double));
    c3x = double_malloc_check(Lp * sizeof(double));
    
    c1y = double_malloc_check(Lp * sizeof(double));
    c2y = double_malloc_check(Lp * sizeof(double));
    c3y = double_malloc_check(Lp * sizeof(double));
    
    spline(Lp, 0, 0, 1, 1, p, x, c1x, c2x, c3x, iflag);
    spline(Lp, 0, 0, 1, 1, p, y, c1y, c2y, c3y, iflag);
    spline(Lp, 0, 0, 1, 1, p, z, c1z, c2z, c3z, iflag);
    
    double dp = 0.1;
    int num_evals = (int) floor((Lp-1) / dp)+ 1;
    
    double *CCx, *CCy, *CCz;
    CCx = double_malloc_check(num_evals * sizeof(double));
    CCy = double_malloc_check(num_evals * sizeof(double));
    CCz = double_malloc_check(num_evals * sizeof(double));
    
    double toeval = 0;
    
    int *last;
    int holder = 0;
    last = &holder;
    
    double *Cpx, *Cpy, *Cp_abs, *Cpz;               /* interpolated curve in p-parameterization */
    Cpx = double_malloc_check(num_evals * sizeof(double));
    Cpy =  double_malloc_check(num_evals * sizeof(double));
    Cpz = double_malloc_check(num_evals * sizeof(double));
    Cp_abs =  double_malloc_check(num_evals * sizeof(double));
    
    for (i = 0; i < num_evals; i++) {
        toeval = (double) i * dp;
        CCx[i] = seval(Lp, toeval, p, x, c1x, c2x, c3x, last);
        CCy[i] = seval(Lp, toeval, p, y, c1y, c2y, c3y, last);
        CCz[i] = seval(Lp, toeval, p, z, c1z, c2z, c3z, last);
        Cpx[i] = deriv(Lp, toeval, p, c1x, c2x, c3x, last);
        Cpy[i] = deriv(Lp, toeval, p, c1y, c2y, c3y, last);
        Cpz[i] = deriv(Lp, toeval, p, c1z, c2z, c3z, last);
        Cp_abs[i] = sqrt(Cpx[i]*Cpx[i] + Cpy[i]*Cpy[i] + Cpz[i]*Cpz[i]);
        
    }
    mxFree(Cpx);    mxFree(Cpy);    mxFree(Cpz);
    
    /* converting to arc-length parameterization from p, using trapezoidal integration */
    
    double *s_of_p;
    s_of_p = double_malloc_check(num_evals * sizeof(double));
    s_of_p[0] = 0;
    
    double sofar = 0;
    
    for (i=0; i < num_evals; i++) {
        
        sofar += (Cp_abs[i]+ Cp_abs[i-1])/2;
        s_of_p[i] =  dp * sofar;
    }
    
    mxFree(Cp_abs);
    
    /* length of the curve */
    double L = s_of_p[num_evals-1];
    
    /* decide ds and compute st for the first point */
    double stt0 = gamma*smax;   /* always assumes first point is max slew */
    double st0 = (stt0*dt)/2;   /* start at half the gradient for accuracy close to g=0 */
    double s0 = st0*dt;

    if (ds_empty == 1) {        /* if a ds value was not specified */
        ds = s0/1.5;     /* smaller step size for numerical accuracy */
    }
    
    int length_of_s =  (int) floor(L/ds);
    int half_ls = (int) floor(L/(ds/2));
    
    *size_sdot = half_ls;
    
    double *s;
    s = double_malloc_check(length_of_s * sizeof(double));
    sta[0] =  double_malloc_check(length_of_s * sizeof(double));
    stb[0] =  double_malloc_check(length_of_s * sizeof(double));
    
    *size_st= length_of_s;
    
    for (i = 0; i< length_of_s; i++) {
        s[i] = i*ds;
        sta[0][i] = 0;
        stb[0][i] = 0;
    }
    double *s_half =  double_malloc_check(half_ls*sizeof(double));
    
    for (i=0; i < half_ls; i++) {
        s_half[i] = (double)i*(ds/2);
    }
    
    double *p_of_s_half;
    p_of_s_half = double_malloc_check(half_ls * sizeof(double));
    
    /* Convert from s(p) to p(s) and interpolate for accuracy */
    double *a1x, *a2x, *a3x;
    a1x = double_malloc_check(num_evals * sizeof(double));
    a2x = double_malloc_check(num_evals * sizeof(double));
    a3x = double_malloc_check(num_evals * sizeof(double));
    
    double sop_num[num_evals];
    for (i=0; i<num_evals; i++) {
        sop_num[i] = i * dp;
    }
    
    spline(num_evals, 0, 0, 1, 1, s_of_p, sop_num, a1x, a2x, a3x, iflag);
    
    for (i=0; i < half_ls; i++) {
        p_of_s_half[i] = seval(num_evals, s_half[i], s_of_p, sop_num, a1x, a2x, a3x, last);
    }
    
    mxFree(a1x); mxFree(a2x); mxFree(a3x);
    mxFree(s_of_p);
    
    int size_p_of_s = half_ls/2;
    
    double *p_of_s;
    p_of_s = double_malloc_check(size_p_of_s * sizeof(double));
    
    for (i=0; i<size_p_of_s; i++) {
        p_of_s[i] = p_of_s_half[2*i];
    }
    
    printf("Compute geometry dependent constraints \n");
    
    double *k;
    k = double_malloc_check(half_ls*sizeof(double));        /* k is the curvature along the curve */
    
    double *Cspx, *Cspy, *Cspz;
    /* Csp is C(s(p)) = [Cx(p(s)) Cy(p(s)) Cz(p(s))] */
    Cspx =  double_malloc_check(length_of_s*sizeof(double));
    Cspy =  double_malloc_check(length_of_s*sizeof(double));
    Cspz =  double_malloc_check(length_of_s*sizeof(double));
    
    for (i=0; i<length_of_s; i++) {
        Cspx[i] = seval(Lp, p_of_s[i], p, x, c1x, c2x, c3x, last);
        Cspy[i] = seval(Lp, p_of_s[i], p, y, c1y, c2y, c3y, last);
        Cspz[i] = seval(Lp, p_of_s[i], p, z, c1z, c2z, c3z, last);
    }
    
    double *Csp1x, *Csp2x, *Csp3x, *Csp1y, *Csp2y, *Csp3y,  *Csp1z, *Csp2z, *Csp3z; /* arrays used by spline function to store coefficients. */
    Csp1x =  double_malloc_check(length_of_s*sizeof(double));
    Csp2x =  double_malloc_check(length_of_s*sizeof(double));
    Csp3x =  double_malloc_check(length_of_s*sizeof(double));
    Csp1y =  double_malloc_check(length_of_s*sizeof(double));
    Csp2y =  double_malloc_check(length_of_s*sizeof(double));
    Csp3y =  double_malloc_check(length_of_s*sizeof(double));
    Csp1z =  double_malloc_check(length_of_s*sizeof(double));
    Csp2z =  double_malloc_check(length_of_s*sizeof(double));
    Csp3z =  double_malloc_check(length_of_s*sizeof(double));
    spline(length_of_s, 0, 0, 1, 1, s, Cspx, Csp1x, Csp2x, Csp3x, iflag);
    spline(length_of_s, 0, 0, 1, 1, s, Cspy, Csp1y, Csp2y, Csp3y, iflag);
    spline(length_of_s, 0, 0, 1, 1, s, Cspz, Csp1z, Csp2z, Csp3z, iflag);
    
    for (i=0; i<half_ls; i++) {
        double kx = (deriv2(length_of_s, s_half[i], s, Csp1x, Csp2x, Csp3x, last));
        double ky = (deriv2(length_of_s, s_half[i], s, Csp1y, Csp2y, Csp3y, last));
        double kz = (deriv2(length_of_s, s_half[i], s, Csp1z, Csp2z, Csp3z, last));
        k[i] =  sqrt(kx*kx + ky*ky + kz*kz);        /* the curvature, magnitude of the second derivative of the curve in arc-length parameterization */
    }
    
    mxFree(CCx);    mxFree(CCy);    mxFree(CCz);
    mxFree(p_of_s_half);
    
    /* computing geomtry dependent constraints (forbidden line curve) */
    
    double *sdot1, *sdot2;
    sdot1 =  double_malloc_check(half_ls*sizeof(double));
    sdot2 =  double_malloc_check(half_ls*sizeof(double));
    
    sdot[0] =  double_malloc_check(half_ls * sizeof(double));
    *size_sdot = half_ls;
    /* Calculating the upper bound for the time parametrization */
    /* sdot (which is a non scaled max gradient constaint) as a function of s. */
    /* sdot is the minimum of gamma*gmax and sqrt(gamma*gmax / k) */
    
    for (i=0; i< half_ls; i++) {
        sdot1[i] = gamma*gmax;
        sdot2[i] = sqrt((gamma*smax) / (fabs(k[i]+(DBL_EPSILON))));
        if (sdot1[i] < sdot2[i]) {
            sdot[0][i] = sdot1[i];
        }else {
            sdot[0][i] = sdot2[i];
        }
    }
    mxFree(sdot1);
    mxFree(sdot2);
    mxFree(s_half);
    mxFree(Cspx); mxFree(Cspy); mxFree(Cspz);
    mxFree(Csp1x); mxFree(Csp2x); mxFree(Csp3x);
    mxFree(Csp1y); mxFree(Csp2y); mxFree(Csp3y);
    mxFree(Csp1z); mxFree(Csp2z); mxFree(Csp3z);
    
    int size_k2 = half_ls+2;    /* extend of k for RK4 */
    double *k2;
    k2 = double_malloc_check(size_k2*sizeof(double));
    
    for(i=0; i < half_ls; i++) {
        k2[i] = k[i];
    }
    
    k2[size_k2-2] = k2[size_k2-3];
    k2[size_k2-1] = k2[size_k2-3];
    
    double g0gamma = g0*gamma + st0;
    double gammagmax = gamma *gmax;
    
    if (g0gamma < gammagmax) {
        sta[0][0] = g0gamma;
    } else {
        sta[0][0] = gammagmax;
    }
    
    /* Solving ODE Forward */
    printf("Solving ODE Forward\n");
    
    for (i=1; i<length_of_s; i++) {
        double k_rk[3];
        k_rk[0] = k2[2*i-2];
        k_rk[1] = k2[2*i-1];
        k_rk[2] = k2[2*i];
        
        double dstds = RungeKutte_riv(ds, sta[0][i-1], k_rk, smax);
        double tmpst = sta[0][i-1] + dstds;
        if (sdot[0][2*i+1] < tmpst) {
            sta[0][i] = sdot[0][2*i+1];
        } else {
            sta[0][i] = tmpst;
        }
    }
    
    mxFree(k2);
    
    /*Solving ODE Backwards: */
    
    printf("Solving ODE Backwards\n");
    
    double max;
    if(gfin_empty == 1) {
        /*if gfin is not provided */
        stb[0][length_of_s-1] = sta[0][length_of_s - 1];
    } else {
        
        if (gfin * gamma > st0) {
            max = gfin*gamma;
        } else {
            max = st0;
        }
        
        if (gamma * gmax < max) {
            stb[0][length_of_s-1] = gamma * gmax;
        } else {
            stb[0][length_of_s-1] = max;
        }
    }
    
    for (i=length_of_s-2; i>-1; i--) {
        double k_rk[3];
        k_rk[0] = k[2*i+2];
        k_rk[1] = k[2*i+1];
        k_rk[2] = k[2*i];
        
        double dstds = RungeKutte_riv(ds, stb[0][i+1], k_rk, smax);
        double tmpst = stb[0][i+1] + dstds;
        
        if (sdot[0][2*i] < tmpst) {
            stb[0][i] = sdot[0][2*i];
        } else {
            stb[0][i] = tmpst;
        }
    }
    
    /* take st(s) to be the minimum of the curves sta and stb */
    double *st_of_s, *st_ds_i;
    st_of_s = double_malloc_check(length_of_s*sizeof(double));
    st_ds_i = double_malloc_check(length_of_s*sizeof(double));
    
    for (i=0; i<length_of_s; i++) {
        if (sta[0][i] < stb[0][i]) {
            st_of_s[i] = sta[0][i];
        }
        else {
            st_of_s[i] = stb[0][i];
        }
        
        st_ds_i[i] = ds*(1/st_of_s[i]);         /* ds * 1/st(s) used in below calculation of t(s) */
    }
    
    /*Final interpolation */
    printf("Final interpolation\n");
    
    /* Converting to the time parameterization, t(s) using trapezoidal integration. t(s) = integral (1/st) ds */
    double *t_of_s= double_malloc_check(length_of_s*sizeof(double));
    t_of_s[0] = 0;
    for (i=1; i < length_of_s; i++) {
        t_of_s[i] =  t_of_s[i-1] + (st_ds_i[i]+ st_ds_i[i-1])/2;
    }
    
    int l_t =  (int) floor(t_of_s[length_of_s-1]/dt);
    *size_interpolated = l_t;       /* size of the interpolated trajectory */
    
    double t[l_t];
    for (i=0; i<l_t; i++) {
        t[i] = i*dt;                /* time array */
    }
    
    double *t1x, *t2x, *t3x;        /* coefficient arrays for spline interpolation of t(s) to get s(t) */
    
    t1x = double_malloc_check(length_of_s*sizeof(double));
    t2x = double_malloc_check(length_of_s*sizeof(double));
    t3x = double_malloc_check(length_of_s*sizeof(double));
    
    double *s_of_t;
    s_of_t = double_malloc_check(l_t * sizeof(double));
    
    spline(length_of_s, 0, 0, 1, 1, t_of_s, s, t1x, t2x, t3x, iflag);
    
    for (i=0; i < l_t; i++){
        s_of_t[i] = seval(length_of_s, t[i], t_of_s, s, t1x, t2x, t3x, last);
    }
    
    mxFree(t1x);    mxFree(t2x);    mxFree(t3x);
    mxFree(st_ds_i);
    mxFree(st_of_s);
    mxFree(t_of_s);
    
    double *p1x, *p2x, *p3x;        /* coefficient arrays for spline interpolation of p(s) with s(t) to get p(s(t)) = p(t) */
    p1x = double_malloc_check(length_of_s*sizeof(double));
    p2x = double_malloc_check(length_of_s*sizeof(double));
    p3x = double_malloc_check(length_of_s*sizeof(double));
    
    spline(length_of_s, 0, 0, 1, 1, s, p_of_s, p1x, p2x, p3x, iflag);
    
    double *p_of_t;
    p_of_t = double_malloc_check(l_t*sizeof(double));
    
    for (i=0; i < l_t; i++){
        p_of_t[i] = seval(length_of_s, s_of_t[i], s, p_of_s, p1x, p2x, p3x, last);
    }
    
    mxFree(s);
    mxFree(p_of_s);
    mxFree(p1x);    mxFree(p2x);    mxFree(p3x);
    mxFree(s_of_t);
    
    /*  interpolated k-space trajectory */
    
    Cx[0] =  double_malloc_check(l_t * sizeof(double));
    Cy[0] =  double_malloc_check(l_t * sizeof(double));
    Cz[0] =  double_malloc_check(l_t * sizeof(double));
    
    for (i=0; i<l_t; i++) {
        Cx[0][i] = seval(Lp, p_of_t[i], p, x, c1x, c2x, c3x, last);
        Cy[0][i] = seval(Lp, p_of_t[i], p, y, c1y, c2y, c3y, last);
        Cz[0][i] = seval(Lp, p_of_t[i], p, z, c1z, c2z, c3z, last);
    }
    
    mxFree(x);  mxFree(y);  mxFree(z);
    mxFree(p_of_t);
    mxFree(c1x);    mxFree(c2x);    mxFree(c3x);
    mxFree(c1y);    mxFree(c2y);    mxFree(c3y);
    mxFree(c1z);    mxFree(c2z);    mxFree(c3z);
    
    /* Final gradient waveforms to be returned */
    gx[0] =  double_malloc_check(l_t * sizeof(double));
    gy[0] =  double_malloc_check(l_t * sizeof(double));
    gz[0] =  double_malloc_check(l_t * sizeof(double));
    
    for (i=0; i< l_t -1; i++) {
        gx[0][i] = (Cx[0][i+1] - Cx[0][i]) / (gamma * dt);
        gy[0][i] = (Cy[0][i+1] - Cy[0][i]) / (gamma * dt);
        gz[0][i] = (Cz[0][i+1] - Cz[0][i]) / (gamma * dt);
    }
    
    gx[0][l_t-1] = gx[0][l_t-2] + gx[0][l_t-2] - gx[0][l_t-3];
    gy[0][l_t-1] = gy[0][l_t-2] + gy[0][l_t-2] - gy[0][l_t-3];
    gz[0][l_t-1] = gz[0][l_t-2] + gz[0][l_t-2] - gz[0][l_t-3];
    
    
    /* k-space trajecoty to be returned (calculated by integrating gradient waveforms by trapezoidal integration) */
    kx[0] =  double_malloc_check(l_t * sizeof(double));
    ky[0] =  double_malloc_check(l_t * sizeof(double));
    kz[0] =  double_malloc_check(l_t * sizeof(double));
    
    double sofarx = 0;
    double sofary = 0;
    double sofarz = 0;
    
    kx[0][0] = 0;
    ky[0][0] = 0;
    kz[0][0] = 0;
    
    for (i=1; i < l_t; i++) {
        sofarx += (gx[0][i] + gx[0][i-1]) / 2;
        sofary += (gy[0][i] + gy[0][i-1]) / 2;
        sofarz += (gz[0][i] + gz[0][i-1]) / 2;
        kx[0][i] = sofarx * dt * gamma;
        ky[0][i] = sofary * dt * gamma;
        kz[0][i] = sofarz * dt * gamma;
    }
    mxFree(k);
    /* slew waveforms to be returned */
    sx[0] =  double_malloc_check(l_t * sizeof(double));
    sy[0] =  double_malloc_check(l_t * sizeof(double));
    sz[0] =  double_malloc_check(l_t * sizeof(double));
    
    for (i=0; i < l_t-1; i++) {
        sx[0][i] = (gx[0][i+1] - gx[0][i])/dt;
        sy[0][i] = (gy[0][i+1] - gy[0][i])/dt;
        sz[0][i] = (gz[0][i+1] - gz[0][i])/dt;
    }
    sx[0][l_t-1] = sx[0][l_t-2];
    sy[0][l_t-1] = sy[0][l_t-2];
    sz[0][l_t-1] = sz[0][l_t-2];
    
    /* total traversal time */
    *time = t[l_t-1];
    printf("Done\n");
}
