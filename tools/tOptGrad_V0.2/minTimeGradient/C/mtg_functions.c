#include <util.h>
#include <float.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>

double sdotdot(double xpp, double ypp, double zpp, double xp, double yp, double zp, double st, double smax) {
    /* Computes maximum possible value for sdotdot , the second time derivative of s, the arc-length.*/
    
    /* xp, yp, zp are the first derivatives of the curve derivatives of the curve in the arc-length(s) parameterization. */
    /* xpp, ypp, and zpp are the second derivatives */
	
    double gamma = 4.257;
    double sx, sy, sz;		/* constraints for sdotdot in x, y, and z directions */
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

double RungeKutte(double ds, double st, double xpp[], double ypp[], double zpp[],
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

double beta(double k, double st, double smax) {
    /* calculates sqrt (gamma^2 * smax^2 - k^2 * st^4) used in RK4 method for rotationally invariant ODE solver */
    double gamma = 4.257;
    return  1/st * sqrt(sqrt((gamma*gamma*smax*smax - k*k*st*st*st*st)*(gamma*gamma*smax*smax - k*k*st*st*st*st)));
}

double RungeKutte_riv(double ds, double st, double k[], double smax) {
    /*  Solves ODE for rotationally invariant solution using Runge-Kutte method */
    double k1 =  beta(k[0], st, smax);
    double k2 =  beta(k[1], st + k1/2, smax);
    double k3 =  beta(k[1], st + k2/2, smax);
    double k4 =  beta(k[2], st + k3/2, smax);
    double rtn = ds*(k1/6 + k2/3 + k3/3 + k4/6);
    return rtn + st;
}

double* double_malloc_check(int L) {
    /* Double array malloc */
    double *rtn = (double*) malloc(L*sizeof(double));
    if(!rtn) {
        printf("malloc failed");
        exit(1);
    }
    return rtn;
}

void minTimeGradientRV(double *x, double *y, double *z, int Lp, double g0, double gfin, double gmax, double smax, double T, double ds, FILE *fo) {
	
	/* 
	 minTimeGradientRV (file_in, file_out, g0, gfin, gmax, smax, T, ds)
	 minTimeGradientRV (file_in);
	 
	 Finds the time optimal gradient waveforms for the rotationally variant constraints case.
	 
	 The inputes are
	 
	 file_in		-		Input file with curve, Nx3 double with columns [Cx Cy Cz]
	 file_out		-		Output file to print results in
	 g0				-		Initial gradient amplitude.
	 gfin			-		Gradient value at the end of the trajectory.
	 If given value is not possible
	 the result would be the largest possible amplitude.
	 enter -1 for default.
	 gmax			-		Maximum gradient [G/cm] (4 default)
	 smax			-		Maximum slew [G/cm/ms] (15 default)
	 T				-		Sampling time intervale [ms] (4e-3 default)
	 ds				-		Step size for ODE integration. Enter -1 to use default.
	 
	 The following are printed into the output file (each value is in a seperate column)
	 
	 [Cx Cy Cz]	-		reparameterized curve, sampled at T[ms]
	 [gx gy gz]	-		gradient waveforms [G/cm]
	 [sx sy sz]	-		slew rate [G/cm/ms]
	 [kx ky kz]	-		exact k-space corresponding to gradient g
	 sta			-		Solution for the forward ODE
	 stb			-		Solution for the backward ODE
	 phi			-		Geometry constrains on amplitude vs. arclength
	 time		-		time array (total time is time(end)).
	 */
	double dt = T;
	double gamma = 4.257;
    int i;
	
	/* iflag used in spline method to signal error */
    int *iflag;
    int iflagp;
    iflag = &iflagp;
	
	/* Representing the curve with parameter p */
	double p [Lp];
	for (i = 0; i < Lp; i++) {
		p[i] = i;
	}
	
    /* Interpolation of curve for gradient accuracy, using cubic spline interpolation	 */
	
	double *c1x, *c2x, *c3x, *c1y, *c2y, *c3y, *c1z, *c2z, *c3z;	/* arrays used by spline function to store coefficients. */
	
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
	
	double *CCx, *CCy;
	CCx = double_malloc_check(num_evals * sizeof(double));
	CCy = double_malloc_check(num_evals * sizeof(double));
	
	double *CCz;
	CCz = double_malloc_check(num_evals * sizeof(double));
	
	double toeval = 0;
	
	int *last;
	int holder = 0;
	last = &holder;
	
	
	double *Cpx, *Cpy,*Cpz, *Cp_abs;		/* interpolated curve in p-parameterization */
	Cpx = double_malloc_check(num_evals * sizeof(double));
	Cpy =  double_malloc_check(num_evals * sizeof(double));;
	Cp_abs =  double_malloc_check(num_evals * sizeof(double)); 
	
	Cpz = double_malloc_check(num_evals * sizeof(double));
	
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
	
	free(Cpx);	free(Cpy);	free(Cpz);
	
	/* converting to arc-length parameterization from p, using trapezoidal integration */
	
	double *s_of_p;
	s_of_p = double_malloc_check(num_evals * sizeof(double));
	s_of_p[0] = 0;
	
	double sofar = 0;
	
	for (i=0; i < num_evals; i++) {
		
		sofar += (Cp_abs[i]+ Cp_abs[i-1])/2;
		s_of_p[i] =  dp * sofar;
		
	}
	
	free(Cp_abs);
	
	/* length of the curve	*/
	double L = s_of_p[num_evals-1];		
	
	/* decide ds and compute st for the first point */
    double stt0 = gamma*smax;	/* always assumes first point is max slew */
    double st0 = (stt0*dt)/2;	/* start at half the gradient for accuracy close to g=0 */
    double s0 = st0*dt;
    if(ds < 0) {
		ds = s0/1.5;			/* smaller step size for numerical accuracy */
	}
    int length_of_s =  (int) floor(L/ds);
    int half_ls = (int) floor(L/(ds/2));
	
	double *s = (double *) double_malloc_check(length_of_s*sizeof(double));	
	
	for (i = 0; i< length_of_s; i++) {
		s[i] = i*ds;
	}
	
	double *s_half = (double *) double_malloc_check(half_ls*sizeof(double));
	
	for (i=0; i < half_ls; i++) {
		s_half[i] = (double)i*(ds/2);
	}
	
	/* Convert from s(p) to p(s) and interpolate for accuracy */
	double *p_of_s_half;
	p_of_s_half = double_malloc_check(half_ls * sizeof(double));
	
	
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
		p_of_s_half[i] = seval (num_evals, s_half[i], s_of_p, sop_num, a1x, a2x, a3x, last);
	}
	
	free(a1x); free(a2x); free(a3x);
    free(s_of_p);
	
	int size_p_of_s = half_ls/2;
	
	double *p_of_s;
	p_of_s = double_malloc_check(size_p_of_s * sizeof(double)); 
	
	for (i=0; i<size_p_of_s; i++) {
		p_of_s[i] = p_of_s_half[2*i];		
	}
	
	printf("Compute geometry dependent constraints \n");
	
	double *k;
	k = double_malloc_check(half_ls*sizeof(double));	/* k is the curvature along the curve */
	
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
    
    double *Csp1x, *Csp2x, *Csp3x, *Csp1y, *Csp2y, *Csp3y,  *Csp1z, *Csp2z, *Csp3z;
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
    
    double *xpp, *ypp, *zpp;		/* xpp, ypp, zpp are d^2/ds^2 - second derivative in s parameterization */
    xpp = double_malloc_check(half_ls*sizeof(double));
    ypp = double_malloc_check(half_ls*sizeof(double));
    zpp = double_malloc_check(half_ls*sizeof(double));
    
    double *xp, *yp, *zp;		   /* d/ds - first derivative in s parameterization */
    xp = double_malloc_check(half_ls*sizeof(double));
    yp = double_malloc_check(half_ls*sizeof(double));
    zp = double_malloc_check(half_ls*sizeof(double));
    
    /* Compute the curvature along the curve */
    
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
        k[i] =  sqrt(kx*kx + ky*ky + kz*kz);		/* the curvature, magnitude of the second derivative of the curve in arc-length parameterization */
    }
	
	free(CCx);	free(CCy); free(CCz);
	free(p_of_s_half);
	
	double *sdot1, *sdot2, *sdot, *sdot3;
	sdot1 = double_malloc_check(half_ls*sizeof(double));
	sdot2 = double_malloc_check(half_ls*sizeof(double));
	sdot3 = double_malloc_check(half_ls*sizeof(double));
	sdot = double_malloc_check(half_ls*sizeof(double));
	
	for (i=0; i< half_ls; i++) {
        sdot1[i] = gamma*gmax/xp[i];
        sdot2[i] = gamma*gmax/yp[i];
        sdot3[i] = gamma*gmax/zp[i];
        sdot[i] = sdot1[i];
        if (sdot2[i] < sdot[i]) {
            sdot[i] = sdot2[i];
        }
        if (sdot3[i] < sdot[i]){
            sdot[i] = sdot3[i];
        }
    }
	
	free(s_half);
    free(Cspx); free(Cspy); free(Cspz);
    free(Csp1x); free(Csp2x); free(Csp3x);
    free(Csp1y); free(Csp2y); free(Csp3y);
    free(Csp1z); free(Csp2z); free(Csp3z);
    free(sdot1);
    free(sdot2);
    free(sdot3);
	
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
	
	double *sta;		/* Forward solution */
	sta = double_malloc_check(length_of_s*sizeof(double));
	
	double *stb;		/* Backward solution */
	stb = double_malloc_check(length_of_s*sizeof(double));
	
	/* Forward solution */
    /* Initial cond. */
	double g0gamma = g0*gamma + st0;
	double gammagmax = gamma *gmax;	
	if (g0gamma < gammagmax) {
		sta[0] = g0gamma;
	} else {
		sta[0] = gammagmax;
	}
	
	/* Solving ODE Forward */
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
		
		double dstds = RungeKutte(ds, sta[i-1], xpp_rk , ypp_rk, zpp_rk, xp_rk, yp_rk, zp_rk, smax);
		double tmpst = sta[i-1] + dstds;
		
		if (sdot[2*i] < tmpst) {
			sta[i] = sdot[2*i];
		} else {
			sta[i] = tmpst;
		}
	}
	
	free(k2);
    
    printf("Solving ODE Backwards\n");
    
    if (gfin < 0) {  /*if gfin is not provided */
        stb[length_of_s-1] = sta[length_of_s - 1];
    } else {
        double max;
        if (gfin * gamma > st0) {
            max = gfin*gamma;
        } else {
            max = st0;
        }
        if (gamma * gmax < max) {
            stb[length_of_s-1] = gamma * gmax;
        } else {
            stb[length_of_s-1] = max;
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
		
		double dstds = RungeKutte(ds, stb[i+1], xpp_rk, ypp_rk, zpp_rk, xp_rk, yp_rk, zp_rk, smax);
		double tmpst = stb[i+1] + dstds;
		
		if (sdot[2*i] < tmpst) {
			stb[i] = sdot[2*i];
		} else {
			stb[i] = tmpst;
		}
	}
	
	double *st_of_s, *st_ds_i;
	st_of_s = double_malloc_check(length_of_s*sizeof(double));
	st_ds_i = double_malloc_check(length_of_s*sizeof(double));
	
	for (i=0; i<length_of_s; i++) {
		if (sta[i] < stb[i]) {
			st_of_s[i] = sta[i];		/* use the minimum of the two solutions, sta and stb */
		}
		else {
			st_of_s[i] = stb[i];		/* ds * 1/st(s) used in below calculation of t(s) */
		}
		
		st_ds_i[i] = ds*(1/st_of_s[i]);
	}
	
	free(k);
    free(xp); free(yp); free(zp);
    free(xp2); free(yp2); free(zp2);
    free(xpp);    free(ypp);    free(zpp);
    free(xpp2);   free(ypp2);   free(zpp2);
	
	/*Final interpolation */
	printf("Final interpolation\n");
	/* Converting to the time parameterization, t(s) using trapezoidal integration. t(s) = integral (1/st) ds	 */
	double *t_of_s= double_malloc_check(length_of_s*sizeof(double));
	t_of_s[0] = 0;
	for (i=1; i < length_of_s; i++) {
		t_of_s[i] =  t_of_s[i-1] + fabs(st_ds_i[i]+ st_ds_i[i-1])/2;
	}
	
	int l_t =  (int) floor(t_of_s[length_of_s-1]/dt);		/* size of the interpolated trajectory */
	
	double t[l_t];
	for (i=0; i<l_t; i++) {
		t[i] = i*dt;				/* time array */
	}
	
    double *t1x, *t2x, *t3x;	/* coefficient arrays for spline interpolation of t(s) to get s(t)	 */
	t1x = double_malloc_check(length_of_s*sizeof(double));
	t2x = double_malloc_check(length_of_s*sizeof(double));
	t3x = double_malloc_check(length_of_s*sizeof(double));
	
	double *s_of_t;
	s_of_t = double_malloc_check(l_t * sizeof(double));
	
	spline(length_of_s, 0, 0, 1, 1, t_of_s, s, t1x, t2x, t3x, iflag);
	
	for (i=0; i < l_t; i++){
		s_of_t[i] = seval(length_of_s, t[i], t_of_s, s, t1x, t2x, t3x, last);
	} 
	
	free(t1x);	free(t2x);	free(t3x);
	
	free(st_ds_i);
    free(st_of_s);
    free(t_of_s);
	
    double *p1x, *p2x, *p3x;	/* coefficient arrays for spline interpolation of p(s) with s(t) to get p(s(t)) = p(t) */
	p1x = double_malloc_check(length_of_s*sizeof(double));
	p2x = double_malloc_check(length_of_s*sizeof(double));
	p3x = double_malloc_check(length_of_s*sizeof(double));
	
	spline(length_of_s, 0, 0, 1, 1, s, p_of_s, p1x, p2x, p3x, iflag);
	
	double *p_of_t;
	p_of_t = double_malloc_check(l_t*sizeof(double));
	
	for (i=0; i < l_t; i++){
		p_of_t[i] = seval(length_of_s, s_of_t[i], s, p_of_s, p1x, p2x, p3x, last);
		
	}
	
	free(s);
    free(p_of_s);
    free(p1x);	free(p2x);	free(p3x);
    free(s_of_t);
	
	double *Cx, *Cy, *Cz;     /*  Reparameterized curve C, sampled at T, to get C(p(t)) */
	Cx = double_malloc_check(l_t * sizeof(double));
	Cy = double_malloc_check(l_t * sizeof(double));
	Cz = double_malloc_check(l_t * sizeof(double));
	
	for (i=0; i<l_t; i++) {
		
		Cx[i] = seval(Lp, p_of_t[i], p, x, c1x, c2x, c3x, last);
		Cy[i] = seval(Lp, p_of_t[i], p, y, c1y, c2y, c3y, last);
		Cz[i] = seval(Lp, p_of_t[i], p, z, c1z, c2z, c3z, last);
	}
	
	free(p_of_t);
    free(c1x);	free(c2x);	free(c3x);
    free(c1y);	free(c2y);	free(c3y);
    free(c1z);	free(c2z);	free(c3z);;
	
	double *gx, *gy, *gz, *g;     /* Final gradient waveforms to be returned */
	gx = double_malloc_check(l_t * sizeof(double));
	gy = double_malloc_check(l_t * sizeof(double));
	gz = double_malloc_check(l_t * sizeof(double));
	
	for (i=0; i< l_t -1; i++) {
		gx[i] = (Cx[i+1] - Cx[i]) / (gamma * dt);
		gy[i] = (Cy[i+1] - Cy[i]) / (gamma * dt);
		gz[i] = (Cz[i+1] - Cz[i]) / (gamma * dt);
	}
	
	gx[l_t-1] = gx[l_t-2] + gx[l_t-2] - gx[l_t-3];
	gy[l_t-1] = gy[l_t-2] + gy[l_t-2] - gy[l_t-3];
	gz[l_t-1] = gz[l_t-2] + gz[l_t-2] - gz[l_t-3];
	
	
	double *kx, *ky, *kz;     /* exact k-space trajectory corresponding to g  */
	/* (calculated by integrating gradient waveforms) */
	kx = double_malloc_check(l_t * sizeof(double));
	ky = double_malloc_check(l_t * sizeof(double));
	kz = double_malloc_check(l_t * sizeof(double));
	
	double sofarx = 0;
	kx[0] = 0;
	double sofary = 0;
	ky[0] = 0;
	double sofarz = 0;
	kz[0] = 0;
	
	for (i=1; i < l_t; i++) {
		sofarx += (gx[i] - gx[i-1]) / 2;
		sofary += (gy[i] - gy[i-1]) / 2;
		sofarz += (gz[i] - gz[i-1]) / 2;
		kx[i] = sofarx * dt * gamma;
		ky[i] = sofary * dt * gamma;
		kz[i] = sofarz * dt * gamma;
	}
	
	double *sx, *sy, *sz;		    /* slew waveforms to be returned */
	sx = double_malloc_check((l_t) * sizeof(double));
	sy = double_malloc_check((l_t) * sizeof(double));
	sz = double_malloc_check((l_t) * sizeof(double));
	
	for (i=0; i < l_t-1; i++) {
		sx[i] = (gx[i+1] - gx[i])/dt;
		sy[i] = (gy[i+1] - gy[i])/dt;
		sz[i] = (gz[i+1] - gz[i])/dt;
	}
	sx[l_t-1] = sx[l_t-2];
	sy[l_t-1] = sy[l_t-2];
	sz[l_t-1] = sz[l_t-2];
	
	/* total traversal time */
	double time = t[l_t-1];
	
	/*	Printing results to file */
	for (i=0; i < l_t; i++) {
		fprintf(fo, "%f\t%f\t%f	\t%f\t%f\t%f\t %f\t%f\t%f\t %f\t%f\t%f\t %f\t%f\t%f \t%f\n",
				Cx[i], Cy[i], Cz[i],
				gx[i], gy[i], gz[i],
				sx[i], sy[i], sz[i],
				kx[i], ky[i], kz[i],
				sta[i], stb[i], sdot[i], t[i]);
	}
	printf("Time:\t%f\n", time);
	
	free(sdot);	free(sta);	free(stb);
	free(Cx); free(Cy); free(Cz);
	free(kx);free(ky);	free(kz);
	free(sx); free(sy);	free(sz);
	free(gx);	free(gy);	free(gz);
	free(x);	free(y); free(z);

	fclose(fo);
}

void minTimeGradientRIV(double *x, double *y, double *z, int Lp, double g0, double gfin, double gmax, double smax, double T, double ds, FILE *fo) {
    
	/* Finds the time optimal gradient waveforms for the rotationally invariant constraints case.
	 
	 file_in		-		Input file with curve, Nx3 double with columns [Cx Cy Cz]
	 file_out		-		Output file to print results in
	 g0				-		Initial gradient amplitude.
	 gfin			-		Gradient value at the end of the trajectory.
	 If given value is not possible
	 the result would be the largest possible amplitude.
	 enter -1 for default.
	 gmax			-		Maximum gradient [G/cm] (4 default)
	 smax			-		Maximum slew [G/cm/ms] (15 default)
	 T				-		Sampling time intervale [ms] (4e-3 default)
	 ds				-		Step size for ODE integration. Enter -1 to use default.
	 
	 The following are printed into the output file (each value is in a seperate column)
	 
	 [Cx Cy Cz]	-		reparameterized curve, sampled at T[ms]
	 [gx gy gz]	-		gradient waveforms [G/cm]
	 [sx sy sz]	-		slew rate [G/cm/ms]
	 [kx ky kz]	-		exact k-space corresponding to gradient g
	 sta			-		Solution for the forward ODE
	 stb			-		Solution for the backward ODE
	 phi			-		Geometry constrains on amplitude vs. arclength
	 time		-		time array (total time is time(end)).
	 */
	    
    int i=0;
	
	/* iflag used in spline method to signal error */
    int *iflag;
    int iflagp;
    iflag = &iflagp;
	
	/* Representing the curve with parameter p */
	double p [Lp];
	for (i = 0; i < Lp; i++) {
		p[i] = i;
	}
	
	double dt = T;
	double gamma = 4.257;
    
    /* Interpolation of curve for gradient accuracy, using cubic spline interpolation */
    
    double *c1x, *c2x, *c3x,
	*c1y, *c2y, *c3y,
	*c1z, *c2z, *c3z;		/* arrays used by spline function to store coefficients. */
    
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
    
    double *Cpx, *Cpy, *Cp_abs, *Cpz;				/* interpolated curve in p-parameterization */
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
    free(Cpx);	free(Cpy);	free(Cpz);
    
    /* converting to arc-length parameterization from p, using trapezoidal integration */
    
    double *s_of_p;
    s_of_p = double_malloc_check(num_evals * sizeof(double));
    s_of_p[0] = 0;
    
    double sofar = 0;
    
    for (i=0; i < num_evals; i++) {
        
        sofar += (Cp_abs[i]+ Cp_abs[i-1])/2;
        s_of_p[i] =  dp * sofar;
    }
    
    /* length of the curve */
    double L = s_of_p[num_evals-1];
    
    /* decide ds and compute st for the first point */
    double stt0 = gamma*smax;	/* always assumes first point is max slew */
    double st0 = (stt0*dt)/2;	/* start at half the gradient for accuracy close to g=0 */
    double s0 = st0*dt;
	if(ds < 0) {
		ds = s0/1.5;			/* smaller step size for numerical accuracy */
	}
    int length_of_s =  (int) floor(L/ds);
    int half_ls = (int) floor(L/(ds/2));
    
    double *s;
    s = double_malloc_check(length_of_s * sizeof(double));	
	
	double *sta;
	sta = double_malloc_check(length_of_s*sizeof(double));
	
	double *stb;
	stb = double_malloc_check(length_of_s*sizeof(double));
	
	for (i = 0; i< length_of_s; i++) {
		s[i] = i*ds;
	}
	
	
	double *s_half = (double *) double_malloc_check(half_ls*sizeof(double));
	
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
    
    free(a1x); free(a2x); free(a3x);
    
    free(s_of_p);
    
    int size_p_of_s = half_ls/2;
    
    double *p_of_s;
    p_of_s = double_malloc_check(size_p_of_s * sizeof(double));
    
    for (i=0; i<size_p_of_s; i++) {
        p_of_s[i] = p_of_s_half[2*i];
    }
    
    printf("Compute geometry dependent constraints \n");
    
    double *k;
    k = double_malloc_check(half_ls*sizeof(double));		/* k is the curvature along the curve */
    
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
    
    double *Csp1x, *Csp2x, *Csp3x, *Csp1y, *Csp2y, *Csp3y,  *Csp1z, *Csp2z, *Csp3z;
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
        k[i] =  sqrt(kx*kx + ky*ky + kz*kz);		/* the curvature, magnitude of the second derivative of the curve in arc-length parameterization */
    }
    
	free(CCx);	free(CCy);    free(CCz);
    free(p_of_s_half);
    
	
	double *sdot1, *sdot2, *sdot;
	sdot1 = double_malloc_check(half_ls*sizeof(double));
	sdot2 = double_malloc_check(half_ls*sizeof(double));
	sdot = double_malloc_check(half_ls*sizeof(double));
	
	/* Calculating the upper bound for the time parametrization */
    /* sdot (which is a non scaled max gradient constaint) as a function of s. */
    /* sdot is the minimum of gamma*gmax and sqrt(gamma*gmax / k) */
    
	
	for (i=0; i< half_ls; i++) {
		sdot1[i] = gamma*gmax;
		sdot2[i] = sqrt((gamma*smax) / (fabs(k[i]+(DBL_EPSILON))));
		if (sdot1[i] < sdot2[i]) {
			sdot[i] = sdot1[i];
		}else {
			sdot[i] = sdot2[i];
		}
		
	}
	
	free(sdot1);
	free(sdot2);
	free(s_half);
    free(Cspx); free(Cspy); free(Cspz);
    free(Csp1x); free(Csp2x); free(Csp3x);
    free(Csp1y); free(Csp2y); free(Csp3y);
    free(Csp1z); free(Csp2z); free(Csp3z);
	
	int size_k2 = half_ls+2; 	/* extend of k for RK4 */
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
		sta[0] = g0gamma;
	} else {
		sta[0] = gammagmax;
	}
	
	/* Solving ODE Forward */
	printf("Solving ODE Forward\n");
	
	for (i=1; i<length_of_s; i++) {
		double k_rk[3];
		k_rk[0] = k2[2*i-2];
		k_rk[1] = k2[2*i-1];
		k_rk[2] = k2[2*i];
		
		double dstds = RungeKutte_riv(ds, sta[i-1], k_rk, smax);
		
		double tmpst = sta[i-1] + dstds;
		if (sdot[2*i+1] < tmpst) {
			sta[i] = sdot[2*i+1];
		} else {
			sta[i] = tmpst;
		}
	}
	
	free(k2);
	
	/*Solving ODE Backwards: */	
	
	printf("Solving ODE Backwards\n");
	double max;
    if(gfin<0) {
        /*if gfin is not provided */
        stb[length_of_s-1] = sta[length_of_s - 1];
    } else {
        
        if (gfin * gamma > st0) {
            max = gfin*gamma;
        } else {
            max = st0;
        }
        
        if (gamma * gmax < max) {
            stb[length_of_s-1] = gamma * gmax;
        } else {
            stb[length_of_s-1] = max;
        }
    }
	
	for (i=length_of_s-2; i>-1; i--) {
		double k_rk[3];
		k_rk[0] = k[2*i+1];
		k_rk[1] = k[2*i];
		k_rk[2] = k[2*i-1];
		
		double *h0;
		double hstart = ds;
		h0 = &hstart;
		
		double dstds = RungeKutte_riv(ds, stb[i+1], k_rk, smax);
		double tmpst = stb[i+1] + dstds;
		
		if (sdot[2*i] < tmpst) {
			stb[i] = sdot[2*i];
		} else {
			stb[i] = tmpst;
		}
	}
	
	/* take st(s) to be the minimum of the curves sta and stb */
	double *st_of_s, *st_ds_i;
	st_of_s = double_malloc_check(length_of_s*sizeof(double));
	st_ds_i = double_malloc_check(length_of_s*sizeof(double));
	
	for (i=0; i<length_of_s; i++) {
		if (sta[i] < stb[i]) {
			st_of_s[i] = sta[i];
		}
		else {
			st_of_s[i] = stb[i];
		}
		
		st_ds_i[i] = ds*(1/st_of_s[i]);		/* ds * 1/st(s) used in below calculation of t(s) */
		
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
	double t[l_t];
	for (i=0; i<l_t; i++) {
		t[i] = i*dt;		/* time array */
	}
	
	double *t1x, *t2x, *t3x;		/* coefficient arrays for spline interpolation of t(s) to get s(t) */
	
	t1x = double_malloc_check(length_of_s*sizeof(double));
	t2x = double_malloc_check(length_of_s*sizeof(double));
	t3x = double_malloc_check(length_of_s*sizeof(double));
	
	double *s_of_t;
	s_of_t = double_malloc_check(l_t * sizeof(double));
	
	spline(length_of_s, 0, 0, 1, 1, t_of_s, s, t1x, t2x, t3x, iflag);
	
	for (i=0; i < l_t; i++){
		s_of_t[i] = seval(length_of_s, t[i], t_of_s, s, t1x, t2x, t3x, last);
	} 
	
	free(t1x);	free(t2x);	free(t3x);
	
	free(st_ds_i);
	free(st_of_s);
	free(t_of_s);
	
	double *p1x, *p2x, *p3x;		/* coefficient arrays for spline interpolation of p(s) with s(t) to get p(s(t)) = p(t) */
	p1x = double_malloc_check(length_of_s*sizeof(double));
	p2x = double_malloc_check(length_of_s*sizeof(double));
	p3x = double_malloc_check(length_of_s*sizeof(double));
	
	spline(length_of_s, 0, 0, 1, 1, s, p_of_s, p1x, p2x, p3x, iflag);
	
	double *p_of_t;
	p_of_t = double_malloc_check(l_t*sizeof(double));
	
	for (i=0; i < l_t; i++){
		p_of_t[i] = seval(length_of_s, s_of_t[i], s, p_of_s, p1x, p2x, p3x, last);
	}
	
	free(s);
    free(p_of_s);
    free(p1x);	free(p2x);	free(p3x);
    free(s_of_t);
	
	/*  interpolated k-space trajectory */
	
	double *Cx, *Cy, *Cz;
	Cx = double_malloc_check(l_t * sizeof(double));
	Cy = double_malloc_check(l_t * sizeof(double));
	Cz = double_malloc_check(l_t * sizeof(double));
	
	
	for (i=0; i<l_t; i++) {
		
		Cx[i] = seval(Lp, p_of_t[i], p, x, c1x, c2x, c3x, last);
		Cy[i] = seval(Lp, p_of_t[i], p, y, c1y, c2y, c3y, last);
		Cz[i] = seval(Lp, p_of_t[i], p, z, c1z, c2z, c3z, last);
		
	}
	
	free(p_of_t);
    free(c1x);	free(c2x);	free(c3x);
    free(c1y);	free(c2y);	free(c3y);
    free(c1z);	free(c2z);	free(c3z);
	
	double *gx, *gy, *gz;
	gx = double_malloc_check(l_t * sizeof(double));
	gy = double_malloc_check(l_t * sizeof(double));
	gz = double_malloc_check(l_t * sizeof(double));
	
	for (i=0; i< l_t -1; i++) {
		gx[i] = (Cx[i+1] - Cx[i]) / (gamma * dt);
		gy[i] = (Cy[i+1] - Cy[i]) / (gamma * dt);
		gz[i] = (Cz[i+1] - Cz[i]) / (gamma * dt);
	}
	
	gx[l_t-1] = gx[l_t-2] + gx[l_t-2] - gx[l_t-3];
	gy[l_t-1] = gy[l_t-2] + gy[l_t-2] - gy[l_t-3];
	gz[l_t-1] = gz[l_t-2] + gz[l_t-2] - gz[l_t-3];
	
	/* k-space trajecoty to be returned (calculated by integrating gradient waveforms by trapezoidal integration) */
	double *kx, *ky, *kz;
	kx = double_malloc_check(l_t * sizeof(double));
	ky = double_malloc_check(l_t * sizeof(double));
	kz = double_malloc_check(l_t * sizeof(double));
	
	double sofarx = 0;
	kx[0] = 0;
	double sofary = 0;
	ky[0] = 0;
	double sofarz = 0;
	kz[0] = 0;
	
	for (i=1; i < l_t; i++) {
		sofarx += (gx[i] - gx[i-1]) / 2;
		sofary += (gy[i] - gy[i-1]) / 2;
		sofarz += (gz[i] - gz[i-1]) / 2;
		
		kx[i] = sofarx * dt * gamma;
		ky[i] = sofary * dt * gamma;
		kz[i] = sofarz * dt * gamma;
	}
	
	/* slew waveforms to be returned */
	double *sx, *sy, *sz;
	sx = double_malloc_check((l_t - 1) * sizeof(double));
	sy = double_malloc_check((l_t - 1) * sizeof(double));
	sz = double_malloc_check((l_t - 1) * sizeof(double));
	
	for (i=0; i < l_t-1; i++) {
		sx[i] = (gx[i+1] - gx[i])/dt;
		sy[i] = (gy[i+1] - gy[i])/dt;
		sz[i] = (gz[i+1] - gz[i])/dt;
	}
	sx[l_t-1] = sx[l_t-2];
	sy[l_t-1] = sy[l_t-2];
	sz[l_t-1] = sz[l_t-2];
	
	
	/* total traversal time */
	double time = t[l_t-1];
	
	/*	Printing results to file */
	
	for (i=0; i < l_t; i++) {
		fprintf(fo, "%f\t%f\t%f\t %f\t%f\t%f\t %f\t%f\t%f\t %f\t%f\t%f\t %f\t%f\t%f\t %f\n",
				Cx[i], Cy[i], Cz[i],
				gx[i], gy[i], gz[i], 
				sx[i], sy[i], sz[i], 
				kx[i], ky[i], kz[i],
				sta[i], stb[i], sdot[i], t[i]);
	}
	
	printf("Time:\t%f\n", time);
	
	free(sdot);	free(sta);	free(stb);
	free(Cx); free(Cy); free(Cz);
	free(kx);free(ky);	free(kz);
	free(sx); free(sy);	free(sz);
	free(gx);	free(gy);	free(gz);
	free(Cp_abs);
	free(x);	free(y); free(z);	
}
