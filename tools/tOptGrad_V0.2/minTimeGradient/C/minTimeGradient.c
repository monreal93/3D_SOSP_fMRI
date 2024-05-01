 /* 
 minTimeGradient (file_in, file_out, 1, g0, gfin, gmax, smax, T, ds)
 minTimeGradient (file_in;
 minTimeGradient(file_in, file_out, 0, g0, -1, gmax, smax, T, -1)
 minTimeGradient (file_in, file_out, 1);


 Finds the time optimal gradient waveforms for the rotationally variant constraints case.

 file_in		-		Input file with curve, Nx3 double with columns [Cx Cy Cz]
 file_out		-		Output file to print results in
 RV/RIV			-		Type of solution. 0 for rotationally invariant, and 1 for rotationally variant.
						Default is rotationally invariant.
 g0				-		Initial gradient amplitude.
 gfin			-		Gradient value at the end of the trajectory.
						If given value is not possible
						the result would be the largest possible amplitude. Enter -1 to use default.
 gmax			-		Maximum gradient [G/cm] (4 default)
 smax			-		Maximum slew [G/cm/ms] (15 default)
 T				-		Sampling time intervale [ms] (4e-3 default)
 ds				-		step size for integration. Enter -1 for defaults.

 The following are printed into the output file (each value is in a seperate column)

	[Cx Cy Cz]	-		reparameterized curve, sampled at T[ms]
	[gx gy gz]	-		gradient waveforms [G/cm]
	[sx sy sz]	-		slew rate [G/cm/ms]
	[kx ky kz]	-		exact k-space corresponding to gradient g
	sta			-		Solution for the forward ODE
	stb			-		Solution for the backward ODE
	phi			-		Geometry constrains on amplitude vs. arclength
	time		-		total time to traverse trajectory */


#include <util.h>
#include <float.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "spline.c"
#include <time.h>
#include <sys/types.h>
#include "mtg_functions.c"

int main(int argc, char **argv) {
	
	double g0, gmax, smax, T, gfin, ds;
	double rv;
	FILE *fi, *fo;				/* input and output files */
	if(argc < 2){
		printf("You gotta give me a curve!");
		exit(1);
	}
	if(argc < 3){
		fi = fopen(argv[1], "r");
		fo = fopen("out.dat", "w");
	} else {
		fi = fopen(argv[1], "r");
		fo = fopen(argv[2], "w");
	}
	if (argc < 4) {
		rv = 0;
	} else {
		rv = atof(argv[3]);
	}

	if (argc < 5) {
		g0 = 0;
	} else {
		g0 = atof(argv[4]);
	}

	if (argc < 6) {
		gfin = -1;	
	} else {
		gfin = atof(argv[5]);
	}
	if(argc < 7) {
		gmax = 4;
	}else {
		gmax = atof(argv[6]);
	}

	if(argc < 8) {
		smax = 15;
	} else {
		smax = atof(argv[7]);
	}

	if (argc < 9) {
		T = 4e-3;
	} else {
		T = atof(argv[8]);
	}
	
	if (argc < 10) {
		ds = -1;
	}else {
		ds = atof(argv[9]);
	}
	
	clock_t start, end;
	double elapsed;
	start = clock();
	
	int i = 0;
	/* Reading the inputed curve from input file */
	double *x, *y, *z;
	int sxy = 5000;
	x = double_malloc_check(sxy * sizeof(double));
	y = double_malloc_check(sxy * sizeof(double));
	z = double_malloc_check(sxy * sizeof(double));
	while (!feof(fi))	{
		if (i>=sxy) {
			sxy = 2*sxy;
			x = (double*)realloc(x, sxy*sizeof(double));
			y = (double*)realloc(y, sxy*sizeof(double));
			z = (double*)realloc(z, sxy*sizeof(double));
		}
		fscanf(fi, "%lf %lf %lf", &x[i], &y[i], &z[i]);
		i++;
	} 
	int Lp = i-1;
 
	/* Determining the solution */
	
	if(rv == 1) {
		printf("Computing rotationally variant solution\n");
		minTimeGradientRV(x, y, z, Lp, g0, gfin, gmax, smax, T, ds, fo);
	} else {
		printf("Computing rotationally invariant solution\n");
		minTimeGradientRIV(x, y, z,Lp, g0, gfin, gmax, smax, T, ds, fo);
	}
	
	fclose(fi);
	fclose(fo);
}

