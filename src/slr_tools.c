// ************************************************************************
//
//  SLR_C_LIB: Shinnar-Le Roux Pulse Design C Library
//
//
//  Written by Christopher M. Sandino, 4/15/15
//  Modified by Eric Gibbons, 2017.02.03
//  Based on rf_tools by Adam Kerr & John Pauly
//  Magnetic Resonance Engineering Research Laboratory
//  Stanford University
//
// ************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "slr_tools.h"
#include "four1.h" // 1-D Fast Fourier Transform

#define MP (0.0000001)

// squared magnitude of complex number
static inline double sqmag(double *arr, int i) {
    return arr[2*i]*arr[2*i] + arr[2*i+1]*arr[2*i+1];
}


void Beta2Alpha(double *a, double *b, int numpts)
{
    
    int i, j; // iteration vars
    double bMax, bMag, im;
    
    // determine zero-pad length (must be power of two)
    int numpts_zp = 1;
    while (numpts_zp < numpts) {
        numpts_zp *= 2;
    }
    numpts_zp *= 16; // number of complex vals in each array
    
    double *b_ft = calloc(2*(size_t)numpts_zp, sizeof(double));
    double *a_ft = calloc(2*(size_t)numpts_zp, sizeof(double));
    double *a_mag = calloc(2*(size_t)numpts_zp, sizeof(double));
    
    if ( b_ft==NULL || a_ft==NULL || a_mag==NULL )
    {
        printf("Error allocating memory");
        return;
    }
    
    memcpy(b_ft, b, 2*(long unsigned int)numpts*sizeof(double));
    
    four1(b_ft-1, numpts_zp, 1);

    for (i = 0, bMax = 0; i < numpts_zp; i++)
    {
    	bMag = sqmag(b_ft,i);
    	if (bMag > bMax)
    	{
    	    bMax = bMag;
    	}
    }

    if (bMax >= 1.0)
    {
    	for (i = 0; i < numpts_zp; i++)
    	{
    	    b_ft[2*i] /= (sqrt(bMax) + MP);
    	    b_ft[2*i+1] /= (sqrt(bMax) + MP);
    	}
    }


    // find |a|
    for (i=0; i<numpts_zp; i++)
    {
        a_mag[2*i] = sqrt(1 - sqmag(b_ft,i));
    }

    
    // now find angle(a) from min phase constraint using Hilbert Transform
    for (i=0; i<numpts_zp; i++)
    {
        a_ft[2*i] = log(a_mag[2*i]);
	a_ft[2*i+1] = 0;
    }
    
    four1(a_ft-1, numpts_zp, 1);

    
    for (i=1; i<numpts_zp/2-1; i++)
    {
        a_ft[2*i] *= 2;
        a_ft[2*i+1] *= 2;
    }
    
    for (i=numpts_zp/2+1; i<numpts_zp; i++)
    {
        a_ft[2*i] = 0;
        a_ft[2*i+1] = 0;
    }

    
    four1(a_ft-1, numpts_zp, -1);
    
    for (i=0; i<numpts_zp; i++)
    {
        a_ft[2*i] /= numpts_zp;
        a_ft[2*i+1] /= numpts_zp;
    }
    
    // compute a = mag(a)*exp(j*angle(a))
    for (i=0; i<numpts_zp; i++)
    {
        im = a_ft[2*i+1];
        a_ft[2*i] = a_mag[2*i]*cos(-im);
        a_ft[2*i+1] = a_mag[2*i]*sin(-im);
    }
    
    four1(a_ft-1, numpts_zp, -1);
    
    for(i=0; i<numpts; i++)
    {
        j = numpts-i-1;
        a[2*i] = a_ft[2*j]/numpts_zp;
        a[2*i+1] = a_ft[2*j+1]/numpts_zp;
    }

    /* for (i = 0; i < numpts; i++) */
    /* { */
    /* 	printf("a[%i] = %f + j%f\n",i,a[2*i],a[2*i+1]); */
    /* 	/\* printf("a_mag[%i] = %f\n",i,a_mag[2*i]); *\/ */
    /* 	/\* printf("b_mag[%i] = %f\n",i,sqmag(b_ft,i)); *\/ */
    /* } */
  
    
    free(b_ft);
    free(a_ft);
    free(a_mag);
}


void InverseSLR(double *rf1, double *a, double *b, int numpts)
{
    int j, k, kp; // iteration vars
    double cj, phi, theta, a_mag, a_re, a_im, b_re, b_im;
    double sj[2];
    
    for (j=numpts-1; j>-1; j--) {
        a_mag = sqmag(a,j);
        cj = sqrt(1/(1 + sqmag(b,j)/a_mag));
        sj[0] = (b[2*j]*a[2*j] + b[2*j+1]*a[2*j+1])/a_mag;
        sj[1] = -1*(b[2*j+1]*a[2*j] - b[2*j]*a[2*j+1])/a_mag;
        phi = 2*atan2(sqrt(sqmag(sj,0)), cj);
        theta = atan2(sj[1], sj[0]);
        
        // compute rf samples
        rf1[2*j] = phi*(double)cos(theta);
        rf1[2*j+1] = phi*sin(theta);
        
        if (j == 0) break;
        
        // compute a and b of next RF sample
        for (k=0; k<j; k++) {
            kp = k+1;
            a_re = cj*a[2*kp] + sj[0]*b[2*kp] - sj[1]*b[2*kp+1];
            a_im = cj*a[2*kp+1] + sj[0]*b[2*kp+1] + sj[1]*b[2*kp];
            b_re = cj*b[2*k] - sj[0]*a[2*k] - sj[1]*a[2*k+1];
            b_im = cj*b[2*k+1] - sj[0]*a[2*k+1] + sj[1]*a[2*k];
            
            a[2*k] = a_re;
            a[2*k+1] = a_im;
            b[2*k] = b_re;
            b[2*k+1] = b_im;
        }
    }
}

void ab2ex(double *mxy, double *a, double *b, int length)
{
    int ii;
    
    for (ii = 0; ii < length; ii++)
    {
	mxy[2*ii] = 2*(a[2*ii]*b[2*ii] + a[2*ii+1]*b[2*ii+1]);
	mxy[2*ii+1] = 2*(-a[2*ii]*b[2*ii+1] + a[2*ii+1]*b[2*ii]);
    }
}

void ab2inv(double *mxy, double *b, int length)
{
    int ii;
    
    for (ii = 0; ii < length; ii++)
    {
	mxy[2*ii] = 1 - 2*(b[2*ii]*b[2*ii] + b[2*ii+1]*b[2*ii+1]);
	mxy[2*ii+1] = 0;
    }
}

void ab2se(double *mxy, double *b, int length)
{
    int ii;
    
    for (ii = 0; ii < length; ii++)
    {
	mxy[2*ii] = -(b[2*ii]*b[2*ii+1] + b[2*ii+1]*b[2*ii]);
	mxy[2*ii+1] = (b[2*ii]*b[2*ii] - b[2*ii+1]*b[2*ii+1]);
    }
}


void abrx(double *rf, double *gx, int ns, double *x, int nx,
	  double *alpha, double *beta) 
{
    double alf[2], bet[2];
    double xLocation;
    int ix;
    
    for (ix = 0; ix < nx; ix++) 
    {
	alf[0] = 1.0; 
	alf[1] = 0.0; 
	bet[0] = 0.0; 
	bet[1] = 0.0;
	xLocation = x[ix];
	abrot(alf, bet, xLocation, rf, gx, ns);
	alpha[2*ix] = alf[0];
	alpha[2*ix+1] = alf[1];
	beta[2*ix] = bet[0];
	beta[2*ix+1] = bet[1];
    }
}

void abrot(double *a, double *b, double x, 
	   double *rf, double *gx, int ns)
{
    
    double nx, ny, nz, snp, csp, cg, cpr, cpi, phi;
    double al[2], be[2], ap[2], bp[2];
    int k;
    
    for (k = 0; k < ns; k++)
    {
	cg = x*gx[k];
	cpr = rf[2*k];
	cpi = rf[2*k + 1];

	phi = sqrt(cg*cg+cpr*cpr+cpi*cpi);

	if (phi > 0.0)
	{
	    nx = cpr/phi; 
	    ny = cpi/phi; 
	    nz = cg/phi;
	} 
	else {
	    nx = 0.0; 
	    ny = 0.0; 
	    nz = 1.0;   /* doesn't matter, phi=0*/
	}

	csp = cos(phi/2); 
	snp = sin(phi/2);
	al[0] = csp; 
	al[1] = nz*snp;
	be[0] = ny*snp; 
	be[1] = nx*snp;
	
	bp[0] = al[0]*b[0]-al[1]*b[1]+be[0]*a[0]-be[1]*(-a[1]);
	bp[1] = al[0]*b[1]+al[1]*b[0]+be[1]*a[0]+be[0]*(-a[1]);
	    
	ap[0] = -(be[0]*b[0]-(-be[1])*b[1]) 
	    + al[0] *a[0]-(-al[1])*(-a[1]);
	ap[1] = -(-(-(be[1])*b[0]+  be[0] *b[1]) 
		  + (-al[1])*a[0]+  al[0] *(-a[1]));
	
	a[0] = ap[0]; 
	a[1] = ap[1]; 
	b[0] = bp[0]; 
	b[1] = bp[1];
    }
    
    return;
}


