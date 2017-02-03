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

/* float dinf(float d1, float d2) */
/* { */

/*     float a1 = (float)5.309e-3; */
/*     float a2 = (float)7.114e-2; */
/*     float a3 = (float)-4.761e-1; */
/*     float a4 = (float)-2.66e-3; */
/*     float a5 = (float)-5.941e-1; */
/*     float a6 = (float)-4.278e-1; */

/*     float l10d1 = (float)log10(d1); */
/*     float l10d2 = (float)log10(d2); */

/*     float d; */

/*     d=(a1*l10d1*l10d1+a2*l10d1+a3)*l10d2+(a4*l10d1*l10d1+a5*l10d1+a6); */

/*     return d; */

/* } */


/* void gen_slr_rf(float* rf1, int numpts, float tb, int ptype, float in_err, float out_err) */
/* { */


/*     // Variable declarations */
/*     int i; */
/*     double bsf, d_inf, w; */
    
/*     // Cayley-Klein coeffs (arrays to hold real and imag vals) */
/*     double *b = malloc(2*(long unsigned int)numpts*sizeof(double)); */
/*     double *a = malloc(2*(long unsigned int)numpts*sizeof(double)); */
    
/*     if ( a==NULL || b==NULL ) { */
/*         printf("Error allocating memory"); */
/*         return; */
/*     } */
    
/*     // Determine effective in-slice and out-of-slice errors */
/*     switch (ptype) { */
/*         case EXCITATION : */
/*             bsf = (float)sqrt(0.5); */
/*             in_err = (float)(sqrt(in_err)*bsf); */
/*             out_err = (float)(out_err*bsf); */
/*             break; */
            
/*         case SMALLTIP : */
/*             bsf = 1; */
/*             break; */
            
/*         case SPINECHO : */
/*             bsf = 1; */
/*             in_err /= 4; */
/*             out_err = (float)sqrt(out_err); */
/*             break; */
            
/*         case SATURATION: */
/*             bsf = (float)sqrt(0.5); */
/*             in_err /= 2; */
/*             out_err = (float)sqrt(out_err); */
/*             break; */
            
/*         case INVERSION : */
/*             bsf = 1; */
/*             in_err /= 8; */
/*             out_err = (float)sqrt(0.5*out_err); */
/*             break; */
            
/*         default : */
/*             printf("Invalid pulse type selected. Options: 1-EXCITATION, 2-SMALLTIP, 3-SPINECHO, 4-SATURATION, 5-INVERSION\n"); */
/*             return; */
/*     } */
    
/*     // Transform scan/slice parameters into filter design parameters */
/*     d_inf = dinf(in_err, out_err); */
/*     w = d_inf/tb; // transition width */
    
/*     // Parks-McClellan inputs */
/*     int filttype = BANDPASS;    // type of filter */
/*     int numbands = 2;           // number of bands (stop + pass) */
/*     double des[2] = {1, 0};     // intended band amplitudes */
/*     double weight[2] = {1, in_err/out_err}; // error weights */
/*     double bands[4] = {0, (1-w)*tb/numpts/2, (1+w)*tb/numpts/2, 0.5}; // location of band edges [0 w1, w2, 0.5] */
/*     double *h = malloc((long unsigned int)numpts*sizeof(double)); // filter coeffs */
    
/*     // compute filter coeffs using parks-mcclellan */
/*     remez(h, numpts, numbands, bands, des, weight, filttype); */
    
    
/*     if (ptype == SMALLTIP) { */
/*         // use small-tip approximation and assume fourier transform relationship */
/*         for (i=0; i<numpts; i++) { */
/*             rf1[i] = (float) h[i]; */
/*         } */
        
/*         free(a); */
/*         free(b); */
/*         free(h); */
/*         return; */
/*     } */
    
/*     // transform into b, array of cayley-klein polynomial coeffs */
/*     for (i=0; i <numpts; i++){ */
/*         b[2*i] = h[i]*bsf; */
/*         b[2*i+1] = 0; */
/*     } */
    
/*     b2a(a, b, numpts); // compute a */
    
/*     inverse_slr(rf1, a, b, numpts); // inverse SLR - (An, Bn)==>(RF samples) */
 
    
/*     // ************************************************************************ */
/*     //  Debug - print a bunch of outputs to screen */
/*     // ************************************************************************ */
    
/*     // printf("d1: %f, d2: %f, D_inf: %f, TB: %f, Lower: %f, Upper: %f\n", in_err, out_err, d_inf, tb, (1-w)*tb/numpts/2, (1+w)*tb/numpts/2); */
    
/*     // for (i=0; i<100; i++) { */
/*     //     printf("%f\n",h[i]); */
/*     // } */
    
/*     // for (i=0; i<200; i++) { */
/*     //     printf("%f\n",a[i]); */
/*     // } */
    
/*     // for (i=0; i<200; i++) { */
/*     //     printf("%f\n",b[i]); */
/*     // } */
    
/*     // for (i=0; i<100; i++) { */
/*     //     printf("%f\n",rf1[i]); */
/*     // } */
    
/*     free(h); */
/*     free(b); */
/*     free(a); */
    
/*     return; */
/* } */

/* void ge_export(float* rf1) { */
    
/*     return; */
/* } */
