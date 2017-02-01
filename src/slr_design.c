
// ************************************************************************
//
//  SLR_C_LIB: Shinnar-Le Roux Pulse Design C Library
//
//    This library includes functions for the design of optimal slice
//    selective RF pulses using the Shinnar-Le Roux algorithm. As described
//    by Pauly et al [1], scan parameters are transformed into filter design
//    parameters, and the Parks-McClellan algorithm is used to design optimal
//    filter coefficients for the specified transition bands and ripple
//    errors. One can design a variety of different types of RF pulses as well,
//    including PI/2 excitation, small tip, spin-echo, saturation, and
//    inversion pulses.
//
//
//  Usage:
//    Call: gen_slr_rf(rf1, numpts, slice_thick, rf1_dur, a_gzrf1, in_err, out_err, type)
//      rf1         - array in which samples of B1 envelope will be stored (pass in float*,
//                  make sure 2*numpts*sizeof(float) bytes of memory has been allocated!)
//      numpts      - number of points describing RF waveform
//      slice_thick - slice thickness (cm)
//      rf1_dur     - duration of RF pulse (ms)
//      a_gzrf1     - amplitude of slice selective gradient waveorm (Hz/G)
//      type        - type of pulse
//      in_err      - percent error of in-slice ripple (recommended: 0.05)
//      out_err     - percent error of out-of-slice ripple (recommended: 0.05)
//
//    type options (input number or NAME):
//      1 - EXCITATION - pi/2 excitation
//      2 - SMALLTIP - small-tip (<pi/2) excitation
//      3 - SPINECHO - pi spin-echo pulse
//      4 - SATURATION - pi/2 saturation pulse
//      5 - INVERSION - inversion pulse
//
//
//
//  Written by Christopher M. Sandino, 4/15/15
//  Based on rf_tools by Adam Kerr & John Pauly
//  Magnetic Resonance Engineering Laboratory
//  University of Southern California
//
// ************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "slr_design.h"
#include "remez.h" // Parks-McClellan Algorithm
#include "four1.h" // 1-D Fast Fourier Transform


// squared magnitude of complex number
static inline double sqmag(double *arr, int i) {
    return arr[2*i]*arr[2*i] + arr[2*i+1]*arr[2*i+1];
}


void b2a(double *a, double *b, int numpts) {
    
    int i, j; // iteration vars
    double im;
    
    // determine zero-pad length (must be power of two)
    int numpts_zp = 1;
    while (numpts_zp < numpts) {
        numpts_zp *= 2;
    }
    numpts_zp *= 16; // number of complex vals in each array
    
    double *b_ft = calloc(2*numpts_zp, sizeof(double));
    double *a_ft = calloc(2*numpts_zp, sizeof(double));
    double *a_mag = calloc(2*numpts_zp, sizeof(double));
    
    if ( b_ft==NULL || a_ft==NULL || a_mag==NULL ) {
        printf("Error allocating memory");
        return;
    }
    
    memcpy(b_ft, b, 2*numpts*sizeof(double));
    
    four1(b_ft-1, numpts_zp, 1);
    
    // find |a|
    for (i=0; i<numpts_zp; i++) {
        a_mag[2*i] = sqrt(1 - sqmag(b_ft,i));
    }
    
    // now find angle(a) from min phase constraint using Hilbert Transform
    for (i=0; i<numpts_zp; i++) {
        a_ft[2*i] = log(a_mag[2*i]);
    }
    
    four1(a_ft-1, numpts_zp, 1);
    
    for (i=1; i<numpts_zp/2-1; i++) {
        a_ft[2*i] *= 2;
        a_ft[2*i+1] *= 2;
    }
    
    for (i=numpts_zp/2+1; i<numpts_zp; i++) {
        a_ft[2*i] = 0;
        a_ft[2*i+1] = 0;
    }
    
    four1(a_ft-1, numpts_zp, -1);
    
    for (i=0; i<numpts_zp; i++) {
        a_ft[2*i] /= numpts_zp;
        a_ft[2*i+1] /= numpts_zp;
    }
    
    // compute a = mag(a)*exp(j*angle(a))
    for (i=0; i<numpts_zp; i++) {
        im = a_ft[2*i+1];
        a_ft[2*i] = a_mag[2*i]*cos(-im);
        a_ft[2*i+1] = a_mag[2*i]*sin(-im);
    }
    
    four1(a_ft-1, numpts_zp, -1);
    
    for(i=0; i<numpts; i++) {
        j = numpts-i-1;
        a[2*i] = a_ft[2*j]/numpts_zp;
        a[2*i+1] = a_ft[2*j+1]/numpts_zp;
    }
    
    free(b_ft);
    free(a_ft);
    free(a_mag);
}


void inverse_slr(float *rf1, double *a, double *b, int numpts) {
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
        rf1[j] = (float)phi*cos(theta);
        // rf1[2*j+1] = phi*sin(theta); // need to keep imag component?
        
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


void gen_slr_rf(float* rf1, int numpts, float slice_thick, float rf1_dur, float a_gzrf1, float in_err, float out_err, int type) {
    
    // Variable declarations
    int i;
    double bsf, tb, d_inf, w, log_in_err, log_out_err;
    
    // Cayley-Klein coeffs (arrays to hold real and imag vals)
    double *b = malloc(2*numpts*sizeof(double));
    double *a = malloc(2*numpts*sizeof(double));
    
    if ( a==NULL || b==NULL ) {
        printf("Error allocating memory");
        return;
    }
    
    // Determine effective in-slice and out-of-slice errors
    switch (type) {
        case EXCITATION :
            bsf = sqrt(0.5);
            in_err = sqrt(in_err)*bsf;
            out_err = out_err*bsf;
            break;
            
        case SMALLTIP :
            bsf = 1;
            break;
            
        case SPINECHO :
            bsf = 1;
            in_err /= 4;
            out_err = sqrt(out_err);
            break;
            
        case SATURATION:
            bsf = sqrt(0.5);
            in_err /= 2;
            out_err = sqrt(out_err);
            break;
            
        case INVERSION :
            bsf = 1;
            in_err /= 8;
            out_err = sqrt(0.5*out_err);
            break;
            
        default :
            printf("Invalid pulse type selected. Options: 1-EXCITATION, 2-SMALLTIP, 3-SPINECHO, 4-SATURATION, 5-INVERSION\n");
            return;
    }
    
    // Transform scan/slice parameters into filter design parameters
    log_in_err = log10(in_err);
    log_out_err = log10(out_err);
    tb = GAMMA*a_gzrf1*slice_thick*rf1_dur; // time-bandwidth product
    d_inf = (5.309e-3*log_in_err*log_in_err + 7.114e-2*log_in_err + -4.761e-1)*log_out_err + -2.66e-3*log_in_err*log_in_err + -5.941e-1*log_in_err + -4.278e-1; // performance measure
    w = d_inf/tb; // transition width
    
    // Parks-McClellan inputs
    int filttype = BANDPASS;    // type of filter
    int numbands = 2;           // number of bands (stop + pass)
    double des[2] = {1, 0};     // intended band amplitudes
    double weight[2] = {1, in_err/out_err}; // error weights
    double bands[4] = {0, (1-w)*tb/numpts/2, (1+w)*tb/numpts/2, 0.5}; // location of band edges [0 w1, w2, 0.5]
    double *h = malloc(numpts*sizeof(double)); // filter coeffs
    
    // compute filter coeffs using parks-mcclellan
    remez(h, numpts, numbands, bands, des, weight, filttype);
    
    
    if (type == SMALLTIP) {
        // use small-tip approximation and assume fourier transform relationship
        for (i=0; i<numpts; i++) {
            rf1[i] = (float) h[i];
        }
        
        free(a);
        free(b);
        free(h);
        return;
    }
    
    // transform into b, array of cayley-klein polynomial coeffs
    for (i=0; i <numpts; i++){
        b[2*i] = h[i]*bsf;
        b[2*i+1] = 0;
    }
    
    b2a(a, b, numpts); // compute a
    
    inverse_slr(rf1, a, b, numpts); // inverse SLR - (An, Bn)==>(RF samples)
 
    
    // ************************************************************************
    //  Debug - print a bunch of outputs to screen
    // ************************************************************************
    
    // printf("d1: %f, d2: %f, D_inf: %f, TB: %f, Lower: %f, Upper: %f\n", in_err, out_err, d_inf, tb, (1-w)*tb/numpts/2, (1+w)*tb/numpts/2);
    
    // for (i=0; i<100; i++) {
    //     printf("%f\n",h[i]);
    // }
    
    // for (i=0; i<200; i++) {
    //     printf("%f\n",a[i]);
    // }
    
    // for (i=0; i<200; i++) {
    //     printf("%f\n",b[i]);
    // }
    
    // for (i=0; i<100; i++) {
    //     printf("%f\n",rf1[i]);
    // }
    
    free(h);
    free(b);
    free(a);
    
    return;
}

void ge_export(float* rf1) {
    
    return;
}
