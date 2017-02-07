//
//  slr_tools.h
//  
//
//  Created by Christopher Sandino on 4/15/15.
//  Modified (heavily) by Eric Gibbons 2017.02.03
//

#ifndef __SLR_TOOLS_H__
#define __SLR_TOOLS_H__

// Pulse type options
#define EXCITATION 0
#define SMALLTIP   1
#define SPINECHO   2
#define SATURATION 3
#define INVERSION  4

// Constants
#define GAMMA 4.26e3

// Only functions that should be called from external code
void Beta2Alpha(double *a, double *b, int numpts);
void InverseSLR(double *rf1, double *a, double *b, int numpts);
void ab2ex(double *mxy, double *a, double *b, int length);
void ab2inv(double *mxy, double *b, int length);
void ab2se(double *mxy, double *b, int length);
void abrx(double *rf, double *gx, int ns, double *x, int nx,
	  double *alpha, double *beta);
void abrot(double *a, double *b, double x, 
	   double *rf, double *gx, int ns);

#endif /* defined(____slr_tools__) */
