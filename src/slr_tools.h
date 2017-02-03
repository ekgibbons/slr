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

#endif /* defined(____slr_tools__) */
