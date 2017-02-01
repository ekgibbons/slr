//
//  slr_design.h
//  
//
//  Created by Christopher Sandino on 4/15/15.
//
//

#ifndef __SLR_DESIGN_H__
#define __SLR_DESIGN_H__

// Pulse type options
#define EXCITATION 0
#define SMALLTIP   1
#define SPINECHO   2
#define SATURATION 3
#define INVERSION  4

// Constants
#define GAMMA 4.26e3

// Only functions that should be called from external code
void gen_slr_rf(float* rf1, int numpts, float tb, int ptype,
		float in_err, float out_err);
void ge_export(float *rf1);


#endif /* defined(____slr_design__) */
