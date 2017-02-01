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
#define EXCITATION 1
#define SMALLTIP   2
#define SPINECHO   3
#define SATURATION 4
#define INVERSION  5

// Constants
#define GAMMA 4.26e3

// Only functions that should be called from external code
void gen_slr_rf(float *rf1, int numpts, float slice_thick, float rf1_dur,
		    float a_gzrf1, float in_err, float out_err, int type);
    
void ge_export(float *rf1);


#endif /* defined(____slr_design__) */
