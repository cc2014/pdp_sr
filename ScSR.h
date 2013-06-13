#ifndef _SCSR_H_
#define _SCSR_H_

#include "stdio.h"
#include "string.h"
#include "math.h"

typedef struct _ParamScSR
{
       double up_scale, lambda, overlap;
       double *Dh, *Dl, *DlTxDl;
	   int Dhw, Dhh, Dlw, Dlh;
       int patch_size, maxIter; 

}ParamScSR;

bool ScSR( unsigned char *im_l_y, int &nrow, int &ncol, ParamScSR &strParamScSR );

#endif

