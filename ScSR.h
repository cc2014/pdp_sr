#ifndef _SCSR_H_
#define _SCSR_H_

#include "stdio.h"
#include "string.h"
#include "math.h"

#if defined(_DEBUG)
#define dprintf(fmt, ...) printf("%s():%d "fmt,__func__,__LINE__,##__VA_ARGS__)
#else
#define dprintf(fmt, ...)
#endif

typedef struct _ParamScSR
{
       double up_scale, lambda, overlap;
       double *Dh, *Dl, *DlTxDl;
	   int Dhw, Dhh, Dlw, Dlh;
       int patch_size, maxIter; 
	   int DlTxDl_h, DlTxDl_w;	
}ParamScSR;

bool ScSR( unsigned char *im_l_y, int &nrow, int &ncol, ParamScSR &strParamScSR );
void resize_image_d( double *src_data, double *dst_data, const int &src_w, const int &src_h, const int &dst_w, const int &dst_h ); 
void L1QP_FeatureSign_yang(const double &lambda, double* A, double* b, double* x_ret, int &d_size);
double sum_of_product( double *vector_a, double *vector_b, int vector_length );

int clust_invert(
  double **a,      /* input/output matrix */
  int    n,        /* dimension */
  double *det_man, /* determinant mantisa */
  int    *det_exp, /* determinant exponent */
  /* scratch space */
  int    *indx,    /* indx = G_alloc_ivector(n);  */
  double **y,      /* y = G_alloc_matrix(n,n); */
  double *col      /* col = G_alloc_vector(n); */
);
#endif

