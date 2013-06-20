#ifndef _IMG_UTILS_H
#define _IMG_UTILS_H
#include <cstdio>
#include <cstring>
#include <cstdlib>

#if defined(_DEBUG)
#define dprintf(fmt, ...) printf("%s():%d "fmt,__func__,__LINE__,##__VA_ARGS__)
#else
#define dprintf(fmt, ...)
#endif



using namespace std;

typedef unsigned char uint8;

void resize_image_bau( unsigned char *src_data, 
		unsigned char *dst_data, const int &src_w, 
		const int &src_h, const int &dst_w, const int &dst_h );

void resize_image_d( double *src_data, double *dst_data, 
		const int &src_w, const int &src_h, const int &dst_w, 
		const int &dst_h ); 

bool convolve2D(double* in, double* out, int dataSizeX, int dataSizeY,    
                double* kernel, int kernelSizeX, int kernelSizeY)   ;

bool convolve2DSeparable(double* in, double* out_horl, 
		double *out_vert, int dataSizeX, int dataSizeY, 
		double* kernelX, int kSizeX, double* kernelY, int kSizeY);

bool convolve2DSeparableBP(double* in, double* out, int dataSizeX, int dataSizeY,    
                         double* kernelX, int kSizeX, double* kernelY, int kSizeY)   ;



void backprojection (double* im_h, const int &width_h, const int &height_h, 
		uint8* im_l, const int &width_l, const int &height_l, const int &maxIter);

int write_pgm_y(char* filename, int x_dim, int y_dim, unsigned char* image);

#endif //_IMG_UTILS_H
