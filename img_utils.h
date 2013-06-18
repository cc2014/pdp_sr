#ifndef _IMG_UTILS_H
#define _IMG_UTILS_H
#include <cstdio>
#include <cstring>
#include <cstdlib>

using namespace std;

void resize_image_bau( unsigned char *src_data, 
		unsigned char *dst_data, const int &src_w, 
		const int &src_h, const int &dst_w, const int &dst_h );

void resize_image_d( double *src_data, double *dst_data, 
		const int &src_w, const int &src_h, const int &dst_w, 
		const int &dst_h ); 

bool convolve2DSeparable(double* in, double* out_horl, 
		double *out_vert, int dataSizeX, int dataSizeY, 
		double* kernelX, int kSizeX, double* kernelY, int kSizeY);
#endif //_IMG_UTILS_H
