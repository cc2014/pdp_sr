#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>

#if defined(_DEBUG)
#define dprintf(fmt, ...) printf("%s():%d "fmt,__func__,__LINE__,##__VA_ARGS__)
#else
#define dprintf(fmt, ...)
#endif

typedef unsigned char uint8;

using namespace std;

double p[25] = {
 0.0001,  0.0021,  0.0058,  0.0021,  0.0001,  
 0.0021,  0.0431,  0.1171,  0.0431,  0.0021,  
 0.0058,  0.1171,  0.3183,  0.1171,  0.0058,  
 0.0021,  0.0431,  0.1171,  0.0431,  0.0021,  
 0.0001,  0.0021,  0.0058,  0.0021,  0.0001 
};

/* create gaussian filter.
 * ====================================================== *
void createFilter(double* gKernel, int size, double sigma)
{
    // set standard deviation to 1.0
    //double sigma = 1.0;
    double r, s = 2.0 * sigma * sigma;
 
    // sum is for normalization
    double sum = 0.0;
	int half=size/2; 

    // generate size*size kernel
    for (int x = 0; x < size; x++)
    {
        for(int y = 0; y < size; y++)
        {
			int x_i=(x-half);
			int y_i=(y-half);
            r = sqrt(x_i*x_i + y_i*y_i);
            gKernel[x*size+y] = (exp(-(r*r)/s))/(M_PI * s);
            sum += gKernel[x*size+y];
        }
    }
 
    // normalize the Kernel
    for(int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j)
            gKernel[i*size+j] /= sum;
 
}
 * ====================================================== *
*/

bool convolve2DSeparable(double* in, double* out, int dataSizeX, int dataSizeY,    
                         float* kernelX, int kSizeX, float* kernelY, int kSizeY)   
{
    int i, j, k, m, n;   
    float *tmp, *sum;                               // intermediate data buffer   
    double *inPtr, *outPtr;                  // working pointers   
    float *tmpPtr, *tmpPtr2;                        // working pointers   
    int kCenter, kOffset, endIndex;                 // kernel indice   
   
    // check validity of params   
    if(!in || !out || !kernelX || !kernelY) return false;   
    if(dataSizeX <= 0 || kSizeX <= 0) return false;   
   
    // allocate temp storage to keep intermediate result   
    tmp = new float[dataSizeX * dataSizeY];   
    if(!tmp)
	{
		delete [] tmp;
		return false;  // memory allocation error   
	}
   
    // store accumulated sum   
    sum = new float[dataSizeX];   
    if(!sum) 
	{
		delete[] sum;
		return false;  // memory allocation error   
	}
   
    // covolve horizontal direction ///////////////////////   
   
    // find center position of kernel (half of kernel size)   
    kCenter = kSizeX >> 1;                          // center index of kernel array   
    endIndex = dataSizeX - kCenter;                 // index for full kernel convolution   
   
    // init working pointers   
    inPtr = in;   
    tmpPtr = tmp;                                   // store intermediate results from 1D horizontal convolution   
   
    // start horizontal convolution (x-direction)   
    for(i=0; i < dataSizeY; ++i)                    // number of rows   
    {   
   
        kOffset = 0;                                // starting index of partial kernel varies for each sample   
   
        // COLUMN FROM index=0 TO index=kCenter-1   
        for(j=0; j < kCenter; ++j)   
        {   
            *tmpPtr = 0;                            // init to 0 before accumulation   
   
            for(k = kCenter + kOffset, m = 0; k >= 0; --k, ++m) // convolve with partial of kernel   
            {   
                *tmpPtr += *(inPtr + m) * kernelX[k];   
            }   
            ++tmpPtr;                               // next output   
            ++kOffset;                              // increase starting index of kernel   
        }   
   
        // COLUMN FROM index=kCenter TO index=(dataSizeX-kCenter-1)   
        for(j = kCenter; j < endIndex; ++j)   
        {   
            *tmpPtr = 0;                            // init to 0 before accumulate   
   
            for(k = kSizeX-1, m = 0; k >= 0; --k, ++m)  // full kernel   
            {   
                *tmpPtr += *(inPtr + m) * kernelX[k];   
            }   
            ++inPtr;                                // next input   
            ++tmpPtr;                               // next output   
        }   
   
        kOffset = 1;                                // ending index of partial kernel varies for each sample   
   
        // COLUMN FROM index=(dataSizeX-kCenter) TO index=(dataSizeX-1)   
        for(j = endIndex; j < dataSizeX; ++j)   
        {   
            *tmpPtr = 0;                            // init to 0 before accumulation   
   
            for(k = kSizeX-1, m=0; k >= kOffset; --k, ++m)   // convolve with partial of kernel   
            {   
                *tmpPtr += *(inPtr + m) * kernelX[k];   
            }   
            ++inPtr;                                // next input   
            ++tmpPtr;                               // next output   
            ++kOffset;                              // increase ending index of partial kernel   
        }   
   
        inPtr += kCenter;                           // next row   
    }   
    // END OF HORIZONTAL CONVOLUTION //////////////////////   
   
    // start vertical direction ///////////////////////////   
   
    // find center position of kernel (half of kernel size)   
    kCenter = kSizeY >> 1;                          // center index of vertical kernel   
    endIndex = dataSizeY - kCenter;                 // index where full kernel convolution should stop   
   
    // set working pointers   
    tmpPtr = tmpPtr2 = tmp;   
    outPtr = out;   
   
    // clear out array before accumulation   
    for(i = 0; i < dataSizeX; ++i)   
        sum[i] = 0;   
   
    // start to convolve vertical direction (y-direction)   
   
    // ROW FROM index=0 TO index=(kCenter-1)   
    kOffset = 0;                                    // starting index of partial kernel varies for each sample   
    for(i=0; i < kCenter; ++i)   
    {   
        for(k = kCenter + kOffset; k >= 0; --k)     // convolve with partial kernel   
        {   
            for(j=0; j < dataSizeX; ++j)   
            {   
                sum[j] += *tmpPtr * kernelY[k];   
                ++tmpPtr;   
            }   
        }   
   
        for(n = 0; n < dataSizeX; ++n)              // convert and copy from sum to out   
        {   
            // covert negative to positive   
            *outPtr = (unsigned char)((float)fabs(sum[n]) + 0.5f);   
            sum[n] = 0;                             // reset to zero for next summing   
            ++outPtr;                               // next element of output   
        }   
   
        tmpPtr = tmpPtr2;                           // reset input pointer   
        ++kOffset;                                  // increase starting index of kernel   
    }   
   
    // ROW FROM index=kCenter TO index=(dataSizeY-kCenter-1)   
    for(i = kCenter; i < endIndex; ++i)   
    {   
        for(k = kSizeY -1; k >= 0; --k)             // convolve with full kernel   
        {   
            for(j = 0; j < dataSizeX; ++j)   
            {   
                sum[j] += *tmpPtr * kernelY[k];   
                ++tmpPtr;   
            }   
        }   
   
        for(n = 0; n < dataSizeX; ++n)              // convert and copy from sum to out   
        {   
            // covert negative to positive   
            *outPtr = (unsigned char)((float)fabs(sum[n]) + 0.5f);   
            sum[n] = 0;                             // reset for next summing   
            ++outPtr;                               // next output   
        }   
   
        // move to next row   
        tmpPtr2 += dataSizeX;   
        tmpPtr = tmpPtr2;   
    }   
   
    // ROW FROM index=(dataSizeY-kCenter) TO index=(dataSizeY-1)   
    kOffset = 1;                                    // ending index of partial kernel varies for each sample   
    for(i=endIndex; i < dataSizeY; ++i)   
    {   
        for(k = kSizeY-1; k >= kOffset; --k)        // convolve with partial kernel   
        {   
            for(j=0; j < dataSizeX; ++j)   
            {   
                sum[j] += *tmpPtr * kernelY[k];   
                ++tmpPtr;   
            }   
        }   
   
        for(n = 0; n < dataSizeX; ++n)              // convert and copy from sum to out   
        {   
            // covert negative to positive   
            *outPtr = (unsigned char)((float)fabs(sum[n]) + 0.5f);   
            sum[n] = 0;                             // reset for next summing   
            ++outPtr;                               // next output   
        }   
   
        // move to next row   
        tmpPtr2 += dataSizeX;   
        tmpPtr = tmpPtr2;                           // next input   
        ++kOffset;                                  // increase ending index of kernel   
    }   
    // END OF VERTICAL CONVOLUTION ////////////////////////   
   
    // deallocate temp buffers   
    delete [] tmp;   
    delete [] sum;   
    return true;   
} 

void resize_image_bau( double *src_data, double *dst_data, const int src_w, const int src_h, const int dst_w, const int dst_h )
{
 
 
    double scalex, scaley;
    double sr, sc, ratior, ratioc, value1, value2, b1, b2; 
    int isr, isc;
    double *imgp;
    int dd, t, stepr, stepc, ii;//, b1, b2;
 
    //unsigned char *data = (unsigned char*)malloc( dst_w*dst_h );
    //scalex = (double)dst_w/(double)src_w;
    //scaley = (double)dst_h/(double)src_h;
    scalex = (double)src_w/(double)dst_w;
    scaley = (double)src_h/(double)dst_h;
     
    b1 = (double)(src_w-1)/scalex;
    b2 = (double)(src_h-1)/scaley;
    for (int dr=0; dr<b2; dr++)
    {
        dd = dr*dst_w;
        sr = dr*scaley;
        isr = (int)sr;
        ratior = (double)sr-isr;
        ii = isr*src_w;
        for (int dc=0; dc<b1; dc++)
        {         
            sc = dc*scalex;
            isc = (int)sc;
            ratioc = sc-isc;
            imgp = src_data + ii+isc;
            value1 = *imgp*(1-ratioc) + *(imgp+1)*ratioc;
            imgp += src_w;
            value2 = *imgp*(1-ratioc) + *(imgp+1)*ratioc;
            dst_data[ dd + dc ] = (double)(value1*(1-ratior)+value2*ratior);
        }
 
        for (int dc=b1; dc<dst_w; dc++)
        {                 
            sc = dc*scalex;
            isc = (int)sc;
            ratioc = sc-isc;
            imgp = src_data + isr*src_w+isc;
            value1 = *imgp;
            imgp += src_w;
            value2 = *imgp;
            dst_data[ dd + dc ] = (double)(value1*(1-ratior)+value2*ratior);
        }
    }
 
    for (int dr=b2; dr<dst_h; dr++)
    {
        dd = dr*dst_w;
        sr = dr*scaley;
        isr = (int)sr;
        ratior = sr-isr;
        for (int dc=0; dc<b1; dc++)
        {         
            sc = dc*scalex;
            isc = (int)sc;
            ratioc = sc-isc;
            imgp = src_data + isr*src_w+isc;
            value1 = *imgp*(1-ratioc) + *(imgp+1)*ratioc;
            value2=value1;
            dst_data[ dd + dc ] = (double)(value1*(1-ratior)+value2*ratior);
        }
 
        for (int dc=b1; dc<dst_w; dc++)
        {         
            sc = dc*scalex;
            isc = (int)sc;
            ratioc = sc-isc;
            imgp = src_data + isr*src_w+isc;
            value1 = *imgp;
            value2 = value1;
            dst_data[ dd + dc ] = (double)(value1*(1-ratior)+value2*ratior);
        }
 
    }

}

/* back projection function of matlab code
 * ====================================================== *
function [im_h] = backprojection(im_h, im_l, maxIter)

[row_l, col_l] = size(im_l);
[row_h, col_h] = size(im_h);

p = fspecial('gaussian', 5, 1);
p = p.^2;
p = p./sum(p(:));

im_l = double(im_l);
im_h = double(im_h);

for ii = 1:maxIter,
    im_l_s = imresize(im_h, [row_l, col_l], 'bicubic');
    im_diff = im_l - im_l_s;
    
    im_diff = imresize(im_diff, [row_h, col_h], 'bicubic');
    im_h = im_h + conv2(im_diff, p, 'same');
end
 * ====================================================== *
*/   


void backprojection (uint8* im_h, int width_h, int height_h, uint8* im_l, int width_l, int height_l, int maxIter, double* im_ret)
{
	/* create gaussian filter. 
	 * ====================================================== *
	int g_size=5;
	double sigma=1.0;
	double gKernel[g_size*g_size];
	double sum=0;

	createFilter(gKernel, g_size, sigma);

    for(int i = 0; i < g_size; ++i)
        for(int j = 0; j < g_size; ++j)
			gKernel[i*g_size+j]=pow(gKernel[i*g_size+j], 2);

    for(int i = 0; i < g_size; ++i)
        for(int j = 0; j < g_size; ++j)
			sum+=gKernel[i*g_size+j];

    for(int i = 0; i < g_size; ++i)
        for(int j = 0; j < g_size; ++j)
			gKernel[i*g_size+j]=gKernel[i*g_size+j]/sum;

	* ====================================================== *
	*/

	double* im_l_s = (double*)malloc(width_l*height_l*sizeof(double));
	double* im_diff = (double*)malloc(width_l*height_l*sizeof(double));
	double* im_diff_h = (double*)malloc(width_h*height_h*sizeof(double));
	double* im_diff_conv = (double*)malloc(width_h*height_h*sizeof(double));

    float kernelX[5] = { 1/16.0f,  4/16.0f,  6/16.0f,  4/16.0f, 1/16.0f };   
    float kernelY[5] = { 1/16.0f,  4/16.0f,  6/16.0f,  4/16.0f, 1/16.0f };   

	double* im_l_d = (double*)malloc(width_l*height_l*sizeof(double));

	/* convert uint8 to double */
	for(int i=0; i<width_h*height_h; i++)
	{
		im_ret[i]=(double)im_h[i];
	}

	for(int i=0; i<width_l*height_l; i++)
	{
		im_l_d[i]=(double)im_l[i];
	}

	/*
	 * first downscale im_h to im_l_s, then calculate diff between im_l and im_l_s, 
	 * upscale the difference to high-res and convolve with gaussian. merge back to im_h.
	 */
	for(int i=1; i<=maxIter; i++)
	{
		resize_image_bau(im_ret, im_l_s, width_h, height_h, width_l, height_l);

		for(int i = 0; i<width_l; i++)
			for(int j = 0; j<height_l; j++)
				im_diff[i*width_l+j]=(double)im_l_d[i*width_l+j] - im_l_s[i*width_l+j];

		/*
	for(int i=0; i<256; i++)
	{
		for(int j=0; j<1; j++)
		{
			printf("%f ", im_diff[i*256+j]);
		}
	}
	*/

		resize_image_bau(im_diff, im_diff_h, width_l, height_l, width_h, height_h);

		convolve2DSeparable(im_diff_h, im_diff_conv, width_h, height_h,  kernelX, 5, kernelY, 5);

		for(int i = 0; i<width_h; i++)
			for(int j = 0; j<height_h; j++)
				im_ret[i*width_h+j]=im_ret[i*width_h+j] + im_diff_conv[i*width_h+j];
	}

	//free(im_h_d);
	free(im_l_d);
	free(im_l_s);
	free(im_diff);
	free(im_diff_h);
	free(im_diff_conv);
}

/* final p used in backprojection */
/*
 0.0001  0.0021  0.0058  0.0021  0.0001  
 0.0021  0.0431  0.1171  0.0431  0.0021  
 0.0058  0.1171  0.3183  0.1171  0.0058  
 0.0021  0.0431  0.1171  0.0431  0.0021  
 0.0001  0.0021  0.0058  0.0021  0.0001 
*/



int main(int argc, char* argv[])
{
	FILE *im_h_fp = NULL;
	FILE *im_l_fp = NULL;

	im_h_fp = fopen("./im_h.dat", "rb");
	im_l_fp = fopen("./im_l.dat", "rb");

	int size_h=256, size_l=128;
	double* im_final = (double*)malloc(size_h*size_h*sizeof(double));
	uint8* im_h = (uint8*)malloc(size_h*size_h*sizeof(uint8));
	uint8* im_l = (uint8*)malloc(size_l*size_l*sizeof(uint8));
	fread(im_h, sizeof(uint8),  256*256, im_h_fp );
	fread(im_l, 1,  128*128, im_l_fp );

	backprojection(im_h, 256, 256, im_l, 128, 128, 20, im_final);

	FILE* im_out=NULL;
	im_out = fopen("./im_out.dat", "wb");
	fwrite(im_final, 1, 256*256, im_out);

	/*
	for(int i=0; i<256; i++)
	{
		for(int j=0; j<256; j++)
		{
			printf("%f ", im_final[i*256+j]);
		}
		printf("\n");
	}
	*/


	free(im_final);
	free(im_h);
	free(im_l);
	return 0;



}
