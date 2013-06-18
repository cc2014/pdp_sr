#include "img_utils.h"

///////////////////////////////////////////////////////////////////////////////   
// double precision float version   
///////////////////////////////////////////////////////////////////////////////   
bool convolve2DSeparable(double* in, double* out_horl, double *out_vert, int dataSizeX, int dataSizeY,    
                         double* kernelX, int kSizeX, double* kernelY, int kSizeY)   
{   
    int i, j, k, m, n;   
    double *tmp, *sum;                              // intermediate data buffer   
    double *inPtr, *outPtr;                         // working pointers   
    double *tmpPtr, *tmpPtr2;                       // working pointers   
    int kCenter, kOffset, endIndex;                 // kernel indice   
   
    // check validity of params   
    if(!in || !out_horl || !out_vert || !kernelX || !kernelY) return false;   
    if(dataSizeX <= 0 || kSizeX <= 0) return false;   
   
    // allocate temp storage to keep intermediate result   
    tmp = new double[dataSizeX * dataSizeY];   
    if(!tmp) return false;  // memory allocation error   
   
    // store accumulated sum   
    sum = new double[dataSizeX];   
    if(!sum) return false;  // memory allocation error   
   
    // covolve horizontal direction ///////////////////////   
   
    // find center position of kernel (half of kernel size)   
    kCenter = kSizeX >> 1;                          // center index of kernel array   
    endIndex = dataSizeX - kCenter;                 // index for full kernel convolution   
   
    // init working pointers   
    inPtr = in;   
    tmpPtr = out_horl;//tmp;                                   // store intermediate results from 1D horizontal convolution   
   
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
    tmpPtr = tmpPtr2 = in;//tmp;   
    outPtr = out_vert;//out;   
   
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
            *outPtr = sum[n];                       // store final result to output array   
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
            *outPtr = sum[n];                       // store final result to output array   
            sum[n] = 0;                             // reset to zero for next summing   
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
            *outPtr = sum[n];                       // store final result to output array   
            sum[n] = 0;                             // reset to zero for next summing   
            ++outPtr;                               // increase ending index of partial kernel   
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

void resize_image_bau( unsigned char *src_data, unsigned char *dst_data, const int &src_w, const int &src_h, const int &dst_w, const int &dst_h ) 
{ 
 
 
    double scalex, scaley; 
    double sr, sc, ratior, ratioc, value1, value2;  
    int dr, dc, isr, isc; 
    unsigned char *imgp; 
    int dd, t, stepr, stepc, ii, b1, b2; 
 
    //unsigned char *data = (unsigned char*)malloc( dst_w*dst_h ); 
    //scalex = (double)dst_w/(double)src_w; 
    //scaley = (double)dst_h/(double)src_h; 
    scalex = (double)src_w/(double)dst_w; 
    scaley = (double)src_h/(double)dst_h; 
     
    b1 = (src_w-1)/scalex; 
    b2 = (src_h-1)/scaley; 
    for (dr=0; dr<b2; dr++) 
    { 
        dd = dr*dst_w; 
        sr = dr*scaley; 
        isr = (int)sr; 
        ratior = sr-isr; 
        ii = isr*src_w; 
        for (dc=0; dc<b1; dc++) 
        {         
            sc = dc*scalex; 
            isc = (int)sc; 
            ratioc = sc-isc; 
            imgp = src_data + ii+isc; 
            value1 = *imgp*(1-ratioc) + *(imgp+1)*ratioc; 
            imgp += src_w; 
            value2 = *imgp*(1-ratioc) + *(imgp+1)*ratioc; 
            dst_data[ dd + dc ] = (unsigned char)(value1*(1-ratior)+value2*ratior); 
        } 
 
        for (dc=b1; dc<dst_w; dc++) 
        {                 
            sc = dc*scalex; 
            isc = (int)sc; 
            ratioc = sc-isc; 
            imgp = src_data + isr*src_w+isc; 
            value1 = *imgp; 
            imgp += src_w; 
            value2 = *imgp; 
            dst_data[ dd + dc ] = (unsigned char)(value1*(1-ratior)+value2*ratior); 
        } 
    } 
 
    for (dr=b2; dr<dst_h; dr++) 
    { 
        dd = dr*dst_w; 
        sr = dr*scaley; 
        isr = (int)sr; 
        ratior = sr-isr; 
        for (dc=0; dc<b1; dc++) 
        {         
            sc = dc*scalex; 
            isc = (int)sc; 
            ratioc = sc-isc; 
            imgp = src_data + isr*src_w+isc; 
            value1 = *imgp*(1-ratioc) + *(imgp+1)*ratioc; 
            value2=value1; 
            dst_data[ dd + dc ] = (unsigned char)(value1*(1-ratior)+value2*ratior); 
        } 
 
        for (dc=b1; dc<dst_w; dc++) 
        {         
            sc = dc*scalex; 
            isc = (int)sc; 
            ratioc = sc-isc; 
            imgp = src_data + isr*src_w+isc; 
            value1 = *imgp; 
            value2 = value1; 
            dst_data[ dd + dc ] = (unsigned char)(value1*(1-ratior)+value2*ratior); 
        } 
 
    } 
}

void resize_image_d( double *src_data, double *dst_data, const int &src_w, const int &src_h, const int &dst_w, const int &dst_h ) 
{ 
 
 
    double scalex, scaley; 
    double sr, sc, ratior, ratioc, value1, value2;  
    int dr, dc, isr, isc; 
    double *imgp; 
    int dd, t, stepr, stepc, ii;
	int b1, b2; 
 
    //unsigned char *data = (unsigned char*)malloc( dst_w*dst_h ); 
    //scalex = (double)dst_w/(double)src_w; 
    //scaley = (double)dst_h/(double)src_h; 
    scalex = (double)src_w/(double)dst_w; 
    scaley = (double)src_h/(double)dst_h; 
     
    b1 = (src_w-1)/scalex; 
    b2 = (src_h-1)/scaley; 
    for (dr=0; dr<b2; dr++) 
    { 
        dd = dr*dst_w; 
        sr = dr*scaley; 
        isr = (int)sr; 
        ratior = sr-isr; 
        ii = isr*src_w; 
        for (dc=0; dc<b1; dc++) 
        {         
            sc = dc*scalex; 
            isc = (int)sc; 
            ratioc = sc-isc; 
            imgp = src_data + ii+isc; 
            value1 = *imgp*(1.0-ratioc) + *(imgp+1)*ratioc; 
            imgp += src_w; 
            value2 = *imgp*(1.0-ratioc) + *(imgp+1)*ratioc; 
            dst_data[ dd + dc ] = (value1*(1.0-ratior)+value2*ratior); 
        } 
 
        for (dc=b1; dc<dst_w; dc++) 
        {                 
            sc = dc*scalex; 
            isc = (int)sc; 
            ratioc = sc-isc; 
            imgp = src_data + isr*src_w+isc; 
            value1 = *imgp; 
            imgp += src_w; 
            value2 = *imgp; 
            dst_data[ dd + dc ] = (value1*(1-ratior)+value2*ratior); 
        } 
    } 
 
    for (dr=b2; dr<dst_h; dr++) 
    { 
        dd = dr*dst_w; 
        sr = dr*scaley; 
        isr = (int)sr; 
        ratior = sr-isr; 
        for (dc=0; dc<b1; dc++) 
        {         
            sc = dc*scalex; 
            isc = (int)sc; 
            ratioc = sc-isc; 
            imgp = src_data + isr*src_w+isc; 
            value1 = *imgp*(1-ratioc) + *(imgp+1)*ratioc; 
            value2=value1; 
            dst_data[ dd + dc ] = (value1*(1-ratior)+value2*ratior); 
        } 
 
        for (dc=b1; dc<dst_w; dc++) 
        {         
            sc = dc*scalex; 
            isc = (int)sc; 
            ratioc = sc-isc; 
            imgp = src_data + isr*src_w+isc; 
            value1 = *imgp; 
            value2 = value1; 
            dst_data[ dd + dc ] = (value1*(1-ratior)+value2*ratior); 
        } 
 
    } 
}
