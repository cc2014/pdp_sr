/*
function [lImFea] = extr_lIm_fea( lIm )

[nrow, ncol] = size(lIm);

lImFea = zeros([nrow, ncol, 4]);

% first order gradient filters
hf1 = [-1,0,1];
vf1 = [-1,0,1]';
 
lImFea(:, :, 1) = conv2(lIm, hf1, 'same');
lImFea(:, :, 2) = conv2(lIm, vf1, 'same');

% second order gradient filters
hf2 = [1,0,-2,0,1];
vf2 = [1,0,-2,0,1]';
 
lImFea(:, :, 3) = conv2(lIm,hf2,'same');
lImFea(:, :, 4) = conv2(lIm,vf2,'same');
//*/

//% first order gradient filters
float hf1 = [-1,0,1];
float vf1 = [-1,0,1];
 
//% second order gradient filters
float hf2 = [1,0,-2,0,1];
float vf2 = [1,0,-2,0,1];

bool convolve2DSeparable(double* in, double* out_horl, double *out_vert, int dataSizeX, int dataSizeY,    
                         double* kernelX, int kSizeX, float* kernelY, int kSizeY);   


void extr_lIm_fea( double *lIm, double *lImFea, int &nrow, int &ncol  )
{

	//lImFea(:, :, 1) = conv2(lIm, hf1, 'same');
	//lImFea(:, :, 2) = conv2(lIm, vf1, 'same');
	//lImFea(:, :, 3) = conv2(lIm,hf2,'same');
	//lImFea(:, :, 4) = conv2(lIm,vf2,'same');

	double *ptr1=lImFea;
	double *ptr2=lImFea+nrow*ncol;	
	double *ptr3=lImFea+2*nrow*ncol;	
	double *ptr4=lImFea+3*nrow*ncol;
	// need border image
	convolve2DSeparable( lIm, ptr1, ptr2, ncol, nrow, hf1, 3, vf1, 3 );
	convolve2DSeparable( lIm, ptr3, ptr4, ncol, nrow, hf2, 5, vf2, 5 );


}


///////////////////////////////////////////////////////////////////////////////   
// double precision float version   
///////////////////////////////////////////////////////////////////////////////   
static bool convolve2DSeparable(double* in, double* out, int dataSizeX, int dataSizeY,    
                         double* kernelX, int kSizeX, float* kernelY, int kSizeY)   
{   
    int i, j, k, m, n;   
    double *tmp, *sum;                              // intermediate data buffer   
    double *inPtr, *outPtr;                         // working pointers   
    double *tmpPtr, *tmpPtr2;                       // working pointers   
    int kCenter, kOffset, endIndex;                 // kernel indice   
   
    // check validity of params   
    if(!in || !out || !kernelX || !kernelY) return false;   
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
