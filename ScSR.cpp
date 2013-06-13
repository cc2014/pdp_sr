#include "stdafx.h"
#include "ScSR.h"

/*
typedef struct _ParamScSR
{
       double up_scale, lambda, overlap;
       double *Dh, *Dl, *DlTxDl;
	   int Dhw, Dhh, Dlw, Dlh;
       int patch_size, maxIter; 

}ParamScSR;
*/
void copy_gray_image_d( const double *pSrc, int &src_imgw, int &src_imgh, int &start_x, int &start_y, double *pDst, int &dst_imgw, int &dst_imgh );
void resize_image_bau( unsigned char *src_data, unsigned char *dst_data, const int &src_w, const int &src_h, const int &dst_w, const int &dst_h );
bool convolve2DSeparable(double* in, double* out_horl, double *out_vert, int dataSizeX, int dataSizeY,    
                         double* kernelX, int kSizeX, double* kernelY, int kSizeY);

						 
//% first order gradient filters
double hf1[] = {-1,0,1};
double vf1[] = {-1,0,1};
 
//% second order gradient filters
double hf2[] = {1,0,-2,0,1};
double vf2[] = {1,0,-2,0,1};

  
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
int write_pgm_y(char* filename, int x_dim, int y_dim, unsigned char* image)
{

  unsigned char* y = image;
  FILE* filehandle;
  filehandle = fopen(filename, "wb");
  if (filehandle) {
    fprintf(filehandle, "P5\n\n%d %d 255\n", x_dim, y_dim);
    fwrite(y, 1, x_dim * y_dim, filehandle);
    fclose(filehandle);
    return 0;
  } else
    return 1;

}

bool ScSR( unsigned char *im_l_y, int &nrow, int &ncol, ParamScSR &strParamScSR )
{
	// % bicubic interpolation of the low-resolution image
	// mIm = single(imresize(lIm, up_scale, 'bicubic'));
	// hIm = zeros(size(mIm));
	// cntMat = zeros(size(mIm));
	// [h, w] = size(mIm);
	// % extract low-resolution image features
	// lImfea = extr_lIm_fea(mIm);

	int i, nrowx2, ncolx2;
	unsigned char *byte_mIm;
	double *mIm, *lImfea, *mPatch, *mPatchFea;
	
	mPatch = new double[nrow*ncol];
	mPatchFea = new double[nrow*ncol];
	
	nrowx2 = (nrow<<1);
	ncolx2 = (ncol<<1);
	byte_mIm = new unsigned char[nrowx2*ncolx2];
	mIm = new double[nrowx2*ncolx2];
	lImfea = new double[nrowx2*ncolx2*4];
	if(!byte_mIm || !mIm || !lImfea){
		printf("Error: No memory!\n" );
	}
	resize_image_bau( im_l_y, byte_mIm, nrow, ncol, nrowx2, ncolx2 ); 
	
	write_pgm_y( "im_l_y_scale2.pgm", nrowx2, ncolx2, byte_mIm );

	for( i=0; i<(nrowx2*ncolx2); i++ )
	{
		mIm[i] = (double)byte_mIm[i];
	}
	
	// % extract low-resolution image features	
	extr_lIm_fea( mIm, lImfea, nrowx2, ncolx2 );
	
	//% patch indexes for sparse recovery (avoid boundary)
	//gridx = 3:patch_size - overlap : w-patch_size-2;
	//gridx = [gridx, w-patch_size-2];
	//gridy = 3:patch_size - overlap : h-patch_size-2;
	//gridy = [gridy, h-patch_size-2];	
	int *gridx = new int[ncolx2];
	int *gridy = new int[nrowx2];
	memset( gridx, 0, sizeof(int)*ncolx2 );
	memset( gridy, 0, sizeof(int)*nrowx2 );
	int step = strParamScSR.patch_size - strParamScSR.overlap;
	//gridx[0]=3;
	//gridy[0]=3;
	int count_idx=0;
	for( i=3; i<=(ncolx2-strParamScSR.patch_size-2); i+=step ){
		gridx[count_idx] = i;	
		count_idx++;
	}
	int gridx_len = count_idx;
	gridx[count_idx] = ncolx2-strParamScSR.patch_size-2;
	
	count_idx=0;
	for( i=3; i<=(nrowx2-strParamScSR.patch_size-2); i+=step ){
		gridy[count_idx] = i;	
		count_idx++;
	}
	gridy[count_idx] = nrowx2-strParamScSR.patch_size-2;
	int gridy_len = count_idx;

	int ii, jj, kk, cnt = 0;

	for( ii=0; ii<gridx_len; ii++ )
	{
		for( jj=0; jj<gridy_len; jj++ )
		{
			   
			//cnt = cnt+1;
			int xx = gridx[ii];
			int yy = gridy[jj];
			//mPatch = mIm(yy:yy+patch_size-1, xx:xx+patch_size-1);
			//mMean = mean(mPatch(:));
			//mPatch = mPatch(:) - mMean;
			//mNorm = sqrt(sum(mPatch.^2));
			int pxp = strParamScSR.patch_size*strParamScSR.patch_size;

			copy_gray_image_d( mIm, ncolx2, nrowx2, xx, yy, 
				mPatch, strParamScSR.patch_size, strParamScSR.patch_size );
			double mean_val=0, mNorm=0, mfNorm=0;
			for( kk=0; kk<(strParamScSR.patch_size*strParamScSR.patch_size); kk++ ){
				mean_val+=mPatch[kk];
			}
			mean_val/=(strParamScSR.patch_size*strParamScSR.patch_size);
			for( kk=0; kk<(strParamScSR.patch_size*strParamScSR.patch_size); kk++ ){
				mPatch[kk]-=mean_val;
				mNorm+=(mPatch[kk]*mPatch[kk]);
			}
			mNorm = sqrt(mNorm);

			//mPatchFea = lImfea(yy:yy+patch_size-1, xx:xx+patch_size-1, :);   
			//mPatchFea = mPatchFea(:);
			//mfNorm = sqrt(sum(mPatchFea.^2));
			//if mfNorm > 1,
			//	y = mPatchFea./mfNorm;
			//else
			//	y = mPatchFea;
			//end

			for( kk=0; kk<4; kk++ ){
				copy_gray_image_d( lImfea+kk*(ncolx2*nrowx2), ncolx2, nrowx2, xx, yy, 
					mPatchFea+kk*(strParamScSR.patch_size*strParamScSR.patch_size), strParamScSR.patch_size, strParamScSR.patch_size );
			}
			for( kk=0; kk<4*pxp; kk++ ){
				mfNorm += (mPatchFea[kk]*mPatchFea[kk]);	
			}
			mfNorm = sqrt(mfNorm);
			if(mfNorm>1){
				for( kk=0; kk<4*pxp; kk++ ){
					mPatchFea[kk] /= mfNorm;	
				}
				
			}

			//b = -Dl'*y;
	  //    
			//% sparse recovery
			//w = L1QP_FeatureSign_yang(lambda, A, b);
	  //      
			//% generate the high resolution patch and scale the contrast
			//hPatch = Dh*w;
			//hPatch = lin_scale(hPatch, mNorm);
	  //      
			//hPatch = reshape(hPatch, [patch_size, patch_size]);
			//hPatch = hPatch + mMean;
	  //      
			//hIm(yy:yy+patch_size-1, xx:xx+patch_size-1) = hIm(yy:yy+patch_size-1, xx:xx+patch_size-1) + hPatch;
			//cntMat(yy:yy+patch_size-1, xx:xx+patch_size-1) = cntMat(yy:yy+patch_size-1, xx:xx+patch_size-1) + 1;



			cnt++;
		}
	}

	



	delete[]gridx;
	delete[]gridy;
	delete[]byte_mIm;
	delete[]lImfea;
	delete[]mIm;
	delete[]mPatch;
	delete[]mPatchFea;
	
	return true;
}


void copy_gray_image( const unsigned char *pSrc, int &src_imgw, int &src_imgh, int &start_x, int &start_y, unsigned char *pDst, int &dst_imgw, int &dst_imgh )
{
	int i, x, y, src_img_dim, end_x, end_y, count;
	end_x = start_x + dst_imgw;
	end_y = start_y + dst_imgh;

	count = 0;

	/*i = 0;
	for( y=0; y<src_imgh; y++ ) {
	for( x=0; x<src_imgw; x++ , i++) {
	if( x>=start_x && x<end_x && y>=start_y && y<end_y ) {
	pDst[count++] = pSrc[i];
	}
	}
	}//*/
	if( start_x<0 ){ start_x = 0; }
	if( start_y<0 ){ start_y = 0; }
	if( (start_x+dst_imgw)>src_imgw ){ dst_imgw = src_imgw - start_x; }
	if( (start_y+dst_imgh)>src_imgh ){ dst_imgh = src_imgh - start_y; }

	unsigned char *ptraa = (unsigned char*)pSrc + start_y*src_imgw + start_x;
	unsigned char *ptrbb = pDst;
	for( y=0; y<dst_imgh; y++ )
	{
		memcpy( ptrbb, ptraa, dst_imgw ); 
		ptraa+=src_imgw;
		ptrbb+=dst_imgw;
	}


}

void copy_gray_image_d( const double *pSrc, int &src_imgw, int &src_imgh, int &start_x, int &start_y, double *pDst, int &dst_imgw, int &dst_imgh )
{
	int i, x, y, src_img_dim, end_x, end_y, count;
	end_x = start_x + dst_imgw;
	end_y = start_y + dst_imgh;

	count = 0;

	/*i = 0;
	for( y=0; y<src_imgh; y++ ) {
	for( x=0; x<src_imgw; x++ , i++) {
	if( x>=start_x && x<end_x && y>=start_y && y<end_y ) {
	pDst[count++] = pSrc[i];
	}
	}
	}//*/
	if( start_x<0 ){ start_x = 0; }
	if( start_y<0 ){ start_y = 0; }
	if( (start_x+dst_imgw)>src_imgw ){ dst_imgw = src_imgw - start_x; }
	if( (start_y+dst_imgh)>src_imgh ){ dst_imgh = src_imgh - start_y; }

	double *ptraa = (double*)pSrc + start_y*src_imgw + start_x;
	double *ptrbb = pDst;
	for( y=0; y<dst_imgh; y++ )
	{
		memcpy( ptrbb, ptraa, sizeof(double)*dst_imgw ); 
		ptraa+=src_imgw;
		ptrbb+=dst_imgw;
	}


}

static void resize_image_bau( unsigned char *src_data, unsigned char *dst_data, const int &src_w, const int &src_h, const int &dst_w, const int &dst_h ) 
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


///////////////////////////////////////////////////////////////////////////////   
// double precision float version   
///////////////////////////////////////////////////////////////////////////////   
static bool convolve2DSeparable(double* in, double* out_horl, double *out_vert, int dataSizeX, int dataSizeY,    
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