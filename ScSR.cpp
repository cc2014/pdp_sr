//#include "stdafx.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "ScSR.h"
#include "img_utils.h"

using namespace std;

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
double sum_of_product( double *vector_a, double *vector_b, int vector_length );
void set_patch_image( double *pSrc, int *labelMtx, int &src_imgw, int &src_imgh, int &start_x, int &start_y, double *pDst, int &dst_imgw, int &dst_imgh );						 
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

	int i, j, ii, nrowx2, ncolx2;
	unsigned char *byte_mIm;
	double *mIm, *lImfea, *mPatch, *mPatchFea;
	double *b, *DlT, *L1QP_w, *hPatch, *hIm;
	int *cntMat;

	mPatch = new double[nrow*ncol];
	mPatchFea = new double[nrow*ncol];
	b = new double[strParamScSR.Dlw];
	DlT = new double[strParamScSR.Dlw*strParamScSR.Dlh];
	L1QP_w = new double[strParamScSR.Dhh*strParamScSR.Dhw];
	hPatch = new double[strParamScSR.Dhh];
    
	// transpose
	for( i=0; i<strParamScSR.Dlh; i++ )
	{
		ii=i*strParamScSR.Dlw;
		for( j=0; j<strParamScSR.Dlw; j++ ){
			DlT[j*strParamScSR.Dlh+i] = strParamScSR.Dl[ii+j];
		}
	}

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

	hIm = new double[nrowx2*ncolx2];
	cntMat = new int[nrowx2*ncolx2];
	memset( hIm, 0, sizeof(double)*nrowx2*ncolx2 );
	memset( cntMat, 0, sizeof(int)*nrowx2*ncolx2 );
	for( i=0; i<(nrowx2*ncolx2); i++ )
	{
		mIm[i] = (double)byte_mIm[i];
	}


	/* read from raw data !!*/
	FILE* fp=fopen("mIm_256x256.dat", "rb");
	fread(mIm, sizeof(double), 256*256, fp);
	fclose(fp);

	// % extract low-resolution image features	
	extr_lIm_fea( mIm, lImfea, nrowx2, ncolx2 );

	fp = fopen("lImfea_1x256x256x4.dat", "rb");
	fread(lImfea, sizeof(double), 256*256*4, fp);
	fclose(fp);
	
int tmp;	
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
		dprintf("gridy[%d]=%d\n", count_idx, gridx[count_idx]);
		count_idx++;
	}
	dprintf("count_idx=%d\n", count_idx);
	int gridx_len = count_idx+1;
	gridx[count_idx] = ncolx2-strParamScSR.patch_size-2;
	
	count_idx=0;
	for( i=3; i<=(nrowx2-strParamScSR.patch_size-2); i+=step ){
		gridy[count_idx] = i;	
		dprintf("gridy[%d]=%d\n", count_idx, gridy[count_idx]);
		count_idx++;
	}
	dprintf("count_idx=%d\n", count_idx);
	cin>>tmp;
	gridy[count_idx] = nrowx2-strParamScSR.patch_size-2;
	int gridy_len = count_idx;//+1;

	int jj, kk, cnt = 0;

	for( ii=0; ii<gridx_len; ii++ )
	{
		for( jj=0; jj<gridy_len; jj++ )
		{
			   
			dprintf("checkpoint jj=%d\n", jj);
			//cnt = cnt+1;
			int xx = gridx[ii]-1;
			int yy = gridy[jj]-1;
			dprintf("xx=%d, yy=%d\n",xx, yy);
			//mPatch = mIm(yy:yy+patch_size-1, xx:xx+patch_size-1);
			//mMean = mean(mPatch(:));
			//mPatch = mPatch(:) - mMean;
			//mNorm = sqrt(sum(mPatch.^2));
			int pxp = strParamScSR.patch_size*strParamScSR.patch_size;

			copy_gray_image_d( mIm, ncolx2, nrowx2, xx, yy, 
				mPatch, strParamScSR.patch_size, strParamScSR.patch_size );
			/*
			for(i=0; i<strParamScSR.patch_size; i++)
			{
				for(j=0; j<strParamScSR.patch_size; j++)
				{
					mPatch[i*strParamScSR.patch_size+j]=mPatch_tmp[j*strParamScSR.patch_size+i];
				}
//				printf("%.02f \n", mPatch[i]);
			}
//			cin>>tmp;
//			*/

	dprintf("checkpoint\n");

			double mean_val=0, mNorm=0, mfNorm=0;
			for( kk=0; kk<pxp; kk++ ){
				mean_val+=mPatch[kk];
			}


			mean_val/=pxp;
	dprintf("checkpoint mean_val=%f pxp=%d\n", mean_val, pxp);
			for( kk=0; kk<pxp; kk++ ){
				dprintf("mPatch[%d] = %f : ", kk, mPatch[kk]);
				mPatch[kk]-=mean_val;
				dprintf("mPatch[%d] = %f \n", kk, mPatch[kk]);
				mNorm+=pow(mPatch[kk], 2);
			}
			dprintf("checkpoint mNorm=%f\n", mNorm);
			mNorm = sqrt(mNorm);
			dprintf("checkpoint mNorm=%f\n", mNorm);

			for(i=0; i<25; i++)
			{
				printf("%.02f \n", mPatch[i]);
			}


			//mPatchFea = lImfea(yy:yy+patch_size-1, xx:xx+patch_size-1, :);   
			//mPatchFea = mPatchFea(:);
			//mfNorm = sqrt(sum(mPatchFea.^2));
			//if mfNorm > 1,
			//	y = mPatchFea./mfNorm;
			//else
			//	y = mPatchFea;
			//end

	dprintf("checkpoint\n");
			for( kk=0; kk<4; kk++ ){
				copy_gray_image_d( lImfea+kk*(ncolx2*nrowx2), ncolx2, nrowx2, xx, yy, 
					mPatchFea+kk*pxp, strParamScSR.patch_size, strParamScSR.patch_size );
			}

	dprintf("checkpoint\n");

			for( kk=0; kk<4*pxp; kk++ ){
				mfNorm += (mPatchFea[kk]*mPatchFea[kk]);	
			}

	dprintf("checkpoint\n");
			mfNorm = sqrt(mfNorm);
			if(mfNorm>1){
				for( kk=0; kk<4*pxp; kk++ ){
					mPatchFea[kk] /= mfNorm;	
				}
			}
	dprintf("checkpoint\n");

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


				for(int j=0; j<strParamScSR.Dlh; j++ )
				{
				dprintf("checkpoint mPatchFea[%d]=%f\n", j, mPatchFea[j]);
				}
				//cin>>tmp;

			memset(b, 0, sizeof(double)*1024);
			for( i=0; i<strParamScSR.Dlw; i++ ){
				b[i] = -sum_of_product( strParamScSR.Dl+i*strParamScSR.Dlh, mPatchFea, strParamScSR.Dlh ); 
			}

				
	dprintf("checkpoint\n");
			for(i=0; i<1024; i++)
			{
				dprintf("checkpoint b[%d]=%f\n", i, b[i]);
				//cin>>tmp;
			}
			//cin>>tmp;

			dprintf("ii=%d jj=%d\n", ii, jj);
			//% sparse recovery
			//w = L1QP_FeatureSign_yang(lambda, A, b);
			int d_size = 1024;//strParamScSR.Dlw*strParamScSR.Dlw;
			/*
			
			int t1, t2;
			FILE *fid = fopen( "b_1024x1.dat", "rb" );    
			fread( &t1, 1, sizeof(int), fid );
			fread( &t2, 1, sizeof(int), fid );    				
			fread( b, (t1*t2), sizeof(double), fid );
			fclose(fid);

			d_size=t1;
			*/

			L1QP_FeatureSign_yang(strParamScSR.lambda, strParamScSR.DlTxDl, b, L1QP_w, d_size );
			/*
	dprintf("GOGOGOGOGO\n");
	sleep(5);
	*/
			/*
			dprintf("checkpoint I GOT OUT\n");
			for(i=0; i<1024; i++)
			{
				if(L1QP_w[i]!=0)
				{
					dprintf("checkpoint %d =  %f \n ", i, L1QP_w[i]);
				}
			}
			*/

			for( i=0; i<strParamScSR.Dhh; i++ ){
				hPatch[i] = sum_of_product( strParamScSR.Dh+i*strParamScSR.Dhw, L1QP_w, strParamScSR.Dhw );  
			}


	dprintf("checkpoint\n");
			// lin_scale
			double hNorm=0,s;
			for( i=0; i<strParamScSR.Dhh; i++ ){
				hNorm+=(hPatch[i]*hPatch[i]);
			}
			hNorm = sqrt(hNorm);

	dprintf("checkpoint\n");
			if(hNorm!=0){
				s=mNorm*1.2/hNorm;
				for( i=0; i<strParamScSR.Dhh; i++ ){
					hPatch[i] *= s;
				}
			}
			for( i=0; i<strParamScSR.Dhh; i++ ){
				hPatch[i] += mean_val;
			}
			
			for(i=0; i<25; i++)
			{
				dprintf("hPatch[%d]=%f\n", i, hPatch[i]);
			}
			//cin>>tmp;


	dprintf("checkpoint\n");
			set_patch_image( hIm, cntMat, ncolx2, nrowx2, xx, yy, 
				hPatch, strParamScSR.patch_size, strParamScSR.patch_size );
	dprintf("checkpoint\n");

		}
	dprintf("checkpoint\n");
	}
	dprintf("checkpoint\n");

	for( i=0; i<(ncolx2*nrowx2); i++ )
	{
		if(cntMat[i]<1){
			hIm[i] = mIm[i];
			cntMat[i] = 1;
		}
		hIm[i]/=cntMat[i];
	}


	dprintf("checkpoint\n");

	// cast to int8


	delete[]gridx;
	delete[]gridy;
	delete[]byte_mIm;
	delete[]lImfea;
	delete[]mIm;
	delete[]mPatch;
	delete[]mPatchFea;
	delete[]DlT;
	delete[]b;
	delete[]L1QP_w;
	delete[]hPatch;
	delete[]hIm;
	delete[]cntMat;


	return true;
}

double sum_of_product( double *vector_a, double *vector_b, int vector_length )
{
	double sum=0;
	int i;
	double *ptraa=vector_a;
	double *ptrbb=vector_b;

	for( i=0; i<vector_length; i++ )
	{
		sum += (*ptraa * *ptrbb);
		ptraa++;
		ptrbb++;
	}

	return sum;
}

void set_patch_image( double *pSrc, int *labelMtx, int &src_imgw, int &src_imgh, int &start_x, int &start_y, double *pDst, int &dst_imgw, int &dst_imgh )
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
	int *ptrcc = labelMtx + start_y*src_imgw + start_x;
	for( y=0; y<dst_imgh; y++ )
	{
		 
		for( x=0; x<dst_imgw; x++ ){
			ptraa[x]+=ptrbb[x];
			ptrcc[x] += 1; 
		}
		ptraa+=src_imgw;
		ptrbb+=dst_imgw;

		ptrcc += src_imgw;
	}


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
	}*/
	if( start_x<0 ){ start_x = 0; }
	if( start_y<0 ){ start_y = 0; }
	if( (start_x+dst_imgw)>src_imgw ){ dst_imgw = src_imgw - start_x; }
	if( (start_y+dst_imgh)>src_imgh ){ dst_imgh = src_imgh - start_y; }

	/*
	double *ptraa = (double*)pSrc + start_y*src_imgw + start_x;
	double *ptrbb = pDst;
	for( y=0; y<dst_imgh; y++ )
	{
		memcpy( ptrbb, ptraa, sizeof(double)*dst_imgw ); 
		ptraa+=src_imgw;
		ptrbb+=dst_imgw;
	}
	*/

	double *new_dst=new double[dst_imgw*dst_imgh];
	
	double *ptraa = (double*)pSrc + start_y*src_imgw + start_x;
	double *ptrbb = new_dst;
	for( y=0; y<dst_imgh; y++ )
	{
		memcpy( ptrbb, ptraa, sizeof(double)*dst_imgw ); 
		ptraa+=src_imgw;
		ptrbb+=dst_imgw;
	}
	
	
	for( y=0; y<dst_imgh; y++ )
	{
		for( x=0; x<dst_imgw; x++ )
		{	
			pDst[x*dst_imgh+y] = new_dst[y*dst_imgw+x];		
		}
	}
	delete[]new_dst;

}

