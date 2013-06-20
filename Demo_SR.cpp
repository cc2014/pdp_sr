// Demo_SR.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include "ScSR.h"
#include "alloc_util.h"
#include "img_utils.h"
#include <opencv/cv.h>
#include <opencv/highgui.h>

using namespace std; 
using namespace cv; 

void set_ParamScSR( ParamScSR &strParamScSR )
{
	strParamScSR.Dh = NULL;
	strParamScSR.Dl = NULL;
	strParamScSR.DlTxDl = NULL;

}

void free_ParamScSR( ParamScSR &strParamScSR )
{
	if(strParamScSR.Dh!=NULL){delete[]strParamScSR.Dh;strParamScSR.Dh = NULL;}
	if(strParamScSR.Dl!=NULL){delete[]strParamScSR.Dl;strParamScSR.Dl = NULL;}
	if(strParamScSR.DlTxDl!=NULL){delete[]strParamScSR.DlTxDl;strParamScSR.DlTxDl = NULL;}
}

ParamScSR g_ParamScSR;

void ReadDictionary()
{
	set_ParamScSR(g_ParamScSR);
	
    FILE *fid = fopen( "sr.dat", "rb" );        
    fread( &g_ParamScSR.patch_size, 1, sizeof(int), fid );    
	fread( &g_ParamScSR.maxIter, 1, sizeof(int), fid  );    
	// double    
    fread( &g_ParamScSR.up_scale, 1, sizeof(double), fid  );    
    fread( &g_ParamScSR.lambda, 1, sizeof(double), fid  );
    fread( &g_ParamScSR.overlap, 1, sizeof(double), fid );
	// double vector	
	
    fread( &g_ParamScSR.Dhh, 1, sizeof(int), fid );    
	fread( &g_ParamScSR.Dhw, 1, sizeof(int), fid );    	
	g_ParamScSR.Dh = new double[g_ParamScSR.Dhw*g_ParamScSR.Dhh];
	if(!g_ParamScSR.Dh){
		printf( "Error: No memory!!\n");
	}

	fread( g_ParamScSR.Dh, g_ParamScSR.Dhw*g_ParamScSR.Dhh, sizeof(double), fid  );
	
    fread( &g_ParamScSR.Dlh, 1, sizeof(int), fid );    
	fread( &g_ParamScSR.Dlw, 1, sizeof(int), fid );    	
	g_ParamScSR.Dl = new double[g_ParamScSR.Dlw*g_ParamScSR.Dlh];
	if(!g_ParamScSR.Dl){
		printf( "Error: No memory!!\n");
	}

	fread( g_ParamScSR.Dl, g_ParamScSR.Dlw*g_ParamScSR.Dlh, sizeof(double), fid );
	
    fread( &g_ParamScSR.DlTxDl_h, 1, sizeof(int), fid );    
	fread( &g_ParamScSR.DlTxDl_w, 1, sizeof(int), fid );    	

	g_ParamScSR.DlTxDl = new double[g_ParamScSR.DlTxDl_h*g_ParamScSR.DlTxDl_w];
	if(!g_ParamScSR.DlTxDl){
		printf( "Error: No memory!!\n");
	}
    fread( g_ParamScSR.DlTxDl, g_ParamScSR.DlTxDl_h*g_ParamScSR.DlTxDl_w, sizeof(double), fid );
    fclose(fid);
    //%% --------------------
}

void Demo_SR(uchar* im_l_y, int width, int height, uint8* img_result)
{
	int nrow, ncol;
	nrow = width;
	ncol = height;
	
	ReadDictionary();

	/*
    FILE *fid = fopen( "im_l_y.dat", "rb" );    
    fread( &nrow, 1, sizeof(int), fid );
	fread( &ncol, 1, sizeof(int), fid );    	
	im_l_y = new unsigned char[nrow*ncol];
	if(!im_l_y){
		printf("Error: No memory!\n");
	}
	fread( im_l_y, (nrow*ncol), sizeof(unsigned char), fid );
    fclose(fid);
	*/

	int hnrow= (nrow<<1) ;
	int hncol = (ncol<<1);
	double* im_h_y = new double[hnrow*hncol];

	ScSR( im_l_y, nrow, ncol, g_ParamScSR, im_h_y );

	backprojection(im_h_y, hnrow, hncol, im_l_y, nrow, ncol, 20);

	for(int i=0;i<(hnrow*hncol);i++ )
		img_result[i] = (uint8)im_h_y[i];
	
	free_ParamScSR(g_ParamScSR);
//	delete [] im_l_y;
	delete [] im_h_y;
}

//int _tmain(int argc, _TCHAR* argv[])
int main(int argc, char* argv[])
{
	if(argc<2)
	{
		fprintf(stderr, "Usage: %s <path-to-file>\n\n", argv[0]);
		exit(-1);
	}

	Mat img;
	img=imread(argv[1], 1);

	if(!img.data)
	{
		fprintf(stderr, "cannot load image: %s\n\n", argv[1]);
		exit(-1);
	}

	Mat img_yuv;
	vector<Mat> img_y;
	// convert to YUV
	cvtColor(img, img_yuv, CV_RGB2YCrCb);

	split(img_yuv, img_y);

	int hnrow=img.size().width*2;
	int hncol=img.size().height*2;
	uint8* img_result = new uint8[hnrow*hncol];

	Demo_SR(img_y[0].data, img.size().width, img.size().height, img_result);

	Mat Y, Cr, Cb;

	resize(img_y[0], Y, img.size()*2, 0, 0, INTER_CUBIC);
	resize(img_y[1], Cr, img.size()*2, 0, 0, INTER_CUBIC);
	resize(img_y[2], Cb, img.size()*2, 0, 0, INTER_CUBIC);

	Y.data = img_result;
	img_y[0] = Y;
	img_y[1] = Cr;
	img_y[2] = Cb;

	namedWindow("Y", CV_WINDOW_AUTOSIZE);
	imshow("Y", img_y[0]);

	namedWindow("Cb", CV_WINDOW_AUTOSIZE);
	imshow("Cb", img_y[1]);
	
	namedWindow("Cr", CV_WINDOW_AUTOSIZE);
	imshow("Cr", img_y[2]);
	waitKey(0);

	Mat result_y, result;
	merge(img_y, result_y);
	cvtColor(result_y, result, CV_YCrCb2RGB);


	Mat img_cubic;
	resize(img, img_cubic, img.size()*2, 0, 0, INTER_CUBIC);

	imwrite("upscaled_img.jpg", result);

	namedWindow("original", CV_WINDOW_AUTOSIZE);
	imshow("original", img);

	namedWindow("bicubic interpolation", CV_WINDOW_AUTOSIZE);
	imshow("bicubic interpolation", img_cubic);

	namedWindow("upscaled", CV_WINDOW_AUTOSIZE);
	imshow("upscaled", result);

	waitKey(0);

	delete [] img_result;

	return 0;
}



#if 0 // sample for finding invert matrix

	double det_man; /* determinant mantisa */
	int    det_exp; /* determinant exponent */
	int Aa_w, Aa_h;
	double *Aa;


    fid = fopen( "Aa_7x7.dat", "rb" );    
    fread( &Aa_w, 1, sizeof(int), fid );
	fread( &Aa_h, 1, sizeof(int), fid );    	
	Aa = new double[nrow*ncol];
	if(!Aa){
		printf("Error: No memory!\n");
	}
	fread( Aa, (Aa_w*Aa_h), sizeof(double), fid );
    fclose(fid);


	double **a=G_alloc_matrix(Aa_w,Aa_w);
	int    *indx = G_alloc_ivector(Aa_w);
	double **y = G_alloc_matrix(Aa_w,Aa_w);
	double *col = G_alloc_vector(Aa_w);

	for( int i=0; i<Aa_w; i++ ){
		for( int j=0; j<Aa_w; j++ ){
			a[i][j] = Aa[i*Aa_w+j];
		}
	}


	clust_invert(
	  a,      /* input/output matrix */
	  Aa_w,        /* dimension */
	  &det_man, /* determinant mantisa */
	  &det_exp, /* determinant exponent */
	  /* scratch space */
	  indx,    /* indx = G_alloc_ivector(n);  */
	  y,      /* y = G_alloc_matrix(n,n); */
	  col      /* col = G_alloc_vector(n); */
	);

	FILE *fpMtx = fopen( "inv_matrix.txt", "w" );
	for( int i=0; i<Aa_w; i++ ){
		for( int j=0; j<Aa_w; j++ ){
			fprintf( fpMtx, "%lf ", a[i][j] );
		}
		fprintf( fpMtx, "\n" );
	}
	fclose(fpMtx);

	G_free_ivector(indx);
	G_free_vector(col);
	G_free_matrix(a);
	G_free_matrix(y);
	delete[]Aa;

#endif
