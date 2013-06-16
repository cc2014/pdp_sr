// Demo_SR.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "ScSR.h"

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
	
	g_ParamScSR.DlTxDl = new double[g_ParamScSR.Dlw*g_ParamScSR.Dlw];
	if(!g_ParamScSR.DlTxDl){
		printf( "Error: No memory!!\n");
	}
    fread( g_ParamScSR.DlTxDl, g_ParamScSR.Dlw*g_ParamScSR.Dlw, sizeof(double), fid );
    fclose(fid);
    //%% --------------------
}
void Demo_SR()
{
	int nrow, ncol;
	unsigned char *im_l_y;
	im_l_y=NULL;
	
	
	ReadDictionary();

    FILE *fid = fopen( "im_l_y.dat", "rb" );    
    fread( &nrow, 1, sizeof(int), fid );
	fread( &ncol, 1, sizeof(int), fid );    	
	im_l_y = new unsigned char[nrow*ncol];
	if(!im_l_y){
		printf("Error: No memory!\n");
	}
	fread( im_l_y, (nrow*ncol), sizeof(unsigned char), fid );
    fclose(fid);
#if 0
	unsigned char *im_h_y = new unsigned char[nrow*ncol*4];
	fid = fopen( "im_h_y_256x256.dat", "rb" );
	fread( im_h_y, (nrow*ncol*4), sizeof(unsigned char), fid );
    fclose(fid);
	
	double *src_image=new double[nrow*ncol*4];
	double *dst_image=new double[nrow*ncol];
	// test
	int src_nrow = (nrow<<1);
	int src_ncol = (ncol<<1);

	int dst_nrow = (nrow);
	int dst_ncol = (ncol);
	for( int i=0; i<(nrow*ncol*4); i++ )
	{
		src_image[i] = im_h_y[i];
	}
	resize_image_d( src_image, dst_image, src_nrow, src_ncol, dst_nrow, dst_ncol ); 

	delete[]dst_image;
	delete[]src_image;
#endif
	ScSR( im_l_y, nrow, ncol, g_ParamScSR );
	
	if(im_l_y!=NULL){delete[]im_l_y;im_l_y=NULL;}
	
	free_ParamScSR(g_ParamScSR);
}

int _tmain(int argc, _TCHAR* argv[])
{

	Demo_SR();



	return 0;
}

