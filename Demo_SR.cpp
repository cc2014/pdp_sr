// Demo_SR.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
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
void Demo_SR()
{
	ParamScSR m_ParamScSR;
	int nrow, ncol;
	unsigned char *im_l_y;
	
	im_l_y=NULL;
	set_ParamScSR(m_ParamScSR);
	
    FILE *fid = fopen( "sr.dat", "rb" );        
    fread( &m_ParamScSR.patch_size, 1, sizeof(int), fid );    
	fread( &m_ParamScSR.maxIter, 1, sizeof(int), fid  );    
	// double    
    fread( &m_ParamScSR.up_scale, 1, sizeof(double), fid  );    
    fread( &m_ParamScSR.lambda, 1, sizeof(double), fid  );
    fread( &m_ParamScSR.overlap, 1, sizeof(double), fid );
	// double vector	
	
    fread( &m_ParamScSR.Dhw, 1, sizeof(int), fid );    
	fread( &m_ParamScSR.Dhh, 1, sizeof(int), fid );    	
	m_ParamScSR.Dh = new double[m_ParamScSR.Dhw*m_ParamScSR.Dhh];
	if(!m_ParamScSR.Dh){
		printf( "Error: No memory!!\n");
	}

	fread( m_ParamScSR.Dh, m_ParamScSR.Dhw*m_ParamScSR.Dhh, sizeof(double), fid  );
	
    fread( &m_ParamScSR.Dlw, 1, sizeof(int), fid );    
	fread( &m_ParamScSR.Dlh, 1, sizeof(int), fid );    	
	m_ParamScSR.Dl = new double[m_ParamScSR.Dlw*m_ParamScSR.Dlh];
	if(!m_ParamScSR.Dl){
		printf( "Error: No memory!!\n");
	}

	fread( m_ParamScSR.Dl, m_ParamScSR.Dlw*m_ParamScSR.Dlh, sizeof(double), fid );
	
	m_ParamScSR.DlTxDl = new double[m_ParamScSR.Dlh*m_ParamScSR.Dlh];
	if(!m_ParamScSR.DlTxDl){
		printf( "Error: No memory!!\n");
	}
    fread( m_ParamScSR.DlTxDl, m_ParamScSR.Dlh*m_ParamScSR.Dlh, sizeof(double), fid );
    fclose(fid);
    //%% --------------------
	
    fid = fopen( "im_l_y.dat", "rb" );    
    fread( &nrow, 1, sizeof(int), fid );
	fread( &ncol, 1, sizeof(int), fid );    	
	im_l_y = new unsigned char[nrow*ncol];
	if(!im_l_y){
		printf("Error: No memory!\n");
	}
	fread( im_l_y, (nrow*ncol), sizeof(unsigned char), fid );
    fclose(fid);

	ScSR( im_l_y, nrow, ncol, m_ParamScSR );
	
	if(im_l_y!=NULL){delete[]im_l_y;im_l_y=NULL;}
	free_ParamScSR(m_ParamScSR);
}
#if 0
int _tmain(int argc, _TCHAR* argv[])
{

	Demo_SR();



	return 0;
}
#endif
