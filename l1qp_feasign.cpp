//#include "stdafx.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
//#include <armadillo>
#include "ScSR.h"
#include "alloc_util.h"

#define dprintf(fmt, ...) 

using namespace std;
//using namespace arma;

/* 
 * matlab code for solving L1QP 
 * =======================================================================
function [x]=L1QP_FeatureSign_yang(lambda,A,b)

A = double(A);
b = double(b);

EPS = 1e-9;
x=zeros(size(A, 1), 1);           %coeff

grad=A*sparse(x)+b;
[ma mi]=max(abs(grad).*(x==0));

while true,
    
    
  if grad(mi)>lambda+EPS,
    x(mi)=(lambda-grad(mi))/A(mi,mi);
  elseif grad(mi)<-lambda-EPS,
    x(mi)=(-lambda-grad(mi))/A(mi,mi);            
  else
    if all(x==0)
      break;
    end
  end    
  
  while true,
    a=x~=0;   %active set
    Aa=A(a,a);
    ba=b(a);
    xa=x(a);

    %new b based on unchanged sign
    vect = -lambda*sign(xa)-ba;
    x_new= Aa\vect;
    idx = find(x_new);
    o_new=(vect(idx)/2 + ba(idx))'*x_new(idx) + lambda*sum(abs(x_new(idx)));
    
    %cost based on changing sign
    s=find(xa.*x_new<=0);
    if isempty(s)
      x(a)=x_new;
      loss=o_new;
      break;
    end
    x_min=x_new;
    o_min=o_new;
    d=x_new-xa;
    t=d./xa;
    for zd=s',
      x_s=xa-d/t(zd);
      x_s(zd)=0;  %make sure it's zero
%       o_s=L1QP_loss(net,Aa,ba,x_s);
      idx = find(x_s);
      o_s = (Aa(idx, idx)*x_s(idx)/2 + ba(idx))'*x_s(idx)+lambda*sum(abs(x_s(idx)));
      if o_s<o_min,
        x_min=x_s;
        o_min=o_s;
      end
    end
    
    x(a)=x_min;
    loss=o_min;
  end 
    
  grad=A*sparse(x)+b;
  
  [ma mi]=max(abs(grad).*(x==0));
  if ma <= lambda+EPS,
    break;
  end
end
 * =======================================================================
*/

// A: 1024x1024 double array
// b: 1024x1 double array
// x_ret: 1024x1 double array, representing resolved coefficients.
void L1QP_FeatureSign_yang(const double &lambda, double* A, double* b, double* x_ret, int &d_size)
{
	int tmp;
	//int d_size=1024;
	double EPS=1e-9;
	double* grad=new double[d_size];
//	memcpy(grad, b, sizeof(double)*d_size);
	memset(grad, 0, sizeof(double)*d_size);
	memset( x_ret, 0, sizeof(double)*d_size );

	double loss;
	double max=-10000000;
	int max_i;
	int i, j;
	double o_s, o_min=10000000;
	//% grad=A*sparse(x)+b;==> grad=A*x+b
	for( i=0; i<d_size; i++ ){
			for(j=0; j<d_size; j++)
			{
					grad[i] += A[i*d_size+j]*x_ret[j]; 
			}
		grad[i] += b[i];
	}

	/*
	for(i=0;i<d_size; i++)
	{
		printf("gradient[%d]=%f\n",i,  grad[i]);
	}
	*/

	//% [ma mi]=max(abs(grad).*(x==0));
	for( i=0; i<d_size; i++){
		if(x_ret[i]==0)
		{
			if( fabs(grad[i]) > max )
			{
				max=fabs(grad[i]);
				max_i=i;
			}
		}
	}


	while(true)
	{

		if(grad[max_i]>lambda+EPS)
		{
			x_ret[max_i]=(lambda-grad[max_i])/A[max_i*d_size+max_i];
		}
		else if(grad[max_i]<-lambda-EPS)
		{
			x_ret[max_i]=(-lambda-grad[max_i])/A[max_i*d_size+max_i];
		}
		else
		{
			// leave loop if all elements zero.
			int i=0;
			for(i=0; i<d_size; i++)
			{
				if(x_ret[i]>0)
					break;
			}
			if(i==d_size)
				break;
		}



		while(true)
		{
			
			//vector<int> a; 

			/*
			 * x_ret is an array of size $d_size
			 * vector $a stores the index of $x_ret where $w_ret[a.item] != 0
			 * Aa=A(a, a): create a 2-D array that stores value from 2-D array $A, 
			 *				of the indexes popped from vector $a.
			 */ 
		
			int *a=new int[d_size];
			int *idx=new int[d_size];
			int i, j, k, a_size=0, idx_count=0;
			//%a=x~=0;   %active set
			for( i=0; i<d_size; i++){
				if(x_ret[i]!=0)	{
					a[a_size++] = i;
				}
			}
			// attention: a_size>0
			double *Aa = new double[a_size*a_size];
			double *Ba = new double[a_size];
			double *Xa = new double[a_size];
			double *vect=new double[a_size];
			double * x_new = new double[d_size];

			//%Aa=A(a,a);
			//%ba=b(a);
			//%xa=x(a);
			for( i=0; i<a_size; i++){
				for( j=0; j<a_size; j++){
					Aa[i*a_size+j]=A[ a[i]*d_size + a[j] ];
				}
				Ba[i]=b[ a[i] ];
				Xa[i]=x_ret[ a[i] ];
			}
	
			//%new b based on unchanged sign
						
			for( i=0; i<a_size; i++)
			{
				if(Xa[i]>0)
					vect[i] = -lambda-Ba[i];
				else if(Xa[i]==0)
					vect[i] = -Ba[i];
				else
					vect[i] = lambda-Ba[i];
			}


			// x_new= Aa\vect; ==> x_new = inv(Aa)*vect; 

			double det_man; /* determinant mantisa */
			int    det_exp; /* determinant exponent */
			//int		Aa_w=a_size;

			double **invAa=G_alloc_matrix(a_size,a_size);
			int    *indx = G_alloc_ivector(a_size);
			double **y = G_alloc_matrix(a_size,a_size);
			double *col = G_alloc_vector(a_size);

			double small_value=10000000;
			for(i=0;i<a_size;i++){
				if( Aa[i*a_size+i]<small_value )
					small_value = Aa[i*a_size+i]; 
			}
			for( i=0; i<a_size; i++ ){
				for( j=0; j<a_size; j++ ){
					invAa[i][j] = Aa[i*a_size+j];
				}
				// for make non-singular matrix
				if(i==j){
					invAa[i][j] += small_value*0.001;
				}
			}


			clust_invert(
			  invAa,      /* input/output matrix */
			  a_size,        /* dimension */
			  &det_man, /* determinant mantisa */
			  &det_exp, /* determinant exponent */
			  /* scratch space */
			  indx,    /* indx = G_alloc_ivector(n);  */
			  y,      /* y = G_alloc_matrix(n,n); */
			  col      /* col = G_alloc_vector(n); */
			);

			double *invAa_raw = new double[a_size*a_size];


			for( i=0; i<a_size; i++ ){
				for( j=0; j<a_size; j++ ){
					invAa_raw[i*a_size+j] = invAa[i][j];
				}
			}
			// inv(Aa)*vect
			for( i=0; i<a_size; i++){
				x_new[i] = sum_of_product( invAa_raw+i*a_size, vect, a_size );
			//	dprintf("checkpoint x_new[%d]=%f\n", i, x_new[i]);
			}
			//cin>>tmp;
			delete [] invAa_raw;

			//%idx = find(x_new);
			//%o_new=(vect(idx)/2 + ba(idx))'*x_new(idx) + lambda*sum(abs(x_new(idx)));
			idx_count=0;
			double sum_value=0;
			// lambda*sum(abs(x_new(idx)))
			for( i=0; i<a_size; i++){
				if(x_new[i]!=0){
					idx[idx_count++] = i;
					sum_value+=fabs(x_new[i]);
				}
			}

			double o_new=0;
			if(idx_count>0)
			{
				double tt_value_aa;
				
				//%(vect(idx)/2 + ba(idx))'*x_new(idx)
				for(i=0;i<idx_count;i++){
					tt_value_aa = vect[ idx[i] ]/2 + Ba[ idx[i] ];
					tt_value_aa *= x_new[ idx[i] ];

					o_new += tt_value_aa;
				}
				o_new += (lambda*sum_value);
				
			}

			//%cost based on changing sign
			int s_size=0;
			int *s=new int[idx_count];
			for(i=0;i<a_size;i++)
			{
				if((Xa[i]*x_new[i])<=0)
				{
					s[s_size++]=i;	
				}
			}

			if(s_size==0)
			{

				for(i=0; i<a_size; i++ )
				{
					x_ret[ a[i] ] = x_new[i];	
				}
				
				loss = o_new;

				G_free_ivector(indx);
				G_free_vector(col);
				G_free_matrix(invAa);
				G_free_matrix(y);

				delete [] a;
				delete [] idx;

				delete [] Aa;
				delete [] Ba;
				delete [] Xa;
				delete [] vect;
				delete [] x_new;
				
				delete [] s;
				
				break;
			}


			//x_min=x_new;
			//o_min=o_new;
			//d=x_new-xa;
			//t=d./xa;

			double *x_min=new double[a_size];
			double o_min;//=new double[a_size];
			double *temp=new double[a_size];
			double *x_s=new double[a_size];

			for( i=0; i<a_size; i++ ){
				x_min[i] = x_new[i];
				//o_min[i] = o_new[i];
				temp[i] = (x_new[i] - Xa[i])/Xa[i];
				dprintf("x_min[%d]=%f temp[%d]=%f\n", i, x_min[i], i, temp[i]);
			}
			o_min = o_new;

			dprintf("o_min=%f\n", o_min);
				//cin>>tmp;


			for( int zd=0; zd<s_size; zd++ )
			{

				dprintf("s[%d]=%d\n", zd, s[zd]);
				
				for( k=0; k<a_size; k++ ){
					x_s[k] = Xa[k] - (x_new[k] - Xa[k])/temp[ s[zd] ]; 
					dprintf("x_s[%d]=%f\n", k, x_s[k]);
					
				}
				//cin>>tmp;
				x_s[s[zd]] = 0; // %make sure it's zero
				//
				idx_count=0;
				for( k=0; k<a_size; k++ ){
					if(x_s[k]!=0){
						dprintf("idx[%d]=%d\n", idx_count, k);
						idx[idx_count++] = k;
					}
				}
				dprintf("idxcount=%d\n", idx_count);

				//
				//%o_s = (Aa(idx, idx)*x_s(idx)/2 + ba(idx))'*x_s(idx)+lambda*sum(abs(x_s(idx)));
				double o_s=0;
				if(idx_count>0)
				{
					double *Aa_idx=new double[idx_count*idx_count];
					double *ttvec=new double[idx_count];
					double *ttvec2=new double[idx_count];

					for( int m=0; m<idx_count; m++ )
					{
						for( int n=0; n<idx_count; n++ )
						{
							Aa_idx[m*idx_count+n] = Aa[ idx[m]*a_size+idx[n] ];
							dprintf("Aa_idx[%d]=%f\n", m*idx_count+n, Aa_idx[m*idx_count+n]);
						}
						ttvec2[m] = x_s[ idx[m] ];
						dprintf("ttvec2[%d]=%f\n", m, ttvec2[m]);
					}
					//cin>>tmp;

					
					for( int m=0; m<idx_count; m++ ){
						ttvec[m] = sum_of_product( Aa_idx+m*idx_count, ttvec2, idx_count );  
						ttvec[m] /= 2;
						ttvec[m] += Ba[idx[m]];
						dprintf("Ba[%d]=%f\n", idx[m], Ba[idx[m]]);
						dprintf("ttvec[m]=%f\n", ttvec[m]);

						o_s += ( ttvec[m]*x_s[ idx[m] ] + lambda*fabs( x_s[ idx[m] ] ) );
					}
					
				dprintf("o_s=%f\n", o_s);
				//cin>>tmp;


					if(o_s<o_min){
						for( int m=0; m<idx_count; m++ ){
							x_min[m] = x_s[m];
						}
						o_min=o_s;
					}

					delete[]Aa_idx;
					delete[]ttvec;
					delete[]ttvec2;
				}
			}

			//delete[]x_min;
			//delete[]o_min;
			

			
			//delete [] pbuffer;
			
			//delete [] idx_buffer2;
			for( i=0; i<a_size; i++ )
			{
				x_ret[a[i]] = x_min[i];				
			}
			loss=o_min;

			dprintf("loss=%f\n", loss);
			//cin>>tmp;

			G_free_ivector(indx);
			G_free_vector(col);
			G_free_matrix(invAa);
			G_free_matrix(y);

			delete [] a;
			delete [] idx;

			delete [] Aa;
			delete [] Ba;
			delete [] Xa;
			delete [] vect;
			delete [] x_new;
			
			delete [] s;
			
			delete[]x_min;
			delete[]temp;
			delete[]x_s;

		}
		//% grad=A*sparse(x)+b;
		//% [ma mi]=max(abs(grad).*(x==0));
	    //% if ma <= lambda+EPS,
		//%break;
	    //% end

		memset(grad, 0, d_size*sizeof(double));

		for( i=0; i<d_size; i++ )
		{
			for(j=0; j<d_size; j++)
			{
				grad[i] += A[i*d_size+j]*x_ret[j]; 
			}
			grad[i] += b[i];
		}


		max=-10000000;
		for( i=0; i<d_size; i++)
		{
			if(x_ret[i]==0){
				if( fabs(grad[i]) > max )
				{
					max=fabs(grad[i]);
					max_i=i;
				}
			}
		}

		if(max<=(lambda+EPS)){
			break;
		}

	}

	delete[]grad;

}

//int main(int argc, char* argv[])
//{
////	L1QP_FeatureSign_yang();
//	return 0;
//}
