//massive parallel cox (marginal) in c
//by msb
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

//for now we hard code no. of subs
#define N 36 //just like Lin's paper
#define COVNO 2 //number of covariates
#define TK 3 //number of time starta
#define MAX_ITER 200 //max iteration for newtown raohson
#define TOLERANCE 1e-16

typedef struct {
  double time[N];
  int status[N];
} DATA;


double Z[N][CONVO] //data matrix covars are const. at the moment

void ccox(DATA *dat, DATA_RES *res){
  int covN = COVNO; //number of covrs
  gsl_vector *beta = gsl_vector_alloc(covN);
  gsl_vector *U = gsl_vector_alloc(covN);
  gsl_vector *delta = gsl_vector_alloc(covN);
  gsl_matrix *I = gsl_matrix_alloc(covN, covN);
    
  
}

int main(){

double T1[] = {};
double T2[] = {};
double T3[] = {};

int D1[] = {};
int D2[] = {};
int D3[] = {};

int TREAT[] ={};


 DATA S1, S2, S3; 


  return 0;

}
