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
#define N 36 //just like Lin's paper -- no. of subjects
#define COVNO 2 //number of covariates
#define TK 3 //number of time starta
#define MAX_ITER 200 //max iteration for newtown raohson
#define TOLERANCE 1e-16

typedef struct {
  double time[N];
  int status[N];
} DATA;

typedef struct{
  double betavals[COVNO];
} DATA_RES;


double Z[N][COVNO]; //data matrix covars are const. at the moment

void ccox(DATA *dat, DATA_RES *res){
  int covN = COVNO; //number of covrs
  gsl_vector *beta = gsl_vector_alloc(covN);
  gsl_vector *U = gsl_vector_alloc(covN);
  gsl_vector *delta = gsl_vector_alloc(covN);
  gsl_matrix *I = gsl_matrix_alloc(covN, covN);
    
  
}

void U_I_Calc(DATA *data, double beta[COVNO], double U[COVNO], double I[COVNO][COVNO]){

  for(int i=0; i < COVNO; i++) U[i] = 0.0;

  for(int i=0; i < COVNO; i++){
    for(int j=0; j < COVNO; j++){
      I[i][j] = 0.0;
    }
  }

  double TiE1[N]; //time at event 1 (so called occured)
  int E1 = 0; //occured status (event = 1)

  for (int i=0; i < N; i++){
    if(data->status[i] == 1){
      int E1_here = 0;
      for (int j=0; j < E1; j++){
        if(TiE1[j] == data->time[i]) E1_here = 1;
      }
      if(!E1_here) TiE1[E1++] = data->time[i]; 
    }
  }

  //main loop for U and I calcs
  for (int m=0; m < E1; m++) {
    double t = TiE1[m];
    double sum_ekb = 0.0;
    double s1[COVNO];
    double info_mat[COVNO][COVNO];

    for (int i=0; i < COVNO; i++) s1[i] = 0.0;
    for (int i=0; i < COVNO; i++){
      for (int j= 0; j < COVNO; j++){
        info_mat[i][j] = 0.0;
      }
    }

    int event_no = 0;

    //risk set sums
    for (int i=0; i < N; i++){
      if(data->time[i] >= t){
        double xb = 0.0;
        for(int k=0; k < COVNO; k++) xb += data->Z[i][k] * beta[k];
        double w = exp(xb);

        sum_ekb += w;
        for(int k=0; k < COVNO; k++){
          s1[k] += w * data->Z[i][k];
          for(int l=0; l < COVNO; l++){
            info_mat[k][l] += w * data->Z[i][k] * data->Z[i][l];
          }
        }
      }
      if(data->time[i] == t && data->status[i] == 1) event_no++;
    }

    double z_bar[COVNO];
    for(int k=0; k < COVNO; k++) z_bar[k] = s1[k] / sum_ekb;

    //update U
    for(int i=0; i < N; i++){
      if(data->time[i] == t && data->status[i] == 1){
        for(int k=0; k < COVNO; k++) U[k] += (data->Z[i][k] - z_bar[k]);
      }
    }

    //update I
    for(int k=0; k < COVNO; k++){
      for(int l=0; l < COVNO; l++){
        I[k][l] += event_no * (info_mat[k][l]/sum_ekb - z_bar[k]*z_bar[l]);
      }
    }
  }
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
