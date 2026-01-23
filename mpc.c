//massive parallel cox (marginal) in c
//by msb
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

//for now we hard code no. of subs
#define N 432 //just like Lin's paper -- no. of subjects
#define COVNO 2 //number of covariates
//#define TK 3 //number of time starta
#define MAX_ITER 200 //max iteration for newtown raohson
#define TOLERANCE 1e-15

typedef struct {
  double time[N];
  int status[N];
} DATA;

typedef struct{
  double betavals[COVNO];
  double inv_hessian[COVNO][COVNO];
} DATA_RES;


double Z[N][COVNO]; //data matrix covars are const. at the moment


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
        for(int k=0; k < COVNO; k++) xb += Z[i][k] * beta[k];
        double w = exp(xb);

        sum_ekb += w;
        for(int k=0; k < COVNO; k++){
          s1[k] += w * Z[i][k];
          for(int l=0; l < COVNO; l++){
	     info_mat[k][l] += w * Z[i][k] * Z[i][l];

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
        for(int k=0; k < COVNO; k++) U[k] += (Z[i][k] - z_bar[k]);
      }
    }

    //update I
    for(int k=0; k < COVNO; k++){
      for(int l=0; l < COVNO; l++){
	I[k][l] += event_no * (info_mat[k][l]/sum_ekb - z_bar[k]*z_bar[l]);
	//I[k][l] += (info_mat[k][l]/sum_ekb - z_bar[k]*z_bar[l]);
      }
    }
  }
}



void ccox(DATA *dat, DATA_RES *res){
  int covN = COVNO; //number of covrs
  gsl_vector *beta = gsl_vector_calloc(covN);
  gsl_vector *U = gsl_vector_alloc(covN);
  gsl_vector *delta = gsl_vector_alloc(covN);
  gsl_matrix *I = gsl_matrix_alloc(covN, covN);
  gsl_permutation *prmute = gsl_permutation_alloc(covN); 

  for (int iter = 0; iter < MAX_ITER; iter++){
    double U_arr[COVNO];
      double I_arr[COVNO][COVNO];
      
      U_I_Calc(dat, beta->data, U_arr,I_arr);

      for(int i=0; i < covN; i++){
	gsl_vector_set(U, i, U_arr[i]);
	for(int j=0;j<covN;j++) gsl_matrix_set(I, i,j, I_arr[i][j]);  }

      int signval;

      gsl_matrix *I_copy = gsl_matrix_alloc(covN, covN);

      gsl_matrix_memcpy(I_copy, I); ///copy of the info matrix for beta and delta calc
      gsl_linalg_LU_decomp(I_copy, prmute, &signval);
      gsl_linalg_LU_solve(I_copy, prmute, U, delta);
      gsl_matrix_free(I_copy);

      gsl_vector_add(beta, delta); //here beta = beta+delta
     
      //printf("abs sum value %d %e\n",iter, gsl_blas_dasum(delta));
      if (gsl_blas_dasum(delta) < TOLERANCE) break; //gsl dasum absolute sum

  }

  int s;
  gsl_matrix_view invH = gsl_matrix_view_array(&res->inv_hessian[0][0], covN, covN);

  gsl_matrix_memcpy(I, I);
  gsl_linalg_LU_decomp(I, prmute, &s);

  gsl_linalg_LU_invert(I, prmute, &invH.matrix);

  memcpy(res->betavals, beta->data, sizeof(double)*covN);

  //free mem
  gsl_vector_free(beta);
  gsl_vector_free(U);
  gsl_vector_free(delta);
  gsl_matrix_free(I);
 gsl_permutation_free(prmute);
  
  
}



//main prog
int main() {

  FILE *file = fopen("rossi.csv", "r");
   if (!file) {
        perror("Error opening file");
        return 1;}

    DATA S1;
    DATA_RES result;
    char line[1024];
    
    if (!fgets(line, sizeof(line), file)) {
        printf("Empty file\n");
        fclose(file);
        return 1;
    }


    int i = 0;
    int count;
    while(fgets(line, sizeof(line), file) && i < N){
      count = sscanf(line,"%lf,%d,%lf,%lf",
			 &S1.time[i],
			 &S1.status[i],
			 &Z[i][0], // fin
			 &Z[i][1]); // age
      
      i++;
    }
    
    fclose(file);

    //  printf("line count, %d ",i);
    //printf("line count, %d ",count);
  
  ccox(&S1, &result);

printf("\n%-10s %-10s %-10s %-10s %-10s %-20s\n", "Variable", "Coef", "SE", "p-val", "HR", "95% CI");
  printf("------------------------------------------------------------------------------------\n");

  char *names[] = {"fin", "age"};
  double z_crit = 1.960; //for 95% confidence interval

  for (int k = 0; k < COVNO; k++) {
    double beta = result.betavals[k];
    double se = sqrt(result.inv_hessian[k][k]);
    double hr = exp(beta);
    
    //wald test and p value
    double z_stat = beta / se;
    double p_val = erfc(fabs(z_stat) / sqrt(2.0));

    //confidence int for HR
    double ci_low = exp(beta - z_crit * se);
    double ci_high = exp(beta + z_crit * se);

    printf("%-10s %-10.4f %-10.4f %-10.4f %-10.4f (%0.4f, %0.4f)\n", 
           names[k], beta, se, p_val, hr, ci_low, ci_high);
  }

  printf("------------------------------------------------------------------------------------\n");
}
