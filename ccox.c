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
//#define N 432 //just like Lin's paper -- no. of subjects
//#define COVNO 2 //number of covariates
//#define TK 3 //number of time starta
//#define MAX_ITER 200 //max iteration for newtown raohson
//#define TOLERANCE 1e-15

//int N;
//int COVNO;

typedef struct {
  double *time;
  int *status;
} DATA;

typedef struct{
  double *betavals;
  double *inv_hessian;
} DATA_RES;


//double **Z; //[N][COVNO]; //data matrix covars are const. at the moment


void U_I_Calc(DATA *data, int N,int COVNO,double beta[COVNO], double U[COVNO], double I[COVNO][COVNO], double **Z){

  //  printf("N and COVNO %d %d \n", N,COVNO);

  for(int i=0; i < COVNO; i++) U[i] = 0.0;

  for(int i=0; i < COVNO; i++){
    for(int j=0; j < COVNO; j++){
      I[i][j] = 0.0;
    }
  }

  //noticed stack overflow when running datasets larger than 1mil
  //double TiE1[N]; //time at event 1 (so called occured)
  double *TiE1 = (double *)malloc(N*sizeof(double));
  
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

	// printf("here....\n");
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

  free(TiE1);

}



void ccox(DATA *dat, DATA_RES *res, int N, int COVNO, double **Z, int MAX_ITER, double TOLERANCE){
  int covN = COVNO; //number of covrs
  gsl_vector *beta = gsl_vector_calloc(covN);
  gsl_vector *U = gsl_vector_alloc(covN);
  gsl_vector *delta = gsl_vector_alloc(covN);
  gsl_matrix *I = gsl_matrix_alloc(covN, covN);
  gsl_permutation *prmute = gsl_permutation_alloc(covN); 

  for (int iter = 0; iter < MAX_ITER; iter++){
    double U_arr[COVNO];
      double I_arr[COVNO][COVNO];
      
      U_I_Calc(dat, N,COVNO,beta->data, U_arr,I_arr,Z);

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
  //  gsl_matrix_view invH = gsl_matrix_view_array(&res->inv_hessian[0][0], covN, covN);

    gsl_matrix_view invH = gsl_matrix_view_array(res->inv_hessian, covN, covN);
  //  gsl_vector_view invH = gsl_vector_view_array(&res->inv_hessian[0], covN);

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
int main(int argc, char *argv[]) {
  if (argc < 2) {
        fprintf(stderr, "Usage: %s <config_file>\n", argv[0]);
        return 1;
    }

    char csv_file[256*8];
    int N = 0; 
    int COVNO = 0;
    char cov_line[512*4];
    int MAX_ITER = 0;
    double TOLERANCE = 0.0;
    
    FILE *config = fopen(argv[1], "r");
    if (!config) {
        perror("cannot open config file"); 
        return 1;
    }

    if (fscanf(config, "MAX_ITER=%d\nTOLERANCE=%le\nfile=%s\nn=%d\ncovno=%d\ncovariates=%s", &MAX_ITER, &TOLERANCE,csv_file, &N, &COVNO, cov_line) != 6) {
        fprintf(stderr, "Error: Config file format incorrect.\n");
        fclose(config);
        return 1;
    }
    fclose(config);

    char **target_names = malloc(COVNO * sizeof(char *));
    char *token = strtok(cov_line, ",");
    for (int i = 0; i < COVNO && token != NULL; i++) {
        target_names[i] = strdup(token);
        token = strtok(NULL, ",");
    }

    FILE *file = fopen(csv_file, "r");
       if (!file) {
        perror("cannot open csv file");
        return 1;
     }

    char line[1024 * 2];
    if (!fgets(line, sizeof(line), file)) {
        fprintf(stderr, "Empty CSV file\n");
        fclose(file);
        return 1;
    }

    
    line[strcspn(line, "\r\n")] = 0;

    
    char *header_cols[10000]; 
    int total_csv_cols = 0;
    char *hdr_token = strtok(line, ",");
    while (hdr_token != NULL && total_csv_cols < 100) {
        header_cols[total_csv_cols++] = strdup(hdr_token);
	//	printf("header stuff %s\n", hdr_token);
	hdr_token = strtok(NULL, ",");
	
    }

    int *col_indices = malloc(COVNO * sizeof(int));
    int time_idx = -1, status_idx = -1;

    // Find "time" and "status" (required for S1)
    for (int j = 0; j < total_csv_cols; j++) {
      //  printf("headers ben %d %s\n",j,header_cols[j]);
        if (strcmp(header_cols[j], "week") == 0) time_idx = j;
        if (strcmp(header_cols[j], "arrest") == 0) status_idx = j;
    }

    //  printf("time %d\n", time_idx);

    for (int i = 0; i < COVNO; i++) {
        col_indices[i] = -1;
        for (int j = 0; j < total_csv_cols; j++) {
            if (strcmp(target_names[i], header_cols[j]) == 0) {
                col_indices[i] = j;
                break;
            }
        }
        if (col_indices[i] == -1) {
            fprintf(stderr, "Error: Covariate '%s' not found in CSV header!\n", target_names[i]);
            return 1;
        }
    }


    DATA S1;
    S1.time = malloc(N * sizeof(double));
    S1.status = malloc(N * sizeof(int));

    double **Z = malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        Z[i] = malloc(COVNO * sizeof(double));
    }

    DATA_RES result;
    result.betavals = malloc(COVNO * sizeof(double));
    result.inv_hessian = malloc(COVNO *COVNO* sizeof(double));
    // for (int i = 0; i < COVNO; i++) {
    //    result.inv_hessian[i] = malloc(COVNO * sizeof(double));
    // }

    // printf("here");
    int row_count = 0;
    while (fgets(line, sizeof(line), file) && row_count < N) {
        int current_col = 0;
        char *val_token = strtok(line, ",");
        
        while (val_token != NULL) {
            // Assign to S1.time or S1.status if indices match
            if (current_col == time_idx) S1.time[row_count] = atof(val_token);
            if (current_col == status_idx) S1.status[row_count] = atoi(val_token);

            // Assign to Covariate Matrix Z
            for (int k = 0; k < COVNO; k++) {
                if (current_col == col_indices[k]) {
                    Z[row_count][k] = atof(val_token);
                }
            }
            val_token = strtok(NULL, ",");
            current_col++;
        }
        row_count++;
    }
    fclose(file);

    //call ccox
    ccox(&S1, &result, N, COVNO, Z, MAX_ITER, TOLERANCE);

    printf("\n%-10s %-10s %-10s %-10s %-10s %-20s\n", "Variable", "Coef", "SE", "p-val", "HR", "95% CI");
    printf("------------------------------------------------------------------------------------\n");

    double z_crit = 1.960;
    //  double se = 0.0;
    //     for (int k = 0; k < COVNO *COVNO; k++) {
    // printf("inv hessian %d  %f\n",k,result.inv_hessian[k]); //need diagonals aka k=0, k=3 for 2x2 mat
    // }
    for (int k = 0; k < COVNO; k++) {
        double beta = result.betavals[k];
	//        double se = sqrt(result.inv_hessian[k][k]);

	int diag_index = (k * COVNO) + k;
      
      double variance = result.inv_hessian[diag_index];

      double se = (variance > 0) ? sqrt(variance) : 0.0; //ensure + vals before sqrt
      
        double hr = exp(beta);
        double z_stat = beta / se;
        double p_val = erfc(fabs(z_stat) / sqrt(2.0));
        double ci_low = exp(beta - z_crit * se);
        double ci_high = exp(beta + z_crit * se);

        printf("%-10s %-10.4f %-10.4f %-10.4f %-10.4f (%0.4f, %0.4f)\n", 
               target_names[k], beta, se, p_val, hr, ci_low, ci_high);
    }
    printf("------------------------------------------------------------------------------------\n");

    //free mem
   
    for (int i = 0; i < total_csv_cols; i++) free(header_cols[i]);
    for (int i = 0; i < COVNO; i++) free(target_names[i]);
    free(target_names);
    free(col_indices);
    free(S1.time);
    free(S1.status);
    for (int i = 0; i < N; i++) free(Z[i]);
    free(Z);
    free(result.betavals);
    //for (int i = 0; i < COVNO; i++) free(result.inv_hessian[i]);
       free(result.inv_hessian);
    // for (int i = 0; i < COVNO; i++) free(result.inv_hessian[i]);

   

    return 0;
}
