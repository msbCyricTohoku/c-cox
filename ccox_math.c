#include "ccox.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

int cmp_double(const void *a, const void *b) {
  
  double da = *(const double *)a;

  double db = *(const double *)b;

  return (da > db) - (da < db);
}

void U_I_Calc(DATA *data, int N, int COVNO, double beta[COVNO], double U[COVNO],
              double I[COVNO][COVNO], double **Z, double *TiE1, int E1, int event_code) {

  for (int i = 0; i < COVNO; i++) {
    U[i] = 0.0;
  }
    

  for (int i = 0; i < COVNO; i++) {

   for (int j = 0; j < COVNO; j++) {

      I[i][j] = 0.0;

   }
 }

#pragma omp parallel
  {

  double core_U[COVNO];

  double core_I[COVNO][COVNO];

  for(int i=0; i < COVNO; i++){

    core_U[i] = 0.0;

    for(int j=0; j < COVNO; j++){

      core_I[i][j] = 0.0;

    }

  }

#pragma omp for
  // main loop for U and I calcs
  for (int m = 0; m < E1; m++) {

    double t = TiE1[m];

    double sum_ekb = 0.0;

    double s1[COVNO];

    double info_mat[COVNO][COVNO];

    for (int i = 0; i < COVNO; i++)

      s1[i] = 0.0;

    for (int i = 0; i < COVNO; i++) {

      for (int j = 0; j < COVNO; j++) {

	info_mat[i][j] = 0.0;

      }
  }

    int event_no = 0;

    // risk set sums
    for (int i = 0; i < N; i++) {

      if (data->start[i] < t && data->stop[i] >= t) {

	double xb = 0.0;

	for (int k = 0; k < COVNO; k++)

	  xb += Z[i][k] * beta[k];

	double w = exp(xb);

        sum_ekb += w;

	for (int k = 0; k < COVNO; k++) {

	  s1[k] += w * Z[i][k];

	  for (int l = 0; l < COVNO; l++) {

	    info_mat[k][l] += w * Z[i][k] * Z[i][l];

	  }

	}

      }

      if (data->stop[i] == t && data->status[i] == event_code)

	event_no++;
  }

    if (sum_ekb > 0.0) {

      double z_bar[COVNO];

      for (int k = 0; k < COVNO; k++)

	z_bar[k] = s1[k] / sum_ekb;

      // update U
      for (int i = 0; i < N; i++) {

	if (data->stop[i] == t && data->status[i] == event_code) {

	  for (int k = 0; k < COVNO; k++)

	    core_U[k] += (Z[i][k] - z_bar[k]);

	}

      }

      // update I
      for (int k = 0; k < COVNO; k++) {

	for (int l = 0; l < COVNO; l++) {

	  core_I[k][l] += event_no * (info_mat[k][l] / sum_ekb - z_bar[k] * z_bar[l]);

     }
   }
   }
  }

#pragma omp critical
  {
    for (int k = 0; k < COVNO; k++) {
       U[k] += core_U[k];
       for (int l = 0; l < COVNO; l++) {
         I[k][l] += core_I[k][l];
       }
    }
  }
  }
}

void compute_robust_variance(DATA *dat, DATA_RES *res, int N, int COVNO, double **Z, 
                             double *TiE1, int E1, int event_code) {

  double *dLambda = (double *)calloc(E1, sizeof(double));

  double *Z_bar_mat = (double *)calloc(E1 * COVNO, sizeof(double));
  
#pragma omp parallel
  {
  double current_beta[COVNO];

  for (int k = 0; k < COVNO; k++) {
    
    current_beta[k] = res->betavals[k];

  }

  double S1[COVNO];
  
#pragma omp for
  for (int m = 0; m < E1; m++) {

    double t = TiE1[m];

    double S0 = 0.0;

    for (int i = 0; i < COVNO; i++)
      {

	S1[i] = 0.0;
    }

    int event_no = 0;

    for (int i = 0; i < N; i++) {

      if (dat->start[i] < t && dat->stop[i] >= t) {

	double xb = 0.0;

	for (int k = 0; k < COVNO; k++)
	  {
	    xb += Z[i][k] * current_beta[k];
	  }

	double w = exp(xb);

	S0 += w;

	for (int k = 0; k < COVNO; k++)
	  {
	    S1[k] += w * Z[i][k];
	  }
      }
      
      if (dat->stop[i] == t && dat->status[i] == event_code) {
	event_no++;
      }
      
    }

    if (S0 > 0.0) {

      dLambda[m] = event_no / S0;

      for (int k = 0; k < COVNO; k++)
	{
	  Z_bar_mat[m * COVNO + k] = S1[k] / S0;
	}

    }
  }

  }

  int max_cluster = 0;

  for (int i = 0; i < N; i++) {
    
    if (dat->cluster[i] > max_cluster) max_cluster = dat->cluster[i];

  }

  
  int num_clusters = max_cluster + 1;

  double *cluster_L = (double *)calloc(num_clusters * COVNO, sizeof(double));

#pragma omp parallel
  {

    double current_beta[COVNO];

    for (int k = 0; k < COVNO; k++) {
      
      current_beta[k] = res->betavals[k];
      
    }

  double L_i[COVNO];

#pragma omp for
  for (int i = 0; i < N; i++) {

    for (int k = 0; k < COVNO; k++) {

      L_i[k] = 0.0;
    }

    double xb = 0.0;

    for (int k = 0; k < COVNO; k++) {
      xb += Z[i][k] * current_beta[k];
    }

    double w = exp(xb);

    int start_idx = 0, low = 0, high = E1 - 1;

    while (low <= high) {

      int mid = low + (high - low) / 2;

      if (TiE1[mid] > dat->start[i]) {
	
	start_idx = mid; high = mid - 1;

      } 

      else {

	low = mid + 1;
     }

    }

    if (low == E1) start_idx = E1;

    for (int m = start_idx; m < E1; m++) {

      if (TiE1[m] > dat->stop[i]) break; 

      for (int k = 0; k < COVNO; k++) {

	L_i[k] -= w * dLambda[m] * (Z[i][k] - Z_bar_mat[m * COVNO + k]);

      }

    }

    if (dat->status[i] == event_code) {

      low = 0; high = E1 - 1;

      int m_idx = -1;

      while (low <= high) {

	int mid = low + (high - low) / 2;

	if (TiE1[mid] == dat->stop[i]) { m_idx = mid; break; }

	else if (TiE1[mid] < dat->stop[i]) low = mid + 1;

	else high = mid - 1;

      }

      if (m_idx != -1) {

	for (int k = 0; k < COVNO; k++)
	  {

	  L_i[k] += (Z[i][k] - Z_bar_mat[m_idx * COVNO + k]);

	  }

      }

    }

    int c = dat->cluster[i];

    for (int k = 0; k < COVNO; k++) {

#pragma omp atomic

      cluster_L[c * COVNO + k] += L_i[k];

   }
  }

  }

  double *B = (double *)calloc(COVNO * COVNO, sizeof(double));

  for (int c = 0; c < num_clusters; c++) {

    for (int k = 0; k < COVNO; k++) {

      for (int l = 0; l < COVNO; l++) {

	B[k * COVNO + l] += cluster_L[c * COVNO + k] * cluster_L[c * COVNO + l];

      }

    }

  }

  gsl_matrix_view B_mat = gsl_matrix_view_array(B, COVNO, COVNO);

  gsl_matrix_view invH = gsl_matrix_view_array(res->inv_hessian, COVNO, COVNO);

  gsl_matrix *temp = gsl_matrix_alloc(COVNO, COVNO);

  gsl_matrix_view rob_var = gsl_matrix_view_array(res->robust_var, COVNO, COVNO);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &invH.matrix, &B_mat.matrix, 0.0, temp);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp, &invH.matrix, 0.0, &rob_var.matrix);

  gsl_matrix_free(temp); free(cluster_L); free(B); free(dLambda); free(Z_bar_mat);

}

void ccox(DATA *dat, DATA_RES *res, int N, int COVNO, double **Z, int MAX_ITER,
          double TOLERANCE, int event_code, int robust) {

  int covN = COVNO; // number of covrs

  gsl_vector *beta = gsl_vector_calloc(covN);

  gsl_vector *U = gsl_vector_alloc(covN);

  gsl_vector *delta = gsl_vector_alloc(covN);

  gsl_matrix *I = gsl_matrix_alloc(covN, covN);

  gsl_permutation *prmute = gsl_permutation_alloc(covN);

  double *temp_times = (double *)malloc(N * sizeof(double));

  int count = 0;


  for (int i = 0; i < N; i++) {

    if (dat->status[i] == event_code) temp_times[count++] = dat->stop[i];

  }


  qsort(temp_times, count, sizeof(double), cmp_double);
  
  double *TiE1 = NULL;

  int E1 = 0;


  if (count > 0) {

    TiE1 = (double *)malloc(count * sizeof(double));

    TiE1[E1++] = temp_times[0];

    for (int i = 1; i < count; i++) {

      if (temp_times[i] != temp_times[i-1]) TiE1[E1++] = temp_times[i];

    }

  }

  free(temp_times);

  if (E1 == 0) {

    fprintf(stderr, "No events found!\n");

    return;
  }

  for (int iter = 0; iter < MAX_ITER; iter++) {

    double U_arr[COVNO];

    double I_arr[COVNO][COVNO];

    U_I_Calc(dat, N, COVNO, beta->data, U_arr, I_arr, Z, TiE1, E1, event_code);

    for (int i = 0; i < covN; i++) {

      gsl_vector_set(U, i, U_arr[i]);

      for (int j = 0; j < covN; j++)

	gsl_matrix_set(I, i, j, I_arr[i][j]);

    }

    int signval;

    gsl_matrix *I_copy = gsl_matrix_alloc(covN, covN);

    gsl_matrix_memcpy(I_copy, I); /// copy of the info matrix for beta and delta calc

    gsl_linalg_LU_decomp(I_copy, prmute, &signval);

    gsl_linalg_LU_solve(I_copy, prmute, U, delta);

    gsl_matrix_free(I_copy);

    gsl_vector_add(beta, delta); // here beta = beta+delta

    if (gsl_blas_dasum(delta) < TOLERANCE)

      break; // gsl dasum absolute sum
  }

  int s;

  gsl_matrix_view invH = gsl_matrix_view_array(res->inv_hessian, covN, covN);

  gsl_matrix_memcpy(I, I);

  gsl_linalg_LU_decomp(I, prmute, &s);

  gsl_linalg_LU_invert(I, prmute, &invH.matrix);

  memcpy(res->betavals, beta->data, sizeof(double) * covN);

  if (robust) {

    compute_robust_variance(dat, res, N, COVNO, Z, TiE1, E1, event_code);

  }

  // free mem
  free(TiE1);

  gsl_vector_free(beta);
  gsl_vector_free(U);
  gsl_vector_free(delta);
  gsl_matrix_free(I);
  gsl_permutation_free(prmute);
}

void compute_cif(DATA *dat, int N, int event_code, int num_predict_times, double predict_times[]) {

  if (num_predict_times == 0) return;

  double *event_times = (double *)malloc(N * sizeof(double)); 

  int count = 0;

  for (int i = 0; i < N; i++) {

    if (dat->status[i] > 0) event_times[count++] = dat->stop[i];
  }

  if (count == 0) { free(event_times); return; }
  
  qsort(event_times, count, sizeof(double), cmp_double);

  int unique_events = 0;

  for (int i = 0; i < count; i++) {

    if (i == 0 || event_times[i] != event_times[i-1]) event_times[unique_events++] = event_times[i];

  }
  
  double *cif1_history = (double *)calloc(unique_events, sizeof(double));

  double *S_prev_history = (double *)calloc(unique_events, sizeof(double));

  double *d1_history = (double *)calloc(unique_events, sizeof(double));

  double *d2_history = (double *)calloc(unique_events, sizeof(double));

  double *n_history = (double *)calloc(unique_events, sizeof(double));

  double S_prev = 1.0, cif1 = 0.0;

  for (int m = 0; m < unique_events; m++) {

    double t = event_times[m];

    int n_j = 0, d1_j = 0, d2_j = 0; 

    for (int i = 0; i < N; i++) {

      if (dat->start[i] < t && dat->stop[i] >= t) n_j++;

      if (dat->stop[i] == t && dat->status[i] > 0) {

	if (dat->status[i] == event_code) d1_j++; else d2_j++; 

     }
    }

    n_history[m] = n_j; d1_history[m] = d1_j; d2_history[m] = d2_j; S_prev_history[m] = S_prev;

    if (n_j > 0) {

      cif1 += S_prev * ((double)d1_j / n_j);

      S_prev = S_prev * (1.0 - (double)(d1_j + d2_j) / n_j);
    }

    cif1_history[m] = cif1;
  }

  printf("\n--- Unadjusted Cumulative Incidence (A-J) ---\n");

     printf("%-10s %-15s %-20s\n", "Time", "CIF", "95% CI");

 printf("-------------------------------------------------\n");

  for (int p = 0; p < num_predict_times; p++) {

    double target_time = predict_times[p];

    int best_m = -1;

    for (int m = 0; m < unique_events; m++) {

      if (event_times[m] <= target_time) best_m = m;

      else

	break;

    }

    if (best_m == -1) {
      
      printf("%-10.2f %-15.4f (%0.4f, %0.4f)\n", target_time, 0.0, 0.0, 0.0); 

      continue; 
    }

    double c_val = cif1_history[best_m];

    double var_cif = 0.0;

    for (int j = 0; j <= best_m; j++) {

      double n_j = n_history[j];

      double d1_j = d1_history[j];
      
      double d2_j = d2_history[j];

        double d_j = d1_j + d2_j;

      if (n_j - d_j <= 0)
	continue; 
      
      double term1 = pow(c_val - cif1_history[j], 2) * d_j / (n_j * (n_j - d_j));
      
      double term2 = pow(S_prev_history[j], 2) * d1_j * (n_j - d1_j) / pow(n_j, 3);

      double term3 = 2.0 * (c_val - cif1_history[j]) * S_prev_history[j] * d1_j / pow(n_j, 2);

      var_cif += term1 + term2 - term3; //final variance term CI for CIF

    }

    double se = sqrt(fmax(0.0, var_cif));

    printf("%-10.2f %-15.4f (%0.4f, %0.4f)\n", target_time, c_val, fmax(0.0, c_val - 1.96 * se), fmin(1.0, c_val + 1.96 * se));

  }
  
  free(event_times); free(cif1_history); free(S_prev_history); free(d1_history); free(d2_history); free(n_history);
  
}
