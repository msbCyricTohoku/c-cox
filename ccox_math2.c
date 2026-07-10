/*==============================================================================
 c-cox: high performance cox regression
 Version 2.0.0
 ==============================================================================
 Description: A high performance C implementation of Cox regression.

 Authors:     msb
 
 License:     Distributed under the GNU General Public License (GPL)
==============================================================================*/
#include "ccox2.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/*-----------------------comparison functions------------------------------*/
int cmp_double(const void *a, const void *b) {
  
  double da = *(const double *)a; /*------outer * derefrenecing and cast------*/
  double db = *(const double *)b;

  return (da > db) - (da < db); /*------1 if da > db, 0 if equal, -1 if da < db-----*/
}

/*-----------------------ccox sorting funcs------------------------------*/

int ccox_cmp_desc(const void *a, const void *b) {
  double t1 = ((ccox_SortItem*)a)->time;
  double t2 = ((ccox_SortItem*)b)->time;
  return (t1 < t2) - (t1 > t2); /* descending */
}

int ccox_cmp_asc(const void *a, const void *b) {
  double t1 = ((ccox_SortItem*)a)->time;
  double t2 = ((ccox_SortItem*)b)->time;
  return (t1 > t2) - (t1 < t2); /* ascending */
}


/*-----------------------U and I matrix calc with running sum algorithm------------------*/
void U_I_Calc(DATA *data, int N, int tie_handling, int COVNO, double beta[COVNO], double U[COVNO],
              double I[COVNO][COVNO], double *Z, double *TiE1, int E1, int event_code,
              ccox_SortItem *stop_items, ccox_SortItem *start_items, ccox_SortItem *event_items, int num_events) {

  /*------------init U and I------------------*/
  for (int i = 0; i < COVNO; i++) {
    U[i] = 0.0;
    for (int j = 0; j < COVNO; j++) {
      I[i][j] = 0.0;
    }
  }

  /*-------------precompute exponents--------------------- */
  double *w = (double *)malloc(N * sizeof(double));

/*----------omp parallel here for openmp-------------*/
#pragma omp parallel for
  for (int i = 0; i < N; i++) {
    double xb = 0.0;

    for (int k = 0; k < COVNO; k++) {
      //xb += Z[i][k] * beta[k];
      xb += Z[i * COVNO + k] * beta[k];
      
    }

    w[i] = exp(xb);
  }

  int ptr_stop = 0, ptr_start = 0, ptr_event = 0;
  double S0 = 0.0;
  double S1[COVNO];
  double S2[COVNO][COVNO];
    
  for (int i = 0; i < COVNO; i++) {

    S1[i] = 0.0;

    for (int j = 0; j < COVNO; j++) {

      S2[i][j] = 0.0;
    }
  }

  /*-------backward sweep (m--) ------*/
  for (int m = E1 - 1; m >= 0; m--) {
    
    double t = TiE1[m];

    while (ptr_stop < N && stop_items[ptr_stop].time >= t) {

      int i = stop_items[ptr_stop].id;

      S0 += w[i];

      for (int k = 0; k < COVNO; k++) {

	// S1[k] += w[i] * Z[i][k];  
      S1[k] += w[i] * Z[i * COVNO + k];
	
	/**************calc upper triangle****************/
        for (int l = k; l < COVNO; l++) {

	  // S2[k][l] += w[i] * Z[i][k] * Z[i][l];
       S2[k][l] += w[i] * Z[i * COVNO + k] * Z[i * COVNO + l];
	  

     }
    }
   ptr_stop++;
  }

    while (ptr_start < N && start_items[ptr_start].time >= t) {

      int i = start_items[ptr_start].id;

      S0 -= w[i];

      for (int k = 0; k < COVNO; k++) {

        //S1[k] -= w[i] * Z[i][k];
	S1[k] -= w[i] * Z[i * COVNO + k];
	

        /**************calc upper triangle****************/
        for (int l = k; l < COVNO; l++) {

          //S2[k][l] -= w[i] * Z[i][k] * Z[i][l];
	 S2[k][l] -= w[i] * Z[i * COVNO + k] * Z[i * COVNO + l];
	  
     }
    }
   ptr_start++;
  }

    double tied_w = 0.0;
    double tied_S1[COVNO];
    double tied_S2[COVNO][COVNO];
    
    for (int k = 0; k < COVNO; k++) {
      tied_S1[k] = 0.0;
      for (int l = 0; l < COVNO; l++) tied_S2[k][l] = 0.0;
    }
    
    int event_no = 0;

    while (ptr_event < num_events && event_items[ptr_event].time > t) {
      ptr_event++;
    }
        
    int temp_ptr = ptr_event;
    while (temp_ptr < num_events && event_items[temp_ptr].time == t) {

      int i = event_items[temp_ptr].id;

      tied_w += w[i];
      
      for (int k = 0; k < COVNO; k++) {

        //tied_S1[k] += w[i] * Z[i][k];
	tied_S1[k] += w[i] * Z[i * COVNO + k];
	
        //U[k] += Z[i][k]; /*----covariate-----*/
	U[k] += Z[i * COVNO + k]; /*----flatten version-----*/
        
	  
        for (int l = k; l < COVNO; l++) {
          //tied_S2[k][l] += w[i] * Z[i][k] * Z[i][l];
	 tied_S2[k][l] += w[i] * Z[i * COVNO + k] * Z[i * COVNO + l];
     
      }
      }

      event_no++;
      temp_ptr++;
    }

    /*----breslow/Efron tied events here------ */
    if (tie_handling == 1) { 
      /*---- Breslow -----*/
      if (S0 > 0.0 && event_no > 0) {
        for (int k = 0; k < COVNO; k++) {
      
          double z_bar_k = S1[k] / S0;
          U[k] -= event_no * z_bar_k;

          for (int l = 0; l < COVNO; l++) {
	    /**************calc upper triangle****************/
            double s2_val = (k <= l) ? S2[k][l] : S2[l][k];
	    
            I[k][l] += event_no * ((s2_val / S0) - (z_bar_k * S1[l] / S0));
       }
      }
      }
    }
    else { 
      /* ------------ Efron method ------------ */
      if (S0 > 0.0 && event_no > 0) {
	
 for (int j = 0; j < event_no; j++) {
     double weight = (double)j / event_no;
      double denom = S0 - weight * tied_w;

  for (int k = 0; k < COVNO; k++) {

      double s1_adj_k = S1[k] - weight * tied_S1[k];
      double z_bar_k = s1_adj_k / denom;

       U[k] -= z_bar_k; 

    for (int l = 0; l < COVNO; l++) {


      /**************calc upper triangle****************/
     double s2_val = (k <= l) ? S2[k][l] : S2[l][k];

   double tied_s2_val = (k <= l) ? tied_S2[k][l] : tied_S2[l][k];
             
    double s2_adj_kl = s2_val - weight * tied_s2_val;

   I[k][l] += (s2_adj_kl / denom) - (z_bar_k * (S1[l] - weight * tied_S1[l]) / denom);
          }
        }
        }
      }
    }
   ptr_event = temp_ptr;
  }

  /*--------free mem-------*/
  free(w);
  
} /*-------end of U_I_Calc function----------*/


/*-------------------------compute robust variance with prefix sum optimization------------------*/
void compute_robust_variance(DATA *dat, DATA_RES *res, int N, int COVNO, double *Z, 
                             double *TiE1, int E1, int event_code) {

  double *dLambda = (double *)calloc(E1, sizeof(double));
  double *Z_bar_mat = (double *)calloc(E1 * COVNO, sizeof(double));

  double *w = (double *)malloc(N * sizeof(double));

  ccox_SortItem *stop_items =  (ccox_SortItem *)malloc(N * sizeof(ccox_SortItem));
  
  ccox_SortItem *start_items = (ccox_SortItem *)malloc(N * sizeof(ccox_SortItem));

  ccox_SortItem *event_items = (ccox_SortItem *)malloc(N * sizeof(ccox_SortItem));

  int num_events = 0;

  for (int i = 0; i < N; i++) {

    double xb = 0.0;

    for (int k = 0; k < COVNO; k++)

      //xb += Z[i][k] * res->betavals[k];

     xb += Z[i * COVNO + k] * res->betavals[k];
    

    w[i] = exp(xb);

    stop_items[i].id = i; 
    stop_items[i].time = dat->stop[i];
    start_items[i].id = i; 
    start_items[i].time = dat->start[i];

    if (dat->status[i] == event_code) {

      event_items[num_events].id = i;

      event_items[num_events].time = dat->stop[i];

      num_events++;
  }
  }

  qsort(stop_items, N, sizeof(ccox_SortItem), ccox_cmp_desc);
  qsort(start_items, N, sizeof(ccox_SortItem), ccox_cmp_desc);
  qsort(event_items, num_events, sizeof(ccox_SortItem), ccox_cmp_desc);

  int ptr_stop = 0,ptr_start = 0, ptr_event = 0;
  double S0 = 0.0;
  double S1[COVNO];
  for (int k = 0; k < COVNO; k++) S1[k] = 0.0;

  /*-------dLambda and Z_bar_mat calcs-----------*/
  for (int m = E1 - 1; m >= 0; m--) {

    double t = TiE1[m];

    while (ptr_stop < N && stop_items[ptr_stop].time >= t) {

      int i = stop_items[ptr_stop].id;

      S0 += w[i];

      for (int k = 0; k < COVNO; k++)

        //S1[k] += w[i] * Z[i][k];
      S1[k] += w[i] * Z[i * COVNO + k];
      

      ptr_stop++;
    }

    while (ptr_start < N && start_items[ptr_start].time >= t) {

      int i = start_items[ptr_start].id;

      S0 -= w[i];

      for (int k = 0; k < COVNO; k++)
        //S1[k] -= w[i] * Z[i][k];

       S1[k] -= w[i] * Z[i * COVNO + k];
      

      ptr_start++;

    }

    int event_no = 0;

    while (ptr_event < num_events && event_items[ptr_event].time > t) ptr_event++;
    
    int temp_ptr = ptr_event;
    while (temp_ptr < num_events && event_items[temp_ptr].time == t) {
      event_no++; 
      temp_ptr++;
    }
    ptr_event = temp_ptr;

    if (S0 > 0.0) {
      dLambda[m] = event_no / S0;
      for (int k = 0; k < COVNO; k++) {
        Z_bar_mat[m * COVNO + k] = S1[k] / S0;
      }
    }
  }

  double *cum_dLambda = (double *)calloc(E1, sizeof(double)); /*------cumulative dLambda-------*/

  double *cum_Z_dLambda = (double *)calloc(E1 * COVNO, sizeof(double));

  cum_dLambda[0] = dLambda[0];
  for (int k = 0; k < COVNO; k++) {
    
    cum_Z_dLambda[0 * COVNO + k] = dLambda[0] * Z_bar_mat[0 * COVNO + k];

  }

  for (int m = 1; m < E1; m++) {

   cum_dLambda[m] = cum_dLambda[m-1] + dLambda[m];

  for (int k = 0; k < COVNO; k++) {

    cum_Z_dLambda[m * COVNO + k] = cum_Z_dLambda[(m-1) * COVNO + k] + (dLambda[m] * Z_bar_mat[m * COVNO + k]);

   }
  }

  int max_cluster = 0;
  for (int i = 0; i < N; i++) {

    if (dat->cluster[i] > max_cluster) max_cluster = dat->cluster[i];

  }
  
  int num_clusters = max_cluster + 1;

  double *cluster_L = (double *)calloc(num_clusters * COVNO, sizeof(double));

/*----------omp parallel here for openmp-------------*/
#pragma omp parallel
  {
    double *local_cluster_L = (double *)calloc(num_clusters * COVNO, sizeof(double));
    double L_i[COVNO];

#pragma omp for
    for (int i = 0; i < N; i++) {
      for (int k = 0; k < COVNO; k++) L_i[k] = 0.0;

      int start_m = E1, low = 0, high = E1 - 1;
      while (low <= high) {
    
        int mid = low + (high - low) / 2;

        if (TiE1[mid] > dat->start[i]) {
      
          start_m = mid; high = mid - 1; 
        }

        else { 
          low = mid + 1;


      }
    }

      int stop_m = -1; low = 0; high = E1 - 1;
      while (low <= high) {

        int mid = low + (high - low) / 2;

        if (TiE1[mid] <= dat->stop[i]) {

          stop_m = mid; low = mid + 1; 
        }

        else { 
          high = mid - 1; 
        }
      }

      /*-------------------------------------------*/
      if (start_m <= stop_m && stop_m != -1 && start_m != E1) {
        double sum_dL = cum_dLambda[stop_m] - (start_m > 0 ? cum_dLambda[start_m - 1] : 0.0);

        for (int k = 0; k < COVNO; k++) {

          double sum_Z_dL = cum_Z_dLambda[stop_m * COVNO + k] - (start_m > 0 ? cum_Z_dLambda[(start_m - 1) * COVNO + k] : 0.0);

          //L_i[k] -= w[i] * (Z[i][k] * sum_dL - sum_Z_dL);
	  L_i[k] -= w[i] * (Z[i * COVNO + k] * sum_dL - sum_Z_dL);
	  
      }
     }

      /*-------------------score contrib aka residual----------------------*/
      if (dat->status[i] == event_code) {

        int m_idx = -1; low = 0; high = E1 - 1;

        while (low <= high) {

          int mid = low + (high - low) / 2;

          if (TiE1[mid] == dat->stop[i]) {

            m_idx = mid; break; 
          }

          else if (TiE1[mid] < dat->stop[i]) { 
            low = mid + 1; 
          }

          else {

            high = mid - 1; 
          }
        }

        if (m_idx != -1) {

    for (int k = 0; k < COVNO; k++)
      // L_i[k] += (Z[i][k] - Z_bar_mat[m_idx * COVNO + k]);
    L_i[k] += (Z[i * COVNO + k] - Z_bar_mat[m_idx * COVNO + k]);
	  
   }
  }

      int c = dat->cluster[i];

      for (int k = 0; k < COVNO; k++) {

        local_cluster_L[c * COVNO + k] += L_i[k];
    }
    }

/*----------------thread safe addition critical--------------------*/ 
#pragma omp critical
    {
      for (int c = 0; c < num_clusters; c++) {
        for (int k = 0; k < COVNO; k++) {
          cluster_L[c * COVNO + k] += local_cluster_L[c * COVNO + k];
        }
      }
    }
    free(local_cluster_L);
  }

  /*--------gsl mat B calc--------*/
  double *B = (double *)calloc(COVNO * COVNO, sizeof(double));
  for (int c = 0; c < num_clusters; c++) {
  for (int k = 0; k < COVNO; k++) {
     for (int l = 0; l < COVNO; l++) {

        B[k * COVNO + l] += cluster_L[c * COVNO + k] * cluster_L[c * COVNO + l];

    }
  }
  }

  /*------matrix view of the existing B------*/
  gsl_matrix_view B_mat = gsl_matrix_view_array(B, COVNO, COVNO);
  gsl_matrix_view invH = gsl_matrix_view_array(res->inv_hessian, COVNO, COVNO);
  gsl_matrix *temp = gsl_matrix_alloc(COVNO, COVNO);
  gsl_matrix_view rob_var = gsl_matrix_view_array(res->robust_var, COVNO, COVNO);

  /*--------matrix multiplication of invH * B_mat------------*/
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &invH.matrix, &B_mat.matrix, 0.0, temp);
  
  /*-----------matrix multiplication of invH * rob_var--------*/
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp, &invH.matrix, 0.0, &rob_var.matrix);

  /*------free temp matrix and other from memory-----------*/
  gsl_matrix_free(temp); 
  free(cluster_L); 
  free(B); 
  free(dLambda); 
  free(Z_bar_mat);
  free(stop_items); 
  free(start_items); 
  free(event_items); 
  free(w);
  free(cum_dLambda); 
  free(cum_Z_dLambda);
  
} /*-------end of compute_robust_variance---------*/


/*-----------------------ccox main function called in main------------------------------*/
void ccox(DATA *dat, DATA_RES *res, int N, int tie_handling, int COVNO, double *Z, int MAX_ITER,
          double TOLERANCE, int event_code, int robust) {

  int covN = COVNO; /*----number of covrs-----*/

  /*------allocating vector for beta coeffs. and initializing it------*/
  gsl_vector *beta = gsl_vector_calloc(covN);
  
  /*---------memory allocation for U -- no init as zero-----------*/
  gsl_vector *U = gsl_vector_alloc(covN);

  /*--------memory allocation for dellta -- no init as zero---------*/
  gsl_vector *delta = gsl_vector_alloc(covN);

  /*------memory allocation for I, note square covN x covN matrix -- no init as zero----*/
  gsl_matrix *I = gsl_matrix_alloc(covN, covN);

  /*------allocate permutation object for LU decomp. gsl is needed---------*/
  gsl_permutation *prmute = gsl_permutation_alloc(covN);
  
  /*--------allocate memory for temp times---------*/
  double *temp_times = (double *)malloc(N * sizeof(double));

  int count = 0;

  for (int i = 0; i < N; i++) {
    if (dat->status[i] == event_code) temp_times[count++] = dat->stop[i];
  }

  /*----------qsort function call--------------*/
  qsort(temp_times, count, sizeof(double), cmp_double); /*-----sort ascending order------*/
  
  double *TiE1 = NULL; /*----since pointer set NULL cannot use 0.0 -- i made mistake here before----*/
  int E1 = 0; /*----le event counter----*/

  if (count > 0) {

    TiE1 = (double *)malloc(count * sizeof(double)); /*---pointer mem alloc---*/

    TiE1[E1++] = temp_times[0]; /*--init--*/

    for (int i = 1; i < count; i++) {
      if (temp_times[i] != temp_times[i-1]) TiE1[E1++] = temp_times[i];
    }
  }

  /*--free mem--*/
  free(temp_times);

  /*--no event passed--*/
  if (E1 == 0) {
    fprintf(stderr, "No events found!\n");
    return;
  }

  /*----------pre sort arrays befoere iteration---------*/
  ccox_SortItem *stop_items = (ccox_SortItem *)malloc(N * sizeof(ccox_SortItem));
  ccox_SortItem *start_items = (ccox_SortItem *)malloc(N * sizeof(ccox_SortItem));
  ccox_SortItem *event_items = (ccox_SortItem *)malloc(N * sizeof(ccox_SortItem));

  int num_events = 0;

  for (int i = 0; i < N; i++) {
    stop_items[i].id = i; 
    stop_items[i].time = dat->stop[i];
    start_items[i].id = i; 
    start_items[i].time = dat->start[i];
    
    if (dat->status[i] == event_code) {
      event_items[num_events].id = i; 
      event_items[num_events].time = dat->stop[i]; 
      num_events++;
    }
  }

  qsort(stop_items, N, sizeof(ccox_SortItem), ccox_cmp_desc);
  qsort(start_items, N, sizeof(ccox_SortItem), ccox_cmp_desc);
  qsort(event_items, num_events, sizeof(ccox_SortItem), ccox_cmp_desc);


     /*----------------------scaling bug fix starts---------------------------*/
  /*----I have found issue when continuous variables and binary are mixed----*/
  /*----issue was due to scaling, as exp() values are calculated, large values
   * cause issues----*/
  /*----here I added scaling, I read that python and R does the same-----*/
  double *means = (double *)calloc(COVNO, sizeof(double));
  double *stds = (double *)calloc(COVNO, sizeof(double));

  for (int k = 0; k < COVNO; k++) {
    for (int i = 0; i < N; i++) means[k] += Z[i * COVNO + k];
    means[k] /= N;

    for (int i = 0; i < N; i++) stds[k] += pow(Z[i * COVNO + k] - means[k], 2);
    stds[k] = sqrt(stds[k] / (N > 1 ? (N - 1) : 1));

    /*----prevent division by zero-------*/
    if (stds[k] < 1e-8) stds[k] = 1.0; 

    /*----standardization to Z----*/
    for (int i = 0; i < N; i++) {
      Z[i * COVNO + k] = (Z[i * COVNO + k] - means[k]) / stds[k];
    }
  }

  /*----------------------scaling bug fix end---------------------------*/

  

  for (int iter = 0; iter < MAX_ITER; iter++) {

    double U_arr[COVNO]; /*---U matrix---*/
    double I_arr[COVNO][COVNO]; /*---inverse hessian, I mat---*/

    /*---call U_I_Calc passing in the globally pre-sorted items----*/
    U_I_Calc(dat, N, tie_handling, COVNO, beta->data, U_arr, I_arr, Z, TiE1, E1, event_code,
             stop_items, start_items, event_items, num_events);

    for (int i = 0; i < covN; i++) {

      /*---gsl vector set for U for LU solver---*/
      gsl_vector_set(U, i, U_arr[i]);

      for (int j = 0; j < covN; j++) {
        /*-----set I matrix for gsl LU solve-----*/
        gsl_matrix_set(I, i, j, I_arr[i][j]);
      }
    }

    int signval;
    gsl_matrix *I_copy = gsl_matrix_alloc(covN, covN); /*-----gsl allocate matrix for I-----*/

    gsl_matrix_memcpy(I_copy, I); /*------copy of the info matrix for beta and delta calc-------*/

    gsl_linalg_LU_decomp(I_copy, prmute, &signval); /*---LU decomposition from gsl------*/

    gsl_linalg_LU_solve(I_copy, prmute, U, delta); /*------LU solve gsl here-------*/

    gsl_matrix_free(I_copy); /*----free matrix------*/

    gsl_vector_add(beta, delta); /*--------here beta = beta+delta-------*/

    if (gsl_blas_dasum(delta) < TOLERANCE) /*-------gsl dasum absolute sum------*/
      break; 
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


      /*----------------------scaling bug fix starts---------------------------*/
  /*-------here I un-scale for beta vals------*/
    for (int k = 0; k < covN; k++) {
    /*----unscale beta directly on the output struct----*/
    res->betavals[k] /= stds[k];

    for (int l = 0; l < covN; l++) {
      /*---unscale inverse Hessian mat----*/
      res->inv_hessian[k * covN + l] /= (stds[k] * stds[l]);
      
      /*-----unscale the robust variance mat----*/
      if (robust) {
        res->robust_var[k * covN + l] /= (stds[k] * stds[l]);
      }
    }

    /*----restore Z-----*/
    for (int i = 0; i < N; i++) {
      Z[i * COVNO + k] = (Z[i * COVNO + k] * stds[k]) + means[k];
    }
  }
  free(means);
  free(stds);

  /*--------free and gls free-------*/
  free(stop_items); 
  free(start_items); 
  free(event_items);
  free(TiE1);
  gsl_vector_free(beta);
  gsl_vector_free(U);
  gsl_vector_free(delta);
  gsl_matrix_free(I);
  gsl_permutation_free(prmute);
}

/*---------------------compute cif cumulative incidence function using forward sweep optimization----------------*/
void compute_cif(DATA *dat, int N, int event_code, int num_predict_times, double predict_times[]) {

  if (num_predict_times == 0) return;

  double *event_times = (double *)malloc(N * sizeof(double)); 

  int count = 0;

  for (int i = 0; i < N; i++) {
    if (dat->status[i] > 0) event_times[count++] = dat->stop[i];
  }

  if (count == 0) {
    free(event_times);
    return;
  }

  /*------------qsort ascending----------------------*/
  qsort(event_times, count, sizeof(double), cmp_double);

  int unique_events = 0;

  for (int i = 0; i < count; i++) {
    if (i == 0 || event_times[i] != event_times[i-1]) event_times[unique_events++] = event_times[i];
  }

  /*---------alloc and zero mem for cif calc----------*/
  double *cif1_history = (double *)calloc(unique_events, sizeof(double));

  double *S_prev_history = (double *)calloc(unique_events, sizeof(double));

  double *d1_history = (double *)calloc(unique_events, sizeof(double));

  double *d2_history = (double *)calloc(unique_events, sizeof(double));

  double *n_history = (double *)calloc(unique_events, sizeof(double));

  ccox_SortItem *stop_items = (ccox_SortItem *)malloc(N * sizeof(ccox_SortItem));

  ccox_SortItem *start_items = (ccox_SortItem *)malloc(N * sizeof(ccox_SortItem));

  ccox_SortItem *event_items = (ccox_SortItem *)malloc(N * sizeof(ccox_SortItem));

  int num_events = 0;

  for (int i = 0; i < N; i++) {
    stop_items[i].id = i; 
    stop_items[i].time = dat->stop[i];
    start_items[i].id = i; 
    start_items[i].time = dat->start[i];
    
    if (dat->status[i] > 0) {
      event_items[num_events].id = i; 
      event_items[num_events].time = dat->stop[i]; 
      num_events++;
    }
  }

  qsort(stop_items, N, sizeof(ccox_SortItem), ccox_cmp_asc);
  qsort(start_items, N, sizeof(ccox_SortItem), ccox_cmp_asc);
  qsort(event_items, num_events, sizeof(ccox_SortItem), ccox_cmp_asc);

  int ptr_stop = 0, ptr_start = 0, ptr_event = 0;
  int current_risk = 0;
  double S_prev = 1.0;
  double cif1 = 0.0;

  for (int m = 0; m < unique_events; m++) {
    double t = event_times[m];

    while (ptr_start < N && start_items[ptr_start].time < t) { 
      current_risk++; 
      ptr_start++; 
    }
    while (ptr_stop < N && stop_items[ptr_stop].time < t) { 
      current_risk--; 
      ptr_stop++; 
    }

    int n_j = current_risk;
    int d1_j = 0, d2_j = 0;

    while(ptr_event < num_events && event_items[ptr_event].time < t) {
      ptr_event++;
    }
        
    int temp_ptr = ptr_event;
    while (temp_ptr < num_events && event_items[temp_ptr].time == t) {
      
     int id = event_items[temp_ptr].id;
     if (dat->status[id] == event_code) d1_j++; else d2_j++;

      temp_ptr++;
    }
    ptr_event = temp_ptr;

    n_history[m] = n_j; 
    d1_history[m] = d1_j; 
    d2_history[m] = d2_j; 
    S_prev_history[m] = S_prev;

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
      else break;
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

      if (n_j - d_j <= 0) continue; 
      
      double term1 = pow(c_val - cif1_history[j], 2) * d_j / (n_j * (n_j - d_j));
      double term2 = pow(S_prev_history[j], 2) * d1_j * (n_j - d1_j) / pow(n_j, 3);
      double term3 = 2.0 * (c_val - cif1_history[j]) * S_prev_history[j] * d1_j / pow(n_j, 2);

      var_cif += term1 + term2 - term3; /*--------final variance term CI for CIF------------*/

    }

    double se = sqrt(fmax(0.0, var_cif)); /*------final se------*/

    printf("%-10.2f %-15.4f (%0.4f, %0.4f)\n", target_time, c_val, fmax(0.0, c_val - 1.96 * se), fmin(1.0, c_val + 1.96 * se));

  }
  
  free(event_times); 
  free(cif1_history); 
  free(S_prev_history); 
  free(d1_history); 
  free(d2_history); 
  free(n_history);
  free(stop_items); 
  free(start_items); 
  free(event_items);
}
