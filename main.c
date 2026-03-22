#include "ccox.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct { int id; int index; } ClusterMap;
int cmp_cluster(const void *a, const void *b) { return ((ClusterMap*)a)->id - ((ClusterMap*)b)->id; }

// main prog
int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <config_file>\n", argv[0]);
    return 1;
 }

  char csv_file[256 * 8];
  int N = 0;
  int COVNO = 0;
  char cov_line[512 * 4];
  int MAX_ITER = 0;
  double TOLERANCE = 0.0;
  char start_col[256] = "";
  char stop_col[256] = "";
  char status_col[256] = "";
  char cluster_col[256] = "";
  int event_code = 1;
  int robust = 0;
  double cif_times[100];
  int num_cif_times = 0;

  FILE *config = fopen(argv[1], "r");

  
  if (!config) {

    perror("cannot open config file");

    return 1;
  }

  char line[1024 * 4]; // line read buffer


  while (fgets(line, sizeof(line), config)) {

    line[strcspn(line, "\r\n")] = 0; // remove trailing chars

    char *key = strtok(line, "=");
    
    char *val = strtok(NULL, "=");

    // printf("%s\n", val);

    if (!key || !val) continue;

    if (strcmp(key, "MAX_ITER") == 0) MAX_ITER = atoi(val);
    else if (strcmp(key, "TOLERANCE") == 0) TOLERANCE = atof(val);
    else if (strcmp(key, "file") == 0) strcpy(csv_file, val);
    else if (strcmp(key, "n") == 0) N = atoi(val);
    else if (strcmp(key, "covno") == 0) COVNO = atoi(val);
    else if (strcmp(key, "covariates") == 0) strcpy(cov_line, val);
    else if (strcmp(key, "start_col") == 0) strcpy(start_col, val);
    else if (strcmp(key, "stop_col") == 0) strcpy(stop_col, val);
    else if (strcmp(key, "status_col") == 0) strcpy(status_col, val);
    else if (strcmp(key, "cluster_col") == 0) strcpy(cluster_col, val);
    else if (strcmp(key, "event_code") == 0) event_code = atoi(val);
    else if (strcmp(key, "rob_se") == 0) robust = atoi(val);
    else if (strcmp(key, "cif_times") == 0) {

  char *t_tok = strtok(val, ",");

  while (t_tok != NULL && num_cif_times < 100) {

    cif_times[num_cif_times++] = atof(t_tok); // ascii to float

    t_tok = strtok(NULL, ",");

    //  printf("%s\n",t_tok);

  }
    }
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

  if (!fgets(line, sizeof(line), file)) {
    fprintf(stderr, "Empty CSV file!\n");
    fclose(file);
    return 1;
  }

  line[strcspn(line, "\r\n")] = 0;

  char *header_cols[10000];
  int total_csv_cols = 0;
  
  // simple strtok loop pattern
  char *hdr_token = strtok(line, ",");

  while (hdr_token != NULL && total_csv_cols < 10000) {

    header_cols[total_csv_cols++] = strdup(hdr_token);

    hdr_token = strtok(NULL, ",");

  }

  int start_idx = -1, stop_idx = -1, status_idx = -1, cluster_idx = -1;

  int *col_indices = malloc(COVNO * sizeof(int));

  for (int j = 0; j < total_csv_cols; j++) {

    if (strlen(start_col) > 0 && strcmp(header_cols[j], start_col) == 0) start_idx = j;

    if (strcmp(header_cols[j], stop_col) == 0) stop_idx = j;

    if (strcmp(header_cols[j], status_col) == 0) status_idx = j;

    if (strlen(cluster_col) > 0 && strcmp(header_cols[j], cluster_col) == 0) cluster_idx = j;

  }

  for (int i = 0; i < COVNO; i++) {

    col_indices[i] = -1;

    for (int j = 0; j < total_csv_cols; j++) {

      if (strcmp(target_names[i], header_cols[j]) == 0) {

	col_indices[i] = j;

	break;
      }
    }

    if (col_indices[i] == -1) {

      fprintf(stderr, "Error: Covariate '%s' not found in CSV header!\n",

	      target_names[i]);
      return 1;
      
    }
  }

  DATA S1;

  S1.start = malloc(N * sizeof(double));

  S1.stop = malloc(N * sizeof(double));

  S1.status = malloc(N * sizeof(int));

  S1.cluster = malloc(N * sizeof(int));

  // 2D array Z
  double **Z = malloc(N * sizeof(double *));
  for (int i = 0; i < N; i++) {

    Z[i] = malloc(COVNO * sizeof(double));

  }

  

  int *raw_clusters = malloc(N * sizeof(int));
  int row_count = 0;

  while (fgets(line, sizeof(line), file) && row_count < N) {

    int current_col = 0;
    
    // empty trailing chars
    line[strcspn(line, "\r\n")] = 0;
    
    char *val_token = strtok(line, ",");
    
    // printf("lines in csv: %s\n", val_token);

    S1.start[row_count] = 0.0;
    raw_clusters[row_count] = row_count;

    while (val_token != NULL) {

      // S1 structures
      if (current_col == start_idx) S1.start[row_count] = atof(val_token);

      if (current_col == stop_idx) S1.stop[row_count] = atof(val_token);

      if (current_col == status_idx) S1.status[row_count] = atoi(val_token);

      if (current_col == cluster_idx) raw_clusters[row_count] = atoi(val_token);

      // asign to covar matrix Z
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

  N = row_count;

  ClusterMap *cmap = malloc(N * sizeof(ClusterMap));

  for (int i = 0; i < N; i++) {
    
    cmap[i].id = raw_clusters[i]; cmap[i].index = i;
    
  }

  qsort(cmap, N, sizeof(ClusterMap), cmp_cluster);
  
  int current_id = 0;

  if (N > 0) S1.cluster[cmap[0].index] = current_id;

  for (int i = 1; i < N; i++) {

    if (cmap[i].id != cmap[i-1].id) current_id++;

    S1.cluster[cmap[i].index] = current_id;
  }


  free(cmap);
  free(raw_clusters);

  DATA_RES result;

  result.betavals = malloc(COVNO * sizeof(double));

  result.inv_hessian = malloc(COVNO * COVNO * sizeof(double));

  result.robust_var = malloc(COVNO * COVNO * sizeof(double));

  // ccox call
  ccox(&S1, &result, N, COVNO, Z, MAX_ITER, TOLERANCE, event_code, robust);

  printf("\n%-10s %-10s %-10s %-10s %-10s %-20s\n", "Variable", "Coef", robust ? "Rob_SE" : "SE",
         "p-val", "HR", "95% CI");
  printf("------------------------------------------------------------------------------------\n");

  double z_crit = 1.960;

  double *var_mat = robust ? result.robust_var : result.inv_hessian;

  for (int k = 0; k < COVNO; k++) {

    double beta = result.betavals[k];
    int diag_index = (k * COVNO) + k;
    double variance = var_mat[diag_index];
    double se = (variance > 0) ? sqrt(variance) : 0.0; // no negative vals before sqrt
    double hr = exp(beta);
    double z_stat = (se > 0) ? beta / se : 0.0;
    double p_val = erfc(fabs(z_stat) / sqrt(2.0));
    double ci_low = exp(beta - z_crit * se);
    double ci_high = exp(beta + z_crit * se);

    printf("%-10s %-10.4f %-10.4f %-10.4f %-10.4f (%0.4f, %0.4f)\n",
           target_names[k], beta, se, p_val, hr, ci_low, ci_high);
  }
  
  printf("------------------------------------------------------------------------------------\n");

  //call compute cif
  compute_cif(&S1, N, event_code, num_cif_times, cif_times); 

  // free mem
  for (int i = 0; i < total_csv_cols; i++) {
    free(header_cols[i]);
  }

  for (int i = 0; i < COVNO; i++) {
    free(target_names[i]);
  }
  
  free(target_names);
  free(col_indices);
  free(S1.start);
  free(S1.stop);
  free(S1.status);
  free(S1.cluster);

  for (int i = 0; i < N; i++) {
    free(Z[i]);
  }
  
  free(Z);
  free(result.betavals);
  free(result.inv_hessian);
  free(result.robust_var);

  return 0;
}
