/*==============================================================================
 c-cox: high performance cox regression
 Version 2.0.0
 ==============================================================================
 Description: A high performance C implementation of Cox regression.

 Authors:     msb
 
 License:     Distributed under the GNU General Public License (GPL)
==============================================================================*/
#ifndef CCOX_H
#define CCOX_H

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

typedef struct {
  double *start;
  double *stop;
  int *status;
  int *cluster;
} DATA;

typedef struct {
  double *betavals;
  double *inv_hessian;
  double *robust_var;
} DATA_RES;

void U_I_Calc(DATA *data, int N, int COVNO, double beta[COVNO], double U[COVNO],double I[COVNO][COVNO], double **Z, double *TiE1, int E1, int event_code);

void ccox(DATA *dat, DATA_RES *res, int N, int COVNO, double **Z, int MAX_ITER, double TOLERANCE, int event_code, int robust);

void compute_robust_variance(DATA *dat, DATA_RES *res, int N, int COVNO, double **Z, double *TiE1, int E1, int event_code);

void compute_cif(DATA *dat, int N, int event_code, int num_predict_times, double predict_times[]);

#endif
