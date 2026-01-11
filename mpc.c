#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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


double Z[N][P] //data matrix covars are const. at the moment

int main(){

  double T1[] = {};
  double T2[] = {};
  double T3[] = {};

  double D1[] = {};
  double D2[] = {};
  double D3[] = {};


  DATA S1, S2, S3; 


  return 0;

}
