#include <RcppArmadillo.h>
#include <cmath>
#include <math.h>
#include "sort_knn.h"
// #include <boost/math/special_functions/digamma.hpp>
#include <Rmath.h>
#include <stdlib.h>
// #include <random>
// #include <time.h>

#ifdef R_BUILD
#include "R.h"
#include "Rinternals.h"
#else

#include "time.h"

#endif

using namespace std;

// [[Rcpp::export]]
double knn_mi(arma::mat datax,
              arma::mat datay,
              int k) {
  // Suppose sample size is N, then datax and datay are N*N matrix, which are the distance matrices of X and Y
  int K     = k+1;
  int N     = datax.n_rows;
  int vars  = 2;

  // Find out the distance matrix for Z(X,Y)
  arma::mat temp_dist(N,N);
  arma::uvec temp_inds(N);
  arma::mat  Z_dist(N,K);
  arma::imat Z_inds(N,K);
  arma::uvec ini_inds(N);
  arma::rowvec temp_dist_i(N);

  for(int i = 0; i < N; i++){
    for(int j = i; j < N; j++){
      double dist_z;
      if(i == j){
        dist_z = 0.0;
      } else{
        dist_z = datax(i,j);
        if(datay(i,j) > dist_z){
          dist_z = datay(i,j);
        }
      }
      temp_dist(i,j) = dist_z;
      temp_dist(j,i) = dist_z;
    }
  }

  for (int i = 0; i < N; i++) {
    ini_inds(i) = i;
  }

  for (int i = 0; i < N; i++) {
    // temp_inds        = arma::sort_index(temp_dist.row(i));
    // temp_dist.row(i) = arma::sort(temp_dist.row(i));

    temp_inds = ini_inds;
    temp_dist_i = temp_dist.row(i);

    FindKMin(temp_dist_i, 0, N - 1, K, temp_inds);
    QuickSort(temp_dist_i, 0, K - 1, temp_inds);

    for (int j = 0; j < K; j++) {
      Z_dist(i,j) = temp_dist(i,j);
      Z_inds(i,j) = temp_inds(j);
    }
  }

  double digamma_x = 0;
  // double mi = boost::math::digamma(N) + boost::math::digamma(k);
  double mi = R::digamma(N) + R::digamma(k);
  int N_x;
  double epsilon;
  double dist;

  for (int i = 0; i < N; i++) { // for point i
    for (int j = 1; j <= vars; j++) { // count marginal sum for jth block
      N_x = 0;
      epsilon = 0;

      for (int m = 1; m < K; m++) { // for all these k points, we should find out the maximum marginal dist
        dist = 0;
        if(j == 1){
          dist = datax(i,Z_inds(i,m));
        } else {
          dist = datay(i,Z_inds(i,m));
        }
        if (dist > epsilon) {
          epsilon = dist;
        }
      }

      for (int m = 0; m < N; m++) { // iterate over all points
        dist = 0;
        if(j == 1){
          dist = datax(i,m);
        } else {
          dist = datay(i,m);
        }
        if (dist <= epsilon) N_x++;
      }
      // digamma_x = digamma_x + boost::math::digamma(N_x-1);
      digamma_x = digamma_x + R::digamma(N_x-1);
    }
  }
  mi = mi - digamma_x/(double)N - (vars - 1)/(double)k ;


  return mi;
}

void resample(arma::ivec& i_perm, int N) {
#ifdef R_BUILD
  GetRNGstate();
  for(int i = N - 1; i > 0; --i) {
    int j = ((int) round(RAND_MAX * unif_rand())) % (i + 1);
    int temp = i_perm(j);
    i_perm(j) = i_perm(i);
    i_perm(i) = temp;
  }
  PutRNGstate();
#else
  srand((unsigned) time(NULL));
  for(int i = N - 1; i > 0; --i) {
    int j = rand() % (i + 1);
    int temp = i_perm(j);
    i_perm(j) = i_perm(i);
    i_perm(i) = temp;
  }
#endif

  // std::default_random_engine e;
}

// [[Rcpp::export]]
double mi_test(arma::mat datax, arma::mat datay, 
               int k, int perm_num, double ini_mi) {
  int N = datax.n_rows;
  arma::ivec i_perm(N);
  arma::mat datay_pert(N, N);
  arma::vec mi_list(perm_num);
  
  for(int i = 0; i < N; i++) {
    i_perm(i) = i;
  }
  
  for(int p = 0; p < perm_num; p++) {
    resample(i_perm, N);

    for(int i = 0; i < N; i++) {
      for(int j = 0; j < N; j++) {
        datay_pert(i, j) = datay(i_perm(i), i_perm(j));
      }
    }
    
    mi_list(p) = knn_mi(datax, datay_pert, k);
  }
  
  int outlier_num = 0;
  for(int i = 0; i < perm_num; i++) {
    if(fabs(mi_list(i)) >= fabs(ini_mi)) {
      outlier_num++;
    }
  }
  double p_value = ((double)(1 + outlier_num)) / ((double)(1 + perm_num));
  
  return p_value;
}
