#include <RcppArmadillo.h>
using namespace std;

int partition(arma::rowvec&, int, int, int&, arma::uvec&);

int FindKMin(arma::rowvec&, int, int, int, arma::uvec&);

void QuickSort(arma::rowvec&, int, int, arma::uvec&);
