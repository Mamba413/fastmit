#include <RcppArmadillo.h>
using namespace std;



// Use QuickSort to find the K smallest elements in the array

int partition(arma::rowvec& a, int s,int t,int &k, arma::uvec& indexa)
{
  int i,j,y;
  double x;
  x = a(s);
  y = indexa(s); // Take the partition element
  i = s;        // Set the initial value for ptr
  j = t;
  do
  {
    while((a(j) > x) && i < j) j--;   // Scan from right to left
    if(i < j) {
      a(i) = a(j);
      indexa(i) = indexa(j);
      i++;
    }           // Move the small elements to the left
    while((a(i) <= x) && i < j) i++;      // Scan from left to right
    if(i < j) {
      a(j) = a(i);
      indexa(j) = indexa(i);
      j--;
    }            // Move the large elements to the right
  }while(i < j); // Till i == j
  a(i) = x;
  indexa(i) = y; // get the partition element
  k = i;
  return k;
}

// Find K smallest elements in the array
// index: return the index of the K_th smallest elements (start from 0), and high represents the largest index in the array
int FindKMin(arma::rowvec& a, int low, int high, int k, arma::uvec& indexa)
{
  int q;
  int index = -1;
  if(low < high)
  {
    partition(a, low, high, q, indexa);
    int len = q - low + 1; // the position
    if(len == k)
      index = q; // return the k_th position
    else if(len < k)
      index = FindKMin(a, q + 1, high, k-len, indexa);
    else
      index = FindKMin(a, low, q - 1, k, indexa);
  }
  return index;
}

void QuickSort(arma::rowvec& a, int low, int high, arma::uvec& indexa)
{
  int q;
  if(low < high)
  {
    int pivot = partition(a, low, high, q, indexa);
    QuickSort(a, low, pivot-1, indexa);
    QuickSort(a, pivot+1, high, indexa);
  }
}
