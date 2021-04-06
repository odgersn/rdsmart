// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <iostream>
#include <algorithm>
#include <math.h>
using namespace Rcpp;
using namespace std;

//' C++ Sorting algorithm
//' 
//' @param x A numeric vector to be sorted
//' @param decreasing set to TRUE or FALSE to sort in decreasing order
//' @param nalast Place NA values at the end of the string?
//' @export
// [[Rcpp::export]]
NumericVector sort_cpp(NumericVector x, bool decreasing = false, bool nalast = false) {
  NumericVector y = clone(x);
  if(nalast) {
    if(decreasing) {
      std::sort(std::begin(y), std::end(y),
                [](double d0, double d1) {
                  if( isnan(d0) ) return false;
                  if( isnan(d1) ) return true;
                  return d0 > d1;
                });
    } else {
      std::sort(std::begin(y), std::end(y),
                [](double d0, double d1) {
                  if( isnan(d0) ) return false;
                  if( isnan(d1) ) return true;
                  return d0 < d1;
                });
    }
  } else {
    y.sort(decreasing);
  }
  return y;
}