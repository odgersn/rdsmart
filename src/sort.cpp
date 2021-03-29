// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <iostream>
#include <algorithm>
#include <math.h>
using namespace Rcpp;
using namespace std;

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