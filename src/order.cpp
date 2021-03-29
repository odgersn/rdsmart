// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

template <int RTYPE>
IntegerVector order_impl(const Vector<RTYPE>& x, bool decreasing) {
  auto n = x.size();
  IntegerVector idx = no_init(n);
  std::iota(idx.begin(), idx.end(), static_cast<size_t>(1));
  if (decreasing) {
    auto comparator = [&x](size_t a, size_t b){ return x[a - 1] > x[b - 1]; };
    std::stable_sort(idx.begin(), idx.end(), comparator);
  } else {
    auto comparator = [&x](size_t a, size_t b){ return x[a - 1] < x[b - 1]; };
    std::stable_sort(idx.begin(), idx.end(), comparator);
    // simulate na.last
    size_t nas = 0;
    for (size_t i = 0; i < n; ++i, ++nas)
      if (!Vector<RTYPE>::is_na(x[idx[i] - 1])) break;
      std::rotate(idx.begin(), idx.begin() + nas, idx.end());
  }
  return idx;
}

// [[Rcpp::export]]
IntegerVector order_cpp(SEXP x, bool decreasing = false) {
  switch(TYPEOF(x)) {
  case INTSXP: return order_impl<INTSXP>(x, decreasing);
  case REALSXP: return order_impl<REALSXP>(x, decreasing);
  case STRSXP: return order_impl<STRSXP>(x, decreasing);
  default: stop("Unsupported type.");
  }
}