#include "Rcpp.h"
#include "scuttle/downsample_vector.h"
#include <stdexcept>

//[[Rcpp::export]]
Rcpp::IntegerVector downsample_run(Rcpp::IntegerVector reads, double prop) {
    Rcpp::IntegerVector output(reads.size());
    scuttle::downsample_vector(reads.begin(), reads.end(), output.begin(), prop);
    return output;
}
