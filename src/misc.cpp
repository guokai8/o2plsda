#include <assert.h>
#include <math.h>
#include <iostream>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;


//' @title sort string
//' @keywords internal 
//' @useDynLib o2plsda
//' @return A vector of string
// [[Rcpp::export]]
CharacterVector sort_str( std::vector< std::string > strings ) {
    
    std::sort( strings.begin(), strings.end() );
    return Rcpp::wrap(strings);
}

//' @title sampling a vector
//' @keywords internal 
//' @useDynLib o2plsda
//' @return a vector of length ‘size’ with element drawn from ‘x’
// [[Rcpp::export]]
IntegerVector sample_cpp(IntegerVector x, int n) {
    return sample(x, n, false); 
}

//' @title calculate RMSE
//' @keywords internal 
//' @useDynLib o2plsda
//' @return root-mean-square error value
// [[Rcpp::export]]
double rcpp_rmse(Rcpp::NumericVector& y, Rcpp::NumericVector& y_hat) {
    Rcpp::NumericVector diff = y - y_hat;
    return sqrt( Rcpp::sum(diff*diff) / y.size() );
}

//' @title order a vector of sting
//' @keywords internal 
//' @useDynLib o2plsda
//' @return An character vector 
// [[Rcpp::export]]
IntegerVector order_str(CharacterVector& x) {
    // Order the elements of x by sorting y
    // First create a vector of indices
    IntegerVector idx = seq_along(x)-1;
    // Then sort that vector by the values of y
    std::sort(idx.begin(), idx.end(), [&](int i, int j){return x[i] < x[j];});
    // And return x in that order
    return idx;
}
//' @title order a vector
//' @keywords internal 
//' @useDynLib o2plsda
//' @return a integer vector
// [[Rcpp::export]]
IntegerVector order_cpp(IntegerVector& x) {
    // Order the elements of x by sorting y
    // First create a vector of indices
    IntegerVector idx = seq_along(x)-1;
    // Then sort that vector by the values of y
    std::sort(idx.begin(), idx.end(), [&](int i, int j){return x[i] < x[j];});
    // And return x in that order
    return idx;
}

//' @title split a vector
//' @keywords internal 
//' @useDynLib o2plsda
//' @return a list of vectors containing the values for the groups
// [[Rcpp::export]]
List split_str(CharacterVector x){
    Function sp("split");
    IntegerVector y = seq_along(x);
    return sp(y,x);
} 

//' @title unlist a list
//' @keywords internal 
//' @useDynLib o2plsda
//' @return a vector of an appropriate mode to hold the list components.
// [[Rcpp::export]]
SEXP unlist_cpp(const List& list)
{
    // https://stackoverflow.com/questions/30175104/how-to-effectively-combine-a-list-of-numericvectors-into-one-large-numericvector
    std::size_t n = list.size();
    
    // Figure out the length of the output vector
    std::size_t total_length = 0;
    for (std::size_t i = 0; i < n; ++i)
        total_length += Rf_length(list[i]);
    
    // Allocate the vector
    NumericVector output = no_init(total_length);
    
    // Loop and fill
    std::size_t index = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        NumericVector el = list[i];
        std::copy(el.begin(), el.end(), output.begin() + index);
        
        // Update the index
        index += el.size();
    }
    
    return output;
    
}

//' @title lapply sampling
//' @keywords internal 
//' @useDynLib o2plsda
//' @return a list
// [[Rcpp::export]]
List sample_lapply(List X, int n){
    int lenx = X.size();
    Rcpp::List res(X);
    IntegerVector x;
    IntegerVector nn = seq_len(n);
    IntegerVector xrep;
    IntegerVector xin;
    for(int i=0;i<lenx;i++){
        x = X[i];
        xrep = rep(nn,std::ceil(x.size()/(double)n));
        xin = xrep[Rcpp::Range(0,x.size()-1)];
        res[i]=sample_cpp(xin,xin.size());
    }
    return(res);
}
//' @title MCCV sampling
//' @keywords internal 
//' @useDynLib o2plsda
//' @return a vector with random sampling based on group-balanced Monte Carlo cross-validation
// [[Rcpp::export]]
IntegerVector getMCCV_cpp(CharacterVector x,int n){
    IntegerVector o = order_str(x);
    CharacterVector groups = sort_str(Rcpp::as<std::vector<std::string>>(x));
    List g1 = split_str(groups);
    List g2 = sample_lapply(g1,n);
    IntegerVector gr = unlist_cpp(g2);
    IntegerVector oo = order_cpp(o);
    IntegerVector res(x.size());
    for(int i=0; i< oo.size();i++){
        res[i]=gr[oo[i]];
    }
    return(res);
}

