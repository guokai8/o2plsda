#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' @title get the sd for a matrix
//' @keywords internal 
//' @useDynLib o2plsda
//' @return standard deviation of the matrix 
// [[Rcpp::export]]
NumericVector colsds(arma::mat &X){
    arma::mat sd = arma::stddev(X, 0, 0);
    return Rcpp::wrap(sd);
}
//' @title sum square of a matrix
//' @keywords internal 
//' @useDynLib o2plsda
//' @return sum square value of the vector
// [[Rcpp::export]]
double s2(arma::mat x){
    return(arma::accu(pow(x,2)));
}

//' @title calcualte the Q value
//' @keywords internal 
//' @useDynLib o2plsda
//' @return a numeric
// [[Rcpp::export]]
double Q(arma::mat y, arma::mat y_hat) {
    arma::mat diff = y - y_hat;
    return(1-s2(diff)/s2(y));
}

//' @title get the total values for a matrix
//' @keywords internal 
//' @useDynLib o2plsda
//' @return sum value of the column
// [[Rcpp::export]]
NumericVector column_sums(arma::mat X){
    arma::mat cols = arma::sum(X);
    return Rcpp::wrap(cols);
}
// modified from https://github.com/l-ramirez-lopez/resemble/blob/master/src/regression_methods.cpp
// since we just have a y vector


//' @title orthogonal scores algorithn of partial leat squares (opls) projection
//' @keywords internal 
//' @useDynLib o2plsda
//' @return a list of opls result
// [[Rcpp::export]]
List opls(arma::mat &X, 
            arma::mat &Y, 
            int ncomp,
            double maxiter,
            double tol){
        
        int ny = Y.n_cols;
        arma::mat weights = arma::zeros(ncomp, X.n_cols);
        arma::mat scores = arma::zeros(X.n_rows, ncomp);
        arma::mat Xloadings = arma::zeros(ncomp, X.n_cols);
        arma::mat Yloadings = arma::zeros(ncomp, ny);
        arma::mat explained_var = arma::zeros(3, ncomp);
        arma::mat yex = arma::zeros(ny, ncomp);
        arma::mat Xt = X;
        arma::mat Ypls = Y;
        arma::mat yp;
        arma::mat Yplsb;
        
        //variance of Xt
        double xvar = sum(pow(colsds(Xt), 2));  
        
        // matrices to declare
        int j;
        bool step;  
        arma::mat previous_ts = arma::zeros(Xt.n_rows, 1);
        arma::mat lb;
        arma::mat cr;
        arma::mat ts;
        arma::mat w;
        arma::mat p;
        arma::mat q;
        arma::mat cx;
        arma::mat cy;
        arma::mat projection_matrix;
        arma::mat sratio = arma::zeros(weights.n_rows, weights.n_cols);
        for (int i = 0; i < ncomp; i++) {
            Yplsb = Ypls;
            // Select the Y variable with the largest standard deviation
            yp = Ypls.col(0);
            previous_ts.fill(0);
            
            j = 0;
            step = true;
            while (step) {
                if (j > 0) {
                    previous_ts = ts;
                }
                //Step 1: Compute a vector of loading weights...
                // 1.1 Compute the 'scaling factor'
                cr = sqrt(trans(yp) * Xt * trans(Xt) * yp);
                // 1.2 The weights are computed as the cross product of
                // X0 and Y0 divided by the 'scaling factor'...
                w = (trans(Xt) * yp) / repmat(cr, Xt.n_cols, 1);
                // Step 2: Compute the scores...
                ts = Xt * w;
                // Step 3: Compute the X-loadings (p) and the Y-loadings (q)...
                p = (trans(Xt) * ts) / repmat((trans(ts) * ts), Xt.n_cols, 1);
                q = (trans(Yplsb) * ts) / repmat((trans(ts) * ts), Yplsb.n_cols, 1);
                yp = (Yplsb * q) / repmat((trans(q) * q), Xt.n_rows, 1) ;
                lb = abs(sum((ts - previous_ts) / ts));
                step = lb[0] > tol;
                j = j + 1;
                if(maxiter <= j) {
                    step = false;
               }
            }
            // Step 4: The residual matrix
            // of X is finally computed and...
            cx = ts * trans(p);
            Xt = Xt - cx;
            // ... the vector of residuals of Y is also computed
            cy = ts * trans(q);
            Ypls = Ypls - cy;
            // save the matrices corresponding to the loadings
            // and scores..
            weights.row(i) = trans(w);
            scores.col(i) = ts;
            Xloadings.row(i) = trans(p);
            Yloadings.row(i) = trans(q);
            explained_var(0,i) = arma::var(scores.col(i));
            explained_var(1,i) = explained_var(0,i)/xvar;
            explained_var(2,i) = sum(explained_var.row(0)) / xvar;
        }
        
        // convert this to standard deviation
        explained_var.row(0) = sqrt(explained_var.row(0));
        projection_matrix = trans(weights) * arma::solve(Xloadings * trans(weights), arma::eye(Xloadings.n_rows, Xloadings.n_rows));
        
        arma::mat yexi;
        arma::mat cop;
        for (int i = 0; i < ny; i++) {
            yexi = scores % arma::repmat(trans(Yloadings.col(i)), scores.n_rows, 1) ;
            cop = pow(arma::cor(Y.col(i), yexi.col(0)), 2);
            yex(i,0) = cop(0,0);
            for(int j = 1; j < ncomp; j++){
                yexi.col(j) = yexi.col(j-1) + yexi.col(j);
                cop = arma::cor(Y.col(i), yexi.col(j));
                yex(i,j) = pow(cop(0,0), 2);
            }
        }
        
        arma::mat ss = arma::zeros(Yloadings.n_rows, Yloadings.n_cols);
        arma::mat wss = arma::zeros(Yloadings.n_rows, Yloadings.n_cols);
        arma::mat ssw = arma::zeros(scores.n_rows, Yloadings.n_rows);
        arma::mat vip = arma::zeros(weights.n_rows, weights.n_cols);
        arma::mat ss1 = pow(Yloadings, 2);
        arma::mat ss2 = column_sums(pow(scores, 2));
        ss = ss1 % ss2;
        
        arma::mat sqweights = trans(pow(weights, 2));
        wss = column_sums(sqweights);
        ssw = sqweights % trans(arma::repmat(ss/wss, 1, weights.n_cols));
        arma::mat cssw = arma::zeros(weights.n_cols, weights.n_rows);
        double sclr = ssw.n_rows;
        int it = ssw.n_rows;
        for(int i = 0; i < it; i++) {
            cssw.row(i) = cumsum(ssw.row(i));
        }
        cssw = cssw * sclr;
        vip = pow(cssw / arma::repmat(trans(cumsum(ss)), weights.n_cols, 1), 0.5);
        
        return Rcpp::List::create(
            Rcpp::Named("ncomp") = ncomp,
            Rcpp::Named("scores") = scores,
            Rcpp::Named("X_loadings") = Xloadings,
            Rcpp::Named("Y_loadings") = Yloadings,
            Rcpp::Named("projection_mat") = projection_matrix,
            Rcpp::Named("vip") = vip,
            Rcpp::Named("Y") = Y,
            Rcpp::Named("variance") = Rcpp::List::create(
                Rcpp::Named("x_var") = explained_var,
                Rcpp::Named("y_var") = yex
            ),
            _["weights"] = weights
        );
    }
    
