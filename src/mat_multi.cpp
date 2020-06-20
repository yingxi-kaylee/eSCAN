#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// Reference: SCANG and GENESIS package
#include <math.h>

using namespace Rcpp;

double Saddle(double q, arma::vec egvalues)
{
  double K2(double x, arma::vec egvalues);
  double Bisection(arma::vec egvalues, double q, double xmin, double xmax);
  double K(double x, arma::vec egvalues);
  
  bool lower = false;
  bool logp = false;
  double xmin = 0.0;
  double xmax = 0.0;
  double v = 0.0;
  double res = 0.0;
  
  int i;
  const int en = egvalues.size();
  double lambdamax = max(egvalues);
  q = q/lambdamax;
  
  for(i = 0; i < en; i++)
  {
    egvalues[i] = egvalues[i]/lambdamax;
  }
  lambdamax = 1.0;
  
  if (q > arma::sum(egvalues))
  {
    xmin = -0.01;
  }else
  {
    xmin = -en/(2 * q);
  }
  
  xmax = 1/(2*lambdamax) * 0.99999;
  
  double xhat = Bisection(egvalues,q,xmin,xmax);
  double w = sqrt(2*(xhat*q-K(xhat,egvalues)));
  if (xhat < 0){
    w = -w;
  }
  
  v = xhat*sqrt(K2(xhat,egvalues));
  if(fabs(xhat)<1e-04)
  {
    res = 2;
  }else
  {
    res =  R::pnorm(w+log(v/w)/w,0.0,1.0,lower,logp);
  }
  return res;
  
  
}

double K(double x, arma::vec egvalues)
{
  double res = 0.0;
  const int en = egvalues.size();
  
  for(int i = 0; i < en; i++)
  {
    res = res + log(1-2*egvalues(i)*x);
  }
  
  res = res*(-0.5);
  
  return res;
}

double K1(double x, arma::vec egvalues, double q)
{
  double res = 0.0;
  const int en = egvalues.size();
  
  for(int i = 0; i < en; i++)
  {
    res = res + egvalues(i)/(1-2*egvalues(i)*x);
  }
  
  res = res - q;
  
  return res;
}

double K2(double x, arma::vec egvalues)
{
  double res = 0.0;
  const int en = egvalues.size();
  
  for(int i = 0; i < en; i++)
  {
    res = res + pow(egvalues(i),2)/pow(1-2*egvalues(i)*x,2.0);
  }
  
  res = res*2;
  
  return res;
}

double Bisection(arma::vec egvalues, double q, double xmin, double xmax)
{
  double K1(double x, arma::vec egvalues, double q);
  
  // the range of x to search
  double xupper = xmax;
  double xlower = xmin;
  
  double x0 = 0.0;
  double K1x0 = 1.0;
  
  
  while (fabs(xupper-xlower) > 1e-08)
  {
    x0 = (xupper + xlower)/2.0;
    K1x0 = K1(x0,egvalues,q);
    
    if(K1x0 == 0){
      break;
    }else if (K1x0 > 0){
      xupper = x0;
    }else{
      xlower = x0;
    }
  }
  
  return x0;
}

// Integer log2
int ilog2(int x)
{
  int l2 = 0;
  for (; x; x >>=1) ++l2;
  return l2;
}

// (a,b) -> (a+b,a-b) without overflow
arma::mat rotate(arma::mat x, int a, int b)
{
  static double t;
  t = x(a);
  x(a) = x(a) + x(b);
  x(b) = t - x(b);
  return x;
}


// Fast Walsh-Hadamard transform

arma::mat fwht(arma::mat data, int *size )
{
  const int l2 = ilog2(*size) - 1;
  arma::mat res;
  res = data;
  
  for (int i = 0; i < l2; ++i)
  {
    for (int j = 0; j < (1 << l2); j += 1 << (i+1))
      for (int k = 0; k < (1 << i ); ++k)
        res = rotate(res, j + k, j + k + (1<<i) );
  }
  return res;
}

// [[Rcpp::export]]
SEXP mfwht(arma::mat data, int nrow, int ncol){
  arma::mat res;
  res = data;
  
  for(int i=0; i < ncol; i++){ 
    res.col(i) = fwht(res.col(i), &nrow);
  }
  
  return Rcpp::wrap(res);
}

// [[Rcpp::export]]
SEXP big_mfwht(arma::mat data){

  int nrow = data.n_rows;
  int ncol = data.n_cols;

  for(int i=0; i < ncol; i++){
    data.col(i) = fwht(data.col(i), &nrow);
  }

  return Rcpp::wrap(data);
}


// [[Rcpp::export]]
SEXP compPval(double Q, arma::mat Covw_sub){
  
  //declare Saddle 
  double Saddle(double, arma::vec);
  int rr;
  double pval_saddle;
  // variants number in the subset
  int p = Covw_sub.n_cols;
  
  // eigenvalue_vec
  arma::vec eigenvals;
  eigenvals.zeros(p);
  
  eigenvals = arma::eig_sym(Covw_sub);
  for(rr = 0; rr < p; rr++)
  {
    if(eigenvals[rr] < 1e-8)
    {
      eigenvals[rr] = 0.0;
    }
  }
  
  // saddlepoint method
  pval_saddle = Saddle(Q,eigenvals);
  
  return Rcpp::wrap(pval_saddle);
}

// [[Rcpp::export]]
SEXP compCov(arma::mat G, arma::mat X, arma::vec working, double sigma, int fam){
  // variants number
  int p = G.n_cols;
  // covariates number
  int q = X.n_cols;
  
  
  // t(X)*G
  arma::mat tX_G;
  tX_G.zeros(q,p);
  
  // Cov: the covariance matrix of pesudo-score
  arma::mat Cov;
  Cov.zeros(p,p);
  
  
  // compute pseudo-score matrix x(dim[1]=B, dim[2]=p)
  // compute covariance matrix of the pseudo-score matrix
  if(fam == 0)
  {
    tX_G = trans(X)*G;
    Cov = trans(G)*G - trans(tX_G)*inv(trans(X)*X)*tX_G;
  }else
  {
    tX_G = trans(X)*(arma::diagmat(working))*G;
    Cov = trans(G)*arma::diagmat(working)*G - trans(tX_G)*inv(trans(X)*arma::diagmat(working)*X)*tX_G;
  }
  
  Cov = Cov*pow(sigma,2);
  
  return Rcpp::wrap(Cov);
}

// [[Rcpp::export]]
SEXP compx(arma::mat G, arma::mat X, arma::vec working, double sigma, int fam, int times){
  // variants number
  int p = G.n_cols;
  // sample size
  int n = G.n_rows;
  // covariates number
  int q = X.n_cols;
  
  // pesudo-residual
  arma::mat y;
  y = arma::randn<arma::mat>(times,n);

  // pesudo-score
  arma::mat x;
  x.zeros(times,p);
  
  // t(X)*G
  arma::mat tX_G;
  tX_G.zeros(q,p);
  
  
  // compute pseudo-score matrix x(dim[1]=B, dim[2]=p)
  // compute covariance matrix of the pseudo-score matrix
  if(fam == 0)
  {
    tX_G = trans(X)*G;
    x = (y*G - y*X*inv(trans(X)*X)*tX_G)*sigma;
  }else
  {
    tX_G = trans(X)*(arma::diagmat(working))*G;
    y = y*arma::diagmat(sqrt(working));
    x = (y*G - y*X*inv(trans(X)*arma::diagmat(working)*X)*tX_G)*sigma;
    
  }
  
  return Rcpp::wrap(x);
}

// [[Rcpp::export]]
SEXP compCovw(int p, int w_num, arma::mat Cov, arma::mat weights){
  
  // Cov matrix with weights
  arma::mat Covw;
  Covw.zeros(w_num*p,p);
  
  // Weights Matrix(intermediate value)
  arma::mat W;
  W.zeros(p,p);
  
  //intermediate parameters
  int k, kk;
  arma::vec weights_vec;
  weights_vec.zeros(p);
  
  //compute Covw
  for(k = 0; k < w_num; k++)
  {
    for(kk = 0; kk < p; kk++)
    {
      weights_vec[kk] = weights(kk,k);
    }
    
    // 'each_col()' is to apply a vector operation to each column or row of a matrix
    W.each_col() = weights_vec;
    Covw(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1)) = W%Cov; 
    // % is element-wise multiplication
    // multiple kth row by weights_vec[k]
    
    W.each_row() = trans(weights_vec);
    Covw(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1)) = Covw(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1))%W;
    // multiple kth col by weights_vec[k]
    
  }
  
  return Rcpp::wrap(Covw);
}


// [[Rcpp::export]]
SEXP matrix_flip(arma::mat G) {
  
  int i,j;
  int n = G.n_rows;
  int p = G.n_cols;
  
  arma::vec AF;
  AF.zeros(p);
  
  arma::vec MAF;
  MAF.zeros(p);
  
  double num = 0;
  
  List res(3);
  
  // Calculate AF
  for(i = 0; i < p; i++)
  {
    num = 0;
    
    for(j = 0; j < n; j++)
    {
      if(G(j,i) > -1)
      {
        AF(i) = AF(i) + G(j,i);
        num = num + 1;
      }
    }
    
    AF(i) = AF(i)/2/num;
  }
  
  // Genotype Imputation
  for(i = 0; i < p; i++)
  {
    if(AF(i) <= 0.5)
    {
      for(j = 0; j < n; j++)
      {
        if(!(G(j,i) > -1))
        {
          G(j,i) = 0;
        }
      }
    }else
    {
      for(j = 0; j < n; j++)
      {
        if(!(G(j,i) > -1))
        {
          G(j,i) = 2;
        }
      }
    }
  }
  
  // Genotype Flip
  for(i = 0; i < p; i++)
  {
    if(AF(i) <= 0.5)
    {
      MAF(i) = AF(i);
    }else
    {
      MAF(i) = 1 - AF(i);
      
      for(j = 0; j < n; j++)
      {
        G(j,i) = 2 - G(j,i);
      }
    }
  }
  
  res = List::create(Named("Geno") = G, Named("AF") = AF, Named("MAF") = MAF);
  
  return Rcpp::wrap(res);
  
}

// [[Rcpp::export]]
SEXP CCT_pval(arma::vec x, arma::vec weights)
{
  double cct_stat = 0.0;
  double pval = 0.0;
  
  bool lower = false;
  bool logp = false;
  
  //after normalizing, weights=1/(size of A)
  weights = weights/sum(weights);
  
  //the size of A in paper  
  int n = x.size();
  int k ;
  
  
  for(k = 0;k < n;k++)
  {
    if(x(k) < 1e-16)
    {
      cct_stat = cct_stat + weights(k)/x(k)/PI;
    }else
    {
      cct_stat = cct_stat + weights(k)*tan((0.5-x(k))*PI);
    }
  }
  
  if (cct_stat > 1e+15)
  {
    pval = (1/cct_stat)/PI;
  }else
  {
    pval = R::pcauchy(cct_stat,0.0,1.0,lower,logp);
  }
  
  return Rcpp::wrap(pval);
  
}