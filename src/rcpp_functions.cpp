#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
  #include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]


// Utility Functions -----------------------------------------------------------
inline arma::uvec getAntiIndices(arma::uword N, arma::uword index)
{
	std::vector< arma::uword > temp;
	for ( arma::uword i = 0; i < N; i++ ){
		if ( i == index )
			continue;
		temp.push_back(i);
	}
	return arma::uvec( &temp[0], temp.size() );
}

// Gaussian-Kernel Multivariate Nonparametric Conditional Expectation-----------

// Multivariate at one point
inline double GaussianMultivarAtOnePoint
( const arma::vec& Y, const arma::mat& X, const arma::rowvec& arg, const arma::rowvec& H )
{
	arma::mat D = X.each_row() - arg;
	arma::mat Z = D.each_row() / H;
	arma::mat Kerh = arma::normpdf( Z ).each_row() / H;
	arma::vec w = arma::prod(Kerh, 1);
	return arma::sum( Y%w ) / arma::sum( w );
}

// Multivariate nonparametric estimator
Rcpp::NumericVector GaussianMultivar
(const arma::vec& Y, const arma::mat& X, const arma::mat& args, arma::rowvec H, int n_cores = 3)
{
	// validity check
	if ( Y.n_rows != X.n_rows ){
		Rcpp::Rcerr << "Unconformable sizes: \"Y\" and \"X\"" << std::endl;
		return Rcpp::NumericVector::create(NA_REAL);
	}
	if ( X.n_cols != args.n_cols || X.n_cols != H.n_cols ){
		Rcpp::Rcerr << "Unconformable sizes: \"X\",  \"args\", and \"H\"" << std::endl;
		return Rcpp::NumericVector::create(NA_REAL);
	}

	// multi-threading
	arma::uword Nargs = args.n_rows, i;
	std::vector< double > ce(Nargs); // output bin
	omp_set_num_threads(n_cores);
	#pragma omp parallel for default(none) shared(Nargs, Y, X, args, H, ce) private(i)
	for ( i = 0; i < Nargs; i++ ){
		ce[i] = GaussianMultivarAtOnePoint(Y, X, args.row(i), H);
	}

	return Rcpp::wrap(ce);
}

// Multivariate nonparametric estimator, leave-one-out version
Rcpp::NumericVector GaussianMultivarLeaveOneOut
(const arma::vec& Y, const arma::mat& X, arma::rowvec H, int n_cores = 3)
// This implements leave-one-out estimation.
{
	// validity check
	if ( Y.n_rows != X.n_rows ){
		Rcpp::Rcerr << "Unconformable sizes: \"Y\" and \"X\"" << std::endl;
		return Rcpp::NumericVector::create(NA_REAL);
	}
	if ( X.n_cols != H.n_cols ){
		Rcpp::Rcerr << "Unconformable sizes: \"X\" and \"H\"" << std::endl;
		return Rcpp::NumericVector::create(NA_REAL);
	}

	// multi-threading
	arma::uword Nargs = X.n_rows, i;
	std::vector< double > ce(Nargs); // output bin
	omp_set_num_threads(n_cores);
	#pragma omp parallel for default(none) shared(Nargs, Y, X, H, ce) private(i)
	for ( i  = 0; i < Nargs; i++ ){
		arma::uvec j = getAntiIndices(Nargs, i);
		ce[i] = GaussianMultivarAtOnePoint(Y.rows(j), X.rows(j), X.row(i), H);
	}

	return Rcpp::wrap(ce);
}

// Gaussian-Kernel Univariate Nonparametric Conditional Expectation-------------

// Univariate at one point
inline double GaussianUnivarAtOnePoint
(const arma::vec& Y, const arma::vec& X, const double& arg, const double& H)
{
	arma::vec Z = ( X - arg ) / H;
	arma::vec w = arma::normpdf( Z ) / H;
	return arma::sum( Y%w ) / arma::sum( w );
}

// Univariate nonparametric estimator
Rcpp::NumericVector GaussianUnivar
(const arma::vec& Y, const arma::vec& X, const arma::vec& args, const double& H, int n_cores = 3)
{
	// validity check
	if ( Y.n_rows != X.n_rows ){
		Rcpp::Rcerr << "Unconformable sizes: \"Y\" and \"X\"" << std::endl;
		return Rcpp::NumericVector::create(NA_REAL);
	}
	if ( X.n_cols != args.n_cols || X.n_cols != 1 ){
		Rcpp::Rcerr << "Unconformable sizes: \"X\" and  \"args\"" << std::endl;
		return Rcpp::NumericVector::create(NA_REAL);
	}

	// multi-threading
	arma::uword Nargs = args.n_rows, i;
	std::vector< double > ce(Nargs); // output bin
	omp_set_num_threads(n_cores);
	#pragma omp parallel for default(none) shared(Nargs, Y, X, args, H, ce) private(i)
	for ( i  = 0; i < Nargs; i++ ){
		ce[i] = GaussianUnivarAtOnePoint(Y, X, args[i], H);
	}

	return Rcpp::wrap(ce);
}

// Univariate nonparametric estimator, leave-one-out version
Rcpp::NumericVector GaussianUnivarLeaveOneOut
(const arma::vec& Y, const arma::vec& X, const double& H, int n_cores = 3)
{
	// validity check
	if ( Y.n_rows != X.n_rows ){
		Rcpp::Rcerr << "Unconformable sizes: \"Y\" and \"X\"" << std::endl;
		return Rcpp::NumericVector::create(NA_REAL);
	}
	if ( X.n_cols != 1 ){
		Rcpp::Rcerr << "Unconformable sizes: \"X\"" << std::endl;
		return Rcpp::NumericVector::create(NA_REAL);
	}

	// multi-threading
	arma::uword Nargs = X.n_rows, i;
	std::vector< double > ce(Nargs); // output bin
	omp_set_num_threads(n_cores);
	#pragma omp parallel for default(none) shared(Nargs, Y, X, H, ce) private(i)
	for ( i = 0; i < Nargs; i++ ){
		arma::uvec j = getAntiIndices(Nargs, i);
		ce[i] = GaussianUnivarAtOnePoint(Y.rows(j), X.rows(j), X[i], H);
	}

	return Rcpp::wrap(ce);
}

// Wrapper functions------------------------------------------------------------

//' Gaussian Kernel Nonparametric Conditional Expectation
//'
//' This calculates the Gaussian kernel nonparametric conditional expectation
//' for given data points using \code{Rcpp} with \code{OpenMP}.
//'
//' @param Y A numeric vector of outcome variable
//' @param X A numeric vector or matrix of explanatory variable(s)
//' @param args A numeric vector or matrix of arguments at which nonparametric conditional
//' expectations are calculated
//' @param H A numeric vector of bandwidth(s)
//' @param n_cores the number of threads to parallize the computation of nonparametric
//' conditional expectation over \code{args}
//' @export
//' @return A numeric vector of nonparametric conditional estimates evaluated at \code{args}
//' @useDynLib semiBRM, .registration = TRUE
//' @importFrom Rcpp evalCpp
// [[Rcpp::export]]
Rcpp::NumericVector GaussianKerNonpar
(SEXP Y, SEXP X, SEXP args, SEXP H, int n_cores = 3)
{
	if ( Rf_isMatrix(X) == TRUE ){
		arma::vec Y_ = Rcpp::as< arma::vec >(Y);
		arma::mat X_ = Rcpp::as< arma::mat >(X);
		arma::mat args_ = Rcpp::as< arma::mat >(args);
		arma::rowvec H_ =  Rcpp::as< arma::rowvec >(H);
		return GaussianMultivar(Y_, X_, args_, H_, n_cores);
	}
	else{
		arma::vec Y_ = Rcpp::as< arma::vec >(Y);
		arma::vec X_ = Rcpp::as< arma::vec >(X);
		arma::vec args_ = Rcpp::as< arma::vec >(args);
		double H_ =  Rcpp::as< double >(H);
		return GaussianUnivar(Y_, X_, args_, H_, n_cores);
	}
}

//' Leave-One-Out Gaussian Kernel Nonparametric Conditional Expectation
//'
//' This calculates the leave-one-out Gaussian kernel nonparametric conditional
//' expectation over data points of explanatory variables using \code{Rcpp} with \code{OpenMP}.
//'
//' @param Y A numeric vector of outcome variable
//' @param X A numeric vector or matrix of explanatory variable(s)
//' @param H A numeric vector of bandwidth(s)
//' @param n_cores the number of threads to parallize the computation of nonparametric
//' conditional expectation over \code{X}
//' @export
//' @return A numeric vector of leave-one-out nonparametric conditional estimates evaluated at \code{X}
//' @useDynLib semiBRM, .registration = TRUE
//' @importFrom Rcpp evalCpp
// [[Rcpp::export]]
Rcpp::NumericVector GaussianKerNonparLeaveOneOut
(SEXP Y, SEXP X, SEXP H, int n_cores = 3)
{
	if ( Rf_isMatrix(X) == TRUE ){
		arma::vec Y_ = Rcpp::as< arma::vec >(Y);
		arma::mat X_ = Rcpp::as< arma::mat >(X);
		arma::rowvec H_ =  Rcpp::as< arma::rowvec >(H);
		return GaussianMultivarLeaveOneOut(Y_, X_, H_, n_cores);
	}
	else{
		arma::vec Y_ = Rcpp::as< arma::vec >(Y);
		arma::vec X_ = Rcpp::as< arma::vec >(X);
		double H_ =  Rcpp::as< double >(H);
		return GaussianUnivarLeaveOneOut(Y_, X_, H_, n_cores);
	}
}














