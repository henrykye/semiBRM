#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
  #include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]


// Utility Functions -------------------------------------------------------------------------------
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

// Gaussian-Kernel Multivariate Nonparametric Conditional Expectation-------------------------------

// Multivariate at one point
inline double GaussianMultivarAtOnePoint
(const arma::vec& Y, const arma::mat& X, const arma::rowvec& arg, const arma::rowvec& H)
{
	arma::mat D = X.each_row() - arg;
	arma::mat Z = D.each_row() / H;
	arma::mat Kerh = arma::normpdf( Z ).each_row() / H;
	arma::vec w = arma::prod(Kerh, 1);
	return arma::sum( Y%w ) / arma::sum( w );
}

// Multivariate nonparametric estimator
// [[Rcpp::export(.GaussianMultivar)]]
Rcpp::NumericVector GaussianMultivar
(const arma::vec& Y, const arma::mat& X, const arma::mat& args, arma::rowvec H)
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

	// loop over args
	arma::uword Nargs = args.n_rows, i;
	std::vector< double > ce(Nargs); // output bin
	for ( i = 0; i < Nargs; i++ ){
		ce[i] = GaussianMultivarAtOnePoint(Y, X, args.row(i), H);
	}

	return Rcpp::wrap(ce);
}

// Multivariate nonparametric estimator with omp
// [[Rcpp::export(.GaussianMultivarOMP)]]
Rcpp::NumericVector GaussianMultivarOMP
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

	arma::uword Nargs = args.n_rows, i;
	std::vector< double > ce(Nargs); // output bin

    #ifdef _OPENMP
	    omp_set_num_threads(n_cores);
        #pragma omp parallel for default(none) shared(Nargs, Y, X, args, H, ce) private(i)
	    for ( i = 0; i < Nargs; i++ ){
	        ce[i] = GaussianMultivarAtOnePoint(Y, X, args.row(i), H);
	    }
    #else
	    for ( i = 0; i < Nargs; i++ ){
	        ce[i] = GaussianMultivarAtOnePoint(Y, X, args.row(i), H);
	    }
    #endif

	return Rcpp::wrap(ce);
}

// Multivariate nonparametric estimator, leave-one-out version
// [[Rcpp::export(.GaussianMultivarLeaveOneOut)]]
Rcpp::NumericVector GaussianMultivarLeaveOneOut
(const arma::vec& Y, const arma::mat& X, arma::rowvec H)
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

	// loop over data points
	arma::uword Nargs = X.n_rows, i;
	std::vector< double > ce(Nargs); // output bin
	for ( i  = 0; i < Nargs; i++ ){
		arma::uvec j = getAntiIndices(Nargs, i);
		ce[i] = GaussianMultivarAtOnePoint(Y.rows(j), X.rows(j), X.row(i), H);
	}

	return Rcpp::wrap(ce);
}

// Multivariate nonparametric estimator with omp, leave-one-out version
// [[Rcpp::export(.GaussianMultivarLeaveOneOutOMP)]]
Rcpp::NumericVector GaussianMultivarLeaveOneOutOMP
(const arma::vec& Y, const arma::mat& X, arma::rowvec H, int n_cores = 3)
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

	arma::uword Nargs = X.n_rows, i;
	std::vector< double > ce(Nargs); // output bin

    #ifdef _OPENMP
	    omp_set_num_threads(n_cores);
        #pragma omp parallel for default(none) shared(Nargs, Y, X, H, ce) private(i)
        for ( i  = 0; i < Nargs; i++ ){
            arma::uvec j = getAntiIndices(Nargs, i);
            ce[i] = GaussianMultivarAtOnePoint(Y.rows(j), X.rows(j), X.row(i), H);
        }
    #else
        for ( i  = 0; i < Nargs; i++ ){
            arma::uvec j = getAntiIndices(Nargs, i);
            ce[i] = GaussianMultivarAtOnePoint(Y.rows(j), X.rows(j), X.row(i), H);
        }
    #endif

	return Rcpp::wrap(ce);
}

// Gaussian-Kernel Univariate Nonparametric Conditional Expectation---------------------------------

// Univariate at one point
inline double GaussianUnivarAtOnePoint
(const arma::vec& Y, const arma::vec& X, const double& arg, const double& H)
{
	arma::vec Z = ( X - arg ) / H;
	arma::vec w = arma::normpdf( Z ) / H;
	return arma::sum( Y%w ) / arma::sum( w );
}

// Univariate nonparametric estimator
// [[Rcpp::export(.GaussianUnivar)]]
Rcpp::NumericVector GaussianUnivar
(const arma::vec& Y, const arma::vec& X, const arma::vec& args, const double& H)
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

	// loop over args
	arma::uword Nargs = args.n_rows, i;
	std::vector< double > ce(Nargs); // output bin
	for ( i  = 0; i < Nargs; i++ ){
		ce[i] = GaussianUnivarAtOnePoint(Y, X, args[i], H);
	}

	return Rcpp::wrap(ce);
}

// Univariate nonparametric estimator with omp
// [[Rcpp::export(.GaussianUnivarOMP)]]
Rcpp::NumericVector GaussianUnivarOMP
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

	arma::uword Nargs = args.n_rows, i;
	std::vector< double > ce(Nargs); // output bin

    #ifdef _OPENMP
	    omp_set_num_threads(n_cores);
        #pragma omp parallel for default(none) shared(Nargs, Y, X, args, H, ce) private(i)
	    for ( i  = 0; i < Nargs; i++ ){
	        ce[i] = GaussianUnivarAtOnePoint(Y, X, args[i], H);
	    }
    #else
    	for ( i  = 0; i < Nargs; i++ ){
    	    ce[i] = GaussianUnivarAtOnePoint(Y, X, args[i], H);
    	}
    #endif

	return Rcpp::wrap(ce);
}

// Univariate nonparametric estimator, leave-one-out version
// [[Rcpp::export(.GaussianUnivarLeaveOneOut)]]
Rcpp::NumericVector GaussianUnivarLeaveOneOut
(const arma::vec& Y, const arma::vec& X, const double& H)
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

	// loop over data point
	arma::uword Nargs = X.n_rows, i;
	std::vector< double > ce(Nargs); // output bin
	for ( i = 0; i < Nargs; i++ ){
		arma::uvec j = getAntiIndices(Nargs, i);
		ce[i] = GaussianUnivarAtOnePoint(Y.rows(j), X.rows(j), X[i], H);
	}

	return Rcpp::wrap(ce);
}

// Univariate nonparametric estimator with omp, leave-one-out version
// [[Rcpp::export(.GaussianUnivarLeaveOneOutOMP)]]
Rcpp::NumericVector GaussianUnivarLeaveOneOutOMP
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

	arma::uword Nargs = X.n_rows, i;
	std::vector< double > ce(Nargs); // output bin

    #ifdef _OPENMP
	    omp_set_num_threads(n_cores);
        #pragma omp parallel for default(none) shared(Nargs, Y, X, H, ce) private(i)
	    for ( i = 0; i < Nargs; i++ ){
	        arma::uvec j = getAntiIndices(Nargs, i);
	        ce[i] = GaussianUnivarAtOnePoint(Y.rows(j), X.rows(j), X[i], H);
	    }
    #else
	    for ( i = 0; i < Nargs; i++ ){
	        arma::uvec j = getAntiIndices(Nargs, i);
	        ce[i] = GaussianUnivarAtOnePoint(Y.rows(j), X.rows(j), X[i], H);
	    }
    #endif

	return Rcpp::wrap(ce);
}
