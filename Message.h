/*!
 *	\file Message.h
 *	\details Implementation of Message class
 *	\author Parthan Kasarapu
 *	\date Modified: Thu 5 Jul 2012
 */

#ifndef MESSAGE_H
#define MESSAGE_H

#include "Data.h"
#include "RandomDataGenerator.h"
#include "Error.h"
#include "Gaussian.h"
#include <liblcb/Vector.h>
#include <liblcb/Matrix.h>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

//! bound on the number of data samples
#define MAX_SAMPLES 100000

//! bound on the maximum number of orthogonal basis functions
#define MAX_FUNCTIONS 200

using namespace boost::numeric::ublas ;
using namespace std ;

/*!
 *	\class Message
 *	\brief This is the Message class abstraction.
 *
 *	The class acts as an interface to compute the minimum message length.
 */
class Message
{
	private:
		struct Parameters parameters ;
		lcb::Vector<long double> weights ;
		Data<long double> xVals, yVals, predictions ;
	public:
    //! null constructor
    Message() ;
		//! constructor
		template <class T>
		Message (struct Parameters, lcb::Matrix<T> &, Data<T> &, Data<T> &, 
							Data<T> &) ;
		//! instantiates a normal distribution
		template <class T>
		Gaussian normalDistribution (lcb::Vector<T> &) ;
		//! encodes the number of data samples and number of basis functions
		long double encodeIntegers () ;
		//! computes the message length using Wallace-Freeman formulation
		long double encodeUsingFreeman (int, long double, long double, long double, long double) ;
		//! measures the information content of data
		template <class T>
		long double encodeX (Data<T>, struct Parameters) ;
		//! encodes the coefficients of the basis functions
		long double encodeWeights() ;
		//! encodes the y values
		long double encodeOutput() ;
		//! computes the message length
		long double messageLength() ;
} ;

matrix<long double> make_matrix(lcb::Matrix<long double> &M, int dimension)
{
	matrix<long double> m(dimension,dimension) ;
	for (unsigned i=0; i<m.size1(); i++)
	{
		matrix_row<matrix<long double> > mr(m,i) ;
		for(unsigned j=0; j<mr.size(); j++)
			mr(j) = M[i][j] ;
	}
	return m ;
}

/* Matrix inversion routine.
 Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T>
bool InvertMatrix(const matrix<T>& input, matrix<T>& inverse)
{
	typedef permutation_matrix<std::size_t> pmatrix;

	// create a working copy of the input
	matrix<T> A(input);

	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());

	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;

	// create identity matrix of "inverse"
	inverse.assign(identity_matrix<T> (A.size1()));

	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);

	return true;
}

lcb::Matrix<long double> make_my_matrix(matrix<long double> &m, int dimension)
{
	lcb::Matrix<long double> result(dimension,dimension) ;
	for (unsigned i=0; i<m.size1(); i++)
	{
		matrix_row<matrix<long double> > mr(m,i) ;
		for (unsigned j=0; j<mr.size(); j++)
			result[i][j] = mr(j) ;
	}
	return result ;
}

/*!
 *	\fn lcb::Matrix<T> computeWeights (lcb::Matrix<T> &phi, Data<T> &yValues, int invChoice)
 *	\brief Computes the coefficients of the basis functions
 *	\param phi a reference to a Matrix object
 *	\param yValues a reference to a Data object
 *	\return a matrix object of weights
 */
template <class T>
lcb::Matrix<T> computeWeights (lcb::Matrix<T> &phi, Data<T> &yValues, int invChoice)
{
	lcb::Matrix<T> phiT = phi.transpose() ;
	lcb::Matrix<T> phiTphi = phiT * phi ;
	unsigned dim = phiTphi.rows() ; 
  lcb::Matrix<long double> pseudoInv ;
  matrix<long double> A,Z ;

  cout << "det = " << phiTphi.determinant() << endl ;
  switch(invChoice)
  {
    case 0:
      /* my implementation of matrix inverse */
	    pseudoInv = phiTphi.inverse() ;
      break ;
    case 1:
      /* boost::ublas implementation of matrix inverse */
  	  A = matrix<long double> (dim,dim) ;
  	  Z = matrix<long double> (dim,dim) ;
	    A = make_matrix(phiTphi,dim) ;	
	    InvertMatrix(A,Z) ;
	    pseudoInv = make_my_matrix(Z,dim) ;
      break ;
    default:
      error("Invalid choice of matrix inverse.") ;
      break ;
  }

	lcb::Matrix<T> temp = pseudoInv * phiT ;
	lcb::Matrix<T> y = yValues.convertToMatrix() ;
	lcb::Matrix<T> weights = temp * y ;

	ofstream phiFile ;
	phiFile.open("phi") ;
	for (int i=0; i<phi.rows(); i++)
	{
		for (int j=0; j<phi.columns(); j++)
			phiFile << phi[i][j] << " " ;
		phiFile << endl ;
	}
	phiFile.close() ;

	ofstream invFile ;
	invFile.open("inverse") ;
	for (int i=0; i<pseudoInv.rows(); i++)
	{
		for (int j=0; j<pseudoInv.columns(); j++)
			invFile << pseudoInv[i][j] << " " ;
		invFile << endl ;
	}
	invFile.close() ;

	return weights ;
}

/*!
 *	\fn long double computeRMSE (lcb::Matrix<T> &weights, lcb::Matrix<T> &phi, 
 *	Data<T> &yVals)
 *	Computes the root mean squared error in approximating the data with a
 *	linear combination of given number of basis functions.
 *	\param weights a reference to a Matrix object
 *	\param phi a reference to a Matrix object
 *	\param yVals a reference to a Data object
 *	\return the root mean squared error
 */
template <class T>
long double computeRMSE (lcb::Matrix<T> &weights, lcb::Matrix<T> &phi, 
										Data<T> &yVals)
{
  lcb::Matrix<T> yEst = phi * weights ; // column matrix
	long double diff, error = 0 ;
	int numSamples = phi.rows() ;
	for (int i=0; i<numSamples; i++)
	{
		diff = yEst[i][0] - yVals[i].x() ;
		error += diff * diff ;
	}
	return sqrt(error/numSamples) ;
}

/*!
 *  Null constructor for the Message class
 */
//template <class T>
Message :: Message ()
{
  /*xVals = Data<T>() ;
  yVals = Data<T>() ;
  predictions = Data<T>() ;
  weights = lcb::Matrix<T>() ;*/
}

/*!
 *	\fn Message :: Message (struct Parameters params, Matrix<T> &w, 
 *	Data<T> &xVals, Data<T> &yVals)
 *	\brief This is a constructor function to instantiate a Message object
 *	\param params a struct of Parameters
 *	\param w a reference to a Matrix object of type T
 *	\param xVals a reference to a Data object of type T
 *	\param yVals a reference to a Data object of type T
 */
template <class T>
Message :: Message (struct Parameters params, lcb::Matrix<T> &w, 
										Data<T> &xVals, Data<T> &yVals, Data<T> &yEst) : 
										xVals (xVals), yVals (yVals), parameters(params), 
										predictions(yEst)
{
	weights = w.getColumn(0) ;
}

/*!
 *	\fn long double Message :: encodeIntegers(void)
 *	\brief The function is used to encode the number of data samples
 *	transmitted and the number of orthogonal basis functions used.
 *	\return The message length to transmit the number of data points and
 *	the number of orthogonal basis functions.
 */
long double Message :: encodeIntegers(void)
{
	long double R = (long double) MAX_FUNCTIONS * (long double) MAX_SAMPLES ;
	long double logR = log2(R) ;
	return ceil(logR) ;
}

/*!
 *	\fn Gaussian Message :: normalDistribution (vector<T> &samples)
 *	\brief The function constructs a normal distribution whose parameters 
 *	are given by the mean and standard deviation of the elements in the list.
 *	\param samples a reference to a std::vector of type T
 *	\return An instantiated Gaussian class object
 */
template <class T>
Gaussian Message :: normalDistribution (lcb::Vector<T> &samples)
{
	long double mean = 0, sigmaSq = 0 ;
	int i,numSamples = samples.length() ;
	for (i=0; i<numSamples; i++) {
		mean += samples[i] ;
	}
	mean = mean / numSamples ;
	//cout << "Mean(w): " <<
	for (i=0; i<numSamples; i++)
		sigmaSq += (samples[i]-mean) * (samples[i]-mean) ;
	sigmaSq = sigmaSq/numSamples ;

	long double sigma = sqrt(sigmaSq) ;
	return Gaussian(mean,sigma) ;
}

/*!
 *	\fn long double encodeUsingFreeman (int N, long double sigma, long double rangeMu, 
 *	long double rangeLogSigma, long double Kn)
 *	\brief Computes the message length used to encode data and parameters
 *	using the Wallace-Freeman approach
 *	\param N an integer
 *	\param sigma a long double
 *	\param rangeMu	a long double
 *	\param rangeLogSigma a long double
 *	\param Kn a long double
 *	\return the length of the encoding based on Wallace Freeman approach
 */
long double Message :: encodeUsingFreeman (int N, long double sigma, long double rangeMu, 
																		long double rangeLogSigma, long double Kn)
{
	long double pi = boost::math::constants::pi<long double>() ;
	long double msgLen = 0.5 * (N-1) * log2l ((N * sigma * sigma)/(N-1)) + 0.5 * (N-1) +
									0.5 * N * log2l (2 * pi / (AOM * AOM)) +
									0.5 * log2l (2 * N * N) + log2l (rangeLogSigma) +
									1 + log2l (Kn) ;
	if (rangeMu == 0) {
		return msgLen ;
  }
	else return msgLen + log2l (rangeMu) ;
}

/*!
 *	\fn long double Message :: encodeX (Data<T> data, struct Parameters parameters)
 *	\brief This function is used to encode the random X values generated.
 *	These values are first sorted to arrange in increasing order. The minimum
 *	(the first) value after this arrangement is then subtracted from all the
 *	data values to make the series begin from 0 (shifting). The delta_X are 
 *	computed between any two successive sorted+shifted X values. These delta_Xs
 *	are then sent over the channel assuming these are sampled from a normal
 *	distribution.
 *	\param data a Data object of type T
 *	\param parameters a struct of Parameters
 *	\return message length of Xs
 */
template <class T>
long double Message :: encodeX (Data<T> data, struct Parameters parameters)
{
	int numSamples = data.nPoints() ;
	data.sortElements() ;
	Data<T> sortedX(data.sortedList()) ;
	for (int i=1; i<numSamples; i++) {
		sortedX[i] = sortedX[i] - sortedX[0] ;
  }
	sortedX[0] = 0 ;
	lcb::Vector<T> diff (numSamples-1) ;
	for (int i=0; i<numSamples-1; i++) {
		diff[i] = sortedX[i+1].x() - sortedX[i].x() ;
  }
	
	Gaussian normal = normalDistribution<T>(diff) ;
	//cout << "mu(dx) = " << normal.mean() << endl ;
	//cout << "sigma(dx) = " << normal.standardDeviation() << endl ;
	long N = numSamples - 1 ; 
	long double sigma = normal.standardDeviation() ;
	if (sigma <= 3 * AOM) {
		sigma = 3 * AOM ;
	}
	long double rangeMu = 2 * (parameters.high - parameters.low) ;
	long double logSigmaLowerBound = log2l(3*AOM) ;
	long double logSigmaUpperBound = log2l(rangeMu) ;
	long double rangeLogSigma = logSigmaUpperBound - logSigmaLowerBound ;
	long double K2 = 5 / (36 * sqrt(3)) ;
	
	long double msgLen = encodeUsingFreeman(N,sigma,rangeMu,rangeLogSigma,K2) ;
	//cout << "bits/sample: " << msgLen/N << endl ;
	return msgLen ;
}

/*	uses Wallace-Boulton formulation
long double Message :: encodeWeights (void)
{
	Gaussian normal = normalDistribution<long double>(weights) ;
	//cout << "mu(w): " << normal.mean() << endl ;
	//cout << "sigma(w): " << normal.standardDeviation() << endl ;

	long double rangeMu = 2 ; // mu \in [-1,1]
	long double codeLengthMu = log2l (rangeMu/AOM) ; 
	
	long double rangeSigma = 1 ;		// sigma \in [0,1]
	long double codeLengthSigma = log2l (rangeSigma/AOM) ;

	int N = parameters.numFunctions ;
	long double sigma = normal.standardDeviation() ;
	long double pi = boost::math::constants::pi<long double>() ;
	long double codeLengthWeights = N * log2l (sigma * sqrt(2*pi) / AOM) +
														 N * 0.5 / log(2) ;

										
	long double wt = codeLengthMu + codeLengthSigma + codeLengthWeights ;	
	//cout << "bits/wt: " << wt/N << endl ;
	return wt ;
}
*/

/*!
 *	\fn long double Message :: encodeWeights (void)
 *	\brief Computes the message length of encoding the coefficients of the 
 *	basis functions used to approximate the data using Wallace-Freeman
 *	formulation.
 *	\return length of encoding the weights
 */
long double Message :: encodeWeights (void)
{
	//weights.print() ;
	Gaussian normal = normalDistribution<long double>(weights) ;
	//cout << "mu(w): " << normal.mean() << endl ;
	//cout << "\nsigma(w): " << normal.standardDeviation() << endl ;

  long double mu = normal.mean() ;
  long double sigma = normal.standardDeviation() ;
	if (sigma <= 3 * AOM) {
		sigma = 3 * AOM ;
	}
  size_t N = parameters.numFunctions ;
	long double rangeMu = 2 ; // mu \in [-1,1]
  long double sigma_max = 1 ;
  long double sigma_min = AOM * 3 ;
	if (sigma_min > sigma_max) {
		throw std::domain_error("minimum sigma value exceeds set maximum limit") ;
	}
	long double rangeLogSigma = log2l(sigma_max)-log2l(sigma_min) ;	
	long double K2 = 5 / (36 * sqrt(3)) ;
  
	long double msgLen = encodeUsingFreeman(N,sigma,rangeMu,rangeLogSigma,K2) ;
	//cout << "bits/wt: " << wt/N << endl ;
	return msgLen ;
}

/*!
 *	\fn long double Message :: encodeOutput (void)
 *	\brief Computes the message length of encoding the difference in the output
 *	y values and the approximated versions using Wallace-Freeman formulation.
 *	\return length of encoding the output
 */
long double Message :: encodeOutput (void)
{
	int N = parameters.numSamples ;
	lcb::Vector<long double> diff(N) ;
	for (int i=0; i<N; i++) {
		diff[i] = yVals[i].x() - predictions[i].x() ;
  }
	Gaussian normal = normalDistribution<long double>(diff) ;
	//cout << "mu(dy) = " << normal.mean() << endl ;
	//cout << "\nsigma(dy) = " << normal.standardDeviation() << endl ;
	
	long double rangeMu = 2 ; // mu \in [-1,1]
	long double sigma = normal.standardDeviation() ;
	if (sigma <= 3 * AOM) {
		sigma = 3 * AOM ;
	}
	//cout << "\nsigma(dy) = " << sigma << endl ;
  long double sigma_max = 2;
  long double sigma_min = AOM * 3 ;
	if (sigma_min > sigma_max) {
		throw std::domain_error("minimum sigma value exceeds set maximum limit") ;
	}
	long double rangeLogSigma = log2l(sigma_max)-log2l(sigma_min) ;	
	long double K1 = 1.0 / 12 ;

	// rangeMu = 0 because not sending mu - centre of distribution is f(x)
	long double msgLen = encodeUsingFreeman(N,sigma,0,rangeLogSigma,K1) ;

	//cout << "bits/dy: " << msgDy/N << endl ;
	return msgLen ;
}

/*!
 *	\fn void Message :: messageLength (void)
 *	\brief This function computes the message length (in bits)
 *	MessageLength = length(parameters) + length(data|parameters)
 */
long double Message :: messageLength (void)
{
	/*	encode numFunctions and numSamples	*/
	//cout << "encoding number of functions and number of samples ..." << endl ;
	long double part1 = encodeIntegers() ;
	//cout << "encoding #functions + #samples: " << part1 << endl ;

	/*	encode weights	*/
	//cout << "encoding weights ..." << endl ;
	long double part2 = encodeWeights() ;
	//cout << "encoding weights         : " << part3 << endl ;

	/*	encode x's	*/
	//cout << "encoding X values ..." << endl ; 
	long double part3 = encodeX(xVals,parameters) ;
	//cout << "encoding X: " << part2 << endl ;

	/*	encode delta_y values	*/
	//cout << "encoding difference in output ..." << endl ;
	long double part4 = encodeOutput() ;
	//cout << "encoding Y: " << part4 << endl ;

	cout << "Int: " << part1 << "\tW: " << part2 << "\tX: " << part3 << "\tY: " << part4 << endl ;
	cout << "Message_1: " << part1 + part2 << endl ;
	cout << "Message_2: " << part3 + part4 << endl ;
  cout << "Total msgLen: " << part1 + part2 + part3 + part4 << endl ;
	return part1 + part2 + part3 + part4 ;
}

#endif

