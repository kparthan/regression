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

//! bound on the number of data samples
#define MAX_SAMPLES 100000

//! bound on the maximum number of orthogonal basis functions
#define MAX_FUNCTIONS 200

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
		lcb::Vector<double> weights ;
		Data<double> xVals, yVals, predictions ;
	public:
		//! constructor
		template <class T>
		Message (struct Parameters, lcb::Matrix<T> &, Data<T> &, Data<T> &, 
							Data<T> &) ;
		//! instantiates a normal distribution
		template <class T>
		Gaussian normalDistribution (lcb::Vector<T> &) ;
		//! encodes the number of data samples and number of basis functions
		double encodeIntegers () ;
		//! computes the message length using Wallace-Freeman formulation
		double encodeUsingFreeman (int, double, double, double, double) ;
		//! measures the information content of data
		template <class T>
		double encodeX (Data<T>, struct Parameters) ;
		//! encodes the coefficients of the basis functions
		double encodeWeights() ;
		//! encodes the y values
		double encodeOutput() ;
		//! computes the message length
		double messageLength() ;
} ;

/*!
 *	\fn lcb::Matrix<T> computeWeights (lcb::Matrix<T> &phi, Data<T> &yValues)
 *	\brief Computes the coefficients of the basis functions
 *	\param phi a reference to a Matrix object
 *	\param yValues a reference to a Data object
 *	\return a matrix object of weights
 */
template <class T>
lcb::Matrix<T> computeWeights (lcb::Matrix<T> &phi, Data<T> &yValues)
{
	lcb::Matrix<T> phiT = phi.transpose() ;
	lcb::Matrix<T> phiTphi = phiT * phi ;
	lcb::Matrix<T> pseudoInv = phiTphi.inverse() ;
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
 *	\fn double computeRMSE (lcb::Matrix<T> &weights, lcb::Matrix<T> &phi, 
 *	Data<T> &yVals)
 *	Computes the root mean squared error in approximating the data with a
 *	linear combination of given number of basis functions.
 *	\param weights a reference to a Matrix object
 *	\param phi a reference to a Matrix object
 *	\param yVals a reference to a Data object
 *	\return the root mean squared error
 */
template <class T>
double computeRMSE (lcb::Matrix<T> &weights, lcb::Matrix<T> &phi, 
										Data<T> &yVals)
{
  lcb::Matrix<T> yEst = phi * weights ; // column matrix
	double diff, error = 0 ;
	int numSamples = phi.rows() ;
	for (int i=0; i<numSamples; i++)
	{
		diff = yEst[i][0] - yVals[i].x() ;
		error += diff * diff ;
	}
	return sqrt(error/numSamples) ;
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
 *	\fn double Message :: encodeIntegers(void)
 *	\brief The function is used to encode the number of data samples
 *	transmitted and the number of orthogonal basis functions used.
 *	\return The message length to transmit the number of data points and
 *	the number of orthogonal basis functions.
 */
double Message :: encodeIntegers(void)
{
	double R = (double) MAX_FUNCTIONS * (double) MAX_SAMPLES ;
	double logR = log2(R) ;
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
	double mean = 0, sigmaSq = 0 ;
	int i,numSamples = samples.length() ;
	for (i=0; i<numSamples; i++) {
		mean += samples[i] ;
	}
	mean = mean / numSamples ;
	//cout << "Mean(w): " <<
	for (i=0; i<numSamples; i++)
		sigmaSq += (samples[i]-mean) * (samples[i]-mean) ;
	sigmaSq = sigmaSq/numSamples ;

	double sigma = sqrt(sigmaSq) ;
	return Gaussian(mean,sigma) ;
}

/*!
 *	\fn double encodeUsingFreeman (int N, double sigma, double rangeMu, 
 *	double rangeLogSigma, double Kn)
 *	\brief Computes the message length used to encode data and parameters
 *	using the Wallace-Freeman approach
 *	\param N an integer
 *	\param sigma a double
 *	\param rangeMu	a double
 *	\param rangeLogSigma a double
 *	\param Kn a double
 *	\return the length of the encoding based on Wallace Freeman approach
 */
double Message :: encodeUsingFreeman (int N, double sigma, double rangeMu, 
																		double rangeLogSigma, double Kn)
{
	double pi = boost::math::constants::pi<double>() ;
	double msgLen = 0.5 * (N-1) * log2l ((N * sigma * sigma)/(N-1)) + 0.5 * (N-1) +
									0.5 * N * log2l (2 * pi / (AOM * AOM)) +
									0.5 * log2l (2 * N * N) + log2l (rangeLogSigma) +
									1 + log2l (Kn) ;
	if (rangeMu == 0)
		return msgLen ;
	else return msgLen + log2l (rangeMu) ;
}

/*!
 *	\fn double Message :: encodeX (Data<T> data, struct Parameters parameters)
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
double Message :: encodeX (Data<T> data, struct Parameters parameters)
{
	int numSamples = data.nPoints() ;
	data.sortElements() ;
	Data<T> sortedX(data.sortedList()) ;
	for (int i=1; i<numSamples; i++)
		sortedX[i] = sortedX[i] - sortedX[0] ;
	sortedX[0] = 0 ;
	lcb::Vector<T> diff (numSamples-1) ;
	for (int i=0; i<numSamples-1; i++)
		diff[i] = sortedX[i+1].x() - sortedX[i].x() ;
	
	Gaussian normal = normalDistribution<T>(diff) ;
	//cout << "mu(dx) = " << normal.mean() << endl ;
	//cout << "sigma(dx) = " << normal.standardDeviation() << endl ;
	long N = numSamples - 1 ; 
	double sigma = normal.standardDeviation() ;
	if (sigma <= 3 * AOM)
	{
		sigma = 3 * AOM ;
	}
	double rangeMu = 2 * (parameters.high - parameters.low) ;
	double logSigmaLowerBound = log2l(3*AOM) ;
	double logSigmaUpperBound = log2l(rangeMu) ;
	double rangeLogSigma = logSigmaUpperBound - logSigmaLowerBound ;
	double K2 = 5 / (36 * sqrt(3)) ;
	
	double msgLen = encodeUsingFreeman(N,sigma,rangeMu,rangeLogSigma,K2) ;
	//cout << "bits/sample: " << msgLen/N << endl ;
	return msgLen ;
}

/*	uses Wallace-Boulton formulation
double Message :: encodeWeights (void)
{
	Gaussian normal = normalDistribution<double>(weights) ;
	//cout << "mu(w): " << normal.mean() << endl ;
	//cout << "sigma(w): " << normal.standardDeviation() << endl ;

	double rangeMu = 2 ; // mu \in [-1,1]
	double codeLengthMu = log2l (rangeMu/AOM) ; 
	
	double rangeSigma = 1 ;		// sigma \in [0,1]
	double codeLengthSigma = log2l (rangeSigma/AOM) ;

	int N = parameters.numFunctions ;
	double sigma = normal.standardDeviation() ;
	double pi = boost::math::constants::pi<double>() ;
	double codeLengthWeights = N * log2l (sigma * sqrt(2*pi) / AOM) +
														 N * 0.5 / log(2) ;

										
	double wt = codeLengthMu + codeLengthSigma + codeLengthWeights ;	
	//cout << "bits/wt: " << wt/N << endl ;
	return wt ;
}
*/

/*!
 *	\fn double Message :: encodeWeights (void)
 *	\brief Computes the message length of encoding the coefficients of the 
 *	basis functions used to approximate the data using Wallace-Freeman
 *	formulation.
 *	\return length of encoding the weights
 */
double Message :: encodeWeights (void)
{
	weights.print() ;
	Gaussian normal = normalDistribution<double>(weights) ;
	cout << "mu(w): " << normal.mean() << endl ;
	//cout << "\nsigma(w): " << normal.standardDeviation() << endl ;

  double mu = normal.mean() ;
  double sigma = normal.standardDeviation() ;
	if (sigma <= 3 * AOM) {
		sigma = 3 * AOM ;
	}
	cout << "sigma(w): " << sigma << endl ;
  size_t N = parameters.numFunctions ;
	double rangeMu = 2 ; // mu \in [-1,1]
  double sigma_max = 1 ;
  double sigma_min = AOM * 3 ;
	if (sigma_min > sigma_max) {
		throw std::domain_error("minimum sigma value exceeds set maximum limit") ;
	}
	double rangeLogSigma = log2l(sigma_max)-log2l(sigma_min) ;	
	double K2 = 5 / (36 * sqrt(3)) ;
  
	double msgLen = encodeUsingFreeman(N,sigma,rangeMu,rangeLogSigma,K2) ;
	//cout << "bits/wt: " << wt/N << endl ;
	return msgLen ;
}

/*!
 *	\fn double Message :: encodeOutput (void)
 *	\brief Computes the message length of encoding the difference in the output
 *	y values and the approximated versions using Wallace-Freeman formulation.
 *	\return length of encoding the output
 */
double Message :: encodeOutput (void)
{
	int N = parameters.numSamples ;
	lcb::Vector<double> diff(N) ;
	for (int i=0; i<N; i++)
		diff[i] = yVals[i].x() - predictions[i].x() ;
	Gaussian normal = normalDistribution<double>(diff) ;
	//cout << "mu(dy) = " << normal.mean() << endl ;
	//cout << "\nsigma(dy) = " << normal.standardDeviation() << endl ;
	
	double rangeMu = 2 ; // mu \in [-1,1]
	double sigma = normal.standardDeviation() ;
	if (sigma <= 3 * AOM) {
		sigma = 3 * AOM ;
	}
	cout << "\nsigma(dy) = " << sigma << endl ;
  double sigma_max = 2;
  double sigma_min = AOM * 3 ;
	if (sigma_min > sigma_max) {
		throw std::domain_error("minimum sigma value exceeds set maximum limit") ;
	}
	double rangeLogSigma = log2l(sigma_max)-log2l(sigma_min) ;	
	double K1 = 1.0 / 12 ;
	//double K2 = 5 / (36 * sqrt(3)) ;

	// rangeMu = 0 because not sending mu - centre of distribution is f(x)
	double msgLen = encodeUsingFreeman(N,sigma,0,rangeLogSigma,K1) ;

	//cout << "bits/dy: " << msgDy/N << endl ;
	return msgLen ;
}

/*!
 *	\fn void Message :: messageLength (void)
 *	\brief This function computes the message length (in bits)
 *	MessageLength = length(parameters) + length(data|parameters)
 */
double Message :: messageLength (void)
{
	/*	encode numFunctions and numSamples	*/
	//cout << "encoding number of functions and number of samples ..." << endl ;
	double part1 = encodeIntegers() ;
	//cout << "encoding #functions + #samples: " << part1 << endl ;

	/*	encode x's	*/
	//cout << "encoding X values ..." << endl ; 
	double part2 = encodeX(xVals,parameters) ;
	//cout << "encoding X: " << part2 << endl ;

	/*	encode weights	*/
	//cout << "encoding weights ..." << endl ;
	double part3 = encodeWeights() ;
	//cout << "encoding weights         : " << part3 << endl ;

	/*	encode delta_y values	*/
	//cout << "encoding difference in output ..." << endl ;
	double part4 = encodeOutput() ;
	//cout << "encoding Y: " << part4 << endl ;

	cout << "Int: " << part1 << "\tX: " << part2 << "\tW: " << part3 << "\tY: " << part4 << endl ;
	cout << "Message_1: " << part1 + part3 << endl ;
	cout << "Message_2: " << part2 + part4 << endl ;
	return part1 + part2 + part3 + part4 ;
}

#endif

