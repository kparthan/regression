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
#define MAX_FUNCTIONS 100

//!	setting the accuracy of measurement value for data samples
#define AOM 0.000001

//! number of iterations
#define NUM_ITERATIONS 100

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
		Message (struct Parameters, lcb::Matrix<T> &, Data<T> &, Data<T> &, Data<T> &) ;
		//! instantiates a normal distribution
		template <class T>
		Gaussian normalDistribution (lcb::Vector<T> &) ;
		//! encodes the number of data samples and number of basis functions
		double encodeIntegers () ;
		//! measures the information content of data
		template <class T>
		double encodeX (Data<T> , struct Parameters) ;
		//!
		double encodeWeights() ;
		//!
		double encodeOutput() ;
		//! computes the message length
 		double messageLength() ;
} ;

template <class T>
lcb::Matrix<T> computeWeights (lcb::Matrix<T> &phi, Data<T> &yValues)
{
	lcb::Matrix<T> phiT = phi.transpose() ;
	lcb::Matrix<T> phiTphi = phiT * phi ;
	lcb::Matrix<T> pseudoInv = phiTphi.inverse() ;
	lcb::Matrix<T> temp = pseudoInv * phiT ;

	lcb::Matrix<T> y = yValues.convertToMatrix() ;
	lcb::Matrix<T> weights = temp * y ;
	return weights ;
}

template <class T>
double computeRMSE (lcb::Matrix<T> &weights, lcb::Matrix<T> &phi, Data<T> &yVals)
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
Message :: Message (struct Parameters params, lcb::Matrix<T> &w, Data<T> &xVals, 
           						Data<T> &yVals, Data<T> &yEst) : xVals (xVals), yVals (yVals),
																				parameters(params), predictions(yEst)
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
 *	\brief The function constructs a normal distribution whose parameters are
 *	given by the mean and standard deviation of the elements in the list.
 *	\param samples a reference to a std::vector of type T
 *	\return An instantiated Gaussian class object
 */
template <class T>
Gaussian Message :: normalDistribution (lcb::Vector<T> &samples)
{
	double mean = 0, sigmaSq = 0 ;
	int i,numSamples = samples.length() ;
	for (i=0; i<numSamples; i++)
		mean += samples[i] ;
	mean = mean / numSamples ;
	for (i=0; i<numSamples; i++)
		sigmaSq += (samples[i]-mean) * (samples[i]-mean) ;
	sigmaSq = sigmaSq/numSamples ;

	double sigma = sqrt(sigmaSq) ;
	return Gaussian(mean,sigma) ;
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
	double pi = boost::math::constants::pi<double>() ;
	long N = numSamples - 1 ; 
	double sigma = normal.standardDeviation() ;
	double vSq = N * sigma * sigma ; 
	double rangeMu = 2 * (parameters.high - parameters.low) ;
	double logSigmaLowerBound = log2(AOM) ;
	double logSigmaUpperBound = log2(rangeMu) ;
	double rangeLogSigma = logSigmaUpperBound - logSigmaLowerBound ;
	double K2 = 5 / (36 * sqrt(3)) ;

	/*double p1 =0.5 * (N-1) * log2l(vSq/(N-1)) + 0.5 * (N-1) ;
	double p2 =0.5 * N * log2l (2 * pi/(AOM * AOM));
	double p3 = 	0.5 * log2l (2 * N * N) ; 
	double p4 = log2 (rangeMu * rangeLogSigma) ;
	double p5 =  	1 + log2 (K2) ;
	cout << "p1 = " << p1 << endl;
	cout << "p2 = " << p2 << endl;
	cout << "p3 = " << p3 << endl;
	cout << "p4 = " << p4 << endl;
	cout << "p5 = " << p5 << endl;*/

	double msgLen = 0.5 * (N-1) * log2l(vSq/(N-1)) + 0.5 * (N-1) +
									0.5 * N * log2 (2 * pi/(AOM * AOM)) +
									0.5 * log2 (2 * N * N) + log2 (rangeMu * rangeLogSigma) +
									1 + log2 (K2) ;
	return msgLen ;
}

double Message :: encodeWeights (void)
{
	//weights.print() ;
	Gaussian normal = normalDistribution<double>(weights) ;
	//cout << normal.mean() << " " << normal.standardDeviation() << endl ;

	double rangeMu = 2 ; // mu \in [-1,1]
	double codeLengthMu = log2l (rangeMu/AOM) ; 
	
	double rangeSigma = 1 ; 	// sigma \in [0,1]
	double codeLengthSigma = log2l (rangeSigma/AOM) ;

	int N = parameters.numFunctions ;
	double sigma = normal.standardDeviation() ;
	double pi = boost::math::constants::pi<double>() ;
	double codeLengthWeights = N * log2l (sigma * sqrt(2*pi) / AOM) +
														 N * 0.5 / log(2) ;
													
	return codeLengthMu + codeLengthSigma + codeLengthWeights ;	
}

double Message :: encodeOutput (void)
{
	int N = parameters.numSamples ;
	lcb::Vector<double> diff(N) ;
	for (int i=0; i<N; i++)
		diff[i] = yVals[i].x() - predictions[i].x() ;
	Gaussian normal = normalDistribution<double>(diff) ;
	//cout << "mu(dy) = " << normal.mean() << endl ;
	//cout << "sigma(dy) = " << normal.standardDeviation() << endl ;
	
	double rangeMu = 2 ; // mu \in [-1,1]
	double codeLengthMu = log2l (rangeMu/AOM) ; 
	
	double rangeSigma = 2 ; 	// sigma \in [0,2]
	double codeLengthSigma = log2l (rangeSigma/AOM) ;

	double sigma = normal.standardDeviation() ;
	double pi = boost::math::constants::pi<double>() ;
	double codeLengthDiff = N * log2l (sigma * sqrt(2*pi) / AOM) +
													N * 0.5 / log(2) ;
													
	return codeLengthMu + codeLengthSigma + codeLengthDiff ;	
}

/*!
 *	\fn void Message :: messageLength (void)
 *	\brief This function computes the message length (in bits)
 *	MessageLength = length(parameters) + length(data|parameters)
 */
double Message :: messageLength (void)
{
	// encode numFunctions and numSamples
	//cout << "encoding number of functions and number of samples ..." << endl ;
	double part1 = encodeIntegers() ;
	//cout << "part 1 = " << part1 << endl ;

	// encode x's
	//cout << "encoding X values ..." << endl ; 
	double part2 = encodeX(xVals,parameters) ;
	//cout << "part 2 = " << part2 << endl ;

	// encode weights
	//cout << "encoding weights ..." << endl ;
	double part3 = encodeWeights() ;
	//cout << "part 3 = " << part3 << endl ;

	// encode delta_y values
	//cout << "encoding difference in output ..." << endl ;
	double part4 = encodeOutput() ;
	//cout << "part 4 = " << part4 << endl ;

	return part1 + part2 + part3 + part4 ;
}

#endif

