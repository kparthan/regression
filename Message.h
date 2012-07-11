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
#include "Matrix.h"
#include "Error.h"
#include "Gaussian.h"
#include <vector>
#include <cmath>
#include <boost/math/constants/constants.hpp>

//! bound on the number of data samples
#define MAX_SAMPLES 100000

//! bound on the maximum number of orthogonal basis functions
#define MAX_FUNCTIONS 100

//!	setting the accuracy of measurement value for data samples
#define AOM 0.000001

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
		Vector<double> weights ;
		Data<double> xVals, yVals ;
	public:
		//! constructor
		template <class T>
		Message (struct Parameters, Matrix<T> &, Data<T> &, Data<T> &) ;
		//! instantiates a normal distribution
		template <class T>
		Gaussian normalDistribution (vector<T> &) ;
		//! encodes the number of data samples and number of basis functions
		double encodeIntegers () ;
		//! measures the information content of data
		template <class T>
		double encodeX (Data<T> , struct Parameters) ;
		//! computes the message length
 		void messageLength() ;
} ;

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
Message :: Message (struct Parameters params, Matrix<T> &w, Data<T> &xVals, 
           						Data<T> &yVals) : xVals (xVals), yVals (yVals),
																				parameters(params)
{
	weights = w.convertToVector() ;
}

/*!
 *	\relates Message
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
 *	\relates Message
 *	\brief The function constructs a normal distribution whose parameters are
 *	given by the mean and standard deviation of the elements in the list.
 *	\param samples a reference to a std::vector of type T
 *	\return An instantiated Gaussian class object
 */
template <class T>
Gaussian Message :: normalDistribution (vector<T> &samples)
{
	double mean = 0, sigmaSq = 0 ;
	int i,numSamples = samples.size() ;
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
 *	\relates Message
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
	vector<T> diff (numSamples-1) ;
	for (int i=0; i<numSamples-1; i++)
		diff[i] = sortedX[i+1].x() - sortedX[i].x() ;
	
	Gaussian normal = normalDistribution<T>(diff) ;
	cout << "mu(dx) = " << normal.mean() << endl ;
	cout << "sigma(dx) = " << normal.standardDeviation() << endl ;
	double pi = boost::math::constants::pi<double>() ;
	int N = numSamples - 1 ; 
	double sigma = normal.standardDeviation() ;
	double vSq = N * sigma * sigma ; 
	double rangeMu = 2 * (parameters.high - parameters.low) ;
	double logSigmaLowerBound = log2(AOM) ;
	double logSigmaUpperBound = log2(rangeMu) ;
	double rangeLogSigma = logSigmaUpperBound - logSigmaLowerBound ;
	double K2 = 5 / (36 * sqrt(3)) ;

	double msgLen = 0.5 * (N-1) * log2l(vSq/(N-1)) + 0.5 * (N-1) +
									0.5 * N * log2 (2 * pi/(AOM * AOM)) +
									0.5 * log2 (2 * N * N) + log2 (rangeMu * rangeLogSigma) +
									1 + log2 (K2) ;
	return msgLen ;
}

/*!
 *	\fn void Message :: messageLength (void)
 *	\brief This function computes the message length (in bits)
 *	MessageLength = length(parameters) + length(data|parameters)
 */
void Message :: messageLength (void)
{
	// encode numFunctions and numSamples
	double part1 = encodeIntegers() ;

	// encode x's
	double part2 = encodeX(xVals,parameters) ;
	cout << "part 1 = " << part1 << endl ;
	cout << "part 2 = " << part2 << endl ;
}

#endif

