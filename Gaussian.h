/*!
 *  \file Gaussian.h
 *  \details Implementation of Gaussian class
 *  \author Parthan Kasarapu
 *  \date Modified: Mon 25 Jun 2012
 */

#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
#include <boost/math/constants/constants.hpp>

using namespace std ;

/*!
 *  \class Gaussian
 *  \brief This is the Gaussian class abstraction.
 *
 *  The class acts as an interface to create and work with a 
 *  Gaussian distribution
 */
class Gaussian 
{
	private:
		long double mu ;
		long double sigma ;
	public:
		//! null constructor
		Gaussian() ;
		//! constructor that sets value of parameters
		Gaussian (long double, long double);
		//! value of function
		long double value(long double) ;
		//! generates a sample from the distribution
		vector<long double> generate () ;
		//! returns the mean 
		const long double mean () ;
		//! returns the standard deviation
		const long double standardDeviation () ;
} ;

/*!
 *  \fn Gaussian :: Gaussian ()
 *  \brief Null constructor
 *  sets default values of mean as 0 and standard deviation as 1
 */
Gaussian :: Gaussian () : mu(0), sigma(1)
{
}

/*!
 *  \fn Gaussian :: Gaussian (long double mu, long double sigma)
 *  \brief constructor function which sets the value of mean and 
 *  standard deviation of the distribution
 *  \param mu a long double
 *  \param sigma a long double
 */
Gaussian :: Gaussian (long double mu, long double sigma) : mu(mu), sigma(sigma)
{
}

/*!
 *  \fn long double Gaussian :: value (long double x)
 *  \brief computes the function value of the distribution
 *  \param x a long double
 *  \return value of the function given x
 */
long double Gaussian :: value (long double x)
{
	long double expNumerator = (-1) * (x-mu) * (x-mu) ;
	long double expDenominator = 2 * sigma * sigma ;
	long double pi = boost::math::constants::pi<long double>() ;
	return (exp (expNumerator/expDenominator)) / ((sqrt (2*pi)) * sigma) ;
}

/*!
 *	\fn vector<long double> Gaussian :: generate (void)
 *  \brief This function generates a data point sampled from this 
 *	Gaussian distribution. Uses Box-Muller method to draw samples from the 
 *	standard normal distribution i.e., N(0,1)
 *	X = sqrt(-2 ln U) cos(2*pi*V)
 *	Y = sqrt(-2 ln U) sin(2*pi*V), where
 *	U,V are drawn from a uniform distribution in (0,1). The resultant X,Y
 *	will be sampled from a standard normal distribution
 *	To generate from a general N(mu,sigma), use the transformation:
 *	Z = mu + sigma * X, where X~N(0,1)
 *	\return a sample drawn from the normal distribution
 */
vector<long double> Gaussian :: generate (void)
{
	//srand(time(0)) ;
	vector<long double> samples (2) ;
	long double u = (long double) rand() / RAND_MAX ;
	long double v = (long double) rand() / RAND_MAX ;
	long double sqroot = sqrt(-2 * log(u)) ;
	long double _2piV = 2 * boost::math::constants::pi<long double>() * v ;
	long double r1 = sqroot * cos (_2piV) ;
	long double r2 = sqroot * sin (_2piV) ;
  //cout << sigma << endl ;
	samples[0] = mu + sigma * r1 ;
	samples[1] = mu + sigma * r2 ;
  //cout << samples[0] << " " << samples[1] << endl ;
	return samples ;
}

/*!
 *  \fn const long double Gaussian :: mean (void)
 *  \brief returns the mean of the distribution
 *  \return the mean of the distribution
 */
const long double Gaussian :: mean (void)
{
	return mu ;
}

/*!
 *  \fn const long double Gaussian :: standardDeviation (void)
 *  \brief returns the standard deviation of the distribution
 *  \return the standard deviation of the distribution
 */
const long double Gaussian :: standardDeviation(void)
{
	return sigma ;
}

#endif

