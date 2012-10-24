/*/ 
 *  \file RandomDataGenerator.h
 *  \details Implementation of RandomDataGenerator class
 *  \author Parthan Kasarapu
 *  \date Modified: Mon 25 Jun 2012
 */

#ifndef RANDOM_DATA_GENERATOR_H
#define RANDOM_DATA_GENERATOR_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <liblcb/Vector.h>
#include "Data.h"
#include "Gaussian.h"
#include "Plot.h"
#include "OrthogonalBasis.h"

//!	accuracy of measurement
#define AOM 0.001

//!	minimum number of terms used in data generation using a
//!	finite linear combination model
#define MIN_TERMS 3

//!	maximum number of terms used in data generation using a
//!	finite linear combination model
#define MAX_TERMS 100

using namespace std ;

/*!
 *	\struct Parameters
 *	\brief This structure holds the parameters of the simulation
 *	\var Parameters::mean
 *	Member 'mean' contains the mean of Gaussian noise
 *	\var Parameters::sigma
 * 	Member 'sigma' contains the standard deviation of the distribution
 *	\var Parameters::low
 *	Member 'low' contains the lower bound of data
 *	\var Parameters::high
 *	Member 'high' contains the upper bound of data
 *	\var Parameters::function
 *	Member 'function' contains the function index 
 *	\var Parameters::timePeriod
 *	Member 'timePeriod' contains the time period of the wave
 *	\var Parameters::peak
 *	Member 'peak' contains the maximum amplitude of the wave
 *	\var Parameters::numSamples
 *	Member 'numSamples' contains the number of data samples
 *	\var Parameters::numFunctions
 *	Member 'numFunctions' contains the number of orthogonal functions
 *	\var Parameters::file
 *	Member 'file' contains the input data to run the experiment 
 *  \var Parameters::iterate
 *  Member 'iterate' determines whether to run as a standalone unit or
 *  iterate over a set of parameters
 *  \var Parameters::inverse
 *  Member 'inverse' determines the method to solve the linear system of
 *  equations to compute the weights in the regression fit
 *  \var Parameters::basis
 *  Member 'basis' contains the choice of orthogonal basis set
 *  \var Parameters:m
 *  Member 'm' is the parameter which specifies the position at which the
 *  the peak of a triangle wave is seen
 *  \var Parameters:lambda
 *  Member 'lambda' is the parameter used in the regularized least squares
 *  version
 */
struct Parameters 
{
	long double mean ;
  long double sigma ;
  long double low ;
  long double high ; 
  int function ;
  long double timePeriod ;
  long double peak ;
  int numSamples ;
  int numFunctions ;
	string file ;
  int iterate ;
  int inverse ;
  int basis ;
  long double m ;
  long double lambda ;
} ;

/*!
 *  \class RandomDataGenerator
 *  \brief This is the RandomDataGenerator class abstraction.
 *
 *  The class acts as an interface to generate random data samples.
 */
template <class T>
class RandomDataGenerator
{
	private:
		struct Parameters parameters ;
		Data<T> xVal ;
		Data<T> fxVal ;
		Data<T> yVal ;
		string fname ;
    string finlincomb ; // used only for finlincomb for verification
    string inferred ;
	public:
    //! null constructor
    RandomDataGenerator() ;
		//! constructor
		RandomDataGenerator(struct Parameters) ;
		//! generates data randomly
		void generate() ;		
		//! computes the function value for corresponding x values
		void computeFunctionValues() ;
		//! add noise to the generated data
		void addNoise() ;
		//! predicts the value based on the regression fit
		Data<T> predict (unsigned, lcb::Matrix<T> &, Data<T> &) ;
		//! returns the random X values generated
		Data<T> randomX () ;
		//! returns the corresponding function values
		Data<T> fxValues () ;
		//! returns the measured y values
		Data<T> yValues () ;
		//! plotting x 
		void plotRandomX() ;	
		//! plotting x vx f(x)
		void plotData() ;	
		//! plotting x vs f(x)+e
		void plotDataWithNoise() ;
		//! plots the predicted output values
		void plotPredictions(Data<T> &, Data<T> &, Data<T> &) ;	
		//! prints the private members of the class
		void print() ;
    //! 
    void constructFunctionString(string &, lcb::Vector<T> &, int) ;
    string getFunctionString() ;
} ;

template <class T>
string RandomDataGenerator<T> :: getFunctionString(void)
{
  return finlincomb ;
}

template <class T>
RandomDataGenerator<T> :: RandomDataGenerator()
{
	xVal = Data<T>() ;
	fxVal = Data<T>() ;
	yVal = Data<T>() ;
}

/*!
 *  \fn RandomDataGenerator<T> :: RandomDataGenerator ()
 *  \brief The null constructor sets the default interval from which
 *  the values of x are sampled
 */
template <class T>
RandomDataGenerator<T> :: RandomDataGenerator (struct Parameters parameters) : parameters (parameters)
{
	xVal = Data<T>() ;
	fxVal = Data<T>() ;
	yVal = Data<T>() ;
}

/*!
 *  \fn RandomDataGenerator<T> ::  generate ()
 *  \brief This function is used to generate data samples randomly
 */
template <class T>
void RandomDataGenerator<T> :: generate ()
{
	//srand ((unsigned)time(0)) ;
	srand (100000) ;
	if (!parameters.file.empty())
	{
		ifstream dataFile (parameters.file.c_str()) ;
		// updating number of samples
		dataFile >> parameters.numSamples ;
		if (parameters.numSamples <= 0)
			error ("In input file: # of data samples should be non-negative!") ;
		T *randomValues = new T [parameters.numSamples] ;
		for (int i=0; i<parameters.numSamples; i++)
			dataFile >> randomValues[i] ;
		dataFile.close() ;
		xVal = Data<T> (randomValues,parameters.numSamples) ;
		delete[] randomValues ;
		// updating lower and upper bounds
		parameters.low = xVal.minimum() ;
		parameters.high = xVal.maximum() ;
	}
	else
	{
		T range = parameters.high - parameters.low ;
		T random ;
		T *randomValues = new T [parameters.numSamples] ;
		for (int i=0; i<parameters.numSamples; i++)
		{
			random = rand() / (T) RAND_MAX ;
			randomValues[i] = parameters.low + random * range ;
		}
		/*long double partition = range / parameters.numSamples ;
		randomValues[0] = parameters.low ;
		for (int i=1; i<parameters.numSamples; i++)
			randomValues[i] = partition + randomValues[i-1] ;*/
		xVal = Data<T> (randomValues,parameters.numSamples) ;
		delete[] randomValues ;
	}

	computeFunctionValues() ;
	addNoise() ;
  /*if (parameters.sigma <= AOM)
		yVal = fxVal ;
	else
		addNoise() ;*/
}

/*!
 *	\relates RandomDataGenerator
 *	\brief This module generates the sawtooth function values for the 
 *	randomly generated x's. There are two versions of sawtooth functions
 * 	that can be used. If using (1), the function needs to be approximated
 *	using a cosine series; if using (2), sine series needs to be used.
 *	Sawtooth function is discontinuous and has an inversion centre
 *	about origin.
 *
 *	1. The sawtooth that passes through the origin.
 *									if (x < 0) => f(x) = -f(-x)
 *									if (x = 0) => f(x) = 0
 *									if (x > 0) => f(x) = f(x-nT), where T is the time period
 *				 further, if (0 < x < T/2) => f(x) = mx, where m is the slope
 *						 		  if (T/2 < x < T) => f(x) = -m(T-x) 	
 *							and if (x = T/2) => f(x) = 0
 *
 *	2. The sawtooth that does not pass through the origin.
 *									if (0 < x < T) => f(x) = 2p(x/T - 1/2), where
 *	T is the time period and p is the peak value of the wave
 *									if (x < 0) => f(x) = -f(-x)
 *									if (x = 0) => function discontinuous (randomly assign 
 *																												+/- peak value)
 *
 *	Currently, generating data as per version 2! If changed to version 1, 
 *	make sure that the harmonic is changed from sine to cosine in the file
 *	"OrthogonalBasis.h" in constructing the design matrix.
 *
 *  \param x a long double
 *	\param timePeriod a long double
 *  \param peak a long double
 *	\param slope a long double
 *	\return the sawtooth function value
 */
long double sawtooth (long double x, long double timePeriod, long double peak, long double slope)
{
	//	Sawtooth [1]
	/*
	if (x < 0)
		return (-1) * sawtooth (-x,timePeriod,peak,slope) ;
	else if (x > 0)
	{
		while (x > timePeriod)
			x = x - timePeriod ;
		if ( 2 * x < timePeriod)
			return slope * x ;
		else if (2 * x > timePeriod)
			return slope * (x - timePeriod) ;
		else
			return 0 ;
	}
	else
		return 0 ;
	*/
	//	Sawtooth [2]
	/*if (x < 0)
		return (-1) * sawtooth (-x,timePeriod,peak,slope) ;
	else if (x > 0)
	{
		while (x > timePeriod)
			x = x - timePeriod ;
		return 2 * peak * ((x / timePeriod) - 0.5) ;
	}
	else if (x == 0)
	{
		long double random = (long double)rand() / RAND_MAX ;
		if (random <= 0.5)
			return peak ;
		else return (-1)*peak ;
	}
  */
  //  Sawtooth [3]
  if (x >= 0 && x < timePeriod) {
    return x * peak / timePeriod ;
  }
  else if (x > timePeriod) {
		while (x > timePeriod) {
			x = x - timePeriod ;
    }
    return sawtooth(x,timePeriod,peak,slope) ;
	}
  else if (x < 0) {
    while (x < 0) {
      x = x + timePeriod ;
    }
    return sawtooth(x,timePeriod,peak,slope) ;
  }  
}

/*!
 *	\relates RandomDataGenerator
 *	\brief This module generates the square wave function values for the
 * 	randomly generated x's. There are two versions of square waves that 
 *	could be used. If using (1), the function needs to be approximated
 *	using odd cosine harmonics; if using (2) sine series needs to used.
 *	Square wave is discontinuous.
 *
 *	1. Square wave symmetric about Y-axis and has +ve peak value at the
 *		 origin, hence
 *					if (x < 0) => f(x) = f(-x)
 *				 	if (x = 0) => f(x) = peak
 *				 	if (x > 0) => f(x) = f(x-nT), where T is the time period
 *				 	further, if (0 < x < T/4) => f(x) = peak
 *								if (T/4 < x < 3T/4) => f(x) = -peak
 *								and if (3T/4 < x < T) => f(x) = peak
 *
 *	2. Square wave that is discontinuous at origin and has a centre of
 *		 inversion about the origin.
 *					if (x < 0) => f(x) = -f(-x)
 *				 	if (0 < x < T/2) => f(x) = peak
 *					if (T/2 < x < T) => f(x) = -peak
 *					if (x = 0 or x = T/2) => function discontinuous (randomly assign 
 *																													 +/- peak value)
 *
 *	Currently, generating data as per version 2! If changed to version 1, 
 *	make sure that the harmonic is changed from sine to cosine in the file
 *	"OrthogonalBasis.h" in constructing the design matrix.
 *
 *  \param x a long double
 *	\param timePeriod a long double
 *	\param peak a long double
 *	\return the square wave function value
 */
long double square (long double x, long double timePeriod, long double peak)
{
	// Square [1]
	/*
	if (x < 0)
		return square (-x,timePeriod,peak) ;
	else if (x > 0)
	{
		while (x > timePeriod)
			x = x - timePeriod ;
		if ((timePeriod < 4 * x) && (3 * timePeriod > 4 * x))
			return (-1) * peak ;
		else 
			return peak ;
	}
	else 
		return peak ;
	*/
	//	Square [2]
	if (x < 0) {
		return (-1) * square (-x,timePeriod,peak) ;
  } 
  else if (x >= 0) {
		while (x >= timePeriod) {
			x = x - timePeriod ;
    }
		if (x > 0 && (2 * x < timePeriod)) {
			return peak ;
    }
		else if ((2 * x > timePeriod) && x < timePeriod) {
			return (-1) * peak ;
    } else {
			long double random = (long double)rand() / RAND_MAX ;
			if (random <= 0.5)
				return peak ;
			else return (-1) * peak ;
		}
	}
}

/*!
 *	\relates RandomDataGenerator
 *	\brief This module generates the triangle wave function values for the
 * 	randomly generated x's. 
 *	Triangle wave is continuous.
 *
 *  If T is the time period and p is the peak value, then the functional 
 *  form of the triangle wave passing through the origin is as follows:-
 *      if (0 <= x <= T/2m)       => f(x) = 2mpx/T
 *      if (T/2m <= x < = T-T/2m) => f(x) = (1-2x/T)(mp/m-1)
 *      if (T-T/2m <= x <= T)     => f(x) = (x-T)(2mp/T)
 *
 *  The triangle wave can be symmetric/asymmetric depending on the 
 *  parameter m. If m = 2 (symmetric), the function reaches a peak at 
 *  x = T/4, otherwise the peak value will be at x = T/2m
 *  (asymmetrical).
 */
long double triangle (long double x, long double timePeriod, 
                      long double peak, long double m)
{
  if (x >=0 && x <= timePeriod) {
    long double bound1 = timePeriod / (2*m) ;
    long double bound2 = timePeriod - bound1 ; 
    if (x <= bound1) {
      return 2 * m * peak * x / timePeriod ;
    }
    else if (x > bound1 && x < bound2) {
      return (1 - (2 * x / timePeriod)) * (m * peak / (m - 1)) ;
    } else {
      return (x - timePeriod) * (2 * m * peak / timePeriod) ;
    }
  }
  else if (x > timePeriod) {
    while (x >= timePeriod) {
      x = x - timePeriod ;
    }
    return triangle(x,timePeriod,peak,m) ;
  } 
  else if (x < 0) {
    while (x < 0) {
      x = x + timePeriod ;
    }
    return triangle(x,timePeriod,peak,m) ;
  }
}

/*!
 *	\relates RandomDataGenerator
 *	\brief 
 */
long double finiteLinearCombination (long double x, long double timePeriod,
                                      lcb::Vector<long double> weights)
{
	long double pi = boost::math::constants::pi<long double>() ;
	int M = weights.size(),k ;
	long double arg,yVal = weights[0] ;
  //arg = 2 * pi * x / timePeriod ;
  //arg = x / timePeriod ;
	for (int i=1; i<M; i++)
	{
    arg = x ;
    if (i % 2 == 1)
    {   
      k = i / 2 + 1 ; 
      arg *= k  ;
      //cout << "sin(" << k<< "*" << x/timePeriod << ") = " << sin(arg) << endl;
      yVal += weights[i] * sin (arg) ;
      /*cout << "w = " << weights[i] << endl ;
      cout << "w*sin(" << k<< "*" << x/timePeriod << ") = " << weights[i]*sin(arg) << endl;
      cout << "yVal:" << yVal << endl;*/
    }   
    else
    {   
      k = i / 2 ; 
      arg *= k ;
      //cout << "cos(" << k<< "*" << x/timePeriod << ") = " << cos(arg) << endl;
      yVal += weights[i] * cos (arg) ;
      /*cout << "w = " << weights[i] << endl ;
      cout << "w*cos(" << k<< "*" << x/timePeriod << ") = " << weights[i]*cos(arg) << endl;
      cout << "yVal:" << yVal << endl;*/
    }
	}
	return yVal ;	
}

/*!
 *	\fn void RandomDataGenerator<T> :: computeFunctionValues (void)
 *	\brief This module is used to compute the corresponding function
 *	values for the X's generated randomly
 */
template <class T>
void RandomDataGenerator<T> :: computeFunctionValues (void)
{
	T *y ;
	long double slope ;
	int M ;
	lcb::Vector<long double> weights ;
	switch (parameters.function)
	{
		case 0:		// sawtooth
			y = new T [parameters.numSamples] ;
			slope = (2 * parameters.peak) / parameters.timePeriod ;
			for (int i=0; i<xVal.nPoints(); i++)
			{
				long double randomX = xVal[i].x() ;
				y[i] = sawtooth(randomX,parameters.timePeriod,parameters.peak,
                        slope) ;
        //cout << randomX << " " << y[i] << endl ;
			}	
			fname = "SAWTOOTH" ;
			break ;
		case 1:		// square wave
			y = new T [parameters.numSamples] ;
			for (int i=0; i<xVal.nPoints(); i++)
			{
				long double randomX = xVal[i].x() ;
				y[i] = square (randomX,parameters.timePeriod,parameters.peak) ;
			}
			fname = "SQUARE" ;
			break ;
    case 2:   // triangular wave
			y = new T [parameters.numSamples] ;
			for (int i=0; i<xVal.nPoints(); i++)
			{
				long double randomX = xVal[i].x() ;
				y[i] = triangle (randomX,parameters.timePeriod,parameters.peak,
                         parameters.m) ;
			}
      fname = "TRIANGLE" ;
      break ;
		case 3:		// using a finite linear combination
			//srand (time(NULL)) ;
      srand(100) ;
			y = new T [parameters.numSamples] ;
			// generate the number of terms (between 3 and 100)
			//M = rand() % (MAX_TERMS-MIN_TERMS+1) + MIN_TERMS ;
      M = 4 ;
			weights = lcb::Vector<long double>(M) ;
			for (int i=0; i<M; i++) {
        // weights are generated randomly in [-1,1]
				weights[i] = 2 * (rand() /(long double) RAND_MAX) - 1 ;
			}
      //weights[0] = -10.1;
      weights[2] = 0;
      weights[3] = 0;
			cout << "M_gen: " << M << endl ;
			weights.print() ;
      
      constructFunctionString(finlincomb,weights,M);
      cout << finlincomb << endl ;
			for (int i=0; i<xVal.nPoints(); i++)
			{
				long double randomX = xVal[i].x() ;
        //cout << "X: " << randomX << endl ;
				y[i] = finiteLinearCombination(randomX,parameters.timePeriod,
                                       weights) ;
        //cout << "exiting ..." << endl ; exit(0) ;
			}
			fname = "FINITE LINEAR COMBINATION" ;
			break ;
		default:
			error ("Function index not appropriate!") ;
			break ;
	}
	fxVal = Data<T> (y,parameters.numSamples) ;
	delete[] y ;
}

template <class T>
string convertToString(T number)
{
  ostringstream convert ;
  convert << number ;
  return convert.str() ;
}

template <class T>
void RandomDataGenerator<T> :: constructFunctionString (string &s, lcb::Vector<T> &weights, int M)
{
  string w,coeffx,fn;
  s = "f(x) = ";
  if (fabs(weights[0]) > std::numeric_limits<double>::epsilon()) {
    w = convertToString<long double>(weights[0]);
    s = s + w ;
  }
  int k;
  for (int i=1; i<M; i++) {
    if (fabs(weights[i]) > std::numeric_limits<double>::epsilon()) {
      if (weights[i] > 0) {
        s = s + " + " ;
      } else {
        s = s + " - " ;
      }
      w = convertToString<long double>(fabs(weights[i]));
      s = s + w + "*";
      if (i % 2 == 0) {
        k = i / 2 ;
        fn = "cos ";
      } else {
        k = i / 2 + 1;
        fn = "sin " ;
      }
      if (k != 1) {
        coeffx = convertToString<int>(k);
      } else {
        coeffx = "" ;
      }
      s = s + fn + coeffx + "x" ;
    }
  }  
  //cout << s << endl ;
}
 
/*!
 *	\fn void RandomDataGenerator<T> :: addNoise (void)
 *	\brief This module adds Gaussian noise to the already generated 
 *	data samples
 */
template <class T>
void RandomDataGenerator<T> :: addNoise (void)
{
	Gaussian noise (parameters.mean,parameters.sigma) ;
	T *y = new T [parameters.numSamples] ;
	for (int i=0; i<parameters.numSamples; i++)
	{
		vector<long double> samples = noise.generate() ;
		y[i] = fxVal[i].x()+samples[0] ;
		if (i != parameters.numSamples - 1)
		{
			i++ ;
			y[i] = fxVal[i].x()+samples[1] ;
		}
	}
	yVal = Data<T> (y,parameters.numSamples) ;
	delete[] y ;
}

/*!
 *	\fn Data<T> RandomDataGenerator<T> :: predict (lcb::Matrix<T> &weights, 
 *	Data<T> &xVals)
 *	\brief This function predicts the value of a given x using the weights it 
 *	computed from the regression analysis
 *	\param weights a reference to a Matrix object of type T
 *	\param xVals a reference to a Data object of type T
 *	\return a Data object containing the predicted values
 */
template <class T>
Data<T> RandomDataGenerator<T> :: predict (unsigned M, lcb::Matrix<T> 
                                  &weights, Data<T> &xVals)
{
	parameters.numFunctions = M ;
	OrthogonalBasis orthogonal (parameters.basis,parameters.numFunctions,
                              parameters.timePeriod,parameters.function) ;
	lcb::Matrix<long double> designMatrix ;
	designMatrix = orthogonal.designMatrix(xVals) ;
	lcb::Matrix<T> yEst = designMatrix * weights ; // column matrix
	return Data<T>(yEst) ;
}

/*!
 *	\fn Data<T> RandomDataGenerator<T> :: randomX (void)
 *	\brief This module returns the randomly samples X values.
 *	\return random Xs generated 
 */
template <class T>
Data<T> RandomDataGenerator<T> :: randomX (void)
{
	return xVal ;
}

/*!
 *	\fn Data<T> RandomDataGenerator<T> :: fxValues (void)
 *	\brief This module returns funtion values of the corresponding
 *	input x values.
 *	\return function values
 */
template <class T>
Data<T> RandomDataGenerator<T> :: fxValues (void)
{
	return fxVal ;
}

/*!
 *	\fn Data<T> RandomDataGenerator<T> :: yValues (void)
 *	\brief This module returns the noise added funtion values.
 *	\return y values generated 
 */
template <class T>
Data<T> RandomDataGenerator<T> :: yValues (void)
{
	return yVal ;
}

/*!
 *  \fn void RandomDataGenerator<T> :: plotRandomX (void)
 *  \brief The function plots the random X values sampled
 */
template <class T>
void RandomDataGenerator<T> :: plotRandomX (void)
{
	vector<string> labels ;
	long double min,max ;
	labels.push_back(fname) ;
	labels.push_back("#") ;
	labels.push_back("x") ;

	Plot graph ;
	graph.label(labels) ;
	
	pair<long double,long double> xrange,yrange ;
	xrange = make_pair(1,parameters.numSamples) ;
	yrange = make_pair(xVal.minimum()-0.5,xVal.maximum()+0.5) ;
	graph.setRange(xrange,yrange) ;
	graph.sketch(xVal) ;
}

/*!
 *  \fn void RandomDataGenerator<T> :: plotData (void)
 *  \brief The function plots the data elements and the corresponding 
 *	function values 
 */
template <class T>
void RandomDataGenerator<T> :: plotData (void)
{
	vector<string> labels ;
	long double min,max ;
	labels.push_back(fname) ;
	labels.push_back("x") ;
	labels.push_back("f(x)") ;

	Plot graph ;
	graph.label(labels) ;
	
	pair<long double,long double> xrange,yrange ;
	xrange = make_pair(parameters.low-0.5,parameters.high+0.5) ;
	yrange = make_pair(fxVal.minimum()-0.5,fxVal.maximum()+0.5) ;
	graph.setRange(xrange,yrange) ;
	graph.sketch(xVal,fxVal) ;
}

/*!
 *  \fn void RandomDataGenerator<T> :: plotDataWithNoise (void)
 *  \brief The function plots the data elements generated 
 */
template <class T>
void RandomDataGenerator<T> :: plotDataWithNoise (void)
{
	vector<string> labels ;
	long double min,max ;
	labels.push_back(fname) ;
	labels.push_back("x") ;
	labels.push_back("y=f(x)+e") ;

	Plot graph ;
	graph.label(labels) ;
	
	pair<long double,long double> xrange,yrange ;
	xrange = make_pair(parameters.low-0.5,parameters.high+0.5) ;
	long double extreme[4] ;
	extreme[0] = fxVal.minimum() ;
	extreme[1] = fxVal.maximum() ;
	extreme[2] = yVal.minimum() ;
	extreme[3] = yVal.maximum() ;
	min = extreme[0] ;
	max = extreme[0] ;
	for (int i=1; i<4; i++)
	{
		if (extreme[i] < min)
			min = extreme[i] ;
		if (extreme[i] > max)
			max = extreme[i] ;
	}
	yrange = make_pair(min-0.5,max+0.5) ;
	graph.setRange(xrange,yrange) ;
	graph.sketch(xVal,fxVal,yVal) ;
}

/*!
 *	\fn void RandomDataGenerator<T> :: plotPredictions (Data<T> &xVals, 
 *	Data<T> &yVals, Data<T> &predictions) 
 *	\brief The function is used to plot the Y values and predicted output 
 *	values over one other (multiplot) against the X values. This is to give
 *	a visual idea of how the regression model fits the experimental data.
 *	\param xVals a reference to a Data object of type T
 *	\param yVals a reference to a Data object of type T
 *	\param predictions a reference to a Data object of type T
 */
template <class T>
void RandomDataGenerator<T> :: plotPredictions (Data<T> &xVals, 
                               Data<T> &yVals, Data<T> &predictions) 
{
	vector<string> labels ;
	long double min,max ;
	labels.push_back(fname) ;
	labels.push_back("x") ;
	labels.push_back("predictions") ;

	Plot graph ;
	graph.label(finlincomb,labels) ;
	
	pair<long double,long double> xrange,yrange ;
	xrange = make_pair(parameters.low-0.5,parameters.high+0.5) ;
	long double extreme[4] ;
	extreme[0] = yVals.minimum() ;
	extreme[1] = yVals.maximum() ;
	extreme[2] = predictions.minimum() ;
	extreme[3] = predictions.maximum() ;
	min = extreme[0] ;
	max = extreme[0] ;
	for (int i=1; i<4; i++)
	{
		if (extreme[i] < min)
			min = extreme[i] ;
		if (extreme[i] > max)
			max = extreme[i] ;
	}
	yrange = make_pair(min-0.5,max+0.5) ;
	graph.setRange(xrange,yrange) ;
	graph.sketch(xVals,yVals,predictions) ;
}

/*!
 *  \fn void RandomDataGenerator<T> :: print (void)
 *  \brief The function displays information about the 
 *  private members of the class
 */
template <class T>
void RandomDataGenerator<T> :: print (void)
{
	cout << "Mean of noise = " << parameters.mean << endl ; 
	cout << "Sigma of noise = " << parameters.sigma << endl ;
	cout << "Lower bound of interval = " << parameters.low << endl ;
	cout << "Upper bound of interval = " << parameters.high << endl ;
	cout << "Function index = " << parameters.function << endl ;
	cout << "Time Period = " << parameters.timePeriod << endl ;
	cout << "Peak value = " << parameters.peak << endl ;
	cout << "Number of samples = " << parameters.numSamples << endl ; 
	cout << "# of orthogonal fucntions = " << parameters.numFunctions 
  << endl ;
	if (!parameters.file.empty()) {
		cout << "Input file: " << parameters.file << endl ;
  }

	if (xVal.nPoints() > 0) {
		xVal.print() ;
  }
	if (fxVal.nPoints() > 0) {
		fxVal.print() ;
  }
}

#endif

