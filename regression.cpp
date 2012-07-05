#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include "Error.h"
#include "Vector.h"
#include "Matrix.h"
#include "Data.h"
#include "RandomDataGenerator.h"
#include "OrthogonalBasis.h"
#include "Message.h"

using namespace std ;

struct Parameters parseCommandLine (int argc, char **argv)
{
	double mean = 0 ;								// -gmean
	double sigma = 1 ;							// -gsigma
	double low = -1 ;								// -low
	double high = 1 ;								// -high
	string fname = "sawtooth" ;			// -fn
	int function = 0 ;						  
	double timePeriod = 0.1 ;				// -t
	double peak = 1 ;								// -peak
	int numSamples = 20 ;						// -nsamples
	int numFunctions = 3 ;					// -nof
	string file ;										// -file
	bool paramFlags[10] = {0} ;
	int i = 1 ;

	while (i < argc)
	{
		if (string(argv[i]).compare("-gmean") == 0)
		{
			mean = atof(argv[i+1]) ;
			paramFlags[0] = 1 ;
		}
		else if (string(argv[i]).compare("-gsigma") == 0)
		{
			sigma = atof(argv[i+1]) ;
			paramFlags[1] = 1 ;
		}
		else if (string(argv[i]).compare("-low") == 0)
		{
			low = atof(argv[i+1]) ;
			paramFlags[2] = 1 ;
		}
		else if (string(argv[i]).compare("-high") == 0)
		{
			high = atof(argv[i+1]) ;
			paramFlags[3] = 1 ;
		}
		else if (string(argv[i]).compare("-fn") == 0)
		{
			fname = argv[i+1] ;
			if (fname.compare("sawtooth") == 0)
				function = 0 ;
			else if (fname.compare("square") == 0)
				function = 1 ;
			else
				error ("Function not supported ...") ;
			paramFlags[4] = 1 ;
		}
		else if (string(argv[i]).compare("-t") == 0)
		{
			timePeriod = atof(argv[i+1]) ;
			paramFlags[5] = 1 ;
		}
		else if (string(argv[i]).compare("-peak") == 0)
		{
			peak = atof(argv[i+1]) ;
			paramFlags[6] = 1 ;
		}
		else if (string(argv[i]).compare("-nsamples") == 0)
		{
			numSamples = atoi(argv[i+1]) ;
			if (numSamples <= 0)
				error ("# of data samples should be non-negative ...") ;
			paramFlags[7] = 1 ;
		}
		else if (string(argv[i]).compare("-nof") == 0)
		{
			numFunctions = atoi(argv[i+1]) ;
			if (numFunctions <= 0)
				error ("# of orthogonal functions should be non-negative ...") ;
			paramFlags[8] = 1 ;
		}
		else if (string(argv[i]).compare("-file") == 0)
		{
			file = argv[i+1] ;
			paramFlags[9] = 1 ;
		}
		else
		{
			cout << "Usage: " << argv[0] << " [options]" << endl ;
			cout << "Valid options:" << endl ;
			cout << "\t [-gmean] Mean of Gaussian distribution" << endl ;
			cout << "\t [-gsigma] Sigma of Gaussian distribution" << endl ;
			cout << "\t [-low] lower bound of interval" << endl ;
			cout << "\t [-high] upper bound of interval" << endl ;
			cout << "\t [-fn] function to be used" << endl ;
			cout << "\t [-t] time period of wave function" << endl ;
			cout << "\t [-peak] maximum amplitude of the wave" << endl ;
			cout << "\t [-nsamples] number of data samples to be used" << endl ;
			cout << "\t [-nof] number of orthogonal basis functions to use" << endl ;
			cout << "\t [-file] input file containing data samples" << endl ;
			error ("Invalid command line argument ...") ;
		}
		i += 2 ;
	}	

	if (paramFlags[9] == 1)
	{
		if (paramFlags[2] == 1)
			cout << "Ignoring argument \"-low\"/default value set" << endl ;
		if (paramFlags[3] == 1)
			cout << "Ignoring argument \"-high\"/default value set" << endl ;
		if (paramFlags[7] == 1)
			cout << "Ignoring argument \"-nsamples\"/default value set" << endl ;
	}

	if (paramFlags[9] != 1)
	{
		if (high <= low)
			error ("Interval's upper bound should be less than the lower bound ...") ;

		if (numSamples <= 0)
			error ("Number of data samples should be non-negative ...") ;
	}

	//	printing parameter values to be used in simulation
	if (paramFlags[0] == 0)
		cout << "Using default value for Gaussian Mean: " << mean <<endl ;
	else
			cout << "Gaussian Mean set to: " << mean << endl ;

	if (paramFlags[1] == 0)
		cout << "Using default value for Gaussian sigma: " << sigma << endl ;
	else
			cout << "Gaussian Sigma set to: " << sigma << endl ;

	if (paramFlags[9] != 1)
	{
		if (paramFlags[2] == 0)
			cout << "Using default value for interval's lower bound: " << low << endl ;
		else
			cout << "Interval lower bound set to: " << low << endl ;

		if (paramFlags[3] == 0)
			cout << "Using default value for interval's higher bound: " << high << endl ;
		else
			cout << "Interval higher bound set to: " << high << endl ;
	
		if (paramFlags[7] == 0)
			cout << "Using default number of samples: " << numSamples << endl ;
		else
			cout << "Number of data samples set to: " << numSamples << endl ;
	}

	if (paramFlags[4] == 0)
		cout << "Using default function to generate data: " << fname << endl ;
	else
		cout << "Function to generate data set to: " << fname << endl ;
	
	if (paramFlags[5] == 0)
		cout << "Using default value for time period: " << timePeriod << endl ;
	else
		cout << "Time period set to: " << timePeriod << endl ;

	if (paramFlags[6] == 0)
		cout << "Using default value for maximum amplitude of wave: " << peak << endl ;
	else
		cout << "Peak value set to: " << peak << endl ;

	if (paramFlags[8] == 0)
		cout << "Using default number of orthogonal functions: " << numFunctions << endl ;
	else
		cout << "Number of orthogonal functions set to: " << numFunctions << endl ;

	if (paramFlags[9] == 0)
		cout << "Using data generated randomly" << endl ;
	else
		cout << "Using data from file: " << file << endl ; 

	struct Parameters params ;
	params.mean = mean ;
	params.sigma = sigma ;
	params.low = low ;
	params.high = high ;
	params.function = function ;
	params.timePeriod = timePeriod ;
	params.peak = peak ;
	params.numSamples = numSamples ;
	params.numFunctions = numFunctions ;
	params.file = file ;

	return params ;
}

template <class T>
Matrix<T> computeWeights (Matrix<T> &phi, Data<T> &yValues)
{
	Matrix<T> phiT = phi.transpose() ;
	Matrix<T> phiTphi = phiT * phi ;
	Matrix<T> pseudoInv = phiTphi.inverse() ;
	Matrix<T> temp = pseudoInv * phiT ;

	Matrix<T> y = yValues.convertToMatrix() ;
	Matrix<T> weights = temp * y ;
	return weights ;
}

template <class T>
double computeRMSE (Matrix<T> &weights, Matrix<T> &phi, Data<T> &yVals)
{
	Matrix<T> yEst = phi * weights ; // column matrix
	double diff, error = 0 ;
	int numSamples = phi.rows() ;
	for (int i=0; i<numSamples; i++)
	{
		diff = yEst[i][0] - yVals[i].x() ;
		error += diff * diff ;
	}
	return sqrt(error/numSamples) ;
}

int main(int argc, char **argv)
{
	struct Parameters parameters = parseCommandLine(argc,argv) ;
	RandomDataGenerator<double> dataGenerator (parameters) ;

	dataGenerator.generate() ;
	Data<double> randomX = dataGenerator.randomX() ;
	//Data<double> yValues = dataGenerator.yValues() ;
	Data<double> yValues = dataGenerator.fxValues() ;
	//dataGenerator.plotData() ;
	//dataGenerator.plotDataWithNoise() ;

	Matrix<double> phi ;
	OrthogonalBasis orthogonal (parameters.numFunctions,parameters.timePeriod,
																							parameters.function) ;
	phi = orthogonal.designMatrix(randomX) ;

	Matrix<double> weights ;
	weights = computeWeights<double>(phi,yValues) ;
	weights.print() ;

	Data<double> predictions ;
	predictions = dataGenerator.predict(weights,randomX) ;
	dataGenerator.plotPredictions(randomX,yValues,predictions) ;

	Message msg (parameters.numFunctions,weights,randomX,yValues) ;
	msg.messageLength() ;

	double rmse = computeRMSE<double>(weights,phi,yValues) ;
	cout << "Error in fitting: " << rmse << endl ;

	return 0 ;
}

