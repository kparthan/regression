#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include "Error.h"
#include <liblcb/Matrix.h>
#include <liblcb/Vector.h>
#include "Data.h"
#include "RandomDataGenerator.h"
#include "OrthogonalBasis.h"
#include "Message.h"

using namespace std ;

struct Parameters parseCommandLine (int argc, char **argv)
{
	long double mean = 0 ;              // -gmean
	long double sigma = 0.2 ;						// -gsigma
	long double low = -1 ;							// -low
	long double high = 1 ;							// -high
	string fname = "sawtooth" ;			    // -fn
	int function = 0 ;						  
	long double timePeriod = 0.1 ;		  // -t
	long double peak = 1 ;							// -peak
	int numSamples = 100 ;						  // -nsamples
	int numFunctions = 3 ;					    // -nof
	string file ;										    // -file
  int iterate = 0 ;                   // -iterate 
                                      //  [1 for yes, 0 for no]
  int inverse = 0 ;                   // -inv
                                      // 0 -- my implementation
<<<<<<< HEAD
                                      // 1 -- from boost
                                      // 2 -- LU decomposition
=======
                                      // 1 -- from boost 
>>>>>>> 57d61fedf3014331915abae9f24ac44ddf0b747b
	bool paramFlags[12] = {0} ;
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
			else if(fname.compare("finlincomb") == 0)
				function = 2 ;
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
		else if (string(argv[i]).compare("-iterate") == 0)
		{
			iterate = atoi(argv[i+1]) ;
			paramFlags[10] = 1 ;
		}
		else if (string(argv[i]).compare("-inv") == 0)
		{
			inverse = atoi(argv[i+1]) ;
			paramFlags[11] = 1 ;
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
      cout << "\t [-iterate] choice of iteration over different values" << endl ;
      cout << "\t [-inv] choice of matrix inverse" << endl ;
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
			error ("Interval's upper bound should be less than the " 
      "lower bound ...") ;

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
			cout << "Using default value for interval's lower bound: " << low 
			<< endl ;
		else
			cout << "Interval lower bound set to: " << low << endl ;

		if (paramFlags[3] == 0)
			cout << "Using default value for interval's higher bound: " << high 
			<< endl ;
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
		cout << "Using default value for maximum amplitude of wave: " << peak 
		<< endl ;
	else
		cout << "Peak value set to: " << peak << endl ;

	if (paramFlags[8] == 0)
		cout << "Using default number of orthogonal functions(terms): " << numFunctions 
		<< endl ;
	else
		cout << "Number of orthogonal functions(terms) set to: " << numFunctions << endl ;

	if (paramFlags[9] == 0)
		cout << "Using data generated randomly ..." << endl ;
	else
		cout << "Using data from file: " << file << "..."  << endl ; 

  if (paramFlags[10] == 0)
    cout << "Running an instance of regression fit ..." << endl ;
  else
    cout << "Iterating over different values of parameters ..." << endl ; 

<<<<<<< HEAD
  if(paramFlags[11]) {
    switch(inverse) {
      case 0:
        cout << "Using my implementation of inverse [using " 
        "partial pivoting] ..." << endl ;
        break ;
      case 1:
        cout << "Using BOOST library implmentation of matrix inverse " 
        "..." << endl ;
        break ;
      case 2:
        cout << "Using my implementation of LU Decomposition to " 
        "solve linear system ..." << endl ;
        break ;
      default:
        error("Invalid choice of matrix inverse.") ;
        break ;
    }
  }
  else
    cout << "Using my implementation of inverse [using \
    partial pivoting] ..." << endl ;
=======
  switch(paramFlags[11]) {
    case 0:
      cout << "Using my implementation of inverse [using partial pivoting] ..." << endl ;
      break ;
    case 1:
      cout << "Using BOOST library implmentation of matrix inverse ..." << endl ;
      break ;
    default:
      error("Invalid choice of matrix inverse.") ;
      break ;
  }
>>>>>>> 57d61fedf3014331915abae9f24ac44ddf0b747b
  
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
  params.iterate = iterate ;
  params.inverse = inverse ;

	return params ;
}

void setPrecision(void)
{
	cout.unsetf(ios::floatfield) ;
	int PRECISION = -log10(AOM) ;
	cout.precision(PRECISION) ;
	cout.setf(ios::fixed,ios::floatfield) ;
}

template <class T>
string convertToString(T number)
{
	ostringstream convert ;
	convert << number ;
	return convert.str() ;
}

int main(int argc, char **argv)
{
	//setPrecision() ;
	struct Parameters parameters = parseCommandLine(argc,argv) ;
  RandomDataGenerator<long double> dataGenerator ;
	lcb::Matrix<long double> phi ;
	lcb::Matrix<long double> weights ;
  long double rmse,msgLen ;
	Data<long double> randomX,yValues,predictions ;
  OrthogonalBasis orthogonal ;
  Message msg ;
	string filename ; 
<<<<<<< HEAD
  int sampVals[] = {100} ;
  std::vector<int> Samples (sampVals,sampVals+sizeof(sampVals)/sizeof(int)) ;
  //long double noiseVals[] = {0.1,0.2,0.3,0.4,0.5} ;
  long double noiseVals[] = {0.5} ;
=======
  int sampVals[] = {1000} ;
  std::vector<int> Samples (sampVals,sampVals+sizeof(sampVals)/sizeof(int)) ;
  long double noiseVals[] = {0.75} ;
>>>>>>> 57d61fedf3014331915abae9f24ac44ddf0b747b
  std::vector<long double> Noise (noiseVals,noiseVals+sizeof(noiseVals)/sizeof(long double)) ;
	
  switch(parameters.iterate) 
  {
    case 0:
      dataGenerator = RandomDataGenerator<long double>(parameters) ;
	    dataGenerator.generate() ;
	    randomX = dataGenerator.randomX() ;
	    //yValues = dataGenerator.yValues() ;
	    yValues = dataGenerator.fxValues() ;
	    dataGenerator.plotData() ;
	    //dataGenerator.plotDataWithNoise() ;
			orthogonal = OrthogonalBasis (parameters.numFunctions,
					parameters.timePeriod,parameters.function) ;
			phi = orthogonal.designMatrix(randomX) ;

			weights = computeWeights<long double>(phi,yValues,parameters.inverse) ;
<<<<<<< HEAD
      weights.print() ;
=======
      //weights.print() ;
>>>>>>> 57d61fedf3014331915abae9f24ac44ddf0b747b

			predictions = dataGenerator.predict(parameters.numFunctions,weights,
                                          randomX) ;
			dataGenerator.plotPredictions(randomX,yValues,predictions) ;
	  	rmse = computeRMSE<long double>(weights,phi,yValues) ;
			cout << "Error in fitting: " << rmse << endl ;
			msg = Message (parameters,weights,randomX,yValues,predictions) ;
			msgLen = msg.messageLength() ;
			//cout << "Msg Len = " << msgLen << endl ;
      break ;

    case 1:
	    for (unsigned i=0; i<Samples.size(); i++) 
      {
		    parameters.numSamples = Samples[i] ;
		    for (unsigned j=0; j<Noise.size(); j++)
		    {
<<<<<<< HEAD
			    filename = "temp/results_n" + convertToString<int>(Samples[i]) + "_s" ;
=======
			    filename = "results/results_n" + convertToString<int>(Samples[i]) + "_s" ;
>>>>>>> 57d61fedf3014331915abae9f24ac44ddf0b747b
			    filename = filename + convertToString<long double>(Noise[j]) + 
                      ".txt" ;
			    ofstream results ;
			    results.open(filename.c_str()) ;	
			    parameters.sigma = Noise[j] ;

          dataGenerator = RandomDataGenerator<long double>(parameters) ;
			    dataGenerator.generate() ;
			    randomX = dataGenerator.randomX() ;
			    yValues = dataGenerator.yValues() ;

			    for (unsigned M=1; M<100; M++) 
			    {
				    cout << "N: " << parameters.numSamples << "\t" ;
				    cout << "S: " << parameters.sigma << "\t" ;
				    cout << "M: " << M << endl ;
				    if (Samples[i] > M+15)
				    {
					    parameters.numFunctions = M ;
					    orthogonal = OrthogonalBasis (parameters.numFunctions,
					            parameters.timePeriod,parameters.function) ;
					    phi = orthogonal.designMatrix(randomX) ;
					    weights = computeWeights<long double>(phi,yValues,parameters.inverse) ;
					    predictions = dataGenerator.predict(M,weights,randomX) ;

					    rmse = computeRMSE<long double> (weights,phi,yValues) ;
					    //cout << "Error in fitting: " << rmse << endl ;

					    msg = Message (parameters,weights,randomX,yValues,
                              predictions) ;
					    msgLen = msg.messageLength() ;
<<<<<<< HEAD
=======
					    //cout << "Msg Len = " << msgLen << endl ;
              /*if (M > 48)i
              {
                cout << "det@M = " << M << " is " << 
              }*/
>>>>>>> 57d61fedf3014331915abae9f24ac44ddf0b747b
					    results << parameters.numFunctions << "\t" ;
					    results << rmse << "\t" ;
					    results << msgLen << endl ;
				    }
          }
          results.close() ;
        }
			}
		  break ;
    default:
      error("Wrong choice entered ...") ;
      break ;
  }
	return 0 ;
}

