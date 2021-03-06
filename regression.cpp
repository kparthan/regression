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
	long double timePeriod = 1 ;		    // -t
	long double peak = 1 ;							// -peak
	int numSamples = 100 ;						  // -nsamples
	int numFunctions = 3 ;					    // -nof
	string file ;										    // -file
  int iterate = 0 ;                   // -iterate 
                                      //  [1 for yes, 0 for no]
  int inverse = 0 ;                   // -inv
                                      // 0 -- my implementation
                                      // 1 -- from boost
                                      // 2 -- LU decomposition
  int basis = 0 ;                     // 0 -- sines & cosines
                                      // 1 -- Legendre Polynomials
  long double m_triangle = 2.0 ;      // determines point at which a 
                                      // triangle wave reaches its peak
                                      // value
  long double lambda = 1 ;
  int printComponents = 0;
	bool paramFlags[16] = {0} ;
	int i = 1 ;

	while (i < argc)
	{
		if (string(argv[i]).compare("-gmean") == 0) {
			mean = atof(argv[i+1]) ;
			paramFlags[0] = 1 ;
		}
		else if (string(argv[i]).compare("-gsigma") == 0) {
			sigma = atof(argv[i+1]) ;
			paramFlags[1] = 1 ;
		}
		else if (string(argv[i]).compare("-low") == 0) {
			low = atof(argv[i+1]) ;
			paramFlags[2] = 1 ;
		}
		else if (string(argv[i]).compare("-high") == 0) {
			high = atof(argv[i+1]) ;
			paramFlags[3] = 1 ;
		}
		else if (string(argv[i]).compare("-print") == 0) {
			printComponents = atoi(argv[i+1]) ;
		}
		else if (string(argv[i]).compare("-fn") == 0) {
			fname = argv[i+1] ;
			if (fname.compare("sawtooth") == 0) {
				function = 0 ;
      } 
      else if (fname.compare("square") == 0) {
				function = 1 ;
      }
			else if(fname.compare("triangle") == 0) {
				function = 2 ;
      }
			else if(fname.compare("parabola") == 0) {
				function = 3 ;
      }
			else if(fname.compare("finlincomb") == 0) {
				function = 4 ;
      } else {
				error ("Function not supported ...") ;
      }
			paramFlags[4] = 1 ;
		}
		else if (string(argv[i]).compare("-t") == 0) {
			timePeriod = atof(argv[i+1]) ;
			paramFlags[5] = 1 ;
		}
		else if (string(argv[i]).compare("-peak") == 0) {
			peak = atof(argv[i+1]) ;
			paramFlags[6] = 1 ;
		}
		else if (string(argv[i]).compare("-nsamples") == 0) {
			numSamples = atoi(argv[i+1]) ;
			if (numSamples <= 0) {
				error ("# of data samples should be non-negative ...") ;
      }
			paramFlags[7] = 1 ;
		}
		else if (string(argv[i]).compare("-nof") == 0) {
			numFunctions = atoi(argv[i+1]) ;
			if (numFunctions <= 0) {
				error ("# of orthogonal functions should be non-negative ...") ;
      }
			paramFlags[8] = 1 ;
		}
		else if (string(argv[i]).compare("-file") == 0) {
			file = argv[i+1] ;
			paramFlags[9] = 1 ;
		}
		else if (string(argv[i]).compare("-iterate") == 0) {
			iterate = atoi(argv[i+1]) ;
			paramFlags[10] = 1 ;
		}
		else if (string(argv[i]).compare("-inv") == 0) {
			inverse = atoi(argv[i+1]) ;
			paramFlags[11] = 1 ;
		}
		else if (string(argv[i]).compare("-basis") == 0) {
			basis = atoi(argv[i+1]) ;
			paramFlags[12] = 1 ;
		}
		else if (string(argv[i]).compare("-m") == 0) {
			m_triangle = atof(argv[i+1]) ;
      if (m_triangle < 1) {
        error("parameter m in triangle wave specification should be "
        "greater than 1 ...") ;
      }
			paramFlags[13] = 1 ;
    }
		else if (string(argv[i]).compare("-lambda") == 0) {
      lambda = atof(argv[i+1]);
      paramFlags[14] = 1;
		} else {
			cout << "Usage: " << argv[0] << " [options]" << endl ;
			cout << "Valid options:" << endl ;
			cout << "\t [-gmean] Mean of Gaussian distribution" << endl ;
			cout << "\t [-gsigma] Sigma of Gaussian distribution" << endl ;
			cout << "\t [-low] lower bound of interval" << endl ;
			cout << "\t [-high] upper bound of interval" << endl ;
			cout << "\t [-fn] function to be used" << endl ;
      cout << "\t [-m] >1 point at which a triangle wave reaches " 
      "its peak/maximum value" << endl ;
			cout << "\t [-t] time period of wave function" << endl ;
			cout << "\t [-peak] maximum amplitude of the wave" << endl ;
			cout << "\t [-nsamples] number of data samples to be used" << endl ;
			cout << "\t [-nof] number of orthogonal basis functions to use" 
      << endl ;
			cout << "\t [-file] input file containing data samples" << endl ;
      cout << "\t [-iterate] choice of iteration over different values" 
      << endl ;
      cout << "\t [-inv] choice of matrix inverse" << endl ;
      cout << "\t [-basis] choice of orthogonal basis set" << endl ;
			error ("Invalid command line argument ...") ;
		}
		i += 2 ;
	}	

	if (iterate == 0) {
		if (paramFlags[9] == 1)
		{
			if (paramFlags[2] == 1) {
				cout << "Ignoring argument \"-low\"/default value set" << endl ;
      }
			if (paramFlags[3] == 1) {
				cout << "Ignoring argument \"-high\"/default value set" << endl ;
      }
			if (paramFlags[7] == 1) {
				cout << "Ignoring argument \"-nsamples\"/default value set" 
        << endl ;
      }
		}

		if (paramFlags[9] != 1) {
			if (high <= low) {
				error ("Interval's upper bound should be less than the " 
	      "lower bound ...") ;
      }
			if (numSamples <= 0) {
				error ("Number of data samples should be non-negative ...") ;
      }
		}

		//	printing parameter values to be used in simulation
		if (paramFlags[0] == 0) {
			cout << "Using default value for Gaussian Mean: " << mean <<endl ;
    } else {
				cout << "Gaussian Mean set to: " << mean << endl ;
    }

		if (paramFlags[1] == 0) {
			cout << "Using default value for Gaussian sigma: " << sigma << endl ;    } else {
				cout << "Gaussian Sigma set to: " << sigma << endl ;
    }

		if (paramFlags[9] != 1) {
			if (paramFlags[7] == 0) {
				cout << "Using default number of samples: " << numSamples << endl ;
      } else {
				cout << "Number of data samples set to: " << numSamples << endl ;
      }
		}

		if (paramFlags[8] == 0) {
			cout << "Using default number of orthogonal functions (terms): " 
      << numFunctions 
			<< endl ;
    } else {
			cout << "Number of orthogonal functions (terms) set to: " 
      << numFunctions << endl ;
    }

		if (paramFlags[9] == 0) {
			cout << "Using data generated randomly ..." << endl ;
		} else {
			cout << "Using data from file: " << file << "..."  << endl ; 
    }
	}	

	if (paramFlags[4] == 0) {
	  cout << "Using default function to generate data: " << fname 
    << endl ;
	} else {
	  cout << "Function to generate data set to: " << fname << endl ;
    if (function == 2) {  // triangle wave
      if (paramFlags[13] == 1) {
        cout << "Parameter m in triangle wave specification set to: "
        << m_triangle << endl ;
      } else {
        cout << "Using default value of parameter m in triangle wave "
        "specification: " << m_triangle << endl ;
      }
    } else {
      if (paramFlags[13] == 1) {
        cout << "Ignoring parameter m value as it does not apply to " 
        << fname << "function ..." << endl ;
      }
    } 
  }
	
  if (function != 4) {		
    if (paramFlags[6] == 0) {
	    cout << "Using default value for maximum amplitude of wave: " << peak
      << endl ;
    } else {
		  cout << "Peak value set to: " << peak << endl ;
    }
  }

	if (paramFlags[2] == 0) {
		cout << "Using default value for interval's lower bound: " << low 
		<< endl ;
  } else {
  	cout << "Interval lower bound set to: " << low << endl ;
  }

	if (paramFlags[3] == 0) {
    cout << "Using default value for interval's higher bound: " << high 
		<< endl ;
  } else {
	  cout << "Interval higher bound set to: " << high << endl ;
  }

	if (paramFlags[5] == 0) {
	  cout << "Using default value for time period: " << timePeriod << endl ;
  } else {
	  cout << "Time period set to: " << timePeriod << endl ;
  }	

  if (paramFlags[10] == 0) {
    cout << "Running an instance of regression fit ..." << endl ;
  } else {
    cout << "Iterating over different values of parameters ..." << endl ; 
  }

  if (lambda > std::numeric_limits<long double>::epsilon()) {
    cout << "Using regularized least squares to compute weights ..." <<
    endl ;
  } else {
    cout << "Using normal least squares to compute weights ..." << endl ;
  }

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

  switch(basis) {
    case 0:
      cout << "Using sine/cosine as orthogonal basis functions ..." 
      << endl ;
      break ;
    case 1:
      cout << "Using Legendre Polynomials as orthogonal basis functions " 
      "..." << endl ;
      break ;
    default:
      error("Invalid choice of orthogonal basis set.") ;
      break ;
  }
  
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
  params.funcName = fname;
	params.file = file ;
  params.iterate = iterate ;
  params.inverse = inverse ;
  params.basis = basis ;
  params.m = m_triangle ;
  params.lambda = lambda ;
  params.printComponents = printComponents;

	return params ;
}

void setPrecision(void)
{
	cout.unsetf(ios::floatfield) ;
	int PRECISION = -log10(AOM) ;
	cout.precision(PRECISION) ;
	cout.setf(ios::fixed,ios::floatfield) ;
}

void plot (const char *file, int numSamples, long double noise, long double lambda, 
           string funcOutput, int printComponents)
{
  ofstream script ;
  script.open("temp/plotMsgLen.p") ;
  script << "set term post eps enhanced" << endl ;
  script << "set autoscale\t" ;
  script << "# scale axes automatically" << endl ;
  script << "set xtic auto\t" ;
  script << "# set xtics automatically" << endl ;
  script << "set ytic auto\t" ;
  script << "# set ytics automatically" << endl ;

  string n = convertToString<int>(numSamples) ;
  string s = convertToString<long double>(noise) ;
  string l = convertToString<long double>(lambda) ;
  string title;
  //title = "N = " + n + ", Sigma = " + s + ", Lambda = " + l ;
  script << "set title \"" << title << "\"" << endl ; 
  script << "set label \"" << funcOutput  << "\" at graph 0.005, graph 0.95 "
  "font \",10\"" << endl ;
  script << "set xlabel \"# of terms\"" << endl ;
  script << "set ylabel \"Message Length\"" << endl ;
  script << "set output \"./" << file << ".eps\"" << endl ;
  script << "plot \"./" << file << "\" using 1:5 notitle with lines lt 1 lc rgb "
  "\"blue\"" << endl ;

  script << "set output \"./" << file << ".comp.eps\"" << endl ;
  script << "set key left top" << endl;
  //if (printComponents) {
    script << "set multiplot" << endl ;
    script << "plot \"./" << file << "\" using 1:3 title 'part1' with lines lt 1 lc rgb "
    "\"red\", \\" << endl ;
    script << "\"./" << file << "\" using 1:4 title 'part2' with lines lt 1 lc rgb "
    "\"green\", \\" << endl ;
    script << "\"./" << file << "\" using 1:5 title 'total' with lines lt 1 lc rgb "
    "\"blue\"" << endl ;
  //} else {
    //script << "plot \"./" << file << "\" using 1:5 notitle with lines lt 1 lc rgb "
    //"\"blue\"" << endl ;
  //}  
  system ("gnuplot -persist temp/plotMsgLen.p") ;
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
	string filename,funcOutput ; 
  Components msglen;
  //plot("results_n1000_s0.2_m1.txt",1000,0.2,1,funcOutput,1) ;

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
			orthogonal = OrthogonalBasis(parameters.basis,parameters.numFunctions,
					parameters.timePeriod,parameters.function) ;
			phi = orthogonal.designMatrix(randomX) ;

			weights = computeWeights<long double>(phi,yValues,parameters.inverse,
                                            parameters.lambda) ;
      weights.print() ;

			predictions = dataGenerator.predict(parameters.numFunctions,weights,
                                          randomX) ;
			dataGenerator.plotPredictions(randomX,yValues,predictions) ;
	  	rmse = computeRMSE<long double>(weights,phi,yValues,parameters.lambda) ;
			msg = Message (parameters,weights,randomX,yValues,predictions) ;
			msglen = msg.messageLength() ;
			cout << "Msg Len = " << msglen.total << endl ;
			cout << "Error in fitting: " << rmse << endl ;
      break ;

    case 1:{
      int sampVals[] = {100} ;
      std::vector<int>Samples(sampVals,sampVals+sizeof(sampVals)/sizeof(int)) ;
      //long double noiseVals[] = {0} ;
      long double noiseVals[] = {0,0.1,0.2,0.3,0.4,0.5} ;
      std::vector<long double>Noise(noiseVals,noiseVals+sizeof(noiseVals)/sizeof(long double)) ;
	    for (unsigned i=0; i<Samples.size(); i++) 
      {
		    parameters.numSamples = Samples[i] ;
		    for (unsigned j=0; j<Noise.size(); j++)
		    {
          filename = "results/" + parameters.funcName + "/";
			    filename = filename + "results_n" + convertToString<int>(Samples[i]) + "_s" ;
			    filename = filename + convertToString<long double>(Noise[j]) + "_m";
			    filename = filename + convertToString<long double>(parameters.lambda) + 
                      ".txt" ;
			    ofstream results ;
			    results.open(filename.c_str()) ;	
			    parameters.sigma = Noise[j] ;

          dataGenerator = RandomDataGenerator<long double>(parameters) ;
			    dataGenerator.generate() ;
			    randomX = dataGenerator.randomX() ;
			    yValues = dataGenerator.yValues() ;
	        //dataGenerator.plotData() ;
	        //dataGenerator.plotDataWithNoise() ;
          funcOutput = dataGenerator.getFunctionString();

			    for (unsigned M=1; M<=200; M++) 
			    {
				    cout << "N: " << parameters.numSamples << "\t" ;
				    cout << "S: " << parameters.sigma << "\t" ;
				    cout << "M: " << M << endl ;
					  parameters.numFunctions = M ;
					  orthogonal = OrthogonalBasis (parameters.basis,
                          parameters.numFunctions,parameters.timePeriod,
                          parameters.function) ;
					  phi = orthogonal.designMatrix(randomX) ;
					  weights = computeWeights<long double>(phi,yValues,
                                parameters.inverse,parameters.lambda) ;
					  predictions = dataGenerator.predict(M,weights,randomX) ;

					  rmse = computeRMSE<long double> (weights,phi,yValues,
                                             parameters.lambda) ;

					  msg = Message (parameters,weights,randomX,yValues,
                              predictions) ;
					  msglen = msg.messageLength() ;
					  results << parameters.numFunctions << "\t" ;
					  results << rmse << "\t" ;
					  results << msglen.part1 << "\t";
					  results << msglen.part2 << "\t";
					  results << msglen.total << endl;
          }
          results.close() ;
          plot(filename.c_str(),Samples[i],Noise[j],parameters.lambda,funcOutput,
               parameters.printComponents) ;
        }
			}
		  break ;}

    default:
      error("Wrong choice entered ...") ;
      break ;
  }
	return 0 ;
}

