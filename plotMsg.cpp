#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std ;
template <class T>

string convertToString(T number)
{
  ostringstream convert ;
  convert << number ;
  return convert.str() ;
}

void plot (const char *file, int numSamples, long double noise)
{
	ofstream script ;
	script.open("temp/plotMsgLen.p") ;
	script << "set term post eps" << endl ;
  script << "set autoscale\t" ;
  script << "# scale axes automatically" << endl ;
  script << "set xtic auto\t" ;
  script << "# set xtics automatically" << endl ;
  script << "set ytic auto\t" ;
  script << "# set ytics automatically" << endl ;

	string n = convertToString<int>(numSamples) ;
	string s = convertToString<long double>(noise) ;
	string title = "N = " + n + ", Sigma = " + s ;
	script << "set title \"" << title << "\"" << endl ; 
	script << "set xlabel \"# of terms\"" << endl ;
	script << "set ylabel \"Message Length\"" << endl ;
	script << "set output \"./" << file << ".eps\"" << endl ;

	script << "plot \"./" << file << "\" using 1:3 notitle with linespoints lc rgb \"blue\"" << endl ;
	system ("gnuplot -persist temp/plotMsgLen.p") ;
}

main()
{			
	string file ;
<<<<<<< HEAD
  int sampVals[] = {100} ;
  std::vector<int> Samples (sampVals,sampVals+sizeof(sampVals)/sizeof(int)) ;
  //long double noiseVals[] = {0.1,0.2,0.3,0.4,0.5} ;
  long double noiseVals[] = {0.5} ;
  std::vector<long double> Noise (noiseVals,noiseVals+sizeof(noiseVals)/sizeof(long double)) ;
=======
  int Samples[1] = {1000} ;
  long double Noise[1] = {0.75} ;
>>>>>>> 57d61fedf3014331915abae9f24ac44ddf0b747b
	
	for (int i=0; i<Samples.size(); i++)
	{
		for (int j=0; j<Noise.size(); j++)
		{
<<<<<<< HEAD
			file = "temp/results_n" + convertToString<int>(Samples[i]) + "_s" ;
=======
			file = "results/results_n" + convertToString<int>(Samples[i]) + "_s" ;
>>>>>>> 57d61fedf3014331915abae9f24ac44ddf0b747b
      file = file + convertToString<long double>(Noise[j]) + ".txt" ;
			plot(file.c_str(),Samples[i],Noise[j]) ;
		}
	}
}

