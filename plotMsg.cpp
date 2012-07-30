#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

using namespace std ;
template <class T>

string convertToString(T number)
{
  ostringstream convert ;
  convert << number ;
  return convert.str() ;
}

void plot (const char *file, int numSamples, double noise)
{
	ofstream script ;
	script.open("plotMsgLen.p") ;
	script << "set term post eps" << endl ;
  script << "set autoscale\t" ;
  script << "# scale axes automatically" << endl ;
  script << "set xtic auto\t" ;
  script << "# set xtics automatically" << endl ;
  script << "set ytic auto\t" ;
  script << "# set ytics automatically" << endl ;

	string n = convertToString<int>(numSamples) ;
	string s = convertToString<double>(noise) ;
	string title = "N = " + n + ", Sigma = " + s ;
	script << "set title \"" << title << "\"" << endl ; 
	script << "set xlabel \"# of terms\"" << endl ;
	script << "set ylabel \"Message Length\"" << endl ;
	script << "set output \"./" << file << ".eps\"" << endl ;

	script << "plot \"./" << file << "\" using 1:3 notitle with linespoints lc rgb \"blue\"" << endl ;
	system ("gnuplot -persist plotMsgLen.p") ;
}

main()
{			
	string file ;
  int Samples[1] = {100} ;
  double Noise[1] = {0} ;
  //double Noise[1] = {0.25} ;
	
	for (int i=0; i<1; i++)
	{
		for (int j=0; j<1; j++)
		{
			file = "Results/results_n" + convertToString<int>(Samples[i]) + "_s" ;
			//file = "Results_Pi_extra/results_n" + convertToString<int>(Samples[i]) + "_s" ;
      file = file + convertToString<double>(Noise[j]) + ".txt" ;
			plot(file.c_str(),Samples[i],Noise[j]) ;
		}
	}
}

