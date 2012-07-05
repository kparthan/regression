/************************************************************************
This program generates random data values and the corresponding function values with added noise.
************************************************************************/

#include <iostream>
#include <stdlib.h>
#include <time.h>

#include <liblcb/Vector.h>
using namespace std ;

class randomNumberGenerator
{
	private: 
		int numDataPoints ;
		float xmin,xmax ;
		float *dataPoints ;
	public:
		void setValues(int,float,float) ;
		float *generate(void) ;
} ;

void randomNumberGenerator::setValues(int num, float min, float max)
{
	srand((unsigned)(time(0))) ;
	numDataPoints = num ;
	xmin = min;
	xmax = max ;
	dataPoints = new float [numDataPoints] ;
}

float *randomNumberGenerator::generate(void)
{
	int i ;
	for (i=0; i<numDataPoints; i++)
		dataPoints[i] = xmin + (xmax-xmin) * (rand() / (RAND_MAX + 1.0)) ;
	return dataPoints ;
}

int main(int argc, char **argv)
{
	int i,numDataPoints ;
	float xmin,xmax ;

	/* Enter data parameters */
	cout << "Enter number of points to generate: " ;
	cin >> numDataPoints ;
	cout << "Enter minimum x value: " ;
	cin >> xmin ;
	cout << "Enter maximum x vaue: " ;
	cin >> xmax ;

	randomNumberGenerator *generator = new randomNumberGenerator() ;
	generator->setValues(numDataPoints,xmin,xmax) ;	
	float *dataPoints = generator->generate() ;
	for (i=0;i<numDataPoints;i++)
		cout << dataPoints[i] << " " ;
	cout << endl ;
//	float *functionValues = 

	return 0 ;
}
