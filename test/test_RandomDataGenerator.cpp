#include "RandomDataGenerator.h"

#define MEAN 0.0
#define SIGMA 1.0
#define LOW -1.0
#define HIGH 1.0
#define FUNCTION 0
#define TIME_PERIOD 0.5
#define PEAK 1.0
#define NUMSAMPLES 100

using namespace std ;

struct Parameters initialize (void)
{
	struct Parameters example ;
	example.mean = MEAN ;
	example.sigma = SIGMA ;
	example.low = LOW ;
	example.high = HIGH ;
	example.function = FUNCTION ;
	example.timePeriod = TIME_PERIOD ;
	example.peak = PEAK ;
	example.numSamples = NUMSAMPLES ;	
	return example ;
}

struct Parameters example = initialize() ;

class TESTING_RandomDataGenerator
{
	public:
		void TESTING_nullConstructor() ;
		void TESTING_generateData() ;
		void TESTING_plotRandomX() ;
		void TESTING_plotData() ;
		void TESTING_plotDataWithNoise() ;
} ;

void TESTING_RandomDataGenerator :: TESTING_nullConstructor (void)
{
	RandomDataGenerator<double> r(example) ;
	r.print() ;
}

void TESTING_RandomDataGenerator :: TESTING_generateData (void)
{
	RandomDataGenerator<double> r(example) ;
	r.generate() ;
}

void TESTING_RandomDataGenerator :: TESTING_plotRandomX (void)
{
	RandomDataGenerator<double> r(example) ;
	r.generate() ;
	r.plotRandomX() ;
}

void TESTING_RandomDataGenerator :: TESTING_plotData (void)
{
	RandomDataGenerator<double> r(example) ;
	r.generate() ;
	r.plotData() ;
}

void TESTING_RandomDataGenerator :: TESTING_plotDataWithNoise(void)
{
	RandomDataGenerator<double> r(example) ;
	r.generate() ;
	r.plotDataWithNoise() ;
}

main()
{
	TESTING_RandomDataGenerator test ;

	cout << "Testing null constructor ..." << endl ;
	test.TESTING_nullConstructor() ;		

	cout << "Testing random data generation method ..." << endl ;
	test.TESTING_generateData() ;		

	cout << "Testing plotting random Xs generated ..." << endl ;
	test.TESTING_plotRandomX() ;		
	
	cout << "Testing plotting function values (without noise) ..." << endl ;
	test.TESTING_plotData() ;		

	cout << "Testing plotting noise added values ..." << endl ;
	test.TESTING_plotDataWithNoise() ;		

}

