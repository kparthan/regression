#include "Gaussian.h"

using namespace std ;

class TESTING_Gaussian
{
	public:
		void TESTING_nullConstructor() ;
		void TESTING_explicitConstructor() ;
		void TESTING_functionValue() ;
} ;

void TESTING_Gaussian :: TESTING_nullConstructor (void)
{
	Gaussian d ;
	cout << d.mean() << " " << d.standardDeviation() << endl ;
}

void TESTING_Gaussian :: TESTING_explicitConstructor (void)
{
	Gaussian d ;
	cout << d.mean() << " " << d.standardDeviation() << endl ;
	d = Gaussian(3.6,7) ;
	cout << d.mean() << " " << d.standardDeviation() << endl ;
}

void TESTING_Gaussian :: TESTING_functionValue(void)
{
	Gaussian d ;
	d = Gaussian(3.6,7) ;
	cout << d.mean() << " " << d.standardDeviation() << endl ;
	cout << d.value(3.6) << endl ;
	cout << d.value(0) << endl ;
}

main()
{
	TESTING_Gaussian test ;

	cout << "Testing null constructor ..." << endl ;
	test.TESTING_nullConstructor() ;

	cout << "Testing explicit constructor ..." << endl ;
	test.TESTING_explicitConstructor() ;

	cout << "Testing function value ..." << endl ;
	test.TESTING_functionValue() ;

}

