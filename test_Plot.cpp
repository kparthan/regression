#include "Plot.h"

using namespace std ;

class TESTING_Plot
{
	public:
		void TESTING_nullConstructor() ;
		void TESTING_labeling() ;
		void TESTING_settingAxesRange() ;
		void TESTING_sketchRandom() ;
		void TESTING_sketchFunction() ;
} ;

void TESTING_Plot :: TESTING_nullConstructor (void)
{
	Plot p ;
}

void TESTING_Plot :: TESTING_labeling (void)
{
	Plot p ;
	vector<string> labels ;
	labels.push_back("Test plot") ;
  labels.push_back("Number") ;
  labels.push_back("X values") ;
	p.label(labels) ;
}

void TESTING_Plot :: TESTING_settingAxesRange(void)
{
	Plot p ;
	pair<int,int> xr ;
	pair<double,double> yr ;
	xr = make_pair (1,20) ;
	yr = make_pair (-5.6,8.88) ;
	p.setRange(xr,yr) ;
}

void TESTING_Plot :: TESTING_sketchRandom(void)
{
	int a[10] = {1,6,-9,0,5,9,-4,88,5,90} ;
	Data<int> d(a,10) ;
	d.print() ;
	Plot p ;
	p.sketch(d) ;
}

void TESTING_Plot :: TESTING_sketchFunction(void)
{
	int a[10],b[10],i ;
	for (i=0; i<10; i++)
	{
		a[i] = i ;
		b[i] = i * i ;
	}
	Data<int> d1(a,10) ;
	d1.print() ;
	Data<int> d2(b,10) ;
	d2.print() ;
	Plot p ;
	p.sketch(d1,d2) ;
}

main()
{
	TESTING_Plot test ;
		
	cout << "Testing null constructor ..." << endl ;
	test.TESTING_nullConstructor() ;

	cout << "Testing graph labeling ..." << endl ;
	test.TESTING_labeling() ;

	cout << "Testing setting axes range ..." << endl ;
	test.TESTING_settingAxesRange() ;

	cout << "Testing sketch random x values only ..." << endl ;
	test.TESTING_sketchRandom() ;

	cout << "Testing sketch function ..." << endl ;
	test.TESTING_sketchFunction() ;

}

