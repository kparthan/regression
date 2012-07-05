#include "Point.h"

using namespace std ;

class TESTING_Point
{
	public:
		void TESTING_nullPoint() ;
		void TESTING_constructor() ;
		void TESTING_copyConstructor() ;
		void TESTING_overloadEqual() ;
		void TESTING_overloadPlusEqual() ;
		void TESTING_overloadPlus() ;
		void TESTING_overloadMinusEqual() ;
		void TESTING_overloadMinus() ;
		void TESTING_overloadLess() ;
} ;

void TESTING_Point :: TESTING_nullPoint(void)
{
	Point<int> p ;
	cout << p.x() << endl ;
}

void TESTING_Point :: TESTING_constructor (void)
{
	Point<int> p (5) ;
	cout << p.x() << endl ;
	Point<float> q(5.6) ;
	cout << q.x() << endl ;
}

void TESTING_Point :: TESTING_copyConstructor (void)
{
	Point<int> p (100) ;
	cout << p.x() << endl ;
	Point<int> q(p) ;
	cout << q.x() << endl ;
}

void TESTING_Point :: TESTING_overloadEqual (void)
{
	Point<double> p,q(10) ;
	cout << p.x() << endl ;
	cout << q.x() << endl ;
	p = q ;
	cout << p.x() << endl ;
	Point<double> r = q ;
	cout << r.x() << endl ;
}

void TESTING_Point :: TESTING_overloadPlusEqual (void)
{
	Point<int> p(-78) ;
	cout << p.x() << endl ;
	Point<int> q(55) ;
	cout << q.x() << endl ;
	q += p ;
	cout << q.x() << endl ;
}

void TESTING_Point :: TESTING_overloadPlus (void)
{
	Point<int> p(10) ;
	cout << p.x() << endl ;
	Point<int> q(55) ;
	cout << q.x() << endl ;
	Point<int> r = p + q ;
	cout << r.x() << endl ;
}

void TESTING_Point :: TESTING_overloadMinusEqual (void)
{
	Point<int> p(-78) ;
	cout << p.x() << endl ;
	Point<int> q(55) ;
	cout << q.x() << endl ;
	q -= p ;
	cout << q.x() << endl ;
}

void TESTING_Point :: TESTING_overloadMinus (void)
{
	Point<int> p(10) ;
	cout << p.x() << endl ;
	Point<int> q(55) ;
	cout << q.x() << endl ;
	Point<int> r = p - q ;
	cout << r.x() << endl ;
}

void TESTING_Point :: TESTING_overloadLess (void)
{
	Point<int> p(1001) ;
	cout << p.x() << endl ;
	Point<int> q(155) ;
	cout << q.x() << endl ;
	//cout << p < q << endl ;
	bool x  ;
	//if (p < q) x = 1 ;
	x = p < q ;
	cout << (p<q) << endl ;
}

main()
{
	TESTING_Point test ;

	cout << "Testing null constructor ..." << endl ;
	test.TESTING_nullPoint() ;

	cout << "Testing constructor ..." << endl ;
	test.TESTING_constructor() ;

	cout << "Testing copy constructor ..." << endl ;
	test.TESTING_copyConstructor() ;

	cout << "Testing overload = operator ..." << endl ;
	test.TESTING_overloadEqual() ;

	cout << "Testing overload += operator ..." << endl ;
	test.TESTING_overloadPlusEqual() ;

	cout << "Testing overload + operator ..." << endl ;
	test.TESTING_overloadPlus() ;

	cout << "Testing overload -= operator ..." << endl ;
	test.TESTING_overloadMinusEqual() ;

	cout << "Testing overload - operator ..." << endl ;
	test.TESTING_overloadMinus() ;

	cout << "Testing overload < operator ..." << endl ;
	test.TESTING_overloadLess() ;

}

