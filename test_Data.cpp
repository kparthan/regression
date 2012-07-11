#include "Data.h"
#include "Matrix.h"
#include <vector>

#define SIZE 10

using namespace std ;

class TESTING_Data
{
	public:
		void TESTING_nullConstructor() ;
		void TESTING_initializeToCppArray() ;
		void TESTING_copyConstructor() ;
		void TESTING_instantiateVecPoints() ;
		void TESTING_convertMatrixToData() ;
		void TESTING_overloadEqual() ;
		void TESTING_overloadPlus() ;
		void TESTING_overloadSqBracket() ;
		void TESTING_sortingElementsList() ;
		void TESTING_minimumElement() ;
		void TESTING_maximumElement() ;
		void TESTING_convertToMatrix() ;
} ;

void TESTING_Data :: TESTING_nullConstructor (void)
{
	Data<int> d ;	
	d.print() ;
}

void TESTING_Data :: TESTING_initializeToCppArray(void)
{
	int *a = new int[SIZE] ;
	for (int i=0; i<SIZE; i++)
		a[i] = 2 * i ;
	Data<int> d(a,SIZE) ;
	d.print() ;
	delete[] a ;
}

void TESTING_Data :: TESTING_copyConstructor(void)
{
	int *a = new int[SIZE] ;
	for (int i=0; i<SIZE; i++)
		a[i] = 2 * i ;
	Data<int> d(a,SIZE) ;
	d.print() ;
	delete[] a ;

	Data<int> d2(d) ;
	d2.print() ;
}

void TESTING_Data :: TESTING_instantiateVecPoints (void)
{
	Point<int> p1(10) ;
	Point<int> p2(20) ;
	Point<int> p3(60) ;
	Point<int> p4(80) ;
	Point<int> p5(-20) ;
	vector<Point<int> > vec ;
	vec.push_back(p1) ;
	vec.push_back(p2) ;
	vec.push_back(p3) ;
	vec.push_back(p4) ;
	vec.push_back(p5) ;
	Data<int> d(vec) ;
	d.print() ;
}

void TESTING_Data :: TESTING_convertMatrixToData(void)
{
	Matrix<int> mat1(3,1) ;
	mat1.print() ;
	mat1[0][0] = 100 ;
	mat1[1][0] = -98 ;
	mat1[2][0] = 5 ;
	mat1.print() ;
	Data<int> d(mat1) ;
	d.print() ;
}

void TESTING_Data :: TESTING_overloadEqual(void)
{
	int *a = new int[SIZE] ;
	for (int i=0; i<SIZE; i++)
		a[i] = 2 * i ;
	Data<int> d(a,SIZE) ;
	d.print() ;
	delete[] a ;

	Data<int> d2 ;
	d2.print() ;
	d2 = d ;
	d2.print() ;

	Data<int> d3 = d2 ;
	d3.print() ;
}

void TESTING_Data :: TESTING_overloadPlus(void)
{
	int *a = new int[SIZE] ;
	for (int i=0; i<SIZE; i++)
		a[i] = 2 * i ;
	Data<int> d1(a,SIZE) ;
	d1.print() ;
	delete[] a ;

	int b[10] = {1,-4,7,9,0,-3,8,5,-33,88} ;
	Data<int> d2(b,SIZE) ;
	d2.print() ;

	Data<int> d3 ;
	d3.print() ;
	d3 = d1 + d2 ;
	d3.print() ;

	Data<int> d4 = d1 + d3 + d2 ;
	d4.print() ;
}

void TESTING_Data :: TESTING_overloadSqBracket(void)
{
	int b[10] = {1,-4,7,9,0,-3,8,5,-33,88} ;
	Data<int> d2(b,SIZE) ;
	d2.print() ;

	cout << d2[0].x() << endl ;

	Point<int> p(10) ;
	d2[0] = p ;
	d2.print() ;
}

void TESTING_Data :: TESTING_sortingElementsList(void)
{
	int b[10] = {1,-4,1,9,0,-3,8,5,-33,88} ;
	Data<int> d2(b,SIZE) ;
	d2.print() ;
	
	d2.sortElements() ;
	d2.print() ;
}

void TESTING_Data :: TESTING_minimumElement(void)
{
	int b[10] = {1,-4,1,9,0,-3,8,5,-33,88} ;
	Data<int> d2(b,SIZE) ;
	d2.print() ;
	
	d2.sortElements() ;
	d2.print() ;
	cout << "Minimum = " << d2.minimum() << endl ;
}

void TESTING_Data :: TESTING_maximumElement(void)
{
	int b[10] = {1,-4,1,9,0,-3,8,5,-33,88} ;
	Data<int> d2(b,SIZE) ;
	d2.print() ;
	
	d2.sortElements() ;
	d2.print() ;
	cout << "Maximum = " << d2.maximum() << endl ;
}

void TESTING_Data :: TESTING_convertToMatrix(void)
{
	int b[10] = {1,-4,1,9,0,-3,8,5,-33,88} ;
	Data<int> d2(b,SIZE) ;
	d2.print() ;
	
	Matrix<int> mat1 = d2.convertToMatrix() ;
	mat1.print() ;
	
	Data<int> d ;	
	d.print() ;
	Matrix<int> mat2 = d.convertToMatrix() ;
	mat2.print() ;
}

main()
{
	TESTING_Data test ;

	cout << "Testing null constructor ..." << endl ;
	test.TESTING_nullConstructor() ;		

	cout << "Testing initialization to C/C++ array ..." << endl ;
	test.TESTING_initializeToCppArray() ;		

	cout << "Testing copy constructor ..." << endl ;
	test.TESTING_copyConstructor() ;		

	cout << "Testing constructor with vector of Points ..." << endl ;
	test.TESTING_instantiateVecPoints() ;		

	cout << "Testing convert Matrix to Data ..." << endl ;
	test.TESTING_convertMatrixToData() ;		

	cout << "Testing overload = ..." << endl ;
	test.TESTING_overloadEqual() ;		

	cout << "Testing overload + ..." << endl ;
	test.TESTING_overloadPlus() ;		

	cout << "Testing overload [] ..." << endl ;
	test.TESTING_overloadSqBracket() ;		

	cout << "Testing sorting data elements ..." << endl ;
	test.TESTING_sortingElementsList() ;		

	cout << "Testing minimum ..." << endl ;
	test.TESTING_minimumElement() ;		

	cout << "Testing maximum ..." << endl ;
	test.TESTING_maximumElement() ;		

	cout << "Testing convertToMatrix() method ..." << endl ;
	test.TESTING_convertToMatrix() ;		

}

