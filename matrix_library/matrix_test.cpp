#include "matrix.h"

main()
{
	int i ;
	Matrix<int> m(6) ;
	m.print() ;
	Matrix<int> mt(2,3,2) ;
	mt.print() ;
	Matrix<int> mt2(mt) ;
	mt2.print() ;
	cout << endl ;
	int c = 5 ;
	Matrix<int> mt1 = identity<int>(5);
	mt1.print() ;

	Matrix<int> mt3(4,3,2) ;
	Matrix<int> mt4 ;
	mt4 = mt3.add(mt2) ;
	mt4.print() ;
	//Matrix<int> mt5 ;
	Matrix<int> mt5 = mt4.transpose() ;
	mt5.print();
	Matrix<int> mt6 = mt2.product(mt5) ;
	mt6.print() ;
}
