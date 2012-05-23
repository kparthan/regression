#include "myVector.h"

#define SIZE 3

int main(int argc, char **argv)
{
	int c = 4,pos ;
	int i,*a = new int[10] ;
	for(i=0; i<SIZE; i++)
		a[i] = 2 * i ;
	MyVector<int> x(a,SIZE) ;
	x.print() ;
	MyVector<int> y,z ;
	y.print() ;
	y = x ;
	y.print() ;
	y = c ;
	y.print() ;
	z = c ;
	z.print() ;
	pos = y[2] ; cout << pos << endl ;
	cout << "y length = " << y.l2Norm() << endl ;
	cout << "dot product (x,y) = " << y.dotProduct(x) << endl ;
	/*int i,*a = new int[10] ;
	for(i=0; i<SIZE; i++)
		a[i] = 2 * i ;
	MyVector<int> x (a,SIZE) ;
	x.print() ;
	
	MyVector<int> y  ;
	y.print() ;
	y = x ;
	y.print() ;
	
	int c = 4 ;
	MyVector <int> z(SIZE)  ;
	z = c ;
	z.print() ;

	cout << z[8]*5 << endl ;*/
	return 0 ;
}
