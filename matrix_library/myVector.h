/**********************************************************************************************************
The vector class provides routines with the following functionalities:
	
Constructors:-
	1. Create a NULL vector	
		MyVector<type> A ;
	2. Create a vector of zeroes 
		MyVector<type> A (SIZE) ;
	3. Create a vector with equal elements	
		MyVector<type> A (constant,SIZE) ;
	4. Create a vector whose elements are initialized to that of an array
		MyVector<type> A (Array,SIZE) ;
	5. Copy constructor
		MyVector<type> A ;
		MyVector<type> B (A) ;

Operator overloading:-
	1. operator =			// copies the contents of an object
		MyVector<type> A,B ;
		B = A ;
	2. operator =			// creates a vector with all elements equal to a constant 'c'
		MyVector<type> A ;
		type c ;
		A = c ;
	3. operator []			// subscripting
		MyVector<type> A ;
		type posVal = A[i] ;

Functions:-
	   MyVector<type> A,B ;
	1. A.print()
		print the elements in the vector
	2. A.length()	
		returns the vector size
	3. A.element(i)
		returns the elements at position i
	4. A.l2Norm()
		returns the magnitude/L2-Norm of the vector
	5. float c = A.dotProduct(B)
		returns the dot product between the vectors A and B

**********************************************************************************************************/

#ifndef _MYVECTOR_H_
#define _MYVECTOR_H_

#include <iostream>
#include <vector>
#include <cmath>
#include "error.h"

using namespace std ;

template <class T>
class MyVector
{
	private:
		int size ;
		vector<T> v ;
	public:
		/* Constructors  */
		MyVector () ;					// creates a null vector
		explicit MyVector (int) ;			// Zero based vector
		MyVector (const T &, int) ;			// Initialize to constant value
		MyVector (const T *, int) ;			// Initialize to C/C++ style array
		MyVector (const MyVector &) ;			// Copy constructor

		/* Overloading = and [] */
		MyVector & operator = (const MyVector &) ;	// assignment to a sourceVector
		MyVector & operator = (const T &) ;		// assignment to a constant
		//inline T & operator [] (const int ) ;		// return i^{th} element

		/* Other sub-routines */
		void print (void) ;				// print the elements of the vector
		inline int length(void) const ;			// return size of vector
		inline T element(const int &) ;			// return i^{th} element
		float l2Norm () ;				// computes the L2-Norm (magnitude) of the vector
		float dotProduct (const MyVector<T> &) ;	// computes the dot product of two vectors

		~MyVector() ;					// destructor
} ;

template <class T>
MyVector<T> :: MyVector() : size(0), v(size)
{
}

template <class T>
MyVector<T> :: MyVector(int n) : size(n)
{
	v.assign(size,0) ;
}

template <class T>
MyVector<T> :: MyVector(const T &a, int n) : size(n)
{
	v.assign(size,a) ;
}

template <class T>
MyVector<T> :: MyVector(const T *a, int n) : size(n)
{
	for (int i=0; i<size; i++)
		v.push_back(a[i]) ;
}

template <class T>
MyVector<T> :: MyVector(const MyVector<T> &sourceVector) : size(sourceVector.size), v(sourceVector.v)
{
}

template <class T>
MyVector<T> & MyVector<T> :: operator = (const MyVector<T> &sourceVector)
{
	if (this != &sourceVector)
	{
		size = sourceVector.size ;
		v = sourceVector.v ;	// the contents are automatically erased
	}
	return *this ;
}

template <class T>
MyVector<T> & MyVector<T> :: operator = (const T &a)
{
	v.clear() ;	// erase the contents before copying 
	for (int i=0; i<size; i++)
		v.push_back(a) ;
	return *this ;
}

/*
template <class T>
inline T & MyVector<T> :: operator [] (const int i)
{
	return v.at(i) ;
}
*/

template <class T>
void MyVector<T> :: print(void)
{
	cout << size << endl ;
	for (int i=0; i<size; i++)
		cout << v.at(i) << " " ;
	cout << endl ;
}

template <class T>
inline int MyVector<T> :: length() const
{
	return size ;
}

template <class T>
inline T MyVector<T> :: element (const int &index)
{
	if (index >= size)
		error ("In accessing vector elements: index out of range ...") ;
	else return v[index] ;
}

template <class T>
float MyVector<T> :: l2Norm(void)
{
	float magnitude = 0 ;
	/*vector<T> :: iterator it ;
	for (it=v.begin(); it<v.end(); it++)
		magnitude += *it * *it ;*/
	for (int i=0; i<size; i++)
		magnitude += v.at(i) * v.at(i) ;
	return sqrt(magnitude) ;
}

template <class T>
float MyVector<T> :: dotProduct (const MyVector<T> &other)
{
	float dotProduct = 0 ;
	if (size != other.size)
		error("In computing dot product: Vector sizes don't match!") ;
	else
	{
		for (int i=0; i<size; i++)
			dotProduct += v.at(i) * other.v.at(i) ;
	}
	return dotProduct ;
}

template <class T>
MyVector<T> :: ~MyVector()
{
	v.clear() ;
}

#endif

