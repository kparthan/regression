/*!
 *  \file myVector.h
 *  \details Implementation of MyVector class
 *  \author Parthan Kasarapu
 *  \version 1.0
 *  \date Modified: Tue 29 May 2012
*/
#ifndef _MYVECTOR_H_
#define _MYVECTOR_H_

#include <iostream>
#include <vector>
#include <cmath>
#include "error.h"

using namespace std ;

/*! 
 *  \class MyVector myVector.h "matrix_library/myVector.h"
 *  \brief This is a MyVector class.
 *
 *  The class acts as an interface to create a n-length vector and store elements in it. Also provides basic functionalities such as norm of the vector, dot product between two vectors.
 */
template <class T>
class MyVector
{
	private:
		int size ;
		vector<T> v ;
	public:
				/* Constructors  */
		//! creates a null vector
		MyVector () ;
		//! Zero based vector					
		explicit MyVector (int) ;			
		//! Initialize to constant value
		MyVector (const T &, int) ;			
		//! Initialize to C/C++ style array
		MyVector (const T *, int) ;			
		//! Copy constructor
		MyVector (const MyVector &) ;			

				/* Overloading = and [] */
		//! assignment to a sourceVector
		MyVector & operator = (const MyVector &) ;	
		//! assignment to a constant
		MyVector & operator = (const T &) ;		
		//! return i^{th} element
		inline T & operator [] (const int ) ;		
		//! computes the dot product between vectors
		float operator * (MyVector &) ;

				/* Other sub-routines */
		//! print the elements of the vector
		void print (void) ;				
		//! return size of vector
		inline int length(void) const ;	
		//! changes the element at position i
		inline void modifyElement(const int, T &) ;
		//! computes the L2-Norm (magnitude) of the vector
		float l2Norm () ;				
		//! destructor
		~MyVector() ;					
		
		/*	OBSOLETE
		// return i^{th} element
		inline T element(const int &) ;	*/
		/*	OBSOLETE - REPLACED BY * OVERLOADING
		// computes the dot product of two vectors
		float dotProduct (const MyVector &) ; */
} ;

/*! 
 *  Null vector Constructor.
 *  \return A new instance of a vector.
 */
template <class T>
MyVector<T> :: MyVector() : size(0), v(size)
{
}

/*! 
 *  \fn MyVector<T> :: MyVector(int n)
 *  \brief This constructor function creates a zero vector of size n.
 *  \param n an integer.
 *  \return A new instance of a zero vector.
 */
template <class T>
MyVector<T> :: MyVector(int n) : size(n)
{
	v.assign(size,0) ;
}

/*! 
 *  \fn MyVector<T> :: MyVector(const T &a, int n)
 *  \brief This constructor function creates a vector of size n with all elements equal.
 *  \param a reference to the constant value
 *  \param n an integer.
 *  \return A new instance of a vector with same elements.
 */
template <class T>
MyVector<T> :: MyVector(const T &a, int n) : size(n)
{
	v.assign(size,a) ;
}

/*! 
 *  \fn MyVector<T> :: MyVector(const T *a, int n) 
 *  \brief This constructor function creates a vector of size n and initializes to a C/C++ style array.
 *  \param a pointer to a C/C++ array
 *  \param n an integer.
 *  \return A new instance of a vector with elements as that of the C/C++ array.
 */
template <class T>
MyVector<T> :: MyVector(const T *a, int n) : size(n)
{
	for (int i=0; i<size; i++)
		v.push_back(a[i]) ;
}

/*! 
 *  \fn MyVector<T> :: MyVector(const MyVector<T> &sourceVector)
 *  \brief This constructor function creates a vector and instantiates it with an existing vector.
 *  \param sourceVector reference to a MyVector object
 *  \return A copy of the source MyVector object.
 */
template <class T>
MyVector<T> :: MyVector(const MyVector<T> &sourceVector) : size(sourceVector.size), v(sourceVector.v)
{
}

/*! 
 *  \fn MyVector<T> & MyVector<T> :: operator = (const MyVector<T> &sourceVector)
 *  \brief The assignment operator is overloaded to assign the sourceVector to the current MyVector object.
 *  \param sourceVector a reference to a MyVector object.
 *  \return An instantiated vector which is a copy of the original sourceVector.
 */
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

/*! 
 *  \fn MyVector<T> & MyVector<T> :: operator = (const T &a)
 *  \brief The assignment operator is overloaded to reassign a vector with equal elements.
 *  \param a a reference to a constant value
 *  \return A modified vector whose all elements are equal.
 */
template <class T>
MyVector<T> & MyVector<T> :: operator = (const T &a)
{
	v.clear() ;	// erase the contents before copying 
	for (int i=0; i<size; i++)
		v.push_back(a) ;
	return *this ;
}

/*! 
 *  \fn inline T & MyVector<T> :: operator [] (const int i)
 *  \brief The [] operator is overloaded which returns the element at position i.
 *  \param i an integer index 
 *  \return The ith vector element.
 */
template <class T>
inline T & MyVector<T> :: operator [] (const int i)
{
	
	if (i >= size)
		error ("In accessing vector elements: index out of range ...") ;
	else return v[i] ;	// alternate vector access
}

/*! 
 *  \fn double MyVector<T> :: operator * (MyVector &vec)
 *  \brief The * operator is overloaded to compute the dot product.
 *  \param vec reference to a MyVector object
 *  \return The dot product of two vectors.
 */
template <class T>
float MyVector<T> :: operator * (MyVector<T> &vec)
{
	float dotProduct = 0 ;
	if (size != vec.size)
		error("In computing dot product: Vector sizes don't match!") ;
	else
	{
		for (int i=0; i<size; i++)
			dotProduct += v[i] * vec[i] ;
	}
	return dotProduct ;
}

/*! 
 *  \fn void MyVector<T> :: print(void)
 *  \brief The method prints the elements in the vector.
 */
template <class T>
void MyVector<T> :: print(void)
{
	cout << size << endl ;
	for (int i=0; i<size; i++)
		cout << v.at(i) << " " ;
	cout << endl ;
}

/*! 
 *  \fn inline int MyVector<T> :: length() const
 *  \brief The method returns the size of the vector.
 *  \return vector size
 */
template <class T>
inline int MyVector<T> :: length() const
{
	return size ;
}

/*! 
 *  \fn inline void MyVector<T> :: modifyElement (const int i)
 *  \brief The method modifies the element at position i.
 *  \param i an integer
 *  \param value reference to a value
 */
template <class T>
inline void MyVector<T> :: modifyElement (const int i, T &value)
{
	if (i >= size)
		error ("In modifying element value: Index out of range ...") ;
	else
		v[i] = value ;
}

/*! 
 *  \fn float MyVector<T> :: l2Norm(void)
 *  \brief The method computes the magnitude of the vector. 
 *  \return vector length (L2-norm)
 */
template <class T>
float MyVector<T> :: l2Norm(void)
{
	float magnitude = 0 ;
	for (int i=0; i<size; i++)
		magnitude += v[i] * v[i] ;
	return sqrt(magnitude) ;
}

/*! 
 *  \fn MyVector<T> :: ~MyVector()
 *  \brief Destructor function of MyVector class
 */
template <class T>
MyVector<T> :: ~MyVector()
{
	v.clear() ;
}

#endif

/* 		OBSOLETE
 *  \fn inline T MyVector<T> :: element (const int &index)
 *  \brief The method returns the element at position 'index'.
 *  \param index an integer reference
 *  \return vector element
 */
/*template <class T>
inline T MyVector<T> :: element (const int &index)
{
	if (index >= size)
		error ("In accessing vector elements: index out of range ...") ;
	else return v[index] ;
}*/

/*	OBSOLETE - REPLACED BY * OVERLOADING
 *  \fn float MyVector<T> :: dotProduct (const MyVector<T> &other)
 *  \brief The method computes the dot product with another vector. 
 *  \param other a reference to another vector
 *  \return dot product
 */
/*template <class T>
float MyVector<T> :: dotProduct (const MyVector<T> &other)
{
	float dotProduct = 0 ;
	if (size != other.size)
		error("In computing dot product: Vector sizes don't match!") ;
	else
	{
		for (int i=0; i<size; i++)
			//dotProduct += v.at(i) * other.v.at(i) ;
			dotProduct += v[i] * other.v[i] ;
	}
	return dotProduct ;
}*/
