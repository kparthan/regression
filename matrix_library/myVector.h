/*!
 *  \file myVector.h
 *  \details Implementation of MyVector class
 *  \author Parthan Kasarapu
 *  \version 1.0
 *  \date Modified: Fri 8 Jun 2012
*/
#ifndef _MYVECTOR_H_
#define _MYVECTOR_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <cassert>
#include "error.h"

#define MINIMUM numeric_limits<double>::epsilon()

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
		MyVector (T , int) ;			
		//! Initialize to C/C++ style array
		MyVector (T *, int) ;			
		//! Copy constructor
		MyVector (const MyVector &) ;			

				/* Overloading = [] + += - -= * operators */
		//! assignment to a sourceVector
		MyVector operator = (MyVector) ;	
		//! assignment to a constant
		MyVector operator = (T) ;		
		//! return i^{th} element
		T & operator [] (int) ;		
		//! computes the dot product between vectors
		float operator * (MyVector &) ;
		//! computes the sum of two vectors
		MyVector operator + (MyVector) ;
		//! adds a vector to the existing vector 
		MyVector & operator += (MyVector &) ;
		//! computes the difference between two vectors
		MyVector operator - (MyVector) ;
		//! subtracts a vector from an existing one 
		MyVector & operator -= (MyVector &) ;
		//! scaling a vector by a constant value
		MyVector operator * (T) ;

				/* Other sub-routines */
		//! print the elements of the vector
		void print (void)  ;				
		//! return size of vector
		 int length(void)  ;	
		//! computes the L2-Norm (magnitude) of the vector
		float l2Norm () ;
		//! normalizes the vector
		MyVector<T> normalize () ;
		
				/* destructor */
		//~MyVector() ;					
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
 *  \fn MyVector<T> :: MyVector(T a, int n)
 *  \brief This constructor function creates a vector of size n with all elements equal.
 *  \param a a constant value of type T
 *  \param n an integer.
 *  \return A new instance of a vector with same elements.
 */
template <class T>
MyVector<T> :: MyVector(T a, int n) : size(n)
{
	v.assign(size,a) ;
}

/*! 
 *  \fn MyVector<T> :: MyVector(T *a, int n) 
 *  \brief This constructor function creates a vector of size n and initializes to a C/C++ style array.
 *  \param a pointer to a C/C++ array
 *  \param n an integer.
 *  \return A new instance of a vector with elements as that of the C/C++ array.
 */
template <class T>
MyVector<T> :: MyVector(T *a, int n) : size(n)
{
	for (int i=0; i<size; i++)
		v.push_back(a[i]) ;
}

/*! 
 *  \fn MyVector<T> :: MyVector(const MyVector<T> &sourceVector)
 *  \brief This is a copy constructor function which creates a vector and instantiates it with an existing vector.
 *  \param sourceVector reference to a MyVector object of type T
 *  \return A copy of the source MyVector object.
 */
template <class T>
MyVector<T> :: MyVector(const MyVector<T> &sourceVector) : size(sourceVector.size), v(sourceVector.v)
{
}

/*! 
 *  \fn MyVector<T> MyVector<T> :: operator = (MyVector<T> sourceVector)
 *  \brief The assignment operator is overloaded to assign the sourceVector to the current MyVector object.
 *  \param sourceVector a MyVector object of type T
 *  \return An instantiated vector which is a copy of the original sourceVector.
 */
template <class T>
MyVector<T> MyVector<T> :: operator = (MyVector<T> sourceVector)
{
	if (this != &sourceVector)
	{
		size = sourceVector.size ;
		v = sourceVector.v ;	// the contents are automatically erased
	}
	return *this ;
}

/*! 
 *  \fn MyVector<T> MyVector<T> :: operator = (T a)
 *  \brief The assignment operator is overloaded to reassign a vector with equal elements.
 *  \param a a reference to a constant value
 *  \return A modified vector whose all elements are equal.
 */
template <class T>
MyVector<T> MyVector<T> :: operator = (T a)
{
	v.clear() ;	// erase the contents before copying 
	for (int i=0; i<size; i++)
		v.push_back(a) ;
	return *this ;
}

/*! 
 *  \fn  T & MyVector<T> :: operator [] ( int i)
 *  \brief The [] operator is overloaded which returns the element at position i.
 *  \param i an integer index 
 *  \return The ith vector element.
 */
template <class T>
T & MyVector<T> :: operator [] (int i) 
{
	if (i >= size)
		error ("In accessing vector elements: index out of range ...") ;
	else return v[i] ;	// alternate vector access
}

/*! 
 *  \fn double MyVector<T> :: operator * (MyVector &vec)
 *  \brief The * operator is overloaded to compute the dot product.
 *  \param vec reference to a MyVector object of type T
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
 *  \fn MyVector<T> MyVector<T> :: operator + (MyVector<T> other)
 *  \brief The + operator is overloaded to compute the sum of two vectors.
 *  \param other a MyVector object of type T
 *  \return The vector sum.
 */
template <class T>
MyVector<T> MyVector<T> :: operator + (MyVector<T> other)
{
	if (size != other.size)
		error("In computing vector sum: dimensions mismatch!") ;
	else
		return MyVector<T>(*this)+=other ;
}

/*! 
 *  \fn MyVector<T> & MyVector<T> :: operator += (MyVector<T> &other)
 *  \brief The += operator is overloaded to add a vector to an existing one.
 *  \param other reference to a MyVector object of type T
 *  \return A reference to the original vector.
 */
template <class T>
MyVector<T> & MyVector<T> :: operator += (MyVector<T> &other)
{
	// check if the vector is empty/null
	// if so, resize to fit the dimensions
	if ((*this).size == 0)
	{
		size = other.size ;
		v.resize(size) ;
	}
	for (int i=0; i<size; i++)
		v[i] += other[i] ;
	return *this ;
}
 
/*! 
 *  \fn MyVector<T> MyVector<T> :: operator - (MyVector other)
 *  \brief The - operator is overloaded to compute the difference of two vectors.
 *  \param other a MyVector object
 *  \return The vector difference.
 */
template <class T>
MyVector<T> MyVector<T> :: operator - (MyVector<T> other)
{
	if (size != other.size)
		error("In computing vector sum: dimensions mismatch!") ;
	else
		return MyVector<T>(*this)-=other ;
}
 
/*! 
 *  \fn MyVector<T> & MyVector<T> :: operator -= (MyVector<T> &other)
 *  \brief The -= operator is overloaded to subtract a vector from an existing one.
 *  \param other reference to a MyVector object of type T
 *  \return A reference to the original vector.
 */
template <class T>
MyVector<T> & MyVector<T> :: operator -= (MyVector<T> &other)
{
	// check if the vector is empty/null
	// if so, resize to fit the dimensions
	if ((*this).size == 0)
	{
		size = other.size ;
		v.resize(size) ;
	}
	for (int i=0; i<size; i++)
		v[i] = v[i] - other[i] ;
	return *this ;
}

/*! 
 *  \fn MyVector<T> MyVector<T> :: operator * (T a)
 *  \brief The * operator is overloaded to scale the vector by a ant value.
 *  \param a constant value of type T
 *  \return The scaled vector.
 */
template <class T>
MyVector<T> MyVector<T> :: operator * (T a)
{
	MyVector<T> result (size) ;
	for (int i=0; i<size; i++)
		result[i] = v[i] * a ;
	return result ;
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
 *  \fn  int MyVector<T> :: length() 
 *  \brief The method returns the size of the vector.
 *  \return vector size
 */
template <class T>
 int MyVector<T> :: length() 
{
	return size ;
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
 *  \fn MyVector<T> MyVector<T> :: normalize(void)
 *  \brief The method normalizes the vector. 
 *  \return Normalized unit vector.
 */
template <class T>
MyVector<T> MyVector<T> :: normalize (void)
{
	if (size == 0)
		error ("In normalize() method: trying to normalize a zero vector!") ;
	else
	{
		MyVector<T> result (*this) ;
		float vectorMagnitude = result.l2Norm() ;
		for (int i=0; i<size; i++)
			result[i] = result[i] / vectorMagnitude ;
		return result ;
	}
}

/*! 
 *  \relates MyVector
 *  \brief This function computes the angle between two vectors.
 *  \param vec1 reference to a MyVector object of type T 
 *  \param vec2 reference to a MyVector object of type T 
 *  \return Angle between the vectors.
 */
template <class T>
float angleBetween (MyVector<T> &vec1, MyVector<T> &vec2)
{
	float dotProduct = vec1 * vec2 ;
	float vector1Magnitude = vec1.l2Norm() ;
	float vector2Magnitude = vec2.l2Norm() ;
	return acos(dotProduct / (vector1Magnitude * vector2Magnitude)) * 180 / M_PI ;
}

/*! 
 *  \relates MyVector
 *  \brief This function determines whether two given vectors are parallel or not 
 *  \param vec1 reference to a MyVector object of type T 
 *  \param vec2 reference to a MyVector object of type T 
 *  \return True if parallel; False if not parallel.
 */
template <class T>
bool isParallel (MyVector<T> &vec1, MyVector<T> &vec2)
{
	float angle = angleBetween (vec1,vec2) ;
	cout << angle << " " ;
	return angle <= numeric_limits<float>::epsilon() ? 1 : 0 ;
}

/*! 
 *  \relates MyVector
 *  \brief This function determines whether two given vectors are mutually perpendicular
 *  \param vec1 reference to a MyVector object of type T 
 *  \param vec2 reference to a MyVector object of type T 
 *  \return True if perpendicular; False if not.
 */
template <class T>
bool isPerpendicular (MyVector<T> &vec1, MyVector<T> &vec2)
{
	float dotProduct = vec1 * vec2 ;
	cout << dotProduct << " " ;
	//return dotProduct <= numeric_limits<float>::epsilon() ? 1 : 0 ;
	return dotProduct <= MINIMUM ? 1 : 0 ;
}

/*! 
 *  \relates MyVector
 *  \brief This function computes the cross product between two vectors (3D)
 *  \param vec1 reference to a MyVector object of type T 
 *  \param vec2 reference to a MyVector object of type T 
 *  \return Cross product 
 */
template <class T>
MyVector<T> crossProduct (MyVector<T> &vec1, MyVector<T> &vec2)
{
	if (vec1.length() != vec2.length())
		error ("In computing cross product: dimensions mismatch!") ;
	else
	{
		assert (vec1.length() == 3) ;
		MyVector<T> result (3) ;
		result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1] ;
		result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2] ;
		result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0] ;
		return result ;
	}
}

/*
 *  \fn MyVector<T> :: ~MyVector()
 *  \brief Destructor function of MyVector class
 *//*
template <class T>
MyVector<T> :: ~MyVector()
{
	v.clear() ;
}*/

#endif

