/*!
 *  \file Vector.h
 *  \details Implementation of MyVector class
 *  \author Parthan Kasarapu
 *  \version 1.0
 *  \date Modified: Fri 8 Jun 2012
 */
#ifndef LIBLCB_VECTOR_H__
#define LIBLCB_VECTOR_H__

#include <stdexcept>
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <cassert>

#include <boost/math/constants/constants.hpp>

//namespace lcb {

/*! 
 *  \class Vector
 *  \brief This is a Vector class.
 *
 *  The class acts as an interface to create a n-length vector and store
 *  elements in it. Also provides basic functionalities such as norm of
 *  the vector, dot product, cross product between two vectors,.
 */
template <class T>
class Vector
{
	private:
		int size ;
		std::vector<T> v ;
	public:
				/* Constructors  */
		//! creates a null vector
		Vector () ;
		//! Zero based vector					
		explicit Vector (int) ;			
		//! Initialize to constant value
		Vector (T , int) ;			
		//! Initialize to C/C++ style array
		Vector (T *, int, int) ;			
		//! Copy constructor
		Vector (const Vector &) ;			

				/* Overloading = [] + += - -= * operators */
		//! assignment to a sourceVector
		Vector operator = (Vector) ;	
		//! assignment to a constant
		Vector operator = (T) ;		
		//! return i^{th} element
		T & operator [] (int) ;		
		//! computes the dot product between vectors
		float operator * (Vector &) ;
		//! computes the sum of two vectors
		Vector operator + (Vector) ;
		//! adds a vector to the existing vector 
		Vector & operator += (Vector &) ;
		//! computes the difference between two vectors
		Vector operator - (Vector) ;
		//! subtracts a vector from an existing one 
		Vector & operator -= (Vector &) ;
		//! scaling a vector by a constant value
		Vector operator * (T) ;

				/* Other sub-routines */
		//! print the elements of the vector
		void print (void)  ;				
		//! return size of vector
		 int length(void)  ;	
		//! computes the L2-Norm (magnitude) of the vector
		float l2Norm () ;
		//! normalizes the vector
		Vector<T> normalize () ;
		
				/* destructor */
		//~Vector() ;					
} ;

/*!
 *  \fn Vector<T> :: Vector() 
 *  \brief Null vector Constructor.
 *  \return A new instance of a vector.
 */
template <class T>
Vector<T> :: Vector() : size(0), v(size)
{
}

/*! 
 *  \fn Vector<T> :: Vector(int n)
 *  \brief This constructor function creates a zero vector of size n.
 *  \param n an integer.
 *  \return A new instance of a zero vector.
 */
template <class T>
Vector<T> :: Vector(int n) : size(n)
{
	v.assign(size,0) ;
}

/*! 
 *  \fn Vector<T> :: Vector(T a, int n)
 *  \brief This constructor function creates a vector of size n with all 
 *  elements equal.
 *  \param a a constant value of type T
 *  \param n an integer.
 *  \return A new instance of a vector with same elements.
 */
template <class T>
Vector<T> :: Vector(T a, int n) : size(n)
{
	v.assign(size,a) ;
}

/*! 
 *  \fn Vector<T> :: Vector(T *a, int vectorSize, int arraySize) 
 *  \brief This constructor function creates a vector of size n and 
 *  initializes to a C/C++ style array.
 *  \param a pointer to a C/C++ array
 *  \param vectorSize an integer.
 *  \param arraySize an integer.
 *  \return A new instance of a vector with elements as that of the C/C++ 
 *  array.
 */
template <class T>
Vector<T> :: Vector(T *a, int vectorSize, int arraySize) : size(vectorSize)
{
	if (vectorSize <= arraySize)
	{
		for (int i=0; i<size; i++)
			v.push_back(a[i]) ;
	}
	else
	{
		int i ;
		for (i=0; i<arraySize; i++)
			v.push_back(a[i]) ;
		for (; i<size; i++)
			v.push_back(0) ;
	}
}

/*! 
 *  \fn Vector<T> :: Vector(const Vector<T> &sourceVector)
 *  \brief This is a copy constructor function which creates a vector 
 *  and instantiates it with an existing vector.
 *  \param sourceVector reference to a Vector object of type T
 *  \return A copy of the source Vector object.
 */
template <class T>
Vector<T> :: Vector(const Vector<T> &sourceVector) : size(sourceVector.size), 
						    v(sourceVector.v)
{
}

/*! 
 *  \fn Vector<T> Vector<T> :: operator = (Vector<T> sourceVector)
 *  \brief The assignment operator is overloaded to assign the sourceVector 
 *  to the current Vector object.
 *  \param sourceVector a Vector object of type T
 *  \return An instantiated vector which is a copy of the original 
 *  sourceVector.
 */
template <class T>
Vector<T> Vector<T> :: operator = (Vector<T> sourceVector)
{
	if (this != &sourceVector)
	{
		size = sourceVector.size ;
		v = sourceVector.v ;	// the contents are 
					// automatically erased
	}
	return *this ;
}

/*! 
 *  \fn Vector<T> Vector<T> :: operator = (T a)
 *  \brief The assignment operator is overloaded to reassign a vector 
 *  with equal elements.
 *  \param a a reference to a constant value
 *  \return A modified vector whose all elements are equal.
 */
template <class T>
Vector<T> Vector<T> :: operator = (T a)
{
	v.clear() ;	// erase the contents before copying 
	for (int i=0; i<size; i++)
		v.push_back(a) ;
	return *this ;
}

/*! 
 *  \fn  T & Vector<T> :: operator [] ( int i)
 *  \brief The [] operator is overloaded which returns the element at 
 *  position i.
 *  \param i an integer index 
 *  \return The ith vector element.
 */
template <class T>
T & Vector<T> :: operator [] (int i) 
{
	if (i >= size)
	  throw std::out_of_range("index out of range ...");
	else return v[i] ;	// alternate vector access
}

/*! 
 *  \fn double Vector<T> :: operator * (Vector &vec)
 *  \brief The * operator is overloaded to compute the dot product.
 *  \param vec reference to a Vector object of type T
 *  \return The dot product of two vectors.
 */
template <class T>
float Vector<T> :: operator * (Vector<T> &vec)
{
	float dotProduct = 0 ;
	if (size != vec.size)
	  throw std::length_error("Vector sizes don't match!") ;
	else
	{
		for (int i=0; i<size; i++)
			dotProduct += v[i] * vec[i] ;
	}
	return dotProduct ;
}

/*! 
 *  \fn Vector<T> Vector<T> :: operator + (Vector<T> other)
 *  \brief The + operator is overloaded to compute the sum of two vectors.
 *  \param other a Vector object of type T
 *  \return The vector sum.
 */
template <class T>
Vector<T> Vector<T> :: operator + (Vector<T> other)
{
	if (size != other.size)
	  throw std::length_error("dimensions mismatch!") ;
	else
		return Vector<T>(*this)+=other ;
}

/*! 
 *  \fn Vector<T> & Vector<T> :: operator += (Vector<T> &other)
 *  \brief The += operator is overloaded to add a vector to an existing one.
 *  \param other reference to a Vector object of type T
 *  \return A reference to the original vector.
 */
template <class T>
Vector<T> & Vector<T> :: operator += (Vector<T> &other)
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
 *  \fn Vector<T> Vector<T> :: operator - (Vector other)
 *  \brief The - operator is overloaded to compute the difference of two 
 *  vectors.
 *  \param other a Vector object
 *  \return The vector difference.
 */
template <class T>
Vector<T> Vector<T> :: operator - (Vector<T> other)
{
	if (size != other.size)
	  throw std::length_error("dimensions mismatch!") ;
	else
		return Vector<T>(*this)-=other ;
}
 
/*! 
 *  \fn Vector<T> & Vector<T> :: operator -= (Vector<T> &other)
 *  \brief The -= operator is overloaded to subtract a vector from an 
 *  existing one.
 *  \param other reference to a Vector object of type T
 *  \return A reference to the original vector.
 */
template <class T>
Vector<T> & Vector<T> :: operator -= (Vector<T> &other)
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
 *  \fn Vector<T> Vector<T> :: operator * (T a)
 *  \brief The * operator is overloaded to scale the vector by a ant value.
 *  \param a constant value of type T
 *  \return The scaled vector.
 */
template <class T>
Vector<T> Vector<T> :: operator * (T a)
{
	Vector<T> result (size) ;
	for (int i=0; i<size; i++)
		result[i] = v[i] * a ;
	return result ;
}

/*! 
 *  \fn void Vector<T> :: print(void)
 *  \brief The method prints the elements in the vector.
 */
template <class T>
void Vector<T> :: print(void) 
{
  std::cout << size << std::endl ;
  for (int i=0; i<size; i++)
    std::cout << v.at(i) << " " ;
  std::cout << std::endl ;
}

/*! 
 *  \fn  int Vector<T> :: length() 
 *  \brief The method returns the size of the vector.
 *  \return vector size
 */
template <class T>
 int Vector<T> :: length() 
{
	return size ;
}

/*! 
 *  \fn float Vector<T> :: l2Norm(void)
 *  \brief The method computes the magnitude of the vector. 
 *  \return vector length (L2-norm)
 */
template <class T>
float Vector<T> :: l2Norm(void)
{
	float magnitude = 0 ;
	for (int i=0; i<size; i++)
		magnitude += v[i] * v[i] ;
	return sqrt(magnitude) ;
}

/*! 
 *  \fn Vector<T> Vector<T> :: normalize(void)
 *  \brief The method normalizes the vector. 
 *  \return Normalized unit vector.
 */
template <class T>
Vector<T> Vector<T> :: normalize (void)
{
	if (size == 0)
	  throw std::length_error("trying to normalize a zero vector!") ;
	else
	{
		Vector<T> result (*this) ;
		float vectorMagnitude = result.l2Norm() ;
		for (int i=0; i<size; i++)
			result[i] = result[i] / vectorMagnitude ;
		return result ;
	}
}

/*! 
 *  \relates Vector
 *  \brief This function computes the angle between two vectors.
 *  \param vec1 reference to a Vector object of type T 
 *  \param vec2 reference to a Vector object of type T 
 *  \return Angle between the vectors in radians.
 */
template <class T>
float angleBetween (Vector<T> &vec1, Vector<T> &vec2)
{
	float dotProduct = vec1 * vec2 ;
	float vector1Magnitude = vec1.l2Norm() ;
	float vector2Magnitude = vec2.l2Norm() ;
	// range of arccos is [0,PI] 
	return acos(dotProduct / (vector1Magnitude * vector2Magnitude));
}

/*! 
 *  \relates Vector
 *  \brief This function determines whether two given vectors are parallel 
 *  or not 
 *  \param vec1 reference to a Vector object of type T 
 *  \param vec2 reference to a Vector object of type T 
 *  \return True if parallel; False if not parallel.
 */
template <class T>
bool isParallel (Vector<T> &vec1, Vector<T> &vec2)
{
	float angle = angleBetween (vec1,vec2) ;
	return angle <= std::numeric_limits<float>::epsilon() ? 1 : 0 ;
}

/*! 
 *  \relates Vector
 *  \brief This function determines whether two given vectors are mutually 
 *  perpendicular
 *  \param vec1 reference to a Vector object of type T 
 *  \param vec2 reference to a Vector object of type T 
 *  \return True if perpendicular; False if not.
 */
template <class T>
bool isPerpendicular (Vector<T> &vec1, Vector<T> &vec2)
{
	float angle = angleBetween (vec1,vec2) ;
	return angle ==  boost::math::constants::pi<float>()/2.0 ? 1 : 0 ;
}

/*! 
 *  \relates Vector
 *  \brief This function computes the cross product between two vectors (3D)
 *  \param vec1 reference to a Vector object of type T 
 *  \param vec2 reference to a Vector object of type T 
 *  \return Cross product 
 */
template <class T>
Vector<T> crossProduct (Vector<T> &vec1, Vector<T> &vec2)
{
	if (vec1.length() != vec2.length())
	  throw std::length_error("dimensions mismatch!") ;
	else if(vec1.length() == 3)
	{
		assert (vec1.length() == 3) ;
		Vector<T> result (3) ;
		result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1] ;
		result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2] ;
		result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0] ;
		return result ;
	}else{
	  throw std::length_error("Cannot handle length other than 3");
	}
}

//} //namespace lcb

#endif /* LIBLCB_VECTOR_H__ */

