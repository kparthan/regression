/*!
 *  \file Data.h
 *  \details Implementation of Data class
 *  \author Parthan Kasarapu
 *  \date Modified: Mon 25 Jun 2012
 */

#ifndef DATA_H
#define DATA_H

#include "Point.h"
#include "Error.h"
#include "Matrix.h"
#include "Data.h"
#include <vector>

using namespace std ;

/*!
 *  \class Data
 *  \brief This is the Data class abstraction.
 *
 *  The class acts as an interface to store data elements and provide
 *  access to them.  
 */
template <class T>
class Data
{
	private:
		int numPoints ;
		vector <Point<T> > elements ;
		vector <Point<T> > sortedElements ;
	public:
				/* Constructors */
		//! null constructor
		Data() ;
		//! initialize to a C/C++ array
		Data (T *, int) ;
		//! copy constructor
		Data (const Data &) ;
		//! convert a Matrix object into a Data object
		Data (Matrix<T> &) ;

				/* Overloading = + operators */
		//! assigns a Data object
		Data operator = (Data) ;
		//! returns the ith data element
		Point<T> & operator [] (int) ;
		//! appends to the elements in the current Data object
		Data operator + (Data) ;

				/* Other sub-routines */
		//! converts the Data object to a Matrix object
		Matrix<T> convertToMatrix () ;
		//! sorts the data elements
		void sortElements() ;
		//! quicksort() algorithm
		template <class U>
		friend void quicksort (vector<Point<T> > &, int, int) ;
		//! partition function to find new pivot in quicksort()
		template <class U>
		friend int partition (vector<Point<T> > &, int, int) ;
		//! retrieves the minimum element
		T minimum() ;
		//! retrieves the maximum element
		T maximum() ;
		//! returns the number of data elements
		int nPoints() ;
		//! prints the data elements
		void print() ;
} ;

/*!
 *  \fn Data<T> :: Data()
 *  \brief null constructor
 *  \return a new instance of Data object
 */
template <class T>
Data<T> :: Data() : numPoints(0), elements(0), sortedElements(0)
{
}

/*!
 *  \fn Data<T> :: Data (T *array, int size)
 *  \brief initializes the elements in the Data object with those of
 *  a C/C++ style array
 *  \param array pointer to the C/C++ array
 *  \param size an integer
 *  \return a new instance of Data object
 */
template <class T>
Data<T> :: Data (T *array, int size) : numPoints(size), sortedElements(0)
{
	for (int i=0; i<numPoints; i++)
	{
		Point<T> point(*array++) ;
		elements.push_back(point) ;
	}
}

/*!
 *  \fn Data<T> :: Data (const Data<T> &sourceData)
 *  \brief This is a copy constructor which creates a Data object and 
 *  instantiates it with an existing object
 *  \param sourceData a reference to a Data object
 *  \return a copy of the Data object
 */
template <class T>
Data<T> :: Data (const Data<T> &sourceData) : numPoints(sourceData.numPoints),
						elements(sourceData.elements),
				     sortedElements(sourceData.sortedElements)
{
}

/*!
 *	\fn Data<T> :: Data (Matrix<T> &source)
 *	\brief This constructor creates a Data object from a Matrix
 *	object.
 *	\param source a reference to a Matrix object of type T
 */
template <class T>
Data<T> :: Data (Matrix<T> &source)
{
	numPoints = source.rows() ;
	for (int i=0; i<numPoints; i++)
	{
		Point<T> point(source[i][0]) ;
		elements.push_back(point) ;
	}
}

/*!
 *  \fn Data<T> Data<T> :: operator = (Data sourceData)
 *  \brief This function is used to assign the sourceData to a Data object
 *  \param sourceData a Data object of type T
 *  \return a copy of the source Data object
 */
template <class T>
Data<T> Data<T> :: operator = (Data sourceData)
{
	if (this != &sourceData)
	{
		numPoints = sourceData.numPoints ;
		elements = sourceData.elements ;
		sortedElements = sourceData.sortedElements ;
	}
	return *this ;
}

/*! 
 *  \fn  Point<T> & Data<T> :: operator [] (int i)
 *  \brief The [] operator is overloaded which returns the element at 
 *  position i.
 *  \param i an integer index 
 *  \return The ith data element.
 */
template <class T>
Point<T> & Data<T> :: operator [] (int i)  
{
        if (i >= numPoints)
		error ("index out of range ...") ;
        else 
		return elements[i] ;      // alternate vector access
}

/*!
 *  \fn Data<T> Data<T> :: operator + (Data<T> sourceData)
 *  \brief This function is used to append the existing data elements
 *  with the elements of source Data object
 *  \param sourceData a Data object of type T
 *  \return Data object whose elements are formed by combining the 
 *  elements from two Data objects
 */
template <class T>
Data<T> Data<T> :: operator + (Data<T> sourceData)
{
	int i ;
	Data<T> result ;
	result.numPoints = numPoints + sourceData.numPoints ;
	for (i=0; i<numPoints; i++)
		result.elements.push_back(elements[i]) ;
	for (; i<result.numPoints; i++)
		result.elements.push_back(sourceData.elements[i-numPoints]) ;
	return result ;
}

/*!
 *	\fn Matrix<T> Data<T> :: convertToMatrix (void)
 *	\brief This function is used to convert a Data object to a Matrix
 *	object to enable matrix operations on data
 *	\return a Matrix object with a single column containing the 
 *	data elements
 */
template <class T>
Matrix<T> Data<T> :: convertToMatrix (void)
{
	if (numPoints == 0)
		return Matrix<T>() ;
	else
	{
		Matrix<T> result (numPoints,1) ;
		for (int i=0; i<numPoints; i++)
			result[i][0] = elements[i].x() ;
		return result ;
	}	
}

/*!
 *  \fn void Data<T> :: sortElements (void)
 *  \brief This module is used to sort the list of data elements
 */
template <class T>
void Data<T> :: sortElements (void)
{
	sortedElements = elements ;
	quicksort (sortedElements,0,numPoints-1) ;
}

/*!
 *  This is an implementation of the classic quicksort() algorithm to sort a
 *  list of data values. The module uses the overloading operator (<) to 
 *  compare two Point<T> objects. 
 *  Pivot is chosen as the right most element in the list (default)
 *  This function is called recursively.
 *  \param list a reference to a std::vector of Point objects of type T
 *  \param left an integer
 *  \param right an integer
 */
template <class T>
void quicksort (vector<Point<T> > &list, int left, int right)
{
	if (left < right)
	{
		int pivotNewIndex = partition (list,left,right) ;
		quicksort (list,left,pivotNewIndex-1) ;
		quicksort (list,pivotNewIndex+1,right) ;
	}
}

/*!
 *  This function is called from the quicksort() routine to compute the new
 *  pivot index.
 *  \param list a reference to a std::vector of Point objects of type T
 *  \param left an integer
 *  \param right an integer
 *  \return the new pivot index
 */
template <class T>
int partition (vector<Point<T> > &list, int left, int right)
{
	Point<T> temp,pivotPoint = list[right] ;
	int storeIndex = left ;
	for (int i=left; i<right; i++)
	{
		if (list[i] < pivotPoint)
		{
			temp = list[i] ;
			list[i] = list[storeIndex] ;
			list[storeIndex] = temp ;
			storeIndex += 1 ;	
		}
	}
	temp = list[storeIndex] ;
	list[storeIndex] = list[right] ;
	list[right] = temp ;
	return storeIndex ;
}

/*!
 *	\fn T Data<T> :: minimum (void)
 *	\brief This function retrieves the minimum data element
 *	\return the minimum of all data elements
 */
template <class T>
T Data<T> :: minimum (void)
{
	if (sortedElements.size() == 0)
		sortElements() ;
	return sortedElements[0].x() ;
}

/*!
 *	\fn T Data<T> :: maximum (void)
 *	\brief This function retrieves the maximum data element
 *	\return the maximum of all data elements
 */
template <class T>
T Data<T> :: maximum (void)
{
	if (sortedElements.size() == 0)
		sortElements() ;
	return sortedElements[numPoints-1].x() ;
}

/*!
 *  \fn const int Data<T> :: nPoints (void)
 *  \brief The function returns the number of data elements
 *  \return number of data points
 */
template <class T>
int Data<T> :: nPoints (void)
{
	return numPoints ;
}
 
/*!
 *  \fn void Data<T> :: print (void)
 *  \brief prints the elements in the Data object
 */
template <class T>
void Data<T> :: print (void)
{
	cout << "#points = " << numPoints << endl ;
	for (int i=0; i<numPoints; i++)
		cout << elements[i].x() << " " ;
	cout << endl ;
	if (sortedElements.size() > 0)
	{
		cout << "sorted list: " ;
		for (int i=0; i<numPoints; i++)
			cout << sortedElements[i].x() << " " ;
		cout << endl ;
	}
}

#endif

