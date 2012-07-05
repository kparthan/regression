/*
 *  \file Point.h
 *  \details Implementation of Point class
 *  \author Parthan Kasarapu
 *  \date Modified: Thu 21 Jun 2012
 */
#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <cstdlib>

/*!
 *  \class Point
 *  \brief This is a Point class.
 *
 *  The class acts as an interface to maintain a point
 *  represented by x coordinate.
 */
template <class T>
class Point
{
	private:
		T xCoordinate ;
	public:
				/* Constructors */
		//! null constructor
		Point() ;
		//! constructor
		Point (T) ;
		//! Copy constructor
		Point (const Point &) ;

				/* Overloading = + += - -= operators */
		//! assignment to a source Point object
		Point operator = (Point) ;
		//! adds the coordinates to the current Point object
		Point & operator += (Point &) ;
		//! sums the coordinates of two Point objects
		Point operator + (Point) ;
		//! subtracts the coordinates from the current Point 
		Point & operator -= (Point &) ;
		//! subtracts the coordinates of two Point objects
		Point operator - (Point) ;
		//! compares which of the Points are closer to the origin
		bool operator < (Point) ;

				/* Other sub-routines  */
		//! retrieve the coordinate value
		T x () ;
} ;

/*!
 *  \fn Point<T> :: Point ()
 *  \brief This constructor creates a null Point object. 
 *  \return a new instance of Point object with coordinate initialized to zero
 */
template <class T>
Point<T> :: Point () : xCoordinate (0)
{
}

/*!
 *  \fn Point<T> :: Point (T value)
 *  \brief This constructor sets the coordinate value of the point.
 *  \param value a scalar of type T
 *  \return a new instance of Point object
 */
template <class T>
Point<T> :: Point (T value) : xCoordinate (value)
{
}

/*!
 *  \fn Point<T> :: Point (const Point<T> &sourcePoint)
 *  \brief The copy constructo creates a new Point object and 
 *   	   instantiates it with a sourcePoint object
 *  \param sourcePoint a reference to a Point object
 *  \return a new instance of Point object
 */
template <class T>
Point<T> :: Point (const Point<T> &sourcePoint) : 
	    xCoordinate (sourcePoint.xCoordinate)
{
}

/*!
 *  \fn Point<T> :: Point (T value)
 *  \brief This constructor sets the coordinate value of the point.
 *  \param value a scalar of type T
 *  \return a new instance of Point object
 */
template <class T>
Point<T> Point<T> :: operator = (Point<T> sourcePoint)
{
	if (this != &sourcePoint)
		xCoordinate = sourcePoint.xCoordinate ;
	return *this ;
}

/*!
 *  \fn Point<T> & Point<T> :: operator += (Point<T> &sourcePoint)
 *  \brief The function adds the sourcePoint to the current Point object 
 *  \param sourcePoint a reference to a Point object 
 *  \return a reference to a Point object of type T 
 */
template <class T>
Point<T> & Point<T> :: operator += (Point<T> &sourcePoint)
{
	xCoordinate += sourcePoint.xCoordinate ;
	return *this ;
}

/*!
 *  \fn Point<T> & Point<T> :: operator + (Point<T> sourcePoint)
 *  \brief The function adds the coordinates of two Point objects 
 *  \param sourcePoint a Point object 
 *  \return sum of coordinates of two Point objects 
 */
template <class T>
Point<T> Point<T> :: operator + (Point<T> sourcePoint)
{
	return Point<T>(*this)+=sourcePoint ;
}

/*!
 *  \fn Point<T> & Point<T> :: operator -= (Point<T> &sourcePoint)
 *  \brief The function subtracts the sourcePoint from the current Point object
 *  \param sourcePoint a reference to a Point object 
 *  \return a reference to a Point object of type T 
 */
template <class T>
Point<T> & Point<T> :: operator -= (Point<T> &sourcePoint)
{
	xCoordinate -= sourcePoint.xCoordinate ;
	return *this ;
}

/*!
 *  \fn Point<T> & Point<T> :: operator - (Point<T> sourcePoint)
 *  \brief The function subtracts the coordinates of two Point objects 
 *  \param sourcePoint a Point object 
 *  \return difference of coordinates of two Point objects 
 */
template <class T>
Point<T> Point<T> :: operator - (Point<T> sourcePoint)
{
	return Point<T>(*this)-=sourcePoint ;
}

/*!
 *  \fn bool Point<T> :: operator < (Point other)
 *  \brief The function compares which of the Point objects are closer
 *  to the origin
 *  \param other a Point object 
 *  \return 1 if current object is closer, 0 otherwise 
 */
template <class T>
bool Point<T> :: operator < (Point other)
{
	return (x() < other.x()) ? 1 : 0 ;
}

/*!
 *  \fn T Point<T> :: x (void) 
 *  \brief This function retrieves the coordinate value of the point
 *  \return the coordinate value of type T 
 */
template <class T>
T Point<T> :: x (void)
{
	return xCoordinate ;
}

#endif

