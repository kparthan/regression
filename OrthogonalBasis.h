/*!
 *	\file OrthogonalBasis.h
 *  \details Implementation of OrthogonalBasis class
 *  \author Parthan Kasarapu
 *  \date Modified: Mon 25 Jun 2012
 */

#ifndef ORTHOGONAL_BASIS_H
#define ORTHOGONAL_BASIS_H

#include <iostream>
#include <cstdlib>
#include "Data.h"
#include "Error.h"
#include "Matrix.h"
#include <boost/math/constants/constants.hpp>

using namespace std ;

/*!
 *  \class OrthogonalBasis
 *  \brief This is the OrthogonalBasis class abstraction.
 *
 *  The class acts as an interface to compute the value of
 *	orthogonal basis functions.
 */

class OrthogonalBasis
{
	private:
		int numFunctions ;
		double timePeriod ;
		int functionIndex ;
	public:
		//! constructor
		OrthogonalBasis (int, double, int) ;
		//! computes the design matrix
		template <class T>
		Matrix<T> designMatrix (Data<T> &) ;
} ;

/*!
 *	\fn OrthogonalBasis :: OrthogonalBasis (int)
 *	This constructor function sets the number of orthogonal 
 *	functions to be used.
 */
OrthogonalBasis :: OrthogonalBasis (int n, double length, int index) : 
								numFunctions(n), timePeriod(length), functionIndex(index)
{
}

/*!
 *	\fn Matrix<T> OrthogonalBasis :: designMatrix (Data<T> &data)
 *	\brief This module computes the design matrix based on the 
 *	number of orthogonal functions - odd sines 
 *	\param data a reference to a data object of class T
 *	\return the design matrix to be used in the computation 
 *	of weights of the corresponding orthogonal functions
 */
template <class T>
Matrix<T> OrthogonalBasis :: designMatrix (Data<T> &data)
{
	int numPoints = data.nPoints() ;
	Matrix<T> phi(numPoints,numFunctions) ;
	T pi = boost::math::constants::pi<T>() ;
	switch (functionIndex)
	{
		case 0:	// sawtooth
			for (int i=0; i<numPoints; i++)
			{
				for (int j=0; j<numFunctions; j++)
				{
					double arg = 2 * pi * (j+1) * data[i].x() / timePeriod ;
					phi[i][j] = sin (arg) ;
				}
			}
			break ;
		case 1:	//square
			for (int i=0; i<numPoints; i++)
			{
				for (int j=0; j<numFunctions; j++)
					{
						double arg = 2.0 * (2*j+1) * pi * data[i].x() / timePeriod ;
						phi[i][j] = sin (arg) ;
					}
			}
			break ;
	}
	return phi ;
}

#endif

