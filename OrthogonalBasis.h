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
#include <liblcb/Matrix.h>
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
    //! null constructor
    OrthogonalBasis() ;
		//! constructor
		OrthogonalBasis (int, double, int) ;
		//! computes the design matrix
		template <class T>
		lcb::Matrix<T> designMatrix (Data<T> &) ;
} ;

/*!
 *  \brief Null constructor for the OrthogonalBasis class
 */
OrthogonalBasis :: OrthogonalBasis () : numFunctions(0),
                   timePeriod(0), functionIndex(-1)
{
}

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
 *	\fn lcb::Matrix<T> OrthogonalBasis :: designMatrix (Data<T> &data)
 *	\brief This module computes the design matrix based on the 
 *	number of orthogonal functions - odd sines 
 *	\param data a reference to a data object of class T
 *	\return the design matrix to be used in the computation 
 *	of weights of the corresponding orthogonal functions
 */
template <class T>
lcb::Matrix<T> OrthogonalBasis :: designMatrix (Data<T> &data)
{
	int numPoints = data.nPoints() ;
	lcb::Matrix<T> phi(numPoints,numFunctions) ;
	T pi = boost::math::constants::pi<T>() ;
	/*switch (functionIndex)
	{
		case 0:	// sawtooth
			for (int i=0; i<numPoints; i++)
			{ 
				int k = 0 ;
				for (int j=0; j<numFunctions; j++)
				{
					double arg = 2 * pi * (j+1) * data[i].x() / timePeriod ;
					phi[i][k++] = sin (arg) ;
					phi[i][k++] = cos (arg) ;
				}
			}
			break ;
		case 1:	//square
			for (int i=0; i<numPoints; i++)
			{
				int k = 0 ;
				for (int j=0; j<numFunctions; j++)
					{
						double arg = 2 * pi * (j+1) * data[i].x() / timePeriod ;
						//double arg = 2.0 * (2*j+1) * pi * data[i].x() / timePeriod ;
						//phi[i][j] = sin (arg) ;
						phi[i][k++] = sin (arg) ;
						phi[i][k++] = cos (arg) ;
					}
			}
			break ;
	}*/
	//double phase = (rand() / RAND_MAX) * 2 * pi ;
	/*for (int i=0; i<numPoints; i++)
	{
		for (int j=0; j<numFunctions; j++)
		{
			int k = j / 2 + 1 ;
			long double arg = 2 * pi * k * data[i].x() / timePeriod ;
			if (j % 2 == 0)
				phi[i][j] = sin (arg) ;
			else
				phi[i][j] = cos (arg) ;
		}
	}*/
  int k ;
  long double arg ;
	/*for (int i=0; i<numPoints; i++)
	{
    phi[i][0] = 1.0 ;
		for (int j=1; j<numFunctions; j++)
		{
			if (j % 2 == 1)
      {
        k = j / 2 + 1 ; 
        arg =  k * data[i].x() / timePeriod ;
				phi[i][j] = sin (arg) ;
      }
			else
      {
        k = j / 2 ; 
        arg = k * data[i].x() / timePeriod ;
				phi[i][j] = cos (arg) ;
      }
		}
	}*/
	for (int i=0; i<numPoints; i++)
	{
    phi[i][0] = 1.0 ;
		for (int j=0; j<numFunctions-1; j++)
		{
      k = j / 2 + 1 ; 
      arg =  k * data[i].x() / timePeriod ;
			if (j % 2 == 0)
      {
				phi[i][j+1] = sin (arg) ;
      }
			else
      {
				phi[i][j+1] = cos (arg) ;
      }
		}
	}
	return phi ;
}

#endif

