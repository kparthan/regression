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
    int basis ;
	public:
    //! null constructor
    OrthogonalBasis() ;
		//! constructor
		OrthogonalBasis (int, int, double, int) ;
		//! computes the design matrix
		template <class T>
		lcb::Matrix<T> designMatrix (Data<T> &) ;
} ;

/*!
 *  \brief Null constructor for the OrthogonalBasis class
 */
OrthogonalBasis :: OrthogonalBasis () : numFunctions(0),
                   timePeriod(0), functionIndex(-1), basis(-1)
{
}

/*!
 *	\fn OrthogonalBasis :: OrthogonalBasis (int)
 *	This constructor function sets the number of orthogonal 
 *	functions to be used.
 */
OrthogonalBasis :: OrthogonalBasis (int set, int n, double timePeriod, 
                   int index) : basis(set), timePeriod(timePeriod), 
                                numFunctions(n), functionIndex(index)
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
  int i,j,k ;
  double temp1,temp2,arg,coefficient ;
  double n = sqrt(2.0 / timePeriod);

  switch(basis) 
  {
    case 0:
      // sines & cosines
      coefficient = 2 * pi / timePeriod ;
      //coefficient = 1 ;
	    for (i=0; i<numPoints; i++)
	    {
        phi[i][0] = 1.0 ;
		    for (j=1; j<numFunctions; j++)
		    {
			    if (j % 2 == 0)
          {
            k = j / 2;
            arg =  coefficient * k * data[i].x() ;
				    phi[i][j] = n * cos (arg) ;
          }
			    else
          {
            k = j / 2 + 1;
            arg =  coefficient * k * data[i].x() ;
				    phi[i][j] = n * sin (arg) ;
          }
		    }
	    }
      break ;

     case 1:
      // legendre polynomials
      for (i=0; i<numPoints; i++)
      {
        phi[i][0] = 1.0 ;
        if (numFunctions >= 2) 
        {
          phi[i][1] = data[i].x() ;
          for (j=2; j<numFunctions; j++)
          {
            temp1 = (2 * j - 1) * data[i].x() * phi[i][j-1] ;
            temp2 = (j-1) * phi[i][j-2];
            phi[i][j] = (temp1 - temp2) / (double)j ;
          }
        }
      }
      break ;
    }
	return phi ;
}

#endif

