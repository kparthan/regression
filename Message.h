/*!
 *	\file Message.h
 *	\details Implementation of Message class
 *	\author Parthan Kasarapu
 *	\date Modified: Thu 5 Jul 2012
 */

#ifndef MESSAGE_H
#define MESSAGE_H

#define PI boost::math::constants::pi<double>()

#include <iomanip>
#include "Data.h"
#include "RandomDataGenerator.h"
#include "Error.h"
#include "Gaussian.h"
#include <liblcb/Vector.h>
#include <liblcb/Matrix.h>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

//! bound on the number of data samples
#define MAX_SAMPLES 100000

//! bound on the maximum number of orthogonal basis functions
#define MAX_FUNCTIONS 200

using namespace boost::numeric::ublas ;
using namespace std ;

struct Components
{
  long double part1;
  long double part2;
  long double total;
};

/*!
 *	\class Message
 *	\brief This is the Message class abstraction.
 *
 *	The class acts as an interface to compute the minimum message length.
 */
class Message
{
	private:
		struct Parameters parameters ;
		lcb::Vector<long double> weights ;
		Data<long double> xVals, yVals, predictions ;
    long double determinant;
	public:
    //! null constructor
    Message() ;
		//! constructor
		template <class T>
		Message (struct Parameters, lcb::Matrix<T> &, Data<T> &, Data<T> &, 
							Data<T> &) ;
		//! computes the message length
		Components messageLength() ;
    //!
    long double estimateSigma();
} ;

long double Message::estimateSigma()
{
  double diff_sq = 0;
  for (int i=0; i<parameters.numSamples; i++) {
    double diff = yVals[i].x() - predictions[i].x();
    diff_sq += diff * diff;
  }
  double diff2 = 0;
  for (int i=0; i<weights.size(); i++) {
    diff2 += weights[i] * weights[i];
  }
  diff_sq += parameters.lambda * diff2;
  return sqrt(diff_sq/parameters.numSamples);
}

long double quantizationConstant(int d)
{
  double c1 = -(1 + log(2*PI)) * d * 0.5;
  double c2 = 0.5 * log(d * PI);
  return c1+c2;
}

Components Message::messageLength()
{
  long double sigma_mml = estimateSigma();
  long double part1=0, part2=0;
  int M = parameters.numFunctions - 1;
  int N = parameters.numSamples;
  double R = 2;
  double m = parameters.lambda;
  
  // PART 1
  part1 += quantizationConstant(M+2);
  part1 += log(R);
  part1 += 0.5 * log(2*(N+M+1));
  part1 += 0.5 * log(N+m) * (M+1);
  part1 += -(M+1) * log(sigma_mml);

  // PART 2
  part2 += 0.5 * (N+M+1) * log(2*PI);
  part2 += -N * log(AOM);
  part2 += -0.5 * (M+1) * log(m);
  part2 += (N+M+1) * log(sigma_mml);
  part2 += N / 2;
  part2 += 0.5 * (M+2);

  Components msglen;
  msglen.part1 = part1/log(2);
  msglen.part2 = part2/log(2);
  msglen.total = (part1+part2)/log(2);
  return msglen;
}

matrix<long double> make_matrix(lcb::Matrix<long double> &M, int dimension)
{
	matrix<long double> m(dimension,dimension) ;
	for (unsigned i=0; i<m.size1(); i++)
	{
		matrix_row<matrix<long double> > mr(m,i) ;
		for(unsigned j=0; j<mr.size(); j++)
			mr(j) = M[i][j] ;
	}
	return m ;
}

/* Matrix inversion routine.
 Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T>
bool InvertMatrix(const matrix<T>& input, matrix<T>& inverse)
{
	typedef permutation_matrix<std::size_t> pmatrix;

	// create a working copy of the input
	matrix<T> A(input);

	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());

	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;

	// create identity matrix of "inverse"
	inverse.assign(identity_matrix<T> (A.size1()));

	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);

	return true;
}

lcb::Matrix<long double> make_my_matrix(matrix<long double> &m, int dimension)
{
	lcb::Matrix<long double> result(dimension,dimension) ;
	for (unsigned i=0; i<m.size1(); i++)
	{
		matrix_row<matrix<long double> > mr(m,i) ;
		for (unsigned j=0; j<mr.size(); j++)
			result[i][j] = mr(j) ;
	}
	return result ;
}

int determinant_sign(const permutation_matrix<size_t> &pm)
{
  int pm_sign=1;
  size_t size = pm.size();
  for (size_t i = 0; i < size; ++i) {
    if (i != pm(i)) {
      pm_sign *= -1.0; // swap_rows would swap a pair of rows here, so we change sign
    }
  }
  return pm_sign;
}

template<class T>
long double determinantBoost(lcb::Matrix<T> &M)
{
  matrix<long double> m = make_matrix(M,M.columns());
  permutation_matrix<size_t> pm(m.size1());
  long double det = 1.0;
  if(lu_factorize(m,pm) ) {
    det = 0.0;
  } else {
    for(int i = 0; i < m.size1(); i++) {
      det *= m(i,i); // multiply by elements on diagonal
    }
    det = det * determinant_sign( pm );
  }
  return det;
}

/*!
 *	\fn lcb::Matrix<T> computeWeights (lcb::Matrix<T> &phi, Data<T> 
 *  &yValues, int invChoice, long double lambda)
 *	\brief Computes the coefficients of the basis functions
 *	\param phi a reference to a Matrix object
 *	\param yValues a reference to a Data object
 *  \param invChoice an integer
 *  \param lambda a long double
 *	\return a matrix object of weights
 */
template <class T>
lcb::Matrix<T> computeWeights (lcb::Matrix<T> &phi, Data<T> &yValues, 
                               int invChoice, long double lambda)
{
	lcb::Matrix<T> phiT = phi.transpose() ;
	lcb::Matrix<T> phiTphi = phiT * phi ;
	unsigned dim = phiTphi.rows() ; 
  lcb::Matrix<long double> pseudoInv,temp,y,weights,constants,net ;
  matrix<long double> A,Z ;
  lcb::Matrix<long double>penalty = lcb::Matrix<long double>::identity(dim);
  /*ofstream detFile("temp/determinant",ios::app);
  detFile << scientific << phiTphi.determinant() << " ";
  detFile << fixed ;*/

  if (lambda > std::numeric_limits<long double>::epsilon()) {
    net = phiTphi + penalty * lambda;
    //cout << "***** IN HERE *****" << endl ;
  } else {
    net = phiTphi ;
  }
  //double determinant = net.determinant();
	y = yValues.convertToMatrix() ;
  switch(invChoice)
  {
    case 0:
      /* my implementation of matrix inverse */
	    //pseudoInv = phiTphi.inverse() ;
	    pseudoInv = net.inverse() ;
      /*detFile << scientific << pseudoInv.determinant() << " ";
      detFile << fixed ;
      detFile << scientific << determinantBoost(phiTphi) << endl;
      detFile << fixed ;
      detFile.close() ;*/
      break ;
    case 1:
      /* boost::ublas implementation of matrix inverse */
  	  A = matrix<long double> (dim,dim) ;
  	  Z = matrix<long double> (dim,dim) ;
	    A = make_matrix(net,dim) ;	
	    InvertMatrix(A,Z) ;
	    pseudoInv = make_my_matrix(Z,dim) ;
      break ;
    case 2:
      constants = phiT * y ;
      weights = net.solveLinearSystem(constants) ;
      break ;
    default:
      error("Invalid choice of matrix inverse.") ;
      break ;
  }

  if (invChoice != 2)
  {
	  temp = pseudoInv * phiT ;
	  weights = temp * y ;
	  /*ofstream phiFile ;
	  phiFile.open("temp/phi") ;
	  for (int i=0; i<phi.rows(); i++)
	  {
		  for (int j=0; j<phi.columns(); j++)
			  phiFile << phi[i][j] << " " ;
		  phiFile << endl ;
	  }
	  phiFile.close() ;

	  ofstream invFile ;
	  invFile.open("temp/inverse") ;
	  for (int i=0; i<pseudoInv.rows(); i++)
	  {
		  for (int j=0; j<pseudoInv.columns(); j++)
			  invFile << std::fixed << std::setprecision(7) << pseudoInv[i][j] << " " ;
		  invFile << endl ;
	  }
	  invFile.close() ; */
  }

	return weights ;
}

/*!
 *	\fn long double computeRMSE (lcb::Matrix<T> &weights, lcb::Matrix<T> &phi, 
 *	Data<T> &yVals)
 *	Computes the root mean squared error in approximating the data with a
 *	linear combination of given number of basis functions.
 *	\param weights a reference to a Matrix object
 *	\param phi a reference to a Matrix object
 *	\param yVals a reference to a Data object
 *	\return the root mean squared error
 */
template <class T>
long double computeRMSE (lcb::Matrix<T> &weights, lcb::Matrix<T> &phi, 
										     Data<T> &yVals, long double lambda)
{
  lcb::Matrix<T> yEst = phi * weights ; // column matrix
	long double diff, error = 0 ;
	int numSamples = phi.rows() ;
	for (int i=0; i<numSamples; i++) {
		diff = yEst[i][0] - yVals[i].x() ;
		error += diff * diff ;
	}
  /*lcb::Matrix<long double> penalty ;
  if (lambda > std::numeric_limits<long double>::epsilon()) {
    penalty = (weights.transpose() * weights) * (lambda/2) ;
	  return error/2+penalty[0][0] ;
  } else {
	  return error/2 ;
  }*/
  return sqrt(error/numSamples);
}

/*!
 *  Null constructor for the Message class
 */
//template <class T>
Message :: Message ()
{
}

/*!
 *	\fn Message :: Message (struct Parameters params, Matrix<T> &w, 
 *	Data<T> &xVals, Data<T> &yVals)
 *	\brief This is a constructor function to instantiate a Message object
 *	\param params a struct of Parameters
 *	\param w a reference to a Matrix object of type T
 *	\param xVals a reference to a Data object of type T
 *	\param yVals a reference to a Data object of type T
 */
template <class T>
Message :: Message (struct Parameters params, lcb::Matrix<T> &w, 
										Data<T> &xVals, Data<T> &yVals, Data<T> &yEst) : 
										xVals (xVals), yVals (yVals), parameters(params), 
										predictions(yEst)
{
	weights = w.getColumn(0) ;
}

#endif

