/*!
 *  \file Matrix.h
 *  \details Implementation of Matrix class
 *  \author Parthan Kasarapu
 *  \version 1.0
 *  \date Modified: Jun 11 2012
 */

#ifndef LIBLCB_MATRIX_H__
#define LIBLCB_MATRIX_H__

#include <stdexcept>
#include <iostream>
#include <limits>

#include "GenericMatrix.h"
#include "Vector.h"
#include "error.h"

//namespace lcb {

/*! 
 *  \class Matrix Matrix.h "Matrix.h"
 *  \brief This is a Matrix class.
 *
 *  The class acts as an interface to do various operations associated 
 *  with linear algebra. Matrix operations such as addition, multiplication,
 *  inverse, and solving a system of linear equations can be done using 
 *  this implementation.
 */
template <class T>
class Matrix : public GenericMatrix
{
	private:
  		//int numRows, numCols ; replaced by rows() and columns()
		std::vector <Vector<T> > m ;
	public:
				/* Constructors */
		//! null matrix	
		Matrix () ;		
		//! rectangular zero matrix		
		Matrix (int, int) ;		
		//! square zero matrix
		Matrix (int) ;	
		//! intialize to constant					
		Matrix (const T, int, int) ;	
		// Copy constructor			
		//Matrix (Matrix &) ;
		//! initialize to C/C++ 1-D array
		Matrix (T *, int, int) ;
		//! initialize to C/C++ 2-D array
		Matrix (T **a, int, int) ;

				/* Overloading = [] + += - -= * operators */
		//! assignment to a source matrix		
		Matrix operator = (Matrix) ;			
		//! assignment to a constant
		Matrix operator = (T) ;			
		//! returns the i^{th} row
		Vector<T> & operator [] (int) ;	
		//! sums two matrices
		Matrix operator + (Matrix) ;	
		//! adds a matrix to the current one
		Matrix & operator += (Matrix &) ;
		//! subtracts two matrices
		Matrix operator - (Matrix) ;	
		//! subtracts a matrix from the current one
		Matrix & operator -= (Matrix &) ;	
		//! multiplies two matrices
		Matrix operator * (Matrix) ;
		//! multiplies a matrix with a constant
		Matrix operator * (T) ;

				/* Other sub-routines */
		//! prints the elements of the matrix		
		void print() ;						
		//! generates a rectangular identity matrix
		template <class U> 
		friend Matrix<U> identity (int, int) ;			
		//! generates a square identity matrix
		template <class U> 
		friend Matrix<U> identity (int) ;			
		//! resizes the matrix
		void changeDimensions (unsigned , unsigned);
		//! constructs a vector from a single row/column matrix
		Vector<T> convertToVector () ;
		//! builds the transpose of the matrix
		Matrix transpose () ;
		//! retrieves a section of a specific column
		Vector<T> getColumnElements(unsigned, unsigned, unsigned) ;
		//! computes the inverse of the matrix using partial pivoting
		Matrix inverse() ;
		//! computes the determinant of a square matrix
		double determinant() ;		
		// void inverse_noPivoting() ;
		// destructor
		//~Matrix() ;
				/* TODO */
		 
		/* 	 LU Decomposition Routine 
		  [dolittle & crout matrix decomposition] */
		// Trivial:- AX = B, where A = LU
		// Partial Pivoting:- AX = B, where PAX = PB and A = LU

		// Eigen decomposition

		// SVD
	
		// Anything else??
} ;

/*! 
 *  \relates Matrix
 *  \brief This function is used to generate a matrix (with unequal 
 *  dimensions) of zeroes.
 *  \param rows an integer.
 *  \param cols an integer.	
 *  \return A zero rectangular matrix.
 */
template <class T>
Matrix<T> zeros (int rows, int cols)
{
	Matrix<T> result (rows,cols) ;
	return result ;
}

/*! 
 *  \relates Matrix
 *  \brief This function is used to generate a matrix (with equal dimensions) 
 *  of zeroes.
 *  \param dimension an integer.
 *  \return A zero square matrix.
 */
template <class T>
Matrix<T> zeros (int dimension)
{
	return zeros<T>(dimension,dimension) ;
}

/*! 
 *  \relates Matrix
 *  \brief This function is used to generate a matrix (with equal dimensions) 
 *  of ones.
 *  \param rows an integer.
 *  \param cols an integer.
 *  \return A rectangular matrix of ones.
 */
template <class T>
Matrix<T> ones (int rows, int cols)
{
	Matrix<T> result (1,rows,cols) ;
	return result ;
}

/*! 
 *  \relates Matrix
 *  \brief This function is used to generate a matrix (with equal dimensions) 
 *  of ones.
 *  \param dimension an integer.
 *  \return A square matrix of ones.
 */
template <class T>
Matrix<T> ones (int dimension)
{
	return ones<T>(dimension,dimension) ;
}

/*! 
 *  This function generates an identity matrix of unequal dimensions.
 *  \param rows an integer.
 *  \param cols an integer.
 *  \return A rectangular identity matrix.
 */
template <class T>
Matrix<T> identity (int rows, int cols)
{
	Matrix<T> result (rows,cols) ;
	for (int i=0; i<rows; i++)
		for (int j=0; j<cols; j++)
			if (i == j)
				result[i][j] = 1 ;
	return result ;
}

/*! 
 *  This function generates a square identity matrix.
 *  \param dimension an integer.
 *  \return A square identity matrix.
 */
template <class T>
Matrix<T> identity (int dimension)
{
	return identity<T> (dimension,dimension) ;
}

/*! 
 *  null martix Constructor.
 *  \return A new instance of a matrix.
 */
template <class T>
Matrix<T> :: Matrix() : GenericMatrix(0,0), m(0)
{
}

/*! 
 *  \fn Matrix<T> :: Matrix(int rows, int cols)
 *  \brief This constructor function creates a rectangular zero matrix.
 *  \param rows an integer.
 *  \param cols an integer.
 *  \return A new instance of a zero matrix.
 */
template <class T>
Matrix<T> :: Matrix(int rows, int cols) : GenericMatrix(rows, cols)
{
	m.resize(rows) ;
  	for (int i=0; i<rows; i++)
	{
		Vector<T> vec (cols) ;
		m[i] = vec ;
		for (unsigned j=0; j<columns(); j++)
			m[i][j] = (T)0 ;
	}
}

/*! 
 *  \fn Matrix<T> :: Matrix(int dimension)
 *  \brief This constructor function creates a square zero matrix.
 *  \param dimension an integer.
 *  \return A new instance of a zero matrix.
 */
template <class T>
Matrix<T> :: Matrix (int dimension) : GenericMatrix(dimension, dimension)
{
	m.resize(dimension) ;
	for (int i=0; i<dimension; i++)
	{
		Vector<T> vec (dimension) ;
		m[i] = vec ;
		for (int j=0; j<dimension; j++)
			m[i][j] = (T)0 ;
	}
}

/*! 
 *  \fn Matrix<T> :: Matrix (const T a, int rows, int cols) 
 *  \brief This constructor function generates a matrix of equal elements.
 *  \param a a constant value of type T.
 *  \param rows an integer.
 *  \param cols an integer.
 *  \return A new instance of a matrix with all elements equal.
 */
template <class T>
Matrix<T> :: Matrix (const T a, int rows, int cols) : 
	     GenericMatrix(rows, cols)
{
	m.resize(rows) ;
	for (int i=0; i<rows; i++)
	{
		Vector<T> vec (cols) ;
		m[i] = vec ;
		for (int j=0; j<cols; j++)
			m[i][j] = (T)a ;
	}
}

/*
 *  \fn Matrix<T> :: Matrix (const Matrix<T> &sourceMatrix)
 *  \brief The copy constructor instantiates a new Matrix object which is a 
 *  copy of the sourceMatrix.
 *  \param sourceMatrix a reference to a matrix object.
 *  \return A new instance of a matrix.
 */
/*template <class T>
Matrix<T> :: Matrix (Matrix<T> &sourceMatrix) : rows()(sourceMatrix.rows()),
columns()(sourceMatrix.columns())
{
	int i,j ;
	m.resize(rows()) ;
	for (i=0; i<rows(); i++)
	{
		Vector<T> vec (columns()) ;
		m[i] = vec ;
		for (j=0; j<columns(); j++)
			m[i][j] = sourceMatrix.m[i][j] ;
	}
}*/

/*! 
 *  \fn Matrix<T> :: Matrix (T *a, int rows, int cols) 
 *  \brief This method assigns the elements in the C/C++ style 1-D array to a *  matrix with appropriate dimensions.
 *  \param a C/C++ style 1-D array of type T
 *  \param rows an integer
 *  \param cols an integer
 *  \return A new instance of a matrix.
 */
template <class T>
Matrix<T> :: Matrix (T *a, int rows, int cols) : GenericMatrix(rows, cols)
{
	m.resize(rows) ;
	for(int i=0; i<rows; i++)
	{
		Vector<T> vec (cols) ;
		m[i] = vec ;
		for(int j=0; j<cols; j++)
		{
		  if (a!=NULL)	// need an elegant way to check the 
				// end of array?
		    m[i][j] = *a++ ;
		  else
		    throw std::invalid_argument("array data insufficient");
		}
	}
}

/*! 
 *  \fn Matrix<T> :: Matrix (T **a, int rows, int cols) 
 *  \brief This method assigns the elements in the C/C++ style 2-D array to a *  matrix with appropriate dimensions.
 *  \param a C/C++ style 2-D array of type T
 *  \param rows an integer
 *  \param cols an integer
 *  \return A new instance of a matrix.
 */
template <class T>
Matrix<T> :: Matrix (T **a, int rows, int cols) : GenericMatrix(rows, cols)
{
	m.resize(rows) ;
	for (int i=0; i<rows; i++)
	{
		Vector<T> vec (cols) ;
		m[i] = vec ;
		for (int j=0; j<cols; j++)
			m[i][j] = a[i][j] ;
	}
}

/*! 
 *  \fn Matrix<T> Matrix<T> :: operator = (const Matrix<T> sourceMatrix)
 *  \brief The assignment operator is overloaded to assign the sourceMatrix 
 *  to the current matrix object. The dimensions of the original matrix are 
 *  altered to that of the assigned matrix.
 *  \param sourceMatrix a matrix object of type T.
 *  \return A copy of the source matrix object.
 */
template <class T>
Matrix<T> Matrix<T> :: operator = (Matrix<T> sourceMatrix)
{
	if (this != &sourceMatrix)
	{
		_number_of_rows = sourceMatrix.rows();
		_number_of_columns = sourceMatrix.columns();
		m.resize(rows()) ;
		for (unsigned i=0; i<rows(); i++)
		{
			Vector<T> vec (columns()) ;
			m[i] = vec ;
			for (unsigned j=0; j<columns(); j++)
				m[i][j] = sourceMatrix[i][j] ;
		}
	}
	return *this ;
}

/*! 
 *  \fn Matrix<T> Matrix<T> :: operator = (const T a)
 *  \brief The assignment operator is overloaded to assign a constant value 
 *  to all the elements in a matrix.
 *  \param a constant.
 *  \return A matrix object whose all elements are initialized to a constant  *  value.
 */
template <class T>
Matrix<T> Matrix<T> :: operator = (const T a)
{
	for (unsigned i=0; i<rows(); i++)
		for (unsigned j=0; j<columns(); j++)
			m[i][j] = a ;
	return *this ;
}

/*! 
 *  \fn Matrix<T> & Matrix<T> :: operator [] (const int row)
 *  \brief The [] operator is overloaded to return the ith row in the matrix.
 *  \param row an integer.
 *  \return A vector.
 */
template <class T>
Vector<T> & Matrix<T> :: operator [] (int row) 
{
	return m[row] ;
}

/*! 
 *  \fn Matrix<T> Matrix<T> :: operator + (const Matrix<T> other)
 *  \brief Adds two matrices 
 *  \param other a Matrix object of type T
 *  \return Sum matrix
 */
template <class T>
Matrix<T> Matrix<T> :: operator + (Matrix<T> other) 
{
	if (rows() != other.rows() || columns() != other.columns())
	  throw std::length_error("dimensions mismatch!") ;
	else
		return Matrix(*this)+=other ;
}

/*! 
 *  \fn Matrix<T> & Matrix<T> :: operator += (Matrix<T> &other)
 *  \brief Adds a matrix to the current matrix
 *  \param other a reference to a Matrix object
 *  \return Sum matrix
 */
template <class T>
Matrix<T> & Matrix<T> :: operator += (Matrix<T> &other)
{
	// check to see if the current matrix is empty/null
	// if so, resize the null matrix to fit the dimensions
	if ((*this).rows() == 0 && (*this).columns() == 0) 
	{
		_number_of_rows = other.rows();
		_number_of_columns = other.columns();
		m.resize(rows()) ;
		for (unsigned i=0; i<rows(); i++)
		{
			Vector<T> vec (columns()) ;
			m[i] = vec ;
		}
	}
	for (unsigned i=0; i<rows(); i++)
		for (unsigned j=0; j<columns(); j++)
			m[i][j] += other[i][j] ;
	return *this ;
}

/*! 
 *  \fn Matrix<T> Matrix<T> :: operator - (const Matrix<T> other)
 *  \brief Subtracts two matrices 
 *  \param other a Matrix object of type T
 *  \return The difference matrix
 */
template <class T>
Matrix<T> Matrix<T> :: operator - (Matrix<T> other) 
{
	if (rows() != other.rows() || columns() != other.columns())
	  throw std::length_error("dimensions mismatch!") ;
	else
		return Matrix(*this)-=other ;
}

/*! 
 *  \fn Matrix<T> & Matrix<T> :: operator -= (Matrix<T> &other)
 *  \brief Subtracts the given (other) matrix from the current one
 *  \param other a reference to a Matrix object
 *  \return The difference matrix
 */
template <class T>
Matrix<T> & Matrix<T> :: operator -= (Matrix<T> &other)
{
	// check to see if the current matrix is empty/null
	// if so, resize the null matrix to fit the dimensions
	if ((*this).rows() == 0 && (*this).columns() == 0) 
	{
		_number_of_rows = other.rows();
		_number_of_columns = other.columns();
		m.resize(rows()) ;
		for (unsigned i=0; i<rows(); i++)
		{
			Vector<T> vec (columns()) ;
			m[i] = vec ;
		}
	}
	for (unsigned i=0; i<rows(); i++)
		for (unsigned j=0; j<columns(); j++)
			m[i][j] -= other[i][j] ;
	return *this ;
}

/*! 
 *  \fn Matrix<T> Matrix<T> :: operator * (const Matrix<T> other)
 *  \brief Multiplies two matrices 
 *  \param other a reference to a Matrix object
 *  \return Product matrix
 */
template <class T>
Matrix<T> Matrix<T> :: operator * (Matrix<T> other)
{
	if (columns() != other.rows())
	  throw std::length_error("dimensions mismatch!") ;
	else
	{
		Matrix<T> result (rows(),other.columns()) ;
		for (unsigned i=0; i<rows(); i++)
			for (unsigned j=0; j<other.columns(); j++)
				for (unsigned k=0; k<columns(); k++)
					result.m[i][j] += m[i][k] * other.m[k][j] ;
		return result ;
	}
}

/*! 
 *  \fn Matrix<T> Matrix<T> :: operator * (T a)
 *  \brief Multiplies all the elements in the matrix with a constant 
 *  \param a a constant of type T
 *  \return A scaled matrix
 */
// keep in mind: operator ordering is important
template <class T>
Matrix<T> Matrix<T> :: operator * (T a)
{
	Matrix<T> result (rows(),columns()) ;
	for (unsigned i=0; i<rows(); i++)
		for (unsigned j=0; j<columns(); j++)
			result.m[i][j] = a * m[i][j] ;
	return result ;
}

/*! 
 *  \fn void Matrix<T> :: changeDimensions (int rows, int cols)
 *  \brief The function is used to resize the matrix.
 *  \param rows an integer
 *  \param cols an integer
 */
template <class T>
void Matrix<T> :: changeDimensions (unsigned rows, unsigned cols)
{
	if (rows != this->rows())
	{
	  _number_of_rows = rows;
	  m.resize(rows) ;
	}
	if (cols != columns())
	{
	  _number_of_columns = cols;
	  for (unsigned i=0; i<rows; i++)
	    {
	      Vector<T> vec (columns()) ;
	      m[i] = vec ;
	    }
	}
}

/*! 
 *  \fn Vector<T> Matrix<T> :: convertToVector ()
 *  \brief The function is used to convert a single row/column matrix
 *	into a vector for easy access.
 *	\return a Vector object of type T
 */
template <class T>
Vector<T> Matrix<T> :: convertToVector (void)
{
	if (rows() == 1 || columns() == 1)
	{
		Vector<T> result ;
		if (rows() == 1)
		{
			result = Vector<T>(columns()) ;
			for (int i=0; i<columns(); i++)
				result[i] = (*this)[0][i] ;
		}
		else if (columns() == 1)
		{
			result = Vector<T>(rows()) ;
			for (int i=0; i<rows(); i++)
				result[i] = (*this)[i][0] ;
		}
		return result ;
	}
	else
		error("In forming a vector: matrix does not have single row/column ...") ;
}

/*! 
 *  \fn Matrix<T> Matrix<T> :: transpose (void)
 *  \brief Creates a transpose of the given matrix
 *  \return Transpose of a matrix
 */
template <class T>
Matrix<T> Matrix<T> :: transpose (void)
{
	Matrix<T> result (columns(),rows()) ;
	for (unsigned i=0; i<columns(); i++)
		for (unsigned j=0; j<rows(); j++)
			result[i][j] = m[j][i] ;
	return result ;
}

/*! 
 *  \fn Matrix<T> Matrix<T> :: getColumnElements (int col, int index1, int 
 *  index2)
 *  \brief gets the column elements within the range [index1,index2]
 *  \param col is the column to get from
 *  \param index1 is the index of the first column element
 *  \param index2 is the index of the last column element
 *  \return column vector within a range 
 */
template <class T>
Vector<T> Matrix<T> :: getColumnElements (unsigned col, unsigned index1, 
					  unsigned index2)
{
	if (col >= columns())
	  throw std::invalid_argument("column index out of range!") ;
	if (index2 >= rows())
	  throw std::invalid_argument("row index out of range!") ;
	Vector<T> columnVector (index2-index1+1) ;
	int j = 0 ;
	for (unsigned i=index1; i<=index2; i++)
		columnVector[j++] = m[i][col] ;
	return columnVector ;
}

/*! 
 *  \relates Matrix
 *  \brief This function finds the row in which the pivot (absolute maximum) 
 *  is present 
 *  \param column a Vector object of type T
 *  \param offset an integer.
 *  \return The row index at which the pivot is present.
 */
template <class T>
int getPivotRow (Vector<T> column, int offset)
{
	int maxIndex=0 ;
	T maxValue = fabs(column[0]) ;
	for (int i=1; i<column.length(); i++)
	{
		if (fabs(column[i]) > maxValue)
		{
			maxValue = fabs(column[i]) ;
			maxIndex = i ;
		}
	}	
	return maxIndex + offset ;
}

/*! 
 *  \fn void Matrix<T> :: print(void)
 *  \brief The method prints the elements in the matrix.
 */
template <class T>
void Matrix<T> :: print(void)
{
  std::cout << "#rows = " << rows() << "; " ;
  std::cout << "#cols = " << columns() << std::endl ;
  for (unsigned i=0; i<rows(); i++)
    {
      for (unsigned j=0; j<columns(); j++)
	std::cout << m[i][j] << " " ;
      std::cout << std::endl ;
    }
}

/*! 
 *  \fn Matrix<T>  Matrix<T> :: inverse(void)
 *  \brief creates the inverse of the given matrix using partial pivoting 
 *  \return Matrix inverse
 */
template <class T>
Matrix<T> Matrix<T> :: inverse (void)
{
	// checking for an empty matrix
	if (rows() == 0 && columns() == 0)
	  throw std::invalid_argument("Null matrix supplied!") ;
	if (rows() != columns())
	  throw std::invalid_argument("Not a square matrix!") ;
	else
	{
		int i,j,pivotPresentAtRow,dimension = rows() ;
		Matrix<double> result = identity<double>(dimension) ;
		Matrix<double> source (*this) ;
		Vector<double> row,splitColumn ;
		double pivot,current ;
		for (i=0; i<dimension; i++)	// iterate over columns
		{
			splitColumn = source.getColumnElements(i,i,dimension-1) ;
			pivotPresentAtRow = getPivotRow (splitColumn,i) ;
			pivot = source[pivotPresentAtRow][i] ;
			if (fabs(pivot) <= std::numeric_limits<float>::epsilon())
			  throw std::domain_error("singular matrix ... ") ;
				// pivot is on a different row
				// need to swap rows
			if (i != pivotPresentAtRow)
			{
				row = source[i] ;
				source[i] = source[pivotPresentAtRow] ;
				source[pivotPresentAtRow] = row ;
				// swap the corresponding rows of the 
				// identity matrix on RHS
				row = result[i] ;
				result[i] = result[pivotPresentAtRow] ;
				result[pivotPresentAtRow] = row ;
			}
			source[i] = source[i] * (1/pivot) ; 
			result[i] = result[i] * (1/pivot) ; 
			for (j=0; j<dimension; j++)	// iterate over rows
			{
				if (j != i)
				{
					current = source[j][i] ;
					row = source[i] * current ;
					source[j] = source[j] - row ;
					row = result[i] * current ;
					result[j] = result[j] - row ; 
				}
			}
		}
		return result ;
	}
}

template <class T>
double Matrix<T> :: determinant (void)
{
	// checking for an empty matrix
	if (rows() == 0 && columns() == 0)
		return 0 ;
	if (rows() != columns())
	  throw std::invalid_argument("Not a square matrix!") ;
	Vector<double> row,splitColumn ;
	double pivot,current,det = 1 ;
	int pivotPresentAtRow,dimension = rows() ;
	Matrix<T> source (*this) ;

	for (int i=0; i<dimension; i++)	// iterate over columns
	{
		
		splitColumn = source.getColumnElements(i,i,dimension-1) ;
		pivotPresentAtRow = getPivotRow (splitColumn,i) ;
		pivot = source[pivotPresentAtRow][i] ;
		if (fabs(pivot) <= std::numeric_limits<float>::epsilon())
			return 0 ;	// singular matrix
			// pivot is on a different row
			// need to swap rows
		if (i != pivotPresentAtRow)
		{
			row = source[i] ;
			source[i] = source[pivotPresentAtRow] ;
			source[pivotPresentAtRow] = row ;
			det *= -1 ;	// change the sign of the determinant
		}
		det *= pivot ;
		for (int j=i+1; j<dimension; j++)	// iterate over rows
		{
			current = source[j][i] ;
			row = source[i] * (current / pivot) ;
			source[j] = source[j] - row ;
		}
	}
	return det ;
}

/* 
 *  \fn Matrix<T> :: ~Matrix<T>
 *  \brief Destructor function of Matrix class
 *//*
template <class T>
Matrix<T> :: ~Matrix()
{
	for (int i=0; i<rows(); i++)
		m[i].~Vector() ;
	m.clear() ;
}*/

/*
 *  \fn Matrix<T>  Matrix<T> :: inverse_noPivoting (void)
 *  \brief creates the inverse of the given matrix, NO PIVOTING [not elegant]
 *  \return Inverse matrix object
 *//*
template <class T>
void Matrix<T> :: inverse_noPivoting (void)
{
	if (rows() != columns())
		error ("In computing inverse: Not a square matrix!") ;
	else
	{
		int i,j,dimension = rows() ;
		Matrix<double> result = identity<double>(dimension) ;
		Matrix<double> source (*this) ;
		Vector<double> row ;
		double current ;
		source.print() ;
		for (i=0; i<dimension; i++)
		{
			current = source[i][i] ;
			source[i] = source[i] * (1/current) ; 
			result[i] = result[i] * (1/current) ; 
			for (j=0; j<dimension; j++)
			{
				if (j != i)
				{
					current = source[j][i] ;
					row = source[i] * current ;
					source[j] = source[j] - row ;
					row = result[i] * current ;
					result[j] = result[j] - row ; 
				}
			}
			source.print() ;
			result.print() ;
		}
	}
}*/

// } //namespace lcb

#endif /* liblcb_MATRIX_H__ */

