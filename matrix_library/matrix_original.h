/*!
 *  \file matrix.h
 *  \details Implementation of Matrix class
 *  \author Parthan Kasarapu
 *  \version 1.0
 *  \date Modified: Wed 30 May 2012
 */

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <cmath>
#include "myVector.h"
#include "error.h"

/*! 
 *  \class Matrix matrix.h "matrix_library/matrix.h"
 *  \brief This is a Matrix class.
 *
 *  The class acts as an interface to do various operations associated with linear algebra. Matrix operations such as addition, multiplication, inverse, and solving a system of linear equations can be done using this implementation.
 */
template <class T>
class Matrix
{
	private:
		int numRows, numCols ;
		vector <vector<T> > m ;
	public:
				/* Constructors */
		//! null matrix	
		Matrix () ;		
		//! rectangular zero matrix		
		Matrix (int, int) ;		
		//! square zero matrix
		Matrix (int) ;	
		//! intialize to constant					
		Matrix (const T &, int, int) ;	
		//! Copy constructor			
		Matrix (const Matrix &) ;
		//! initialize to C/C++ 1-D array
		Matrix (const T *, int, int) ;
		//! initialize to C/C++ 2-D array
		Matrix (T **a, int, int) ;

				/* Overloading = [] = - * operators */
		//! assignment to a source matrix		
		Matrix & operator = (const Matrix &) ;			
		//! assignment to a constant
		Matrix & operator = (const T &) ;			
		//! returns the i^{th} row
		inline MyVector<T> & operator [] (const int) ;	
		//! sums two matrices
		Matrix operator + (const Matrix &) ;	
		//! subtracts two matrices
		Matrix operator - (const Matrix &) ;	
		//! multiplies two matrices
		Matrix operator * (const Matrix &) ;
		//! multiplies a matrix with a constant
		Matrix operator * (const T &) ;

				/* Other sub-routines */
		//! prints the elements of the matrix		
		void print() ;						
		//! gets the number of rows
		inline int nrows() const ;				
		//! gets the number of columns
		inline int ncols() const ;				
		//! generates a rectangular identity matrix
		template <class U> 
		friend Matrix<U> identity (int, int) ;			
		//! generates a square identity matrix
		template <class U> 
		friend Matrix<U> identity (int) ;			
		//! resizes the matrix
		void changeDimensions (const int &, const int &) ;	
		//! gets the (i,j)^{th} matrix element
		inline T element (const int &, const int &) const ; 	
		//! builds the transpose of the matrix
		Matrix transpose () ;					
		//! computes the inverse of the matrix (if it is square)
		Matrix inverse() ;					
		//! destructor
		~Matrix() ;						
} ;

/*! 
 *  \relates Matrix
 *  \brief This function is used to generate a matrix (with unequal dimensions) of zeroes.
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
 *  \brief This function is used to generate a matrix (with equal dimensions) of zeroes.
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
 *  \brief This function is used to generate a matrix (with equal dimensions) of ones.
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
 *  \brief This function is used to generate a matrix (with equal dimensions) of ones.
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
	int i,j ;
	Matrix<T> result (rows,cols) ;
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			if (i == j)
				result.m[i][j] = 1 ;
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
Matrix<T> :: Matrix() : numRows(0), numCols(0), m(0)
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
Matrix<T> :: Matrix(int rows, int cols) : numRows(rows), numCols(cols)
{
	int i,j ;
	m.resize(numRows) ;
	for (i=0; i<numRows; i++)
	{
		m[i].resize(numCols) ;
		for (j=0; j<numCols; j++)
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
Matrix<T> :: Matrix (int dimension) : numRows(dimension), numCols(dimension)
{
	int i,j ;
	m.resize(numRows) ;
	for (i=0; i<numRows; i++)
	{
		m[i].resize(numCols) ;
		for (j=0; j<numCols; j++)
			m[i][j] = (T)0 ;
	}
}

/*! 
 *  \fn Matrix<T> :: Matrix (const T &a, int rows, int cols) 
 *  \brief This constructor function matrix of equal elements.
 *  \param a a reference to a constant value.
 *  \param rows an integer.
 *  \param cols an integer.
 *  \return A new instance of a matrix with all elements equal.
 */
template <class T>
Matrix<T> :: Matrix (const T &a, int rows, int cols) : numRows(rows), numCols(cols)
{
	int i,j ;
	m.resize(numRows) ;
	for (i=0; i<numRows; i++)
	{
		m[i].resize(numCols) ;
		for (j=0; j<numCols; j++)
			m[i][j] = a ;
	}
}

/*! 
 *  \fn Matrix<T> :: Matrix (const Matrix<T> &sourceMatrix)
 *  \brief The copy constructor instantiates a new Matrix object which is a copy of the sourceMatrix.
 *  \param sourceMatrix a reference to a matrix object.
 *  \return A new instance of a matrix.
 */
template <class T>
Matrix<T> :: Matrix (const Matrix<T> &sourceMatrix) : numRows(sourceMatrix.numRows), numCols(sourceMatrix.numCols)
{
	int i,j ;
	m.resize(numRows) ;
	for (i=0; i<numRows; i++)
	{
		m[i].resize(numCols) ;
		for (j=0; j<numCols; j++)
			m[i][j] = sourceMatrix.m[i][j] ;
	}
}

/*! 
 *  \fn Matrix<T> :: Matrix (const T *a, int rows, int cols) 
 *  \brief This method assigns the elements in the C/C++ style 1-D array to a matrix with appropriate dimensions.
 *  \param a C/C++ style 1-D array of type T
 *  \param rows an integer
 *  \param cols an integer
 *  \return A new instance of a matrix.
 */
template <class T>
Matrix<T> :: Matrix (const T *a, int rows, int cols) : numRows(rows), numCols(cols)
{
	int i,j ;
	m.resize(numRows) ;
	for(i=0; i<numRows; i++)
	{
		m[i].resize(numCols) ;
		for(j=0; j<numCols; j++)
		{
			if (a!=NULL)	// need an elegant way to check the end of array?
				m[i][j] = *a++ ;
			else
				error("In initializing to C/C++ array: data insufficient ...") ;
		}
	}
}

/*! 
 *  \fn Matrix<T> :: Matrix (T **a, int rows, int cols) 
 *  \brief This method assigns the elements in the C/C++ style 2-D array to a matrix with appropriate dimensions.
 *  \param a C/C++ style 2-D array of type T
 *  \param rows an integer
 *  \param cols an integer
 *  \return A new instance of a matrix.
 */
template <class T>
Matrix<T> :: Matrix (T **a, int rows, int cols) : numRows(rows), numCols(cols)
{
	int i,j ;
	m.resize(numRows) ;
	for (i=0; i<numRows; i++)
	{
		m[i].resize(numCols) ;
		for (j=0; j<numCols; j++)
			m[i][j] = a[i][j] ;
	}
}

/*! 
 *  \fn Matrix<T> & Matrix<T> :: operator = (const Matrix<T> &sourceMatrix)
 *  \brief The assignment operator is overloaded to assign the sourceMatrix to the current matrix object. The dimensions of the original matrix are altered to that of the assigned matrix.
 *  \param sourceMatrix a reference to a matrix object.
 *  \return A matrix.
 */
template <class T>
Matrix<T> & Matrix<T> :: operator = (const Matrix<T> &sourceMatrix)
{
	if (this != &sourceMatrix)
	{
		int i,j ;
		numRows = sourceMatrix.numRows ;
		numCols = sourceMatrix.numCols ;
		m.resize(numRows) ;
		for (i=0; i<numRows; i++)
		{
			m[i].resize(numCols) ;
			for (j=0; j<numCols; j++)
				m[i][j] = sourceMatrix.m[i][j] ;
		}
	}
	return *this ;
}

/*! 
 *  \fn Matrix<T> & Matrix<T> :: operator = (const T &a)
 *  \brief The assignment operator is overloaded to assign a constant value to all the elements in a matrix.
 *  \param a a reference to a constant.
 *  \return A matrix.
 */
template <class T>
Matrix<T> & Matrix<T> :: operator = (const T &a)
{
	int i,j ;
	for (i=0; i<numRows; i++)
		for (j=0; j<numCols; j++)
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
MyVector<T> & Matrix<T> :: operator [] (const int row)
{
	MyVector<T> *myvec = new MyVector<T> (numCols) ;
	for (int i=0; i<numCols; i++)
		myvec[i] = m[row][i] ;
	return *myvec ; //return m[row] ;
}

/*! 
 *  \fn Matrix<T>  Matrix<T> :: operator + (const Matrix<T> &other)
 *  \brief Adds two matrices 
 *  \param other a reference to a Matrix object
 *  \return Sum matrix
 */
template <class T>
Matrix<T> Matrix<T> :: operator + (const Matrix<T> &other)
{
	if (numRows!=other.numRows || numCols!=other.numCols)
		error ("In adding matrices: dimensions mismatch!") ;
	else
	{
		int i,j ;
		Matrix<T> result (numRows,numCols) ;
		for (i=0; i<numRows; i++)
			for (j=0; j<numCols; j++)
				result.m[i][j] = m[i][j] + other.m[i][j] ;
		return result ;
	}
}

/*! 
 *  \fn Matrix<T> Matrix<T> :: operator - (const Matrix<T> &other)
 *  \brief Subtracts two matrices 
 *  \param other a reference to a Matrix object
 *  \return Difference matrix
 */
template <class T>
Matrix<T> Matrix<T> :: operator - (const Matrix<T> &other)
{
	if (numRows!=other.numRows || numCols!=other.numCols)
		error ("In subtracting matrices: dimensions mismatch!") ;
	else
	{
		int i,j ;
		Matrix<T> result (numRows,numCols) ;
		for (i=0; i<numRows; i++)
			for (j=0; j<numCols; j++)
				result.m[i][j] = m[i][j] - other.m[i][j] ;
		return result ;
	}
}

/*! 
 *  \fn Matrix<T> Matrix<T> :: operator * (const Matrix<T> &other)
 *  \brief Multiplies two matrices 
 *  \param other a reference to a Matrix object
 *  \return Product matrix
 */
template <class T>
Matrix<T> Matrix<T> :: operator * (const Matrix<T> &other)
{
	if (numCols != other.numRows)
		error ("In computing matrix product: dimensions mismatch!") ;
	else
	{
		int i,j,k ;
		Matrix<T> result (numRows,other.numCols) ;
		for (i=0; i<numRows; i++)
			for (j=0; j<other.numCols; j++)
				for (k=0; k<numCols; k++)
					result.m[i][j] += m[i][k] * other.m[k][j] ;
		return result ;
	}
}

/*! 
 *  \fn Matrix<T> Matrix<T> :: operator * (const T &a)
 *  \brief Multiplies all the elements in the matrix with a constant 
 *  \param a a constant of type T
 *  \return A scaled matrix
 */
// keep in mind: operator ordering is important
template <class T>
Matrix<T> Matrix<T> :: operator * (const T &a)
{
	int i,j ;
	Matrix<T> result (numRows,numCols) ; 
	for (i=0; i<numRows; i++)
		for (j=0; j<numCols; j++)
			result.m[i][j] = a * m[i][j] ;
	return result ;
}

/*! 
 *  \fn void Matrix<T> :: print(void)
 *  \brief The method prints the elements in the matrix.
 */
template <class T>
void Matrix<T> :: print(void)
{
	int i,j ;
	cout << "#rows = " << numRows << "; " ;
	cout << "#cols = " << numCols << endl ;
	for (i=0; i<numRows; i++)
	{
		for (j=0; j<numCols; j++)
			cout << m[i][j] << " " ;
		cout << endl ;
	}
}

/*! 
 *  \fn int Matrix<T> :: nrows() const
 *  \brief The function is used to return the number of rows in the matrix.
 *  \return the number of rows in the matrix.
 */
template <class T>
int Matrix<T> :: nrows() const
{
	return numRows ;
}

/*! 
 *  \fn int Matrix<T> :: ncols() const
 *  \brief The function is used to return the number of columns in the matrix.
 *  \return the number of columns in the matrix.
 */
template <class T>
int Matrix<T> :: ncols() const
{
	return numCols ;
}

/*! 
 *  \fn void Matrix<T> :: changeDimensions (const int &rows, const int &cols)
 *  \brief The function is used to resize the matrix.
 *  \param rows an integer reference
 *  \param cols an integer reference
 */
template <class T>
void Matrix<T> :: changeDimensions (const int &rows, const int &cols)
{
	if (rows != numRows)
	{
		numRows = rows ;
		m.resize(rows) ;
	}
	if (cols != numCols)
	{
		numCols = cols ;
		for (int i=0; i<numRows; i++)
			m[i].resize(numCols) ;
	}
}

/*! 
 *  \fn inline T Matrix<T> :: element (const int &row, const int &col) const
 *  \brief The function is used to resize the matrix.
 *  \param row an integer reference
 *  \param col an integer reference
 */
template <class T>
inline T Matrix<T> :: element (const int &row, const int &col) const
{
	if (row >= numRows || col >= numCols)
		error ("In accessing element of the matrix: index out of range ...") ;
	else return m[row][col] ;
}

/*! 
 *  \fn Matrix<T> Matrix<T> :: transpose (void)
 *  \brief Creates a transpose of the given matrix
 *  \return Transpose of a matrix
 */
template <class T>
Matrix<T> Matrix<T> :: transpose (void)
{
	int i,j ;
	Matrix<T> result (numCols,numRows) ;
	for (i=0; i<numCols; i++)
		for (j=0; j<numRows; j++)
			result.m[i][j] = m[j][i] ;
	return result ;
}

/*! 
 *  \fn Matrix<T>  Matrix<T> :: inverse (void)
 *  \brief creates the inverse of the given matrix, NO PIVOTING [not elegant]
 *  \return Inverse matrix object
 */
template <class T>
Matrix<T> Matrix<T> :: inverse (void)
{
	if (numRows != numCols)
		error ("In computing inverse: Not a square matrix!") ;
	else
	{
		int dimension = numRows ;
		Matrix<double> result = identity<double>(dimension) ;
		Matrix<double> source (*this) ;
		MyVector<double> row ;
		double current ;
		source.print() ;
		for (int i=0; i<dimension; i++)
		{
			current = source[i][(const int)i] ;
			//cout << current << endl ;
			//source[i].print() ;
			row = source[i] * 5 ;row.print() ;
			//source[i].copy(row) ;
			source[i] = row ;
			cout << source[i][0] << endl;
			source[i].print() ;
			//source.print() ;
		}
		return result ;
	}
}

/*! 
 *  \fn Matrix<T> :: ~Matrix<T>
 *  \brief Destructor function of Matrix class
 */
template <class T>
Matrix<T> :: ~Matrix()
{
	for (int i=0; i<numRows; i++)
		m[i].clear() ;
	m.clear() ;
}

#endif

/* OBSOLETE Add and Product Methods - replaced by overloading + and * respectively
 *  \fn Matrix<T>  Matrix<T> :: add (const Matrix<T> &other)
 *  \brief Adds two matrices 
 *  \param other a reference to a Matrix object
 *  \return Sum matrix
 *
template <class T>
Matrix<T>  Matrix<T> :: add (const Matrix<T> &other)
{
	if (numRows!=other.numRows || numCols!=other.numCols)
		error ("In adding matrices: dimensions mismatch!") ;
	else
	{
		int i,j ;
		Matrix<T> result (numRows,numCols) ;
		for (i=0; i<numRows; i++)
			for (j=0; j<numCols; j++)
				result.m[i][j] = m[i][j] + other.m[i][j] ;
		return result ;
	}
}
*/
/* 
 *  \fn Matrix<T>  Matrix<T> :: product (const Matrix<T> &other)
 *  \brief Multiplies two matrices 
 *  \param other a reference to a Matrix object
 *  \return Product matrix
 */
/*template <class T>
Matrix<T> Matrix<T> :: product (const Matrix<T> &other)
{
	if (numCols != other.numRows)
		error ("In computing matrix product: dimensions mismatch!") ;
	else
	{
		int i,j,k ;
		Matrix<T> result (numRows,other.numCols) ;
		for (i=0; i<numRows; i++)
			for (j=0; j<other.numCols; j++)
				for (k=0; k<numCols; k++)
					result.m[i][j] += m[i][k] * other.m[k][j] ;
		return result ;
	}
}
*/
