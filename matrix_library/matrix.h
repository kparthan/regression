/****************************************************************************************************
The Matrix class provides routines with the following functionalities:
   
Constructors:-
	1. Create a NULL vector 
		Matrix<type> A ;
        2. Create a rectangular matrix of zeroes 
		Matrix<type> A (ROWS,COLS) ;
	3. Create a square matrix of zeroes
		Matrix<type> A (DIMENSION) ;
        4. Create a rectangular matrix with all elements constant
                Matrix<type> A (constant,ROWS,COLS) ;
        5. Copy constructor
                Matrix<type> B ;
                Matrix<type> A (B) ;
  
Operator overloading:-
        1. operator =   		// copies the contents of an object
                Matrix<type> A,B ;
                B = A ;
        2. operator =   		// assigns all the elements of a matrix to a constant 'c'
                Matrix<type> A ;
                type c ;
                A = c ;
        3. operator []  		// subscripting - returns the i^{th} row in the matrix
		Matrix<type> A ;
                MyVector<type> row = A[i] ;
  
Member Functions:-
	   Matrix<type> A,B,C ;
        1. A.print()
		prints the details of the matrix (dimensions and elements) - for error checking
        2. A.nrows()
		gets the # of rows
	3. A.ncols()
		gets the # of cols
	4. A.changeDimensions (rows,cols)
		resizes the matrix to new dimensions and fills the matrix with zeroes
	5. A.element(i,j)
		to access the (i,j)^{th} element in the matrix
	6. A.transpose()
		generates the matrix transpose
	7. C = A.add(B)  
		creates a matrix C (= A + B)
	8. C = A.product(B)
		creates a matrix C (= A * B)
	9. A.inverse()
		generates the matrix inverse

Other functions:-
	1. A = zeros (ROWS,COLS)
		creates a rectangular matrix of zeros
	2. A = zeros (DIMENSION)
		creates a square matrix of zeros
	3. A = ones (ROWS,COLS)
		creates a rectangular matrix of ones
	4. A = ones (DIMENSION)	
		creates a square matrix of ones
	5. A = identity (ROWS,COLS)
		creates a rectangular identity matrix
	6. A = identity (DIMENSION)
		creates a square identity matrix

****************************************************************************************************/

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <cmath>
#include <vector>
#include "error.h"

template <class T>
class Matrix
{
	private:
		int numRows, numCols ;
		vector <vector<T> > m ;
	public:
		/* Constructors */
		Matrix () ;						// null matrix
		Matrix (int, int) ;					// rectangular zero matrix
		Matrix (int) ;						// square zero matrix
		Matrix (const T &, int , int) ;				// intialize to constant
		Matrix (const Matrix &) ;				// Copy constructor

		/* Overloading = and [] */
		Matrix & operator = (const Matrix &) ;			// assignment to a source matrix
		Matrix & operator = (const T &) ;			// assignment to a constant
		//inline Matrix<T> operator [] (const int) ;		// returns the i^{th} row

		/* Other sub-routines */
		void print() ;						// prints the elements of the matrix
		inline int nrows() const ;				// gets the number of rows
		inline int ncols() const ;				// gets the number if columns
		template <class U> 
		friend Matrix<U> identity (int, int) ;			// generates a rectangular/square identity matrix
		void changeDimensions (const int &, const int &) ;	// resizes the matrix
		inline T element (const int &, const int &) const ; 	// gets the (i,j)^{th} matrix element
		Matrix transpose () ;					// builds the transpose of the matrix
		Matrix add (const Matrix &) ;				// sums two matrices
		Matrix product (const Matrix &) ;			// multiplies two matrices
		Matrix inverse() ;					// computes the inverse of the matrix (if it is square)
		
		~Matrix() ;						// destructor
} ;

template <class T>
Matrix<T> zeros (int rows, int cols)
{
	Matrix<T> result (rows,cols) ;
	return result ;
}

template <class T>
Matrix<T> zeros (int dimension)
{
	return zeros<T>(dimension,dimension) ;
}

template <class T>
Matrix<T> ones (int rows, int cols)
{
	Matrix<T> result (1,rows,cols) ;
	return result ;
}

template <class T>
Matrix<T> ones (int dimension)
{
	return ones<T>(dimension,dimension) ;
}

template <class U>
Matrix<U> identity (int rows, int cols)
{
	int i,j ;
	Matrix<U> result (rows,cols) ;
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			if (i == j)
				result.m[i][j] = 1 ;
	return result ;
}

template <class U>
Matrix<U> identity (int dimension)
{
	return identity<U> (dimension,dimension) ;
}

template <class T>
Matrix<T> :: Matrix() : numRows(0), numCols(0), m(0)
{
}

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

template <class T>
Matrix<T> & Matrix<T> :: operator = (const T &a)
{
	int i,j ;
	for (i=0; i<numRows; i++)
		for (j=0; j<numCols; j++)
			m[i][j] = a ;
	return *this ;
}

/* 
template <class T>
MyVector<T> Matrix<T> :: operator [] (const int row)
{
	MyVector<T> myvec (numCols) ;
	for (int i=0; i<numCols; i++)
		myvec[i] = m[row][i] ;
	return myvec ;
}
*/

template <class T>
void Matrix<T> :: print(void)
{
	int i,j ;
	cout << "#rows = " << numRows << endl ;
	cout << "#cols = " << numCols << endl ;
	for (i=0; i<numRows; i++)
	{
		for (j=0; j<numCols; j++)
			cout << m[i][j] << " " ;
		cout << endl ;
	}
}

template <class T>
int Matrix<T> :: nrows() const
{
	return numRows ;
}

template <class T>
int Matrix<T> :: ncols() const
{
	return numCols ;
}

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

template <class T>
inline T Matrix<T> :: element (const int &row, const int &col) const
{
	if (row >= numRows || col >= numCols)
		error ("In accessing element of the matrix: index out of range ...") ;
	else return m[row][col] ;
}

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

template <class T>
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

template <class T>
Matrix<T> Matrix<T> :: inverse (void)
{
	if (numRows != numCols)
		error ("In computing inverse: Not a square matrix!") ;
	else
	{
	}
}

template <class T>
Matrix<T> :: ~Matrix()
{
	for (int i=0; i<numRows; i++)
		m[i].clear() ;
	m.clear() ;
}

#endif
