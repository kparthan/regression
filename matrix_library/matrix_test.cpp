#include "matrix.h"

class TEST_CLASS_FOR_Matrix
{
	public:
		void TESTING_null() ;
		void TESTING_rectZero() ;
		void TESTING_sqZero() ;
		void TESTING_initializeToConstant() ;
		void TESTING_copyConstructor() ;
		void TESTING_initializeToCppArray() ;
		void TESTING_initializeToCppArray2D() ;
		void TESTING_overloadAssignment() ;
		void TESTING_overloadAssignmentConstant() ;
		void TESTING_overloadSqBracketReturnRow() ;
		void TESTING_overloadPlus() ;
		void TESTING_overloadMinus() ;
		void TESTING_overloadAsterisk() ;
		void TESTING_overloadAsteriskConstant() ;
} ;

void TEST_CLASS_FOR_Matrix :: TESTING_null (void)
{
	Matrix<int> mat ;
	mat.print() ;
}

void TEST_CLASS_FOR_Matrix :: TESTING_rectZero (void)
{
	Matrix<int> mat(5,6) ;
	mat.print() ;
}

void TEST_CLASS_FOR_Matrix :: TESTING_sqZero (void)
{
	Matrix<int> mat(3) ;
	mat.print() ;
}

void TEST_CLASS_FOR_Matrix :: TESTING_initializeToConstant (void)
{
	Matrix<int> mat(7,3,4) ;
	mat.print() ;
}

void TEST_CLASS_FOR_Matrix :: TESTING_copyConstructor (void)
{
	Matrix<float> mat1(4.668,4,5) ;
	mat1.print() ;
	Matrix<float> mat2(mat1) ;
	mat2.print() ;
}

void TEST_CLASS_FOR_Matrix :: TESTING_initializeToCppArray (void)
{
	int a[9] = {1,2,4,6,0,-2,-4,9,0} ;
	Matrix<int> mat(a,4,3) ;
	mat.print() ;
}

void TEST_CLASS_FOR_Matrix :: TESTING_initializeToCppArray2D (void)
{
	int NUM_ROWS = 3, NUM_COLS = 2 ;
	int **a ;
	a = (int **)malloc(NUM_ROWS * sizeof(int *)) ;
	for (int i=0; i<NUM_ROWS; i++)
		a[i] = (int *)malloc(NUM_COLS * sizeof(int)) ;
	a[0][0] = 1 ; a[0][1] = -9 ;
	a[1][0] = -7 ; a[1][1] = 8 ;
	a[2][0] = 7 ; a[2][1] = 400 ;
	Matrix<int> mat1(a,3,2) ;
	mat1.print() ;
}

void TEST_CLASS_FOR_Matrix :: TESTING_overloadAssignment (void)
{
	Matrix<int> mat1 (3,2,2) ;
	mat1.print() ;
	Matrix<int> mat2(2,3) ;
	mat2.print() ;
	mat2 = mat1 ;
	mat2.print() ;
}

void TEST_CLASS_FOR_Matrix :: TESTING_overloadAssignmentConstant (void)
{
	Matrix<int> mat1(3,4) ;
	mat1.print() ;
	mat1 = 6 ;
	mat1.print() ;
	Matrix<float> mat2(3) ;
	mat2.print() ;
	mat2 = 7.5 ;
	mat2.print() ;
}

void TEST_CLASS_FOR_Matrix :: TESTING_overloadSqBracketReturnRow(void)
{
	Matrix<int> mat1(2,3,2) ;
	mat1.print() ;
	MyVector<int> row = mat1[0] ;
	row.print() ;
}

void TEST_CLASS_FOR_Matrix :: TESTING_overloadPlus(void)
{
	Matrix<int> mat1(2,3,2) ;
	mat1.print() ;
	Matrix<int> mat2 (44,3,2) ;
	mat2.print() ;
	Matrix<int> mat3 = mat1 + mat2 ;
	mat3.print() ;
}

void TEST_CLASS_FOR_Matrix :: TESTING_overloadMinus(void)
{
	Matrix<int> mat1(2,3,2) ;
	mat1.print() ;
	Matrix<int> mat2 (44,3,2) ;
	mat2.print() ;
	Matrix<int> mat3 = mat1 - mat2 ;
	mat3.print() ;
}

void TEST_CLASS_FOR_Matrix :: TESTING_overloadAsterisk(void)
{
	Matrix<int> mat1(2,3,2) ;
	mat1.print() ;
	Matrix<int> mat2(-3,2,3) ;
	mat2.print() ;
	Matrix<int> mat3 ;
	mat3 = mat2 * mat1 ;
	mat3.print() ;
}

void TEST_CLASS_FOR_Matrix :: TESTING_overloadAsteriskConstant(void)
{
	Matrix<int> mat1(2,3,4) ;
	mat1.print() ;
	Matrix<int> mat2 ;
	mat2.print() ;
	mat2 = (mat1 * 6)*-4 ;
	mat2.print() ;
	Matrix<int> mat3(2,2,2) ;
	mat3.print() ;
	Matrix<int> mat4(5,2,2) ;
	mat4.print() ;
	Matrix<int> mat5 ;
	mat5 = (mat4 * 5) * mat3 ;
	mat5.print() ;
	Matrix<int> mat6 = mat5 * 5 ;
	mat6.print() ;
}

int main(int argc, char **argv)
{
	TEST_CLASS_FOR_Matrix test ;

	cout << "Testing null matrix ..." << endl ;
	test.TESTING_null() ;

	cout << "Testing rectangular zero matrix ..." << endl ;
	test.TESTING_rectZero() ;

	cout << "Testing square zero matrix ..." << endl ;
	test.TESTING_sqZero() ;

	cout << "Testing initializtion to a constant ..." << endl ;
	test.TESTING_initializeToConstant() ;

	cout << "Testing copy constructor ..." << endl ;
	test.TESTING_copyConstructor() ;

	cout << "Testing initialization with C/C++ 1-D array ..." << endl ;
	test.TESTING_initializeToCppArray() ;
	
	cout << "Testing initialization with C/C++ 2-D array ..." << endl ;
	test.TESTING_initializeToCppArray2D() ;

	cout << "Testing overloading = ..." << endl ;
	test.TESTING_overloadAssignment() ;

	cout << "Testing overloading = constant ..." << endl ;
	test.TESTING_overloadAssignmentConstant() ;

	cout << "Testing overloading [] which returns row i ..." << endl ;
	test.TESTING_overloadSqBracketReturnRow() ;

	cout << "Testing overloading + ..." << endl ;
	test.TESTING_overloadPlus() ;

	cout << "Testing overloading - ..." << endl ;
	test.TESTING_overloadMinus() ;
	
	cout << "Testing overloading * ..." << endl ;
	test.TESTING_overloadAsterisk() ;

	cout << "Testing overloading * with constant ..." << endl ;
	test.TESTING_overloadAsteriskConstant() ;

	return 0 ;
}
