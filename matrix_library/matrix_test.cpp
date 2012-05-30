#include "matrix.h"

class TEST_CLASS_FOR_Matrix
{
	public:
		void TESTING_null() ;
		void TESTING_rectZero() ;
		void TESTING_sqZero() ;
		void TESTING_initializeToConstant() ;
		void TESTING_copyConstructor() ;
		void TESTING_overloadAssignment() ;
		void TESTING_overloadAssignmentConstant() ;
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

	cout << "Testing overloading = ..." << endl ;
	test.TESTING_overloadAssignment() ;

	cout << "Testing overloading = constant ..." << endl ;
	test.TESTING_overloadAssignmentConstant() ;

	return 0 ;
}
