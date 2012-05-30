#include "myVector.h"

#define SIZE 5

class TEST_CLASS_FOR_MyVector
{
	public:
		void TESTING_null() ;
		void TESTING_zero() ;
		void TESTING_initializeToConstant() ;
		void TESTING_initializeToCppArray() ;
		void TESTING_copyConstructor() ;
		void TESTING_overloadAssignmentOperator() ;
		void TESTING_overloadAssignmentOperatorConstant() ;
		void TESTING_overloadSquareBracket() ;
		void TESTING_overloadAsterisk() ;
		void TESTING_modifyElement() ;
		void TESTING_l2Norm() ;
		//void TESTING_dotProduct() ;
} ;

void TEST_CLASS_FOR_MyVector :: TESTING_null (void)
{
	MyVector<int> vec1 ;
	vec1.print() ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_zero ()
{
	MyVector<int> vec1 (SIZE) ;
	vec1.print() ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_initializeToConstant (void)
{
	int c1 = 10 ;
	MyVector<int> vec1 (c1,SIZE) ;
	vec1.print() ;

	double c2 = -100.77 ;
	MyVector<double> vec2 (c2,SIZE) ;
	vec2.print() ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_initializeToCppArray(void)
{
	int i ;
	int *intArray = new int (SIZE) ;
	for (int i=0; i<SIZE; i++)
		intArray[i] = 2 * i ;
	MyVector <int> vec1 (intArray,SIZE) ;
	vec1.print() ;
	
	int intArray2[] = {2,3,4,1,2} ;
	MyVector<int> vec2 (intArray2,SIZE) ;
	vec2.print() ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_copyConstructor(void)
{
	MyVector <int> source (12,SIZE) ;
	source.print() ;

	MyVector <int> vec1 (source) ;
	vec1.print() ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_overloadAssignmentOperator(void)
{
	MyVector<int> source (39,SIZE) ;
	source.print() ;

	MyVector<int> vec1 = source ;
	vec1.print() ;

	MyVector<int> vec2(45,2*SIZE) ;
	vec2.print() ;
	vec2 = source ;
	vec2.print() ;
}
 
void TEST_CLASS_FOR_MyVector :: TESTING_overloadAssignmentOperatorConstant(void)
{
	MyVector<int> vec1(3,SIZE) ;
	vec1.print() ;
	vec1 = 15 ;
	vec1.print() ;
}
 
void TEST_CLASS_FOR_MyVector ::  TESTING_overloadSquareBracket(void)
{
	MyVector<int> vec1 (3,2*SIZE) ;
	vec1.print() ;
	cout << vec1[4] << endl ;
	//cout << vec1[100] << endl ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_overloadAsterisk(void)
{
	MyVector<int> vec1 (SIZE) ;
	vec1.print() ;
	MyVector<int> vec2 (2*SIZE) ;
	vec2.print() ;
	MyVector<int> vec3 (3,SIZE) ;
	vec3.print() ;
		
	vec1 = 10 ;
	vec2 = 23 ;
	cout << "dot product[1,3] = " << vec1 * vec3 << endl ;
	cout << "dot product[3,1] = " << vec3 * vec1 << endl ;
	//cout << "dot product[1,2] = " << vec1 * vec2 << endl ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_modifyElement(void)
{
	MyVector<int> vec (100,SIZE) ;
	vec.print() ;
	//vec[40] = 99 ;
	vec[0] = 576 ; vec[2] = 29 ;
	vec.print() ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_l2Norm(void) 
{
	MyVector<int> vec1 (3*SIZE) ;
	vec1.print() ;
	cout << "Magnitude = " << vec1.l2Norm() << endl ;
	vec1 = 22 ;
	vec1.print() ;
	cout << "Magnitude = " << vec1.l2Norm() << endl ;
}

/*
void TEST_CLASS_FOR_MyVector :: TESTING_dotProduct(void)
{
	MyVector<int> vec1 (SIZE) ;
	vec1.print() ;
	MyVector<int> vec2 (2*SIZE) ;
	vec2.print() ;
	MyVector<int> vec3 (3,SIZE) ;
	vec3.print() ;
		
	vec1 = 10 ;
	vec2 = 23 ;
	cout << "dot product[1,3] = " << vec1.dotProduct(vec3) << endl ;
	cout << "dot product[3,1] = " << vec3.dotProduct(vec1) << endl ;
	//cout << "dot product[1,2] = " << vec1.dotProduct(vec2) << endl ;
}
*/

int main (int argc, char **argv)
{
	TEST_CLASS_FOR_MyVector test ;
	
	cout << "Testing Null instance creation ..." << endl ;
	test.TESTING_null() ;

	cout << "Testing zero vector creation ..." << endl ;
	test.TESTING_zero() ;

	cout << "Testing initialization to constant value ..." << endl ;
	test.TESTING_initializeToConstant() ;

	cout << "Testing initialization to C/C++ array ..." << endl ;
	test.TESTING_initializeToCppArray() ;

	cout << "Testing copy constructor ..." << endl ;
	test.TESTING_copyConstructor() ;

	cout << "Testing overloading = ..." << endl ;
	test.TESTING_overloadAssignmentOperator() ;

	cout << "Testing overloading = constant ..." << endl ;
	test.TESTING_overloadAssignmentOperatorConstant() ;

	cout << "Testing overloading [] ..." << endl ;
	test.TESTING_overloadSquareBracket() ;

	cout << "Testing overloading * ..." << endl ;
	test.TESTING_overloadAsterisk() ;

	cout << "Testing modifyElement() method ..." << endl ;
	test.TESTING_modifyElement() ;

	cout << "Testing l2Norm() method ..." << endl ;
	test.TESTING_l2Norm() ;

	/*cout << "Testing dotProduct() method ..." << endl ;
	test.TESTING_dotProduct() ;*/		//OBSOLETE

	return 0 ;
}


