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
		void TESTING_overloadPlus() ;
		void TESTING_overloadPlusEqual() ;
		void TESTING_overloadMinus() ;
		void TESTING_overloadMinusEqual() ;
		void TESTING_overloadAsteriskConstant() ;
		void TESTING_modifyElement() ;
		void TESTING_l2Norm() ;
		void TESTING_normalize() ;
		void TESTING_angleBetweenVectors() ;
		void TESTING_isParallel() ;
		void TESTING_isPerpendicular() ;
		void TESTING_crossProduct() ;
} ;

void TEST_CLASS_FOR_MyVector :: TESTING_null (void)
{
	MyVector<int> vec1 ;
	vec1.print() ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_zero ()
{
	MyVector<int> vec1 (10) ;
	vec1.print() ;
	MyVector<int> vec2 (SIZE) ;
	vec2.print() ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_initializeToConstant (void)
{
	int c1 = 10 ;
	MyVector<int> vec1 (10,SIZE) ;
	vec1.print() ;

	double c2 = -100.77 ;
	MyVector<double> vec2 (c2,SIZE) ;
	vec2.print() ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_initializeToCppArray(void)
{
	int i ;
	int *intArray = new int [SIZE] ;
	for (i=0; i<SIZE; i++)
		intArray[i] = 2 * i ;
	MyVector <int> vec1 (intArray,SIZE) ;
	vec1.print() ;
	delete[] intArray ;
	
	int intArray2[] = {2,3,4,1,2} ;
	MyVector<int> vec2 (intArray2,SIZE) ;
	vec2.print() ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_copyConstructor(void)
{
	int n = 12 ;
	MyVector <int> source (n,SIZE) ;
	source.print() ;

	MyVector <int> vec1 (source) ;
	vec1.print() ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_overloadAssignmentOperator(void)
{
	int n1=39,n2=45 ;
	MyVector<int> source (n1,SIZE) ;
	source.print() ;

	MyVector<int> vec1 ; 
	vec1 = source ;
	vec1.print() ;

	MyVector<int> vec2(n2,3) ;
	vec2.print() ;
	vec2 = source ;
	vec2.print() ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_overloadAssignmentOperatorConstant(void)
{
	int n1=3,n2=15 ;
	MyVector<int> vec1(3,SIZE) ;
	vec1.print() ;
	vec1 = n2 ;
	vec1.print() ;
}

void TEST_CLASS_FOR_MyVector ::  TESTING_overloadSquareBracket(void)
{
	int n=3;
	MyVector<int> vec1 (n,2*SIZE) ;
	vec1.print() ;
	cout << vec1[4] << endl ;
	//cout << vec1[100] << endl ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_overloadAsterisk(void)
{
	int n1=2,n2=3;
	MyVector<int> vec1 (SIZE) ;
	vec1.print() ;
	MyVector<int> vec2 (n1*SIZE) ;
	vec2.print() ;
	MyVector<int> vec3 (n2,SIZE) ;
	vec3.print() ;
		
	int n3=10,n4=23;
	vec1 = n3 ;
	vec2 = 23 ;
	float dp = vec1 * vec3 ;
	cout << "dot product[1,3] = " << vec1 * vec3 << endl ;
	cout << "dot product[3,1] = " << vec3 * vec1 << endl ;
	cout << "dp = " << dp << endl ;
}

void TEST_CLASS_FOR_MyVector ::  TESTING_overloadPlus(void)
{
	int n1=3,n2=9 ;
	MyVector<int> vec1(3,SIZE) ;
	vec1.print() ;
	MyVector<int> vec2(n2,SIZE) ;
	vec2.print() ;
	MyVector<int> vec3 (vec1 + vec2); 
	vec3.print() ;
	(vec1 + vec2).print() ;
	vec3 = vec1 + (vec1 + vec2) * 5  ;
	vec3.print() ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_overloadPlusEqual(void)
{
	MyVector<int> vec1(3,SIZE) ;
	vec1.print() ;
	MyVector<int> vec2(10,SIZE) ;
	vec2.print() ;
	vec1 += vec2 ;
	vec1.print() ;
	MyVector<int> vec3 ;
	vec3.print() ;
	vec3 += vec2 ;
	vec3.print() ;
}

void TEST_CLASS_FOR_MyVector ::  TESTING_overloadMinus(void)
{
	MyVector<int> vec1(3,SIZE) ;
	vec1.print() ;
	MyVector<int> vec2(9,SIZE) ;
	vec2.print() ;
	MyVector<int> vec3 = vec1 - (vec2 + vec2) ;
	vec3.print() ;
	MyVector<int> vec4 = vec3 * -1 ;
	vec4.print() ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_overloadMinusEqual(void)
{
	MyVector<int> vec1(3,SIZE) ;
	vec1.print() ;
	MyVector<int> vec2(10,SIZE) ;
	vec2.print() ;
	vec1 -= vec2 ;
	vec1.print() ;
	MyVector<int> vec3 ;
	vec3.print() ;
	vec3 -= vec2 ;
	vec3.print() ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_overloadAsteriskConstant(void)
{
	MyVector<int> vec1(3,SIZE) ;
	vec1.print() ;
	MyVector<int> vec2 = vec1 * 5 ;
	vec2.print() ;
	MyVector<int> vec3(8,3) ;
	vec3.print() ;
	vec3 = vec1 * 10.6 ;
	vec3.print() ;
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

void TEST_CLASS_FOR_MyVector :: TESTING_normalize(void) 
{
	MyVector<float> vec1 (3,SIZE) ;
	vec1.print() ;
	MyVector<float> vec2 = vec1.normalize() ;
	vec2.print() ;
	vec1.print() ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_angleBetweenVectors(void) 
{
	MyVector<float> vec1 (5,SIZE) ;
	vec1.print() ;
	float *intArray = new float [SIZE] ;
	intArray[0] = 2 ;
	intArray[1] = 3 ;
	intArray[2] = 4 ;
	intArray[3] = 1 ;
	intArray[4] = 2 ;
	MyVector<float> vec2(intArray,SIZE) ;
	vec2.print() ;
	float angle = angleBetween(vec1,vec2) ;
	cout << "angle = " << angle << " degrees" << endl ;
	delete[] intArray ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_isParallel(void) 
{
	MyVector<float> vec1(5,SIZE) ;
	MyVector<float> vec2(8,SIZE) ;
	bool parallel = isParallel(vec1,vec2) ;
	cout << "parallel(1,2) = " << parallel << endl ;
	float *intArray = new float [SIZE] ;
	intArray[0] = 2 ;
	intArray[1] = 3 ;
	intArray[2] = 4 ;
	intArray[3] = 1 ;
	intArray[4] = 2 ;
	MyVector<float> vec3(intArray,SIZE) ;
	vec3.print() ;
	parallel = isParallel(vec1,vec3) ;
	cout << "parallel(1,3) = " << parallel << endl ;
	delete[] intArray ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_isPerpendicular(void) 
{
	MyVector<float> vec1(5,SIZE) ;
	MyVector<float> vec2(8,SIZE) ;
	bool perpendicular = isPerpendicular(vec1,vec2) ;
	cout << "perpendicular(1,2) = " << perpendicular << endl ;
	float *intArray = new float [SIZE] ;
	intArray[0] = 2 ;
	intArray[1] = 3 ;
	intArray[2] = 4 ;
	intArray[3] = 1 ;
	intArray[4] = 2 ;
	MyVector<float> vec3(intArray,SIZE) ;
	vec3.print() ;
	perpendicular = isPerpendicular(vec1,vec3) ;
	cout << "perpendicular(1,3) = " << perpendicular<< endl ;
	delete[] intArray ;
	MyVector<float> vec4 (SIZE) ;
	perpendicular = isPerpendicular(vec1,vec4) ;
	cout << "perpendicular(1,4) = " << perpendicular<< endl ;
}

void TEST_CLASS_FOR_MyVector :: TESTING_crossProduct(void) 
{
	float *intArray = new float [SIZE] ;
	intArray[0] = 2 ;
	intArray[1] = 3 ;
	intArray[2] = 4 ;
	MyVector<float> vec1(intArray,3) ;
	vec1.print() ;
	delete[] intArray ;
	MyVector<float> vec2(5,3) ;
	vec2.print() ;
	MyVector<float> vec3 = crossProduct(vec1,vec2) ;
	vec3.print() ;
}

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

	cout << "Testing overloading + ..." << endl ;
	test.TESTING_overloadPlus() ;
	
	cout << "Testing overloading += ..." << endl ;
	test.TESTING_overloadPlusEqual() ;

	cout << "Testing overloading - ..." << endl ;
	test.TESTING_overloadMinus() ;

	cout << "Testing overloading -= ..." << endl ;
	test.TESTING_overloadMinusEqual() ;

	cout << "Testing overloading * constant ..." << endl ;
	test.TESTING_overloadAsteriskConstant() ;

	cout << "Testing l2Norm() method ..." << endl ;
	test.TESTING_l2Norm() ;	

	cout << "Testing normalize() method ..." << endl ;
	test.TESTING_normalize() ;	

	cout << "Testing angle computation method ..." << endl ;
	test.TESTING_angleBetweenVectors() ;	
	
	cout << "Testing parallel condition ..." << endl ;
	test.TESTING_isParallel() ;	

	cout << "Testing perpendicular condition ..." << endl ;
	test.TESTING_isPerpendicular() ;	

	cout << "Testing cross product ..." << endl ;
	test.TESTING_crossProduct() ;	

	cout << MINIMUM << endl ;

	return 0 ;
}


