#include "Vector.h"

#define SIZE 5

using namespace std ;

class TEST_CLASS_FOR_Vector
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

void TEST_CLASS_FOR_Vector :: TESTING_null (void)
{
	Vector<int> vec1 ;
	vec1.print() ;
}

void TEST_CLASS_FOR_Vector :: TESTING_zero ()
{
	Vector<int> vec1 (10) ;
	vec1.print() ;
	Vector<int> vec2 (SIZE) ;
	vec2.print() ;
}

void TEST_CLASS_FOR_Vector :: TESTING_initializeToConstant (void)
{
	int c1 = 10 ;
	Vector<int> vec1 (10,SIZE) ;
	vec1.print() ;

	double c2 = -100.77 ;
	Vector<double> vec2 (c2,SIZE) ;
	vec2.print() ;
}

void TEST_CLASS_FOR_Vector :: TESTING_initializeToCppArray(void)
{
	int i ;
	int *intArray = new int [SIZE] ;
	for (i=0; i<SIZE; i++)
		intArray[i] = 2 * i ;
	Vector <int> vec1 (intArray,SIZE,SIZE) ;
	vec1.print() ;
	delete[] intArray ;
	
	int intArray2[] = {2,3,4,1,2} ;
	Vector<int> vec2 (intArray2,SIZE,SIZE) ;
	vec2.print() ;
}

void TEST_CLASS_FOR_Vector :: TESTING_copyConstructor(void)
{
	int n = 12 ;
	Vector <int> source (n,SIZE) ;
	source.print() ;

	Vector <int> vec1 (source) ;
	vec1.print() ;
}

void TEST_CLASS_FOR_Vector :: TESTING_overloadAssignmentOperator(void)
{
	int n1=39,n2=45 ;
	Vector<int> source (n1,SIZE) ;
	source.print() ;

	Vector<int> vec1 ; 
	vec1 = source ;
	vec1.print() ;

	Vector<int> vec2(n2,3) ;
	vec2.print() ;
	vec2 = source ;
	vec2.print() ;
}

void TEST_CLASS_FOR_Vector :: TESTING_overloadAssignmentOperatorConstant(void)
{
	int n1=3,n2=15 ;
	Vector<int> vec1(3,SIZE) ;
	vec1.print() ;
	vec1 = n2 ;
	vec1.print() ;
}

void TEST_CLASS_FOR_Vector ::  TESTING_overloadSquareBracket(void)
{
	int n=3;
	Vector<int> vec1 (n,2*SIZE) ;
	vec1.print() ;
	cout << vec1[4] << endl ;
	//cout << vec1[100] << endl ;
}

void TEST_CLASS_FOR_Vector :: TESTING_overloadAsterisk(void)
{
	int n1=2,n2=3;
	Vector<int> vec1 (SIZE) ;
	vec1.print() ;
	Vector<int> vec2 (n1*SIZE) ;
	vec2.print() ;
	Vector<int> vec3 (n2,SIZE) ;
	vec3.print() ;
		
	int n3=10,n4=23;
	vec1 = n3 ;
	vec2 = 23 ;
	float dp = vec1 * vec3 ;
	cout << "dot product[1,3] = " << vec1 * vec3 << endl ;
	cout << "dot product[3,1] = " << vec3 * vec1 << endl ;
	cout << "dp = " << dp << endl ;
}

void TEST_CLASS_FOR_Vector ::  TESTING_overloadPlus(void)
{
	int n1=3,n2=9 ;
	Vector<int> vec1(3,SIZE) ;
	vec1.print() ;
	Vector<int> vec2(n2,SIZE) ;
	vec2.print() ;
	Vector<int> vec3 (vec1 + vec2); 
	vec3.print() ;
	(vec1 + vec2).print() ;
	vec3 = vec1 + (vec1 + vec2) * 5  ;
	vec3.print() ;
}

void TEST_CLASS_FOR_Vector :: TESTING_overloadPlusEqual(void)
{
	Vector<int> vec1(3,SIZE) ;
	vec1.print() ;
	Vector<int> vec2(10,SIZE) ;
	vec2.print() ;
	vec1 += vec2 ;
	vec1.print() ;
	Vector<int> vec3 ;
	vec3.print() ;
	vec3 += vec2 ;
	vec3.print() ;
}

void TEST_CLASS_FOR_Vector ::  TESTING_overloadMinus(void)
{
	Vector<int> vec1(3,SIZE) ;
	vec1.print() ;
	Vector<int> vec2(9,SIZE) ;
	vec2.print() ;
	Vector<int> vec3 = vec1 - (vec2 + vec2) ;
	vec3.print() ;
	Vector<int> vec4 = vec3 * -1 ;
	vec4.print() ;
}

void TEST_CLASS_FOR_Vector :: TESTING_overloadMinusEqual(void)
{
	Vector<int> vec1(3,SIZE) ;
	vec1.print() ;
	Vector<int> vec2(10,SIZE) ;
	vec2.print() ;
	vec1 -= vec2 ;
	vec1.print() ;
	Vector<int> vec3 ;
	vec3.print() ;
	vec3 -= vec2 ;
	vec3.print() ;
}

void TEST_CLASS_FOR_Vector :: TESTING_overloadAsteriskConstant(void)
{
	Vector<int> vec1(3,SIZE) ;
	vec1.print() ;
	Vector<int> vec2 = vec1 * 5 ;
	vec2.print() ;
	Vector<int> vec3(8,3) ;
	vec3.print() ;
	vec3 = vec1 * 10.6 ;
	vec3.print() ;
}

void TEST_CLASS_FOR_Vector :: TESTING_l2Norm(void) 
{
	Vector<int> vec1 (3*SIZE) ;
	vec1.print() ;
	cout << "Magnitude = " << vec1.l2Norm() << endl ;
	vec1 = 22 ;
	vec1.print() ;
	cout << "Magnitude = " << vec1.l2Norm() << endl ;
}

void TEST_CLASS_FOR_Vector :: TESTING_normalize(void) 
{
	Vector<float> vec1 (3,SIZE) ;
	vec1.print() ;
	Vector<float> vec2 = vec1.normalize() ;
	vec2.print() ;
	vec1.print() ;
}

void TEST_CLASS_FOR_Vector :: TESTING_angleBetweenVectors(void) 
{
	Vector<float> vec1 (5,SIZE) ;
	vec1.print() ;
	float *intArray = new float [SIZE] ;
	intArray[0] = 2 ;
	intArray[1] = 3 ;
	intArray[2] = 4 ;
	intArray[3] = 1 ;
	intArray[4] = 2 ;
	Vector<float> vec2(intArray,SIZE,SIZE) ;
	vec2.print() ;
	float angle = angleBetween(vec1,vec2) ;
	cout << "angle = " << angle << " degrees" << endl ;
	delete[] intArray ;
}

void TEST_CLASS_FOR_Vector :: TESTING_isParallel(void) 
{
	Vector<float> vec1(5,SIZE) ;
	Vector<float> vec2(8,SIZE) ;
	bool parallel = isParallel(vec1,vec2) ;
	cout << "parallel(1,2) = " << parallel << endl ;
	float *intArray = new float [SIZE] ;
	intArray[0] = 2 ;
	intArray[1] = 3 ;
	intArray[2] = 4 ;
	intArray[3] = 1 ;
	intArray[4] = 2 ;
	Vector<float> vec3(intArray,SIZE,SIZE) ;
	vec3.print() ;
	parallel = isParallel(vec1,vec3) ;
	cout << "parallel(1,3) = " << parallel << endl ;
	delete[] intArray ;
}

void TEST_CLASS_FOR_Vector :: TESTING_isPerpendicular(void) 
{
	Vector<float> vec1(5,SIZE) ;
	Vector<float> vec2(8,SIZE) ;
	bool perpendicular = isPerpendicular(vec1,vec2) ;
	cout << "perpendicular(1,2) = " << perpendicular << endl ;
	float *intArray = new float [SIZE] ;
	intArray[0] = 2 ;
	intArray[1] = 3 ;
	intArray[2] = 4 ;
	intArray[3] = 1 ;
	intArray[4] = 2 ;
	Vector<float> vec3(intArray,SIZE,SIZE) ;
	vec3.print() ;
	perpendicular = isPerpendicular(vec1,vec3) ;
	cout << "perpendicular(1,3) = " << perpendicular << endl ;
	delete[] intArray ;
	Vector<float> vec4 (SIZE) ;
	perpendicular = isPerpendicular(vec1,vec4) ;
	cout << "perpendicular(1,4) = " << perpendicular << endl ;
}

void TEST_CLASS_FOR_Vector :: TESTING_crossProduct(void) 
{
	float *intArray = new float [SIZE] ;
	intArray[0] = 2 ;
	intArray[1] = 3 ;
	intArray[2] = 4 ;
	Vector<float> vec1(intArray,3,3) ;
	vec1.print() ;
	delete[] intArray ;
	Vector<float> vec2(5,3) ;
	vec2.print() ;
	Vector<float> vec3 = crossProduct(vec1,vec2) ;
	vec3.print() ;
}

int main (int argc, char **argv)
{
	TEST_CLASS_FOR_Vector test ;
	
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

	return 0 ;
}


