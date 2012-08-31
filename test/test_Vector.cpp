#define BOOST_TEST_MODULE test_Vector
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <liblcb/Vector.h>
#include <liblcb/Matrix.h>

#define SIZE 5

BOOST_AUTO_TEST_CASE( null_vector )
{
  std::cout << "Testing null matrix ..." << std::endl ;
  lcb::Vector<int> vec1 ;
  BOOST_CHECK( vec1.length() == 0 );
}

BOOST_AUTO_TEST_CASE( zero_vector )
{
  std::cout << "Testing zero vector ..." << std::endl ;
  lcb::Vector<int> vec1 (10) ;
  BOOST_CHECK( vec1.length() == 10 );
  for(int i = 0; i < 10; ++i){
    BOOST_CHECK( vec1[i] == 0 );
  }

  lcb::Vector<int> vec2 (SIZE) ;
  BOOST_CHECK( vec2.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec2[i] == 0 );
  }
}

BOOST_AUTO_TEST_CASE( initializeToConstant )
{
  std::cout << "Testing initialization to a constant ..." << std::endl ;
  int c1 = 10 ;
  lcb::Vector<int> vec1 (10,SIZE) ;
  BOOST_CHECK( vec1.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec1[i] == c1 );
  }

  double c2 = -100.77 ;
  lcb::Vector<double> vec2 (c2,SIZE) ;
  BOOST_CHECK( vec2.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec2[i] == c2 );
  }
}

BOOST_AUTO_TEST_CASE( initializeToArray )
{
  std::cout << "Testing initialization to an C/C++ array ..." << std::endl ;
  int i ;
  int *intArray = new int [SIZE] ;
  for(i = 0; i < SIZE; ++i)
    intArray[i] = 2 * i ;
  lcb::Vector<int> vec1 (intArray,SIZE,SIZE) ;
  BOOST_CHECK( vec1.length() == SIZE );
  for(i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec1[i] == (2 * i) );
  }
  delete[] intArray ;
	
  int intArray2[] = {2,3,4,1,2} ;
  lcb::Vector<int> vec2 (intArray2,SIZE,SIZE) ;
  BOOST_CHECK( vec2.length() == SIZE );
  for(i = 0; i < 5; ++i){
    BOOST_CHECK( vec2[i] == intArray2[i] );
  }
}

BOOST_AUTO_TEST_CASE( initializeToArrayOfSmallerSize )
{
  std::cout << "Testing initialization to an array of smaller size ..." << std::endl ;
  int i;
  int *intArray = new int[SIZE];
  for(i = 0; i < SIZE; ++i){
    intArray[i] = i;
  }
  lcb::Vector<int> vec1(intArray, SIZE*2,SIZE);
  BOOST_CHECK( vec1.length() == SIZE*2 );
  for(i = 0; i < SIZE; ++i){
    BOOST_CHECK( intArray[i] == vec1[i] );
  }
  for(; i < SIZE*2; ++i){
    BOOST_CHECK( vec1[i] == 0 );
  }
  delete[] intArray;
}

BOOST_AUTO_TEST_CASE( initializeToArrayOfLargerSize )
{
  std::cout << "Testing initialization to an array of larger size ..." << std::endl ;
  int i;
  int *intArray = new int[SIZE*2];
  for(i = 0; i < SIZE*2; ++i){
    intArray[i] = i;
  }
  lcb::Vector<int> vec1(intArray, SIZE,SIZE*2);
  BOOST_CHECK( vec1.length() == SIZE );
  for(i = 0; i < SIZE; ++i){
    BOOST_CHECK( intArray[i] == vec1[i] );
  }
  delete[] intArray;
}

BOOST_AUTO_TEST_CASE( copyConstructor )
{
  std::cout << "Testing copy constructor ..." << std::endl ;
  int n = 12 ;
  lcb::Vector<int> source (n,SIZE) ;
  BOOST_CHECK( source.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( source[i] == n );
  }

  lcb::Vector<int> vec1 (source) ;
  BOOST_CHECK( vec1.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec1[i] == n );
    BOOST_CHECK( vec1[i] == source[i] );
  }
}

BOOST_AUTO_TEST_CASE( overloadAssignmentOperator )
{
  std::cout << "Testing overload = operator ..." << std::endl ;
  int n1=39,n2=45 ;
  lcb::Vector<int> source (n1,SIZE) ;
  BOOST_CHECK( source.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( source[i] == n1 );
  }

  lcb::Vector<int> vec1 ; 
  vec1 = source ;
  BOOST_CHECK( vec1.length() == SIZE );
  BOOST_CHECK( source.length() == vec1.length() );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec1[i] == source[i] );
  }

  lcb::Vector<int> vec2(n2,SIZE-2) ;
  BOOST_CHECK( vec2.length() == SIZE-2 );
  for(int i = 0; i < SIZE-2; ++i){
    BOOST_CHECK( vec2[i] == n2 );
  }
  vec2 = source ;
  BOOST_CHECK( vec2.length() == source.length() );
  BOOST_REQUIRE( vec2.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec2[i] == source[i] );
  }
}

BOOST_AUTO_TEST_CASE( overloadAssignmentOperatorConstant )
{
  std::cout << "Testing overload = operator with a constant ..." << std::endl ;
  const int n1=3,n2=15 ;
  lcb::Vector<int> vec1(n1,SIZE) ;
  BOOST_CHECK( vec1.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec1[i] == n1 );
  }
  vec1 = n2 ;
  BOOST_CHECK( vec1.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec1[i] == n2 );
  }
}

BOOST_AUTO_TEST_CASE( overloadAsterisk )
{
  std::cout << "Testing overload * operator ..." << std::endl ;
  int n2=3;

  lcb::Vector<int> vec1 (SIZE) ;
  BOOST_CHECK( vec1.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec1[i] == 0 );
  }

  lcb::Vector<int> vec3 (n2,SIZE) ;
  BOOST_CHECK( vec3.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec3[i] == n2 );
  }
		
  vec1 = 10 ;
  BOOST_CHECK( (vec1 * vec3) == 150.0 );
  BOOST_CHECK( (vec1 * vec3) == (vec3 * vec1) );
}

BOOST_AUTO_TEST_CASE( dotProductDifferentSizes )
{
  std::cout << "Testing dot product different size vectors ..." << std::endl ;
  lcb::Vector<int> vec1 (SIZE);
  vec1 = 3;

  lcb::Vector<int> vec2 (SIZE*2);
  vec2 = 5;

  BOOST_CHECK_THROW((vec1 * vec2), std::length_error);
}

BOOST_AUTO_TEST_CASE(  overloadPlus )
{
  std::cout << "Testing overload + operator ..." << std::endl ;
  int n1=3,n2=9 ;
  lcb::Vector<int> vec1(n1,SIZE) ;
  BOOST_CHECK( vec1.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec1[i] == n1 );
  }

  lcb::Vector<int> vec2(n2,SIZE) ;
  BOOST_CHECK( vec2.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec2[i] == n2 );
  }

  lcb::Vector<int> vec3 (vec1 + vec2); 
  BOOST_CHECK( vec3.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec3[i] == (vec1[i] + vec2[i]) );
  }
}

BOOST_AUTO_TEST_CASE( overloadPlusEqual )
{
  std::cout << "Testing overload += operator ..." << std::endl ;
  lcb::Vector<int> vec1(3,SIZE) ;
  BOOST_CHECK( vec1.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec1[i] == 3 );
  }

  lcb::Vector<int> vec2(10,SIZE) ;
  BOOST_CHECK( vec2.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec2[i] == 10 );
  }

  vec1 += vec2 ;
  BOOST_CHECK( vec1.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec1[i] == 13 );
  }
}

BOOST_AUTO_TEST_CASE(  overloadMinus )
{
  std::cout << "Testing overload - operator ..." << std::endl ;
  lcb::Vector<int> vec1(3,SIZE) ;
  BOOST_CHECK( vec1.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec1[i] == 3 );
  }

  lcb::Vector<int> vec2(9,SIZE) ;
  BOOST_CHECK( vec2.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec2[i] == 9 );
  }

  lcb::Vector<int> vec3 = vec1 - vec2 ;
  BOOST_CHECK( vec3.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec3[i] == -6 );
  }
}

BOOST_AUTO_TEST_CASE( overloadMinusEqual )
{
  std::cout << "Testing overload -= operator ..." << std::endl ;
  lcb::Vector<int> vec1(3,SIZE) ;
  BOOST_CHECK( vec1.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec1[i] == 3 );
  }

  lcb::Vector<int> vec2(10,SIZE) ;
  BOOST_CHECK( vec2.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec2[i] == 10 );
  }

  vec1 -= vec2 ;
  BOOST_CHECK( vec1.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec1[i] == -7 );
  }
}

BOOST_AUTO_TEST_CASE( overloadAsteriskConstant )
{
  std::cout << "Testing overload * with a constant ..." << std::endl ;
  lcb::Vector<int> vec1(3,SIZE) ;
  BOOST_CHECK( vec1.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec1[i] == 3 );
  }

  lcb::Vector<int> vec2 = vec1 * 5 ;
  BOOST_CHECK( vec2.length() == SIZE );
  for(int i = 0; i < SIZE; ++i){
    BOOST_CHECK( vec2[i] == 15 );
  }
}

BOOST_AUTO_TEST_CASE( l2Norm ) 
{
  std::cout << "Testing l2Norm() ..." << std::endl ;
  lcb::Vector<int> vec1 (4) ;
  BOOST_CHECK( vec1.l2Norm() == 0.0 );

  vec1 = 2 ;
  BOOST_CHECK( vec1.l2Norm() == 4 );
}

BOOST_AUTO_TEST_CASE( normalize ) 
{
  std::cout << "Testing normalize() ..." << std::endl ;
  lcb::Vector<float> vec1 (2,4) ;
  BOOST_CHECK( vec1.l2Norm() == 4 );

  lcb::Vector<float> vec2 = vec1.normalize() ;
  BOOST_CHECK( vec2.l2Norm() == 1.0 );
  for(int i = 0; i < vec2.length(); ++i){
    BOOST_CHECK( vec2[i] == 0.5 );
  }
}

BOOST_AUTO_TEST_CASE( angleBetweenVectors ) 
{
  std::cout << "Testing angle between vectors ..." << std::endl ;
  float vec1INIT[2] = {0, 1};
  float vec2INIT[2] = {1, 0};
  lcb::Vector<float> vec1 (vec1INIT, 2,2);
  lcb::Vector<float> vec2 (vec2INIT, 2,2);

  BOOST_CHECK( angleBetween(vec1, vec2) ==
	       boost::math::constants::pi<float>()/2.0 );
}

BOOST_AUTO_TEST_CASE( isParallel ) 
{
  std::cout << "Testing whether vectors are parallel ..." << std::endl ;
  float vec1INIT[3] = {1, 1, 1};
  float vec2INIT[3] = {1, 1, 1};
  float vec3INIT[3] = {2, 2, 2};
  float vec4INIT[3] = {1, 0, 0};
  lcb::Vector<float> vec1(vec1INIT, 3, 3);
  lcb::Vector<float> vec2(vec2INIT, 3, 3);
  lcb::Vector<float> vec3(vec3INIT, 3, 3);
  lcb::Vector<float> vec4(vec4INIT, 3, 3);

  BOOST_CHECK( lcb::isParallel(vec1, vec2) );
  BOOST_CHECK( lcb::isParallel(vec2, vec3) );
  BOOST_CHECK( !lcb::isParallel(vec2, vec4) );

  float vec5INIT[2] = {1, 0};
  float vec6INIT[2] = {-1, 0};
  lcb::Vector<float> vec5 (vec5INIT, 2, 2);
  lcb::Vector<float> vec6 (vec6INIT, 2, 2);
  BOOST_CHECK( !lcb::isParallel(vec5, vec6) );
}

BOOST_AUTO_TEST_CASE( isPerpendicular ) 
{
  std::cout << "Testing whether vectors are perpendicular ..." << std::endl ;
  float vec1INIT[2] = {0, 1};
  float vec2INIT[2] = {1, 0};
  float vec3INIT[2] = {0, -1};

  lcb::Vector<float> vec1(vec1INIT, 2, 2);
  lcb::Vector<float> vec2(vec2INIT, 2, 2);
  lcb::Vector<float> vec3(vec3INIT, 2, 2);
  BOOST_CHECK( lcb::isPerpendicular(vec1, vec2) );
  BOOST_CHECK( !lcb::isPerpendicular(vec1, vec3) );

  float vec4INIT[3] = {0, 0, 1};
  float vec5INIT[3] = {0, 0, 2};
  float vec6INIT[3] = {0, 1, 0};
  lcb::Vector<float> vec4(vec4INIT, 3, 3);
  lcb::Vector<float> vec5(vec5INIT, 3, 3);
  lcb::Vector<float> vec6(vec6INIT, 3, 3);
  BOOST_CHECK( lcb::isPerpendicular(vec4, vec6) );
  BOOST_CHECK( !lcb::isPerpendicular(vec4, vec5) );
}


BOOST_AUTO_TEST_CASE( crossProduct ) 
{
  std::cout << "Testing cross product(3D) ..." << std::endl ;
  float *intArray = new float [SIZE] ;
  intArray[0] = 2 ;
  intArray[1] = 3 ;
  intArray[2] = 4 ;
  lcb::Vector<float> vec1(intArray,3,3) ;
  delete[] intArray ;
  lcb::Vector<float> vec2(5,3) ;
  lcb::Vector<float> vec3 = lcb::crossProduct(vec1,vec2) ;
  BOOST_CHECK( vec3[0] == -5 );
  BOOST_CHECK( vec3[1] == 10 );
  BOOST_CHECK( vec3[2] == -5 );
}
