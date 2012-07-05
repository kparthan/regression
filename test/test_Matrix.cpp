#define BOOST_TEST_MODULE test_Matrix
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <liblcb/liblcb.h>

BOOST_AUTO_TEST_CASE( Matrix_null )
{
  lcb::Matrix<int> mat ;
  BOOST_CHECK( mat.rows() == 0 );
  BOOST_CHECK( mat.columns() == 0 );
}

BOOST_AUTO_TEST_CASE( Matrix_rectZero )
{
  lcb::Matrix<int> mat(5,6) ;
  BOOST_CHECK( mat.rows() == 5 );
  BOOST_CHECK( mat.columns() == 6 );
  for(unsigned i = 0; i < mat.rows(); ++i){
    for(unsigned j = 0; j < mat.columns(); ++j){
      BOOST_CHECK( mat[i][j] == 0 );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_sqZero )
{
  lcb::Matrix<int> mat(3) ;
  BOOST_CHECK( mat.rows() == 3 );
  BOOST_CHECK( mat.columns() == 3 );
  for(unsigned i = 0; i < mat.rows(); ++i){
    for(unsigned j = 0; j < mat.columns(); ++j){
      BOOST_CHECK( mat[i][j] == 0 );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_rectZero_nonMember )
{
  lcb::Matrix<int> mat = lcb::zeros<int>(5,6) ;
  BOOST_CHECK( mat.rows() == 5 );
  BOOST_CHECK( mat.columns() == 6 );
  for(unsigned i = 0; i < mat.rows(); ++i){
    for(unsigned j = 0; j < mat.columns(); ++j){
      BOOST_CHECK( mat[i][j] == 0 );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_sqZero_nonMember)
{
  lcb::Matrix<int> mat = lcb::zeros<int>(3) ;
  BOOST_CHECK( mat.rows() == 3 );
  BOOST_CHECK( mat.columns() == 3 );
  for(unsigned i = 0; i < mat.rows(); ++i){
    for(unsigned j = 0; j < mat.columns(); ++j){
      BOOST_CHECK( mat[i][j] == 0 );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_rectOnes )
{
  lcb::Matrix<int> mat = lcb::ones<int>(5,6) ;
  BOOST_CHECK( mat.rows() == 5 );
  BOOST_CHECK( mat.columns() == 6 );
  for(unsigned i = 0; i < mat.rows(); ++i){
    for(unsigned j = 0; j < mat.columns(); ++j){
      BOOST_CHECK( mat[i][j] == 1 );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_sqOnes )
{
  lcb::Matrix<int> mat = lcb::ones<int>(5,5) ;
  BOOST_CHECK( mat.rows() == 5 );
  BOOST_CHECK( mat.columns() == 5 );
  for(unsigned i = 0; i < mat.rows(); ++i){
    for(unsigned j = 0; j < mat.columns(); ++j){
      BOOST_CHECK( mat[i][j] == 1 );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_identityRectangular)
{
  lcb::Matrix<int> mat1 = lcb::identity<int>(4,6) ;
  BOOST_CHECK( mat1.rows() == 4 );
  BOOST_CHECK( mat1.columns() == 6 );
  for(unsigned i = 0; i < mat1.rows(); ++i){
    for(unsigned j = 0; j < mat1.columns(); ++j){
      if(i == j)
	BOOST_CHECK( mat1[i][j] == 1 );
      else
	BOOST_CHECK( mat1[i][j] == 0 );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_identitySquare)
{
  lcb::Matrix<int> mat1 = lcb::identity<int>(4) ;
  BOOST_CHECK( mat1.rows() == 4 );
  BOOST_CHECK( mat1.columns() == 4 );
  for(unsigned i = 0; i < mat1.rows(); ++i){
    for(unsigned j = 0; j < mat1.columns(); ++j){
      if(i == j)
	BOOST_CHECK( mat1[i][j] == 1 );
      else
	BOOST_CHECK( mat1[i][j] == 0 );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_initializeToConstant )
{
  lcb::Matrix<int> mat(7,3,4) ;
  BOOST_CHECK( mat.rows() == 3 );
  BOOST_CHECK( mat.columns() == 4 );
  for(unsigned i = 0; i < mat.rows(); ++i){
    for(unsigned j = 0; j < mat.columns(); ++j){
      BOOST_CHECK( mat[i][j] == 7 );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_copyConstructor )
{
  float val = 4.668;
  lcb::Matrix<float> mat1(val,4,5) ;
  BOOST_CHECK( mat1.rows() == 4 );
  BOOST_CHECK( mat1.columns() == 5 );
  for(unsigned i = 0; i < mat1.rows(); ++i){
    for(unsigned j = 0; j < mat1.columns(); ++j){
      BOOST_CHECK( mat1[i][j] == val );
    }
  }
  lcb::Matrix<float> mat2 = mat1 ;
  BOOST_CHECK( mat2.rows() == 4 );
  BOOST_CHECK( mat2.columns() == 5 );
  for(unsigned i = 0; i < mat2.rows(); ++i){
    for(unsigned j = 0; j < mat2.columns(); ++j){
      BOOST_CHECK( mat2[i][j] == val );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_initializeToCppArray )
{
  int a[9] = {1,2,4,6,0,-2,-4,9,0} ;
  lcb::Matrix<int> mat(a,3,3) ;
  BOOST_CHECK( mat.rows() == 3 );
  BOOST_CHECK( mat.columns() == 3 );
  for(unsigned i = 0; i < mat.rows(); ++i){
    for(unsigned j = 0; j < mat.columns(); ++j){
      BOOST_CHECK( mat[i][j] == a[(i*3)+j] );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_initializeToCppArray2D )
{
  const int NUM_ROWS = 3, NUM_COLS = 2 ;
  int** a = new int * [NUM_ROWS];
  for(int i=0; i<NUM_ROWS; i++)
    a[i] = new int[NUM_COLS];
  a[0][0] = 1 ; a[0][1] = -9 ;
  a[1][0] = -7 ; a[1][1] = 8 ;
  a[2][0] = 7 ; a[2][1] = 400 ;

  lcb::Matrix<int> mat1(a,3,2) ;
  BOOST_CHECK( mat1.rows() == 3 );
  BOOST_CHECK( mat1.columns() == 2 );
  for(unsigned i = 0; i < mat1.rows(); ++i){
    for(unsigned j = 0; j < mat1.columns(); ++j){
      BOOST_CHECK( mat1[i][j] == a[i][j] );
    }
  }

  for (int i=0; i<NUM_ROWS;i++)
    delete[] a[i] ;
  delete[] a ;	
}

BOOST_AUTO_TEST_CASE( Matrix_overloadAssignment )
{
  lcb::Matrix<int> mat1 (3,2,2) ;
  lcb::Matrix<int> mat2 = mat1;
  BOOST_CHECK( mat2.rows() == 2 );
  BOOST_CHECK( mat2.columns() == 2 );
  for(unsigned i = 0; i < mat2.rows(); ++i){
    for(unsigned j = 0; j < mat2.columns(); ++j){
      BOOST_CHECK( mat2[i][j] == 3 );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_overloadAssignmentConstant )
{
  lcb::Matrix<int> mat1(3,4) ;
  mat1 = 6 ;
  BOOST_CHECK( mat1.rows() == 3 );
  BOOST_CHECK( mat1.columns() == 4 );
  for(unsigned i = 0; i < mat1.rows(); ++i){
    for(unsigned j = 0; j < mat1.columns(); ++j){
      BOOST_CHECK( mat1[i][j] == 6 );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_overloadPlus )
{
  lcb::Matrix<int> mat1(2,3,2) ;
  lcb::Matrix<int> mat2 (44,3,2) ;
  lcb::Matrix<int> mat3 = mat1 + mat2 ; 
  BOOST_CHECK( mat3.rows() == 3 );
  BOOST_CHECK( mat3.columns() == 2 );
  for(unsigned i = 0; i < mat3.rows(); ++i){
    for(unsigned j = 0; j < mat3.columns(); ++j){
      BOOST_CHECK( mat3[i][j] == 46 );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_overloadPlusEqual )
{
  lcb::Matrix<int> mat1(2,3,2) ;
  lcb::Matrix<int> mat2(45,3,2) ;
  mat1 += mat2 ;
  BOOST_CHECK( mat1.rows() == 3 );
  BOOST_CHECK( mat1.columns() == 2 );
  for(unsigned i = 0; i < mat1.rows(); ++i){
    for(unsigned j = 0; j < mat1.columns(); ++j){
      BOOST_CHECK( mat1[i][j] == 47 );
    }
  }
}


BOOST_AUTO_TEST_CASE( Matrix_overloadMinus )
{
  lcb::Matrix<int> mat1(2,3,2) ;
  lcb::Matrix<int> mat2 (44,3,2) ;
  lcb::Matrix<int> mat3 = mat1 - mat2 ;
  BOOST_CHECK( mat3.rows() == 3 );
  BOOST_CHECK( mat3.columns() == 2 );
  for(unsigned i = 0; i < mat3.rows(); ++i){
    for(unsigned j = 0; j < mat3.columns(); ++j){
      BOOST_CHECK( mat3[i][j] == -42 );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_overloadMinusEqual )
{
  lcb::Matrix<int> mat1(2,3,2) ;
  lcb::Matrix<int> mat2(45,3,2) ;
  mat1 -= mat2 ;
  BOOST_CHECK( mat1.rows() == 3 );
  BOOST_CHECK( mat1.columns() == 2 );
  for(unsigned i = 0; i < mat1.rows(); ++i){
    for(unsigned j = 0; j < mat1.columns(); ++j){
      BOOST_CHECK( mat1[i][j] == -43 );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_overloadAsterisk)
{
  lcb::Matrix<int> mat1(2,3,2) ;
  lcb::Matrix<int> mat2(-3,2,3) ;
  lcb::Matrix<int> mat3 = mat2 * mat1 ;
  BOOST_REQUIRE( mat3.rows() == 2 );
  BOOST_REQUIRE( mat3.columns() == 2 );
  for(unsigned i = 0; i < mat3.rows(); ++i){
    for(unsigned j = 0; j < mat3.columns(); ++j){
      BOOST_CHECK( mat3[i][j] == -18 );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_overloadAsteriskConstant )
{
  lcb::Matrix<int> mat1(2,3,4) ;
  lcb::Matrix<int> mat2 = mat1*2;
  BOOST_CHECK( mat2.rows() == 3 );
  BOOST_CHECK( mat2.columns() == 4 );
  for(unsigned i = 0; i < mat2.rows(); ++i){
    for(unsigned j = 0; j < mat2.columns(); ++j){
      BOOST_CHECK( mat2[i][j] == 4 );
    }
  }
}

//FAILURE: Should this not keep any previous values from before the resize?
BOOST_AUTO_TEST_CASE( Matrix_changeDimensions )
{
  lcb::Matrix<int> mat1(2,3,4) ;
  mat1.changeDimensions(4,5) ;
  BOOST_CHECK( mat1.rows() == 4 );
  BOOST_CHECK( mat1.columns() == 5 );
  for(unsigned i = 0; i < mat1.rows(); ++i){
    for(unsigned j = 0; j < mat1.columns(); ++j){
      BOOST_CHECK( mat1[i][j] == 0 );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_transpose )
{
  lcb::Matrix<int> mat1(2,3,4) ;
  lcb::Matrix<int> mat2 = mat1.transpose() ;
  BOOST_CHECK( mat2.rows() == 4 );
  BOOST_CHECK( mat2.columns() == 3 );
  for(unsigned i = 0; i < mat2.rows(); ++i){
    for(unsigned j = 0; j < mat2.columns(); ++j){
      BOOST_CHECK( mat2[i][j] == 2 );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_getColumnElements )
{
  lcb::Matrix<int> mat1(2,3,4) ;
  lcb::Vector<int> cvec = mat1.getColumnElements(1,1,2) ;
  BOOST_REQUIRE( cvec.length() == 2 );
  BOOST_CHECK( cvec[0] == 2 );
  BOOST_CHECK( cvec[1] == 2 );
}

BOOST_AUTO_TEST_CASE( Matrix_inverse )
{
  lcb::Matrix<double> inv ;
  int NUM_ROWS = 3, NUM_COLS = 3 ;
  double **a ;
  a = new double*[NUM_ROWS] ;
  for (int i=0; i<NUM_ROWS; i++)
    a[i] = new double[NUM_COLS] ;
  a[0][0] = 1.5 ; a[0][1] = -9.9 ; a[0][2] = 12 ;
  a[1][0] = -7 ; a[1][1] = 8.1 ; a[1][2] = 1.6 ;
  a[2][0] = 70 ; a[2][1] = 400 ; a[2][2] = 2 ;
  lcb::Matrix<double> mat1(a,NUM_ROWS,NUM_COLS) ;
  for (int i=0; i<NUM_ROWS; i++)
    delete[] a[i] ;
  delete[] a ;
  inv = mat1.inverse() ;
  lcb::Matrix<double> I = mat1 * inv;
  BOOST_CHECK( I.rows() == (unsigned)NUM_ROWS );
  BOOST_CHECK( I.columns() == (unsigned)NUM_COLS );
  for(unsigned i = 0; i < I.rows(); ++i){
    for(unsigned j = 0; j < I.columns(); ++j){
      if(i == j){
	BOOST_CHECK( (I[i][j] - 1) <= std::numeric_limits<double>::epsilon());
      }
      else
	BOOST_CHECK( I[i][j] <= std::numeric_limits<double>::epsilon() );
    }
  }
}

BOOST_AUTO_TEST_CASE( Matrix_determinant )
{
  int NUM_ROWS = 3, NUM_COLS = 3 ;
  double **a ;
  a = new double*[NUM_ROWS] ;
  for (int i=0; i<NUM_ROWS; i++)
    a[i] = new double[NUM_COLS] ;
  a[0][0] = 1.5 ; a[0][1] = -9.9 ; a[0][2] = 12 ;
  a[1][0] = -7 ; a[1][1] = 8.1 ; a[1][2] = 1.6 ;
  a[2][0] = 70 ; a[2][1] = 400 ; a[2][2] = 2 ;
  lcb::Matrix<double> mat1(a,NUM_ROWS,NUM_COLS) ;
  for (int i=0; i<NUM_ROWS; i++)
    delete[] a[i] ;
  delete[] a ;
  double det = mat1.determinant();
  BOOST_CHECK( (det + 42587.1) <= std::numeric_limits<double>::epsilon() );
	
  a = new double*[2] ;
  for (int i=0; i<2; i++)
    a[i] = new double[2] ;
  a[0][0] = 1 ; a[0][1] = 1 ;
  a[1][0] = 1 ; a[1][1] = 1 ;
  //determinant = 0
  lcb::Matrix<double> mat2(a,2,2) ;
  for (int i=0; i<2; i++)
    delete[] a[i] ;
  delete[] a ;
  BOOST_CHECK( mat2.determinant() == 0.0 );
}
