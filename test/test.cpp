#include <stdexcept>
#include "Vector.h"
#include "Matrix.h"
#include <iostream>
#include <cstdlib>

using namespace lcb ;

main()
{
	std::cout << "hello world" << std::endl ;
	Vector<int> v1(2,3) ;
	v1.print() ;
	Matrix<double> mat1(4) ;
	mat1.print() ;
	mat1.inverse().print() ;
}
