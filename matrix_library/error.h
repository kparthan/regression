#ifndef _ERROR_H_
#define _ERROR_H_

#include <iostream>
#include <cstdlib>

using namespace std ;

inline void error(const char* errorText)
{
	cerr << "Run-time error ..." << endl ;
	cerr << errorText << endl ;
	cerr << "exiting ..." << endl ;
	exit(1) ;
}

#endif
