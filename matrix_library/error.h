/*!
 *  \file error.h
 *  \details Implementation of error ouput routine
 *  \author Parthan Kasarapu
 *  \version 1.0
 *  \date Wed 23 May 2012
*/

#ifndef _ERROR_H_
#define _ERROR_H_

#include <iostream>
#include <cstdlib>

using namespace std ;
/*#define ASSERT(x) \
                if (! (x)) \
                { \
                   cout << "ERROR!! Assert " << #x << " failed\n"; \
                   cout << " on line " << __LINE__  << "\n"; \
                   cout << " in file " << __FILE__ << "\n";  \
                }*/
/*! 
*  \fn inline void error (const char* errorMessage)
*  \brief This is an error handling function - informs the user of a potential breach in input parameters/variables.
*  \param errorText pointer to a character string
*/

inline void error (const char* errorMessage)
{
	cerr << "Run-time error ..." << endl ;
	cerr << errorMessage << endl ;
	cerr << "exiting ..." << endl ;
	exit(1) ;
}

#endif
