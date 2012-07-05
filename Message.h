/*!
 *	\file Message.h
 *	\details Implementation of Message class
 *	\author Parthan Kasarapu
 *	\date Modified: Thu 5 Jul 2012
 */

#ifndef MESSAGE_H
#define MESSAGE_H

#include "Data.h"
#include "Matrix.h"
#include <vector>

using namespace std ;

/*!
 *	\class Message
 *	\brief This is the Message class abstraction.
 *
 *	The class acts as an interface to compute the minimum message length.
 */
class Message
{
	private:
		int numFunctions ;
		vector<double> weights ;
		Data<double> xVals, yVals ;
	public:
		//! constructor
		template <class T>
		Message (int, Matrix<T>, Data<T>, Data<T>) ;
		//! computes the message length
 		void messageLength() ;
} ;

#endif

