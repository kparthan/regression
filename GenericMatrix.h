/*!
 * Generic Matrix: An abstract Matrix class
 * Copyright (C) 2012  James Collier
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LIBLCB_GENERIC_MATRIX_H__
#define LIBLCB_GENERIC_MATRIX_H__

#include <utility>
#include <cstdlib>

//namespace lcb {
  /*!
   *  \class GenericMatrix
   *  \brief This is a GenericMatrix class
   *  
   *  The class acts as an interface for inheritance for Matrix class.
   */
  class GenericMatrix {
  public:
    /*!
     *  \fn GenericMatrix()
     *  \brief Null constructor
     *  \return a null instance of GenericMatrix object 	
     */
    GenericMatrix(): _number_of_rows(0),
		     _number_of_columns(0)
    {}

    /*!
     *  \fn GenericMatrix(const std::size_t rows, const std::size_t cols)
     *  \brief The constructor creates an object with the specified
     *  number of rows and columns
     *  \param rows an integer
     *  \param cols an integer
     *  \return a new instance of GenericMatrix object
     */
    GenericMatrix(const std::size_t rows,
		  const std::size_t cols)
      : _number_of_rows(rows), _number_of_columns(cols)
    {}

    /*!
     *  \fn std::size_t rows() const
     *  \brief This function is used to returns the number of rows in 
     *  the matrix
     *  \return number of rows
     */
    std::size_t
      rows() const
    { return _number_of_rows; }

    /*!
     *  \fn std::size_t columns() const
     *  \brief This function is used to returns the number of columns in 
     *  the matrix
     *  \return number of columns
     */
    std::size_t
      columns() const
    { return _number_of_columns; }

    /*!
     *  \fn const std::pair<std::size_t, std::size_t> dimensions() const
     *  \brief The function is used to return the dimension of the matrix
     *  as a 'pair'
     *  \return the number of rows and columns
     */
    const std::pair<std::size_t, std::size_t>
      dimensions() const
      {
	return std::pair<std::size_t, std::size_t>
	  (_number_of_rows, _number_of_columns);
      }

  protected:
    std::size_t _number_of_rows;
    std::size_t _number_of_columns;
  };
//}

#endif /* LIBLCB_GENERIC_MATRIX_H__ */

