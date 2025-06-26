/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "GIDI.hpp"

namespace GIDI {

/*! \class Matrix
 * This class stores a mathematical matrix and has methods that perform several matrix operations (e.g., addition, subtraction).
 */

/* *********************************************************************************************************//**
 *
 * @param a_rows            [in]    Number of rows of the matrix.
 * @param a_columns         [in]    Number of columns of the matrix.
 ***********************************************************************************************************/

Matrix::Matrix( std::size_t a_rows, std::size_t a_columns ) {

    m_matrix.resize( a_rows );

    for( std::size_t i1 = 0; i1 < a_rows; ++i1 ) m_matrix[i1].resize( a_columns );
}

/* *********************************************************************************************************//**
 *
 * @param a_matrix          [in]    Matrix to copy.
 ***********************************************************************************************************/

Matrix::Matrix( Matrix const &a_matrix ) {

    std::size_t rows = a_matrix.size( );

    m_matrix.resize( rows );

    for( std::size_t i1 = 0; i1 < rows; ++i1 ) m_matrix[i1] = a_matrix[i1];
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Matrix::~Matrix( ) {

}

/* *********************************************************************************************************//**
 * The assignment operator. This method sets the members of *this* to those of *a_rhs*.
 *
 * @param a_rhs                     [in]    Instance whose member are used to set the members of *this*.
 *
 * @return                                  A reference to the updated Matrix instance.
 ***********************************************************************************************************/

Matrix &Matrix::operator=( Matrix const &a_rhs ) {

    if( this != &a_rhs ) {
        m_matrix = a_rhs.matrix( );
    }

    return( *this );
}

/* *********************************************************************************************************//**
 * Returns a new Matrix whose cells are *this* plus *a_value*.
 *
 * @param a_value       [in]    The value to add to each cell.
 * @return                      New Matrix whose cells are *this* plus *a_value*.
 ***********************************************************************************************************/

Matrix Matrix::operator+( double a_value ) const {

    Matrix gidiMatrix( *this );

    gidiMatrix += a_value;
    return( gidiMatrix );
}

/* *********************************************************************************************************//**
 * Adds *a_value* to each cell of *this*.
 *
 * @param a_value       [in]    The value to add to each cell.
 * @return                      Returns reference to *this*.
 ***********************************************************************************************************/

Matrix &Matrix::operator+=( double a_value ) {

    for( std::vector<Vector>::iterator iter = m_matrix.begin( ); iter < m_matrix.end( ); ++iter ) *iter += a_value;

    return( *this );
}

/* *********************************************************************************************************//**
 * Adds two Matrices.
 *
 * @param a_rhs         [in]    Matrix to add to *this*.
 * @return                      New Matrix that is the matrix sum of *this* and *a_rhs*.
 ***********************************************************************************************************/

Matrix Matrix::operator+( Matrix const &a_rhs ) const {

    Matrix gidiMatrix( *this );

    gidiMatrix += a_rhs;
    return( gidiMatrix );
}

/* *********************************************************************************************************//**
 * Adds *a_rhs* to *this*.
 *
 * @param a_rhs         [in]    Matrix to add to *this*.
 * @return                      Returns reference to *this*.
 ***********************************************************************************************************/

Matrix &Matrix::operator+=( Matrix const &a_rhs ) {

    std::size_t i1, rhs_size = a_rhs.size( );

    if( rhs_size == 0 ) return( *this );                           // Do nothing if rhs is empty.

    if( size( ) == 0 ) {
        for( i1 = 0; i1 < rhs_size; ++i1 ) m_matrix.push_back( Vector( a_rhs.numberOfColumns( ) ) );
    }

    if( size( ) != a_rhs.size( ) ) throw Exception( "matrix sizes differ." );
    if( m_matrix[0].size( ) != a_rhs[0].size( ) ) throw Exception( "matrix colums numbers differ." );

    i1 = 0;
    for( std::vector<Vector>::iterator iter = m_matrix.begin( ); iter < m_matrix.end( ); ++iter, ++i1 ) *iter += a_rhs[i1];

    return( *this );
}

/* *********************************************************************************************************//**
 * Returns a new Matrix whose cells are *this* minus *a_value*.
 *
 * @param a_value       [in]    The value to subtract from each cell.
 * @return                      New Matrix whose cells are *this* plus *a_value*.
 ***********************************************************************************************************/

Matrix Matrix::operator-( double a_value ) const {

    Matrix gidiMatrix( *this );

    gidiMatrix -= a_value;
    return( gidiMatrix );
}

/* *********************************************************************************************************//**
 * Subtracts *a_value* from each cell of *this*.
 *
 * @param a_value       [in]    The value to subtract from each cell.
 * @return                      Returns reference to *this*.
 ***********************************************************************************************************/

Matrix &Matrix::operator-=( double a_value ) {

    for( std::vector<Vector>::iterator iter = m_matrix.begin( ); iter < m_matrix.end( ); ++iter ) *iter -= a_value;

    return( *this );
}

/* *********************************************************************************************************//**
 * Subtracts *a_rhs* from *this*.
 *
 * @param a_rhs         [in]    Matrix to subtract from *this*.
 * @return                      New Matrix that is *this* minus *a_rhs*.
 ***********************************************************************************************************/

Matrix Matrix::operator-( Matrix const &a_rhs ) const {

    Matrix gidiMatrix( *this );

    gidiMatrix -= a_rhs;
    return( gidiMatrix );
}

/* *********************************************************************************************************//**
 * Subtracts *a_rhs* to *this*.
 *
 * @param a_rhs         [in]    Matrix to subtract from *this*.
 * @return                      Returns reference to *this*.
 ***********************************************************************************************************/

Matrix &Matrix::operator-=( Matrix const &a_rhs ) {

    std::size_t i1, rhs_size = a_rhs.size( );

    if( rhs_size == 0 ) return( *this );                           // Do nothing if rhs is empty.

    if( size( ) == 0 ) {
        for( i1 = 0; i1 < rhs_size; ++i1 ) m_matrix.push_back( Vector( a_rhs.numberOfColumns( ) ) );
    }

    if( size( ) != a_rhs.size( ) ) throw Exception( "matrix sizes differ." );
    if( m_matrix[0].size( ) != a_rhs[0].size( ) ) throw Exception( "matrix colums numbers differ." );

    i1 = 0;
    for( std::vector<Vector>::iterator iter = m_matrix.begin( ); iter < m_matrix.end( ); ++iter, ++i1 ) *iter -= a_rhs[i1];

    return( *this );
}

/* *********************************************************************************************************//**
 * Returns a new Matrix whose cells are *this* multiplied by *a_value*.
 *
 * @param a_value       [in]    The value to multiply each cell by.
 * @return                      New Matrix whose cells are *this* multiply by *a_value*.
 ***********************************************************************************************************/

Matrix Matrix::operator*( double a_value ) const {

    Matrix gidiMatrix( *this );

    gidiMatrix *= a_value;
    return( gidiMatrix );
}

/* *********************************************************************************************************//**
 * Multiplies each cell of *this* by *a_value*.
 *
 * @param a_value       [in]    The value to multiply each cell by.
 * @return                      Returns reference to *this*.
 ***********************************************************************************************************/

Matrix &Matrix::operator*=( double a_value ) {

    for( std::vector<Vector>::iterator iter = m_matrix.begin( ); iter < m_matrix.end( ); ++iter ) *iter *= a_value;

    return( *this );
}

/* *********************************************************************************************************//**
 * Returns a new Matrix whose cells are *this* divided by *a_value*.
 *
 * @param a_value       [in]    The value to divide each cell by.
 * @return                      New Matrix whose cells are *this* divided by *a_value*.
 ***********************************************************************************************************/

Matrix Matrix::operator/( double a_value ) const {

    Matrix gidiMatrix( *this );

    gidiMatrix /= a_value;
    return( gidiMatrix );
}

/* *********************************************************************************************************//**
 * Divides each cell of *this* by *a_value*.
 *
 * @param a_value       [in]    The value to divide each cell by.
 * @return                      Returns reference to *this*.
 ***********************************************************************************************************/

Matrix &Matrix::operator/=( double a_value ) {

    if( a_value == 0 ) throw Exception( "divide by zero." );
    for( std::vector<Vector>::iterator iter = m_matrix.begin( ); iter < m_matrix.end( ); ++iter ) *iter /= a_value;

    return( *this );
}

/* *********************************************************************************************************//**
 * Returns the number of columns of *this*.
 *
 * @return                      The number of columns of *this*.
 ***********************************************************************************************************/

std::size_t Matrix::numberOfColumns( ) const {

    if( size( ) == 0 ) return( 0 );
    return( m_matrix[0].size( ) );
}

/* *********************************************************************************************************//**
 * Adds a row to *this*.
 *
 * @param a_vector          [in]    The Vector to add to *this* as another row.
 ***********************************************************************************************************/

void Matrix::push_back( Vector const &a_vector ) {

    if( size( ) > 0 ) {
        if( (*this)[0].size( ) != a_vector.size( ) ) throw Exception( "matrix::push_back: size different" );
    }
    m_matrix.push_back( a_vector );
}

/* *********************************************************************************************************//**
 * Transposes the cells of *this*.
 ***********************************************************************************************************/

Matrix Matrix::transpose( ) {

    std::size_t __numberOfColumns( numberOfColumns( ) );
    Matrix __matrix( __numberOfColumns, size( ) );

    for( std::size_t i1 = 0; i1 < size( ); ++i1 ) {
        for( std::size_t i2 = 0; i2 < __numberOfColumns; ++i2 ) __matrix( i2, i1, (*this)[i1][i2] );
    }
    return( __matrix );
}

/* *********************************************************************************************************//**
 * Reverse the rows of *this*.
 ***********************************************************************************************************/

void Matrix::reverse( ) {

    std::size_t i2 = size( ), n_2 = i2 / 2;

    for( std::size_t i1 = 0; i1 < i2; ++i1 ) m_matrix[i1].reverse( );
    --i2;
    for( std::size_t i1 = 0; i1 < n_2; ++i1, --i2 ) {
        Vector temp = m_matrix[i1];

        m_matrix[i1] = m_matrix[i2];
        m_matrix[i2] = temp;
    }
}

/* *********************************************************************************************************//**
 * Prints the contents of *this* by calling the *print* method of each row with **a_prefixForRow** as an argument.
 *
 * @param a_prefixForRow    [in]    Argument passed to each row's print method.
 ***********************************************************************************************************/

void Matrix::print( std::string const &a_prefixForRow ) const {

    for( std::vector<Vector>::const_iterator iter = m_matrix.begin( ); iter < m_matrix.end( ); ++iter ) iter->print( a_prefixForRow );
}

}
