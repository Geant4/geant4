/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdio.h>

#include "GIDI.hpp"

namespace GIDI {

/*! \class Vector
 * This class stores a mathematical vector and has methods that perform several vector operations (e.g., addition, subtraction).
 */

/* *********************************************************************************************************//**
 *
 * @param a_size            [in]    Number of initial elements of the matrix. All elements are initialized to 0.
 ***********************************************************************************************************/

Vector::Vector( std::size_t a_size ) {

    m_vector.resize( a_size, 0.0 );
}

/* *********************************************************************************************************//**
 *
 * @param a_values          [in]    A list of doubles to initialize *this* with.
 ***********************************************************************************************************/

Vector::Vector( std::vector<double> const &a_values ) {

    m_vector = a_values;
}

/* *********************************************************************************************************//**
 *
 * @param a_number          [in]    This number of element pointed to by *a_values*.
 * @param a_values          [in]    A list of doubles to initialize *this* with.
 ***********************************************************************************************************/

Vector::Vector( std::size_t a_number, double const *a_values ) {

    m_vector.resize( a_number );
    for( std::size_t i1 = 0; i1 < a_number; ++i1 ) m_vector[i1] = a_values[i1];
}

/* *********************************************************************************************************//**
 *
 * @param a_vector          [in]    Vector to copy.
 ***********************************************************************************************************/

Vector::Vector( Vector const &a_vector ) :
        m_vector( a_vector.m_vector ) {

}

/* *********************************************************************************************************//**
 * Returns a new Vector whose elements are *this* plus *a_rhs*.
 *
 * @param a_rhs         [in]    The value to add to each element.
 * @return                      New Vectors whose elements are *this* plus *a_rhs*.
 ***********************************************************************************************************/

Vector &Vector::operator=( Vector const &a_rhs ) {

    if( this != &a_rhs ) {
        m_vector = a_rhs.m_vector;
    }

    return( *this );
}
/*
=========================================================
*/
Vector::~Vector( ) {

}

/* *********************************************************************************************************//**
 * Returns a new Vector whose elements are *this* plus *a_value*.
 *
 * @param a_value       [in]    The value to add to each element.
 * @return                      New Vector whose elements are *this* plus *a_value*.
 ***********************************************************************************************************/

Vector Vector::operator+( double a_value ) const {

    Vector gidiVector( *this );

    gidiVector += a_value;
    return( gidiVector );
}

/* *********************************************************************************************************//**
 * Adds *a_value* to each element of *this*.
 *
 * @param a_value       [in]    The value to add to each element.
 * @return                      Returns reference to *this*.
 ***********************************************************************************************************/

Vector &Vector::operator+=( double a_value ) {

    for( std::vector<double>::iterator iter = m_vector.begin( ); iter < m_vector.end( ); ++iter ) *iter += a_value;

    return( *this );
}

/* *********************************************************************************************************//**
 * Adds two Vectors.
 *
 * @param a_rhs         [in]    Vector to add to *this*.
 * @return                      New Vector that is the vector sum of *this* and *a_rhs*.
 ***********************************************************************************************************/

Vector Vector::operator+( Vector const &a_rhs ) const {

    Vector gidiVector( *this );

    gidiVector += a_rhs;
    return( gidiVector );
}

/* *********************************************************************************************************//**
 * Adds *a_rhs* to *this*.
 *
 * @param a_rhs         [in]    Vector to add to *this*.
 * @return                      Returns reference to *this*.
 ***********************************************************************************************************/

Vector &Vector::operator+=( Vector const &a_rhs ) {

    if( a_rhs.size( ) == 0 ) return( *this );

    if( size( ) == 0 ) resize( a_rhs.size( ) );
    if( size( ) != a_rhs.size( ) ) throw Exception( "vector sizes differ." );

    std::size_t i1 = 0;
    for( std::vector<double>::iterator iter = m_vector.begin( ); iter < m_vector.end( ); ++iter, ++i1 ) *iter += a_rhs[i1];

    return( *this );
}

/* *********************************************************************************************************//**
 * Returns a new Vector whose elements are *this* minus *a_value*.
 *
 * @param a_value       [in]    The value to subtract from each element.
 * @return                      New Vector whose elements are *this* plus *a_value*.
 ***********************************************************************************************************/

Vector Vector::operator-( double a_value ) const {

    Vector gidiVector( *this );

    gidiVector -= a_value;
    return( gidiVector );
}

/* *********************************************************************************************************//**
 * Subtracts *a_value* from each element of *this*.
 *
 * @param a_value       [in]    The value to subtract from each element.
 * @return                      Returns reference to *this*.
 ***********************************************************************************************************/

Vector &Vector::operator-=( double a_value ) {

    for( std::vector<double>::iterator iter = m_vector.begin( ); iter < m_vector.end( ); ++iter ) *iter -= a_value;

    return( *this );
}

/* *********************************************************************************************************//**
 * Subtracts *a_rhs* from *this*.
 *
 * @param a_rhs         [in]    Vector to subtract from *this*.
 * @return                      New Vector that is *this* minus *a_rhs*.
 ***********************************************************************************************************/

Vector Vector::operator-( Vector const &a_rhs ) const {

    Vector gidiVector( *this );

    gidiVector -= a_rhs;
    return( gidiVector );
}

/* *********************************************************************************************************//**
 * Subtracts *a_rhs* to *this*.
 *
 * @param a_rhs         [in]    Vector to subtract from *this*.
 * @return                      Returns reference to *this*.
 ***********************************************************************************************************/

Vector &Vector::operator-=( Vector const &a_rhs ) {

    if( a_rhs.size( ) == 0 ) return( *this );

    if( size( ) == 0 ) resize( a_rhs.size( ) );
    if( size( ) != a_rhs.size( ) ) throw Exception( "vector sizes differ." );

    std::size_t i1 = 0;
    for( std::vector<double>::iterator iter = m_vector.begin( ); iter < m_vector.end( ); ++iter, ++i1 ) *iter -= a_rhs[i1];

    return( *this );
}

/* *********************************************************************************************************//**
 * Returns a new Vector whose elements are *this* multiplied by *a_value*.
 *
 * @param a_value       [in]    The value to multiply each element by.
 * @return                      New Vector whose elements are *this* multiply by *a_value*.
 ***********************************************************************************************************/

Vector Vector::operator*( double a_value ) const {

    Vector gidiVector( *this );

    gidiVector *= a_value;
    return( gidiVector );
}

/* *********************************************************************************************************//**
 * Multiplies each element of *this* by *a_value*.
 *
 * @param a_value       [in]    The value to multiply each element by.
 * @return                      Returns reference to *this*.
 ***********************************************************************************************************/

Vector &Vector::operator*=( double a_value ) {

    for( std::vector<double>::iterator iter = m_vector.begin( ); iter < m_vector.end( ); ++iter ) *iter *= a_value;

    return( *this );
}

/* *********************************************************************************************************//**
 * Returns a new Vector whose elements are *this* divided by *a_value*.
 *
 * @param a_value       [in]    The value to divide each element by.
 * @return                      New Vector whose elements are *this* divided by *a_value*.
 ***********************************************************************************************************/

Vector Vector::operator/( double a_value ) const {

    Vector gidiVector( *this );

    gidiVector /= a_value;
    return( gidiVector );
}

/* *********************************************************************************************************//**
 * Divides each element of *this* by *a_value*.
 *
 * @param a_value       [in]    The value to divide each element by.
 * @return                      Returns reference to *this*.
 ***********************************************************************************************************/

Vector &Vector::operator/=( double a_value ) {

    if( a_value == 0 ) throw Exception( "divide by zero." );
    for( std::vector<double>::iterator iter = m_vector.begin( ); iter < m_vector.end( ); ++iter ) *iter /= a_value;

    return( *this );
}

/* *********************************************************************************************************//**
 * Reverse the elements of *this*.
 ***********************************************************************************************************/

void Vector::reverse( ) {

    std::size_t i2 = size( ), n_2 = i2 / 2;

    --i2;
    for( std::size_t i1 = 0; i1 < n_2; ++i1, --i2 ) {
        double temp = m_vector[i1];

        m_vector[i1] = m_vector[i2];
        m_vector[i2] = temp;
    }
}

/* *********************************************************************************************************//**
 * Sets all elements in the range [*a_start*,*a_end*) to *a_value*.
 *
 * @param a_start       [in]    The starting flat-cell index of *this* to fill with *a_value*.
 * @param a_end         [in]    One after the last flat-cell index of *this* to fill with *a_value*.
 * @param a_value       [in]    The value to set each double in the range to.
 ***********************************************************************************************************/

void Vector::setToValueInFlatRange( std::size_t a_start, std::size_t a_end, double a_value ) {

    a_end = std::min( a_end, m_vector.size( ) );
    for( ; a_start < a_end; ++a_start ) m_vector[a_start] = a_value;
}
/* *********************************************************************************************************//**
 * Returns the sum over the values of *this*.
 ***********************************************************************************************************/

double Vector::sum( ) {

    double sum1 = 0.0;

    for( std::size_t i1 = 0; i1 < size( ); ++i1 ) sum1 += m_vector[i1];

    return( sum1 );
}

/* *********************************************************************************************************//**
 * Prints the contents of *this* to std::cout as one line prefixed with **a_prefix**.
 *
 * @param a_prefix      [in]    Prefix to add to line.
 ***********************************************************************************************************/

void Vector::print( std::string const &a_prefix ) const {

    std::cout << a_prefix;
    for( std::vector<double>::const_iterator iter = m_vector.begin( ); iter < m_vector.end( ); ++iter ) printf( "%19.11e", *iter );
    std::cout << std::endl;
}

/* *********************************************************************************************************//**
 * Writes the contents of *this* to *a_file* as one line prefixed with **a_prefix**.
 *
 * @param a_file        [in]    A pointer to an opened C FILE instance where the data are to be written.
 * @param a_prefix      [in]    Prefix to add to line.
 ***********************************************************************************************************/

void Vector::write( FILE *a_file, std::string const &a_prefix ) const {

    if( a_prefix.size( ) > 0 ) fprintf( a_file, "# %s\n", a_prefix.c_str( ) );
    for( std::vector<double>::const_iterator iter = m_vector.begin( ); iter < m_vector.end( ); ++iter ) fprintf( a_file, "%19.11e\n", *iter );
}

/* *********************************************************************************************************//**
 * This method writes the contents of *this* to output *a_file* with each pair (boundary, value) as one line. Each line is written
 * with the C *printf* style format as specified via the *a_format* arguement. If *a_epsilon* is zero then ( size() + 1 )
 * lines are printed with the last y-value being the last value in *this*. If *a_epsilon* is not zero then ( 2 * size() )
 * lines are printed with each boundary with index 1 to size() being print at first ( boundary[index] * ( 1 - a_epsilon ) )
 * and then ( boundary[index] * ( 1 + a_epsilon ) ) as :
 *
 *      boundary[index] * ( 1 - a_epsilon ) ), this[index-1]
 *      boundary[index] * ( 1 + a_epsilon ) ), this[index]
 * 
 *
 * @param a_file        [in]    A pointer to an opened C FILE instance where the data are to be written.
 * @param a_format      [in]    A C *printf* style format to wrint one line of each (boundary, value) pair. Must include the line feed character.
 * @param a_boundaries  [in]    The list of boundaries. Must contain (size() + 1) values except for the case when size() is 0.
 * @param a_epsilon     [in]    Prefix to add to line.
 ***********************************************************************************************************/

void Vector::writeWithBoundaries( FILE *a_file, char const *a_format, std::vector<double> const &a_boundaries, double a_epsilon ) const {

    if( size( ) == 0 ) {
        if( a_boundaries.size( ) > 0 ) {
            Vector vector( a_boundaries.size( ) - 1 );
            vector.writeWithBoundaries2( a_file, a_format, a_boundaries, a_epsilon ); } }
    else {
        if( size( ) + 1 != a_boundaries.size( ) ) throw Exception( "Vector::writeWithBoundaries: Vector size and number of boundaries are not compatible." );
        writeWithBoundaries2( a_file, a_format, a_boundaries, a_epsilon );
    }
}

/* *********************************************************************************************************//**
 * For internal use only. See method **writeWithBoundaries** for description.
 *
 * @param a_file        [in]    A pointer to an opened C FILE instance where the data are to be written.
 * @param a_format      [in]    A C *printf* style format to wrint one line of each (boundary, value) pair. Do not include the line feed character.
 * @param a_boundaries  [in]    The list of boundaries. Must contain (size() + 1) values except for the case when size() is 0.
 * @param a_epsilon     [in]    Prefix to add to line.
 ***********************************************************************************************************/
void Vector::writeWithBoundaries2( FILE *a_file, char const *a_format, std::vector<double> const &a_boundaries, double a_epsilon ) const {

    std::size_t numberOfValues = size( );

    if( a_epsilon == 0.0 ) {
        for( std::size_t index = 0; index < numberOfValues; ++index ) fprintf( a_file, a_format, a_boundaries[index], m_vector[index] ); }
    else {
        if( numberOfValues > 0 ) fprintf( a_file, a_format, a_boundaries[0], m_vector[0] );
        for( std::size_t index = 1; index < numberOfValues; ++index ) {
            fprintf( a_file, a_format, a_boundaries[index] * ( 1.0 - a_epsilon ), m_vector[index-1] );

            fprintf( a_file, a_format, a_boundaries[index] * ( 1.0 + a_epsilon ), m_vector[index] );
        }
    }
    if( numberOfValues > 0 ) fprintf( a_file, a_format, a_boundaries[numberOfValues], m_vector[numberOfValues-1] );
}

}
