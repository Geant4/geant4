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

/*! \class Data1d
 * Currently not used.
 */

/*
============================================================
 *
 * @param a_number
 * @param a_xs
 * @return
 */
Data1d::Data1d( std::size_t a_number, double const *const a_xs ) {

    m_xs.resize( a_number );
    m_ys.resize( a_number );

    for( std::size_t i1 = 0; i1 < a_number; ++i1 ) {
        m_xs[i1] = a_xs[i1];
        m_ys[i1] = 0.0;
    }
}
/*
============================================================
 *
 * @param a_number
 * @param a_xs
 * @param a_ys
 * @return
 */
Data1d::Data1d( std::size_t a_number, double const *const a_xs, double const *const a_ys ) {

    m_xs.resize( a_number );
    m_ys.resize( a_number );
    for( std::size_t i1 = 0; i1 < a_number; ++i1 ) {
        m_xs[i1] = a_xs[i1];
        m_ys[i1] = a_ys[i1];
    }
}
/*
============================================================
 *
 * @param a_xs
 * @return
 */
Data1d::Data1d( std::vector<double> const &a_xs ) {

    m_xs = a_xs;
    m_ys.resize( a_xs.size( ) );
    for( std::vector<double>::iterator iter = m_ys.begin( ); iter < m_ys.end( ); ++iter ) *iter = 0.0;
}
/*
============================================================
 *
 * @param a_xs
 * @param a_ys
 * @return
 */
Data1d::Data1d( std::vector<double> const &a_xs, std::vector<double> const &a_ys ) {

    if( a_xs.size( ) != a_ys.size( ) ) throw Exception( "sizes not the same" );
    m_xs = a_xs;
    m_ys = a_ys;
}
/*
============================================================
 *
 * @param a_gidi_1dData
 * @return
 */
Data1d::Data1d( Data1d const &a_gidi_1dData ) {

    m_xs = a_gidi_1dData.m_xs;
    m_ys = a_gidi_1dData.m_ys;
}
/*
============================================================
*/
Data1d::~Data1d( ) {

}

/*
============================================================
============================ add ===========================
============================================================
 */
Data1d Data1d::operator+( double a_value ) const {

    Data1d gidi_1dData( *this );

    gidi_1dData += a_value;
    return( gidi_1dData );
}
/*
============================================================
 */
Data1d &Data1d::operator+=( double a_value ) {

    for( std::vector<double>::iterator iter = m_ys.begin( ); iter < m_ys.end( ); ++iter ) *iter += a_value;
    return( *this );
}
/*
============================================================
 */
Data1d Data1d::operator+( Data1d const &a_rhs ) const {

    Data1d gidi_1dData( *this );

    gidi_1dData += a_rhs;
    return( gidi_1dData );
}
/*
============================================================
 */
Data1d &Data1d::operator+=( Data1d const &a_rhs ) {

    std::size_t __size = size( );

    if( __size != a_rhs.size( ) ) throw Exception( "data1d sizes differ." );
    for( std::size_t i1 = 0; i1 < __size; ++i1 ) m_ys[i1] += a_rhs.m_ys[i1];
    return( *this );
}

/*
============================================================
========================== subtract ========================
============================================================
 */
Data1d Data1d::operator-( double a_value ) const {

    Data1d gidi_1dData( *this );

    gidi_1dData -= a_value;
    return( gidi_1dData );
}
/*
============================================================
 */
Data1d &Data1d::operator-=( double a_value ) {

    for( std::vector<double>::iterator iter = m_ys.begin( ); iter < m_ys.end( ); ++iter ) *iter -= a_value;
    return( *this );
}
/*
============================================================
 */
Data1d Data1d::operator-( Data1d const &a_rhs ) const {

    Data1d gidi_1dData( *this );

    gidi_1dData -= a_rhs;
    return( gidi_1dData );
}
/*
============================================================
 */
Data1d &Data1d::operator-=( Data1d const &a_rhs ) {

    std::size_t __size = size( );

    if( __size != a_rhs.size( ) ) throw Exception( "data1d sizes differ." );
    for( std::size_t i1 = 0; i1 < __size; ++i1 ) m_ys[i1] -= a_rhs.m_ys[i1];
    return( *this );
}

/*
============================================================
========================== multiply ========================
============================================================
 */
Data1d Data1d::operator*( double a_value ) const {

    Data1d gidi_1dData( *this );

    gidi_1dData *= a_value;
    return( gidi_1dData );
}
/*
============================================================
 */
Data1d &Data1d::operator*=( double a_value ) {

    for( std::vector<double>::iterator iter = m_ys.begin( ); iter < m_ys.end( ); ++iter ) *iter *= a_value;
    return( *this );
}

/*
============================================================
========================== divide ==========================
============================================================
 */
Data1d Data1d::operator/( double a_value ) const {

    Data1d gidi_1dData( *this );

    gidi_1dData /= a_value;
    return( gidi_1dData );
}
/*
============================================================
 */
Data1d &Data1d::operator/=( double a_value ) {

    if( a_value == 0 ) throw Exception( "divide by zero." );
    for( std::vector<double>::iterator iter = m_ys.begin( ); iter < m_ys.end( ); ++iter ) *iter /= a_value;
    return( *this );
}

/*
============================================================
========================== others ==========================
============================================================
 */
void Data1d::print( std::string const &a_prefix ) const {

    std::size_t __size = size( );

    for( std::size_t i1 = 0; i1 < __size; ++i1 ) {
        std::cout << a_prefix;
        printf( "%18.11e %18.11e", m_xs[i1], m_ys[i1] );
        std::cout << std::endl;
    }
}

}
