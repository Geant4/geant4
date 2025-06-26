/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "GIDI.hpp"

namespace GIDI {

namespace Transporting {

/*! \class MultiGroup
 * Specifies the flux data for a specified Legendre order (see class Flux).
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

MultiGroup::MultiGroup( ) {

}

/* *********************************************************************************************************//**
 * @param a_label           [in]    The label for the MultiGroup.
 * @param a_length          [in]    The number of boundaries values.
 * @param a_boundaries      [in]    The list of boundaries.
 ***********************************************************************************************************/

MultiGroup::MultiGroup( std::string const &a_label, int a_length, double const *a_boundaries ) :
        m_label( a_label ) {

    for( int i1 = 0; i1 < a_length; ++i1 ) m_boundaries.push_back( a_boundaries[i1] );
}

/* *********************************************************************************************************//**
 * @param a_label           [in]    The label for the MultiGroup.
 * @param a_boundaries      [in]    The list of boundaries.
 ***********************************************************************************************************/

MultiGroup::MultiGroup( std::string const &a_label, std::vector<double> const &a_boundaries ) :
        m_label( a_label ),
        m_boundaries( a_boundaries ) {

}

/* *********************************************************************************************************//**
 * @param a_group           [in]    The Group used to set *this*.
 ***********************************************************************************************************/

MultiGroup::MultiGroup( Group const &a_group ) :
        m_label( a_group.label( ) ),
        m_boundaries( a_group.data( ) ) {

}

/* *********************************************************************************************************//**
 * @param a_multiGroup      [in]    The MultiGroup instance to copy.
 ***********************************************************************************************************/

MultiGroup::MultiGroup( MultiGroup const &a_multiGroup ) :
        m_label( a_multiGroup.label( ) ),
        m_boundaries( a_multiGroup.boundaries( ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

MultiGroup::~MultiGroup( ) {

}

/* *********************************************************************************************************//**
 * The assignment operator. This method sets the members of *this* to those of *a_rhs*.
 *
 * @param a_rhs                     [in]    Instance whose member are used to set the members of *this*.
 *
 * @return                                  A reference to the updated MultiGroup instance.
 ***********************************************************************************************************/

MultiGroup &MultiGroup::operator=( MultiGroup const &a_rhs ) {

    if( this != &a_rhs ) {
        m_label = a_rhs.label( );
        m_boundaries = a_rhs.boundaries( );
    }

    return( *this );
}

/* *********************************************************************************************************//**
 * Returns the multi-group index whose boundaries enclose *a_energy*. If *a_encloseOutOfRange* is true and
 * *a_energy* is below the lowest boundary, 0 is returned, otherwise -2 is returned. If *a_encloseOutOfRange* is true and
 * *a_energy* is above the highest boundary, the last multi-group index is returned, otherwise -1 is returned.
 *
 * @param a_energy                  [in]    The energy of the whose index is to be returned.
 * @param a_encloseOutOfRange               Determines the action if energy is below or above the domain of the boundaries.
 * @return                                  The index whose boundaries enclose *a_energy*.
 ***********************************************************************************************************/

int MultiGroup::multiGroupIndexFromEnergy( double a_energy, bool a_encloseOutOfRange ) const {

    int iMin = 0, iMid, iMax = (int) m_boundaries.size( ), iMaxM1 = iMax - 1;

    if( iMax == 0 ) return( -3 );
    if( a_energy < m_boundaries[0] ) {
        if( a_encloseOutOfRange ) return( 0 );
        return( -2 );
    }
    if( a_energy > m_boundaries[iMaxM1] ) {
        if( a_encloseOutOfRange ) return( iMax - 2 );
        return( -1 );
    }
    while( 1 ) {
        iMid = ( iMin + iMax ) >> 1;
        if( iMid == iMin ) break;
        if( a_energy < m_boundaries[iMid] ) {
            iMax = iMid; }
        else {
            iMin = iMid;
        }
    }
    if( iMin == iMaxM1 ) iMin--;
    return( iMin );
}

/* *********************************************************************************************************//**
 * @param a_label           [in]    The label for *this*.
 * @param a_boundaries      [in]    The boundaries to set *this* to.
 ***********************************************************************************************************/

void MultiGroup::set( std::string const &a_label, std::vector<double> const &a_boundaries ) {

    m_label = a_label;
    m_boundaries = a_boundaries;
}

/* *********************************************************************************************************//**
 * Print the MultiGroup to std::cout. Mainly for debugging.
 *
 * @param a_indent                  [in]    The std::string to print at the beginning.
 * @param a_outline                 [in]    If true, does not print the flux values.
 * @param a_valuesPerLine           [in]    The number of points (i.e., energy, flux pairs) to print per line.
 ***********************************************************************************************************/

void MultiGroup::print( std::string const &a_indent, bool a_outline, int a_valuesPerLine ) const {

    int nbs = size( );
    bool printIndent( true );

    std::cout << a_indent << "GROUP: label = '" << m_label << "': length = " << nbs << std::endl;
    if( a_outline ) return;
    for( int ib = 0; ib < nbs; ib++ ) {
        if( printIndent ) std::cout << a_indent;
        printIndent = false;
        std::cout << LUPI::Misc::argumentsToString( "%16.8e", m_boundaries[ib] );
        if( ( ( ib + 1 ) % a_valuesPerLine ) == 0 ) {
            std::cout << std::endl;
            printIndent = true;
        }
    }
    if( nbs % a_valuesPerLine ) std::cout << std::endl;
}

/*! \class Groups_from_bdfls
 * Specifies the data for a specified Legendre order (see class Flux).
 */

/* *********************************************************************************************************//**
 * Reads in multi-group data from a *bdfls* file as a list of MultiGroup instances.
 *
 * @param a_fileName                [in]    The *bdfls* file name.
 ***********************************************************************************************************/

Groups_from_bdfls::Groups_from_bdfls( std::string const &a_fileName ) {

    initialize( a_fileName.c_str( ) );
}

/* *********************************************************************************************************//**
 * Reads in multi-group data from a *bdfls* file as a list of MultiGroup instances.
 *
 * @param a_fileName                [in]    The *bdfls* file name.
 ***********************************************************************************************************/

Groups_from_bdfls::Groups_from_bdfls( char const *a_fileName ) {

    initialize( a_fileName );
}

/* *********************************************************************************************************//**
 * Used by constructors to do most of the work.
 *
 * @param a_fileName                [in]    The *bdfls* file name.
 ***********************************************************************************************************/

void Groups_from_bdfls::initialize( char const *a_fileName ) {

    char buffer[132], *pEnd, cValue[16];
    FILE *fIn = fopen( a_fileName, "r" );
    if( fIn == nullptr ) throw Exception( "Groups_from_bdfls::initialize: Could not open bdfls file." );

    while( true ) {
        int gid( -1 );
        if( fgets( buffer, 132, fIn ) == nullptr ) throw Exception( "Groups_from_bdfls::initialize: fgets failed for gid." );
        if( strlen( buffer ) > 73 ) {
            if( buffer[72] == '1' ) break;
        }
        gid = (int) strtol( buffer, &pEnd, 10 );
        if( gid == -1 ) throw Exception( "Groups_from_bdfls::initialize: converting gid to long failed." );
        std::string label( LLNL_gidToLabel( gid ) );

        long numberOfBoundaries( -1 );
        if( fgets( buffer, 132, fIn ) == nullptr ) throw Exception( "Groups_from_bdfls::initialize: fgets failed for numberOfBoundaries." );
        numberOfBoundaries = strtol( buffer, &pEnd, 10 );
        if( numberOfBoundaries == -1 ) throw Exception( "Groups_from_bdfls::initialize: converting gid to long failed." );

        long index( 0 );
        std::vector<double> boundaries( numberOfBoundaries );
        while( numberOfBoundaries > 0 ) {
            long i1, n1( 6 );
            if( numberOfBoundaries < 6 ) n1 = numberOfBoundaries;
            if( fgets( buffer, 132, fIn ) == nullptr ) throw Exception( "Groups_from_bdfls::initialize: fgets failed for boundaries." );
            for( i1 = 0; i1 < n1; ++i1, ++index ) {
                strncpy( cValue, &buffer[12*i1], 12 );
                cValue[12] = 0;
                boundaries[index] = strtod( cValue, &pEnd );
            }
            numberOfBoundaries -= n1;
        }
        m_multiGroups.push_back( MultiGroup( label, boundaries ) );
    }

    fclose( fIn );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Groups_from_bdfls::~Groups_from_bdfls( ) {

}

/* *********************************************************************************************************//**
 * Returns the MultiGroup whose *label* is *a_label*.
 *
 * @param a_label           [in]    The *label* of the MultiGroup to return.
 * @return                          Returns the MultiGroup whose *label* is *a_label*.
 ***********************************************************************************************************/

MultiGroup Groups_from_bdfls::viaLabel( std::string const &a_label ) const {

    for( int ig = 0; ig < (int) m_multiGroups.size( ); ++ig ) {
        if( m_multiGroups[ig].label( ) ==  a_label ) return( m_multiGroups[ig] );
    }
    throw Exception( "Groups_from_bdfls::viaLabel: label not found." );
}

/* *********************************************************************************************************//**
 * Returns the MultiGroup whose *gid* is *a_gid*.
 *
 * @param a_gid             [in]    The bdfls *gid*.
 * @return                          Returns the MultiGroup whose *label* is *a_label*.
 ***********************************************************************************************************/

MultiGroup Groups_from_bdfls::getViaGID( int a_gid ) const {

    std::string label( LLNL_gidToLabel( a_gid ) );

    return( viaLabel( label ) );
}

/* *********************************************************************************************************//**
 * Returns a list of *label*'s for all the MultiGroup's present in *this*.
 *
 * @return                          Returns the MultiGroup whose *label* is *a_label*.
 ***********************************************************************************************************/

std::vector<std::string> Groups_from_bdfls::labels( ) const {

    int size = (int) m_multiGroups.size( );
    std::vector<std::string> _labels( size );

    for( int if1 = 0; if1 < size; ++if1 ) _labels[if1] = m_multiGroups[if1].label( );
    return( _labels );
}

/* *********************************************************************************************************//**
 * Returns a list of *gid*'s for all the MultiGroup's present in *this*.
 *
 * @return                          The list of *gid*'s.
 ***********************************************************************************************************/

std::vector<int> Groups_from_bdfls::GIDs( ) const {

    int size = (int) m_multiGroups.size( );
    std::vector<int> fids( size );
    char *e;

    for( int if1 = 0; if1 < size; ++if1 ) {
        fids[if1] = (int) strtol( &(m_multiGroups[if1].label( ).c_str( )[9]), &e, 10 );
    }
    return( fids );
}

/* *********************************************************************************************************//**
 * Print each MultiGroup to std::cout in *this*. Mainly for debugging.
 *
 * @param a_outline                 [in]    Passed to each MultiGroup print method.
 * @param a_valuesPerLine           [in]    Passed to each MultiGroup print method.
 ***********************************************************************************************************/

void Groups_from_bdfls::print( bool a_outline, int a_valuesPerLine ) const {

    int ngs = (int) m_multiGroups.size( );

    std::cout << "BDFLS GROUPs: number of groups = " << ngs << std::endl;
    for( int if1 = 0; if1 < ngs ; ++if1 ) m_multiGroups[if1].print( "  ", a_outline, a_valuesPerLine );
}

}

}
