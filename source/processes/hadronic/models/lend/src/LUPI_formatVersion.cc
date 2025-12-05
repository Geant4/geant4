/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <limits.h>

#include "LUPI.hpp"

namespace LUPI {

/*! \class FormatVersion
 * Class to store GNDS format.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

FormatVersion::FormatVersion( ) :
        m_format( "" ),
        m_major( -1 ),
        m_minor( -1 ),
        m_patch( "" ) {

}

/* *********************************************************************************************************//**
 * @param a_formatVersion       [in]    The GNDS format.
 ***********************************************************************************************************/

FormatVersion::FormatVersion( std::string const &a_formatVersion ) :
        m_format( a_formatVersion ),
        m_major( -1 ),
        m_minor( -1 ),
        m_patch( "" ) {

    setFormat( a_formatVersion );
}

/* *********************************************************************************************************//**
 * @param a_formatVersion       [in]    The GNDS format.
 ***********************************************************************************************************/

FormatVersion::FormatVersion( FormatVersion const &a_formatVersion ) :
        m_format( a_formatVersion.format( ) ),
        m_major( a_formatVersion.major( ) ),
        m_minor( a_formatVersion.minor( ) ),
        m_patch( a_formatVersion.patch( ) ) {

}

/* *********************************************************************************************************//**
 * The assignment operator. This method sets the members of *this* to those of *a_rhs* except for those
 * not set by base classes.
 *
 * @param a_rhs                     [in]    Instance whose member are used to set the members of *this*.
 ***********************************************************************************************************/

FormatVersion &FormatVersion::operator=( FormatVersion const &a_rhs ) {

    if( this != &a_rhs ) {
        m_format = a_rhs.format( );
        m_major = a_rhs.major( );
        m_minor = a_rhs.minor( );
        m_patch = a_rhs.patch( );
    }

    return( *this );
}

/* *********************************************************************************************************//**
 * Set the format to *a_formatVersion* and parse its components.
 *
 * @param   a_formatVersion       [in]      The GNDS format.
 *
 * @return                                  true if format of the form "MAJOR.MINOR[.PATCH]" where MAJOR and MINOR are integers. Otherwise returns false.
 ***********************************************************************************************************/

bool FormatVersion::setFormat( std::string const &a_formatVersion ) {

    m_format = a_formatVersion;

    std::vector<std::string> formatItems = Misc::splitString( a_formatVersion, '.' );

    if( ( formatItems.size( ) < 2 ) || ( formatItems.size( ) > 3 ) ) goto err;

    if( !Misc::stringToInt( formatItems[0], m_major ) ) goto err;
    if( !Misc::stringToInt( formatItems[1], m_minor ) ) goto err;

    if( formatItems.size( ) == 3 ) m_patch = formatItems[2];

    return( true );

err:
    m_major = -1;
    m_minor = -1;
    m_patch = "";
    return( false );
}

/* *********************************************************************************************************//**
 * Returns true if m_format is a supported format and false otherwise;
 *  
 * @return                                  true if format is supported and false otherwise.
 ***********************************************************************************************************/

bool FormatVersion::supported( ) const {

    if( m_format == GNDS_formatVersion_1_10Chars ) return( true );
    if( m_format == GNDS_formatVersion_2_0Chars ) return( true );
    if( m_format == GNDS_formatVersion_2_0_LLNL_4Chars ) return( true );
    if( m_format == GNDS_formatVersion_2_1Chars ) return( true );

    return( false );
}

}               // End namespace LUPI.
