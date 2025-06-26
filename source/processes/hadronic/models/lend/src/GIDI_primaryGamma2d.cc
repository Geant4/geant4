/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "GIDI.hpp"
#include <HAPI.hpp>

namespace GIDI {

namespace Functions {

/*! \class PrimaryGamma2d
 * Class for the GNDS <**primaryGamma**> node.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the XYs2d.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

PrimaryGamma2d::PrimaryGamma2d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Function2dForm( a_construction, a_node, a_setupInfo, FormType::primaryGamma2d, a_parent ),
        m_domainMin( a_node.attribute( GIDI_domainMinChars ).as_double( ) ),
        m_domainMax( a_node.attribute( GIDI_domainMaxChars ).as_double( ) ),
        m_value( a_node.attribute( GIDI_valueChars ).as_double( ) ),
        m_finalState( a_node.attribute_as_string( GIDI_finalStateChars ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

PrimaryGamma2d::~PrimaryGamma2d( ) {

}

/* *********************************************************************************************************//**
 * The value of the primary gamma energy at the projectile energy *a_x2*.
 *
 * @param a_x2          [in]    The projectile's energy.
 * @param a_x1          [in]    Unknown.
 * @return                      Fix me.
 ***********************************************************************************************************/

double PrimaryGamma2d::evaluate( LUPI_maybeUnused double a_x2, LUPI_maybeUnused double a_x1 ) const {

// FIXME - Do we need to check domain?
#if !defined(__NVCC__) && !defined(__HIP__)
    throw Exception( "PrimaryGamma2d::evaluate: not implemented." );
#endif
    return( m_value );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 * @param       a_embedded          [in]        If *true*, *this* function is embedded in a higher dimensional function.
 * @param       a_inRegions         [in]        This is not used in this method.
 ***********************************************************************************************************/
 
void PrimaryGamma2d::toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, LUPI_maybeUnused bool a_embedded, LUPI_maybeUnused bool a_inRegions ) const {
 
    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string attributes = a_writeInfo.addAttribute( GIDI_valueChars, LUPI::Misc::doubleToShortestString( value( ) ) );

    attributes += a_writeInfo.addAttribute( GIDI_domainMinChars, LUPI::Misc::doubleToShortestString( domainMin( ) ) );
    attributes += a_writeInfo.addAttribute( GIDI_domainMaxChars, LUPI::Misc::doubleToShortestString( domainMax( ) ) );
    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );

    axes( ).toXMLList( a_writeInfo, indent2 );

    a_writeInfo.addNodeEnder( moniker( ) );
}

}               // End namespace Functions.

}               // End namespace GIDI.
