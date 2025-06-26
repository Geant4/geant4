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

static ResonanceBackgroundRegion1d *resonanceBackground1dParse( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );

/*! \class ResonancesWithBackground1d
 * Class for the GNDS <**resonancesWithBackground**> node.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent          [in]     The parent GIDI::Suite.
 ***********************************************************************************************************/

ResonancesWithBackground1d::ResonancesWithBackground1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, 
                SetupInfo &a_setupInfo, Suite *a_parent ) :
        Function1dForm( a_construction, a_node, a_setupInfo, FormType::resonancesWithBackground1d, a_parent ),
        m_resonances( a_node.child( GIDI_resonancesChars ).attribute_as_string( GIDI_hrefChars ) ),
        m_background( a_construction, a_node.child( GIDI_resonanceBackground1dChars ), a_setupInfo, nullptr ) {

    for( HAPI::Node child = a_node.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
        std::string name( child.name( ) );

        if( name == "resonanceRegion" ) {
            m_resonances = a_node.attribute_as_string( "href" ); }
        else {
          /*
            if( m_background != NULL ) {
                throw new std::runtime_error("Oops");
            }
            m_background = data1dParse( child, NULL );
           */
        }
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

ResonancesWithBackground1d::~ResonancesWithBackground1d( ) {

}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 * @param       a_embedded          [in]        If *true*, *this* function is embedded in a higher dimensional function.
 * @param       a_inRegions         [in]        If *true*, *this* is in a Regions2d container.
 ***********************************************************************************************************/

void ResonancesWithBackground1d::toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, LUPI_maybeUnused bool a_embedded, LUPI_maybeUnused bool a_inRegions ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string attributes;

    a_writeInfo.addNodeStarter( a_indent, moniker( ) );

    attributes += a_writeInfo.addAttribute( GIDI_hrefChars, m_resonances );
    a_writeInfo.addNodeStarter( indent2, GIDI_resonancesChars, attributes );

    m_background.toXMLList_func( a_writeInfo, indent2, false, false );

    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class ResonanceBackground1d
 * Class for the GNDS <**background**> node.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent          [in]     The parent GIDI::Suite.
 ***********************************************************************************************************/

ResonanceBackground1d::ResonanceBackground1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo,   
                Suite *a_parent ) :
        Function1dForm( a_construction, a_node, a_setupInfo, FormType::resonanceBackground1d, a_parent ),
        m_resolvedRegion( resonanceBackground1dParse( a_construction, a_node.child( GIDI_resolvedRegionChars ), a_setupInfo, nullptr ) ),
        m_unresolvedRegion( resonanceBackground1dParse( a_construction, a_node.child( GIDI_unresolvedRegionChars ), a_setupInfo, nullptr ) ),
        m_fastRegion( resonanceBackground1dParse( a_construction, a_node.child( GIDI_fastRegionChars ), a_setupInfo, nullptr ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

ResonanceBackground1d::~ResonanceBackground1d( ) {

    delete m_resolvedRegion;
    delete m_unresolvedRegion;
    delete m_fastRegion;
}

/* *********************************************************************************************************//**
 * Returns the domain minimum for the instance.
 *
 * @return          The domain minimum for the instance.
 ***********************************************************************************************************/

double ResonanceBackground1d::domainMin( ) const {

    throw Exception( "ResonanceBackground1d::domainMin: not implemented" );
}

/* *********************************************************************************************************//**
 * Returns the domain maximum for the instance.
 *
 * @return              The domain maximum for the instance.
 ***********************************************************************************************************/

double ResonanceBackground1d::domainMax( ) const {

    throw Exception( "ResonanceBackground1d::domainMax: not implemented" );
}

/* *********************************************************************************************************//**
 * The value of *y(x1)* at the point *a_x1*.
 * Currently not implemented.
 *
 * @param a_x1          [in]    The point for the *x1* axis.
 * @return                      The value of the function at the point *a_x1*.
 ***********************************************************************************************************/

double ResonanceBackground1d::evaluate( LUPI_maybeUnused double a_x1 ) const {

    throw Exception( "ResonanceBackground1d::evaluate: not implemented" );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 * @param       a_embedded          [in]        If *true*, *this* function is embedded in a higher dimensional function.
 * @param       a_inRegions         [in]        If *true*, *this* is in a Regions2d container.
 ***********************************************************************************************************/

void ResonanceBackground1d::toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, LUPI_maybeUnused bool a_embedded, LUPI_maybeUnused bool a_inRegions ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    a_writeInfo.addNodeStarter( a_indent, moniker( ) );

    m_resolvedRegion->toXMLList_func( a_writeInfo, indent2, false, false );
    m_unresolvedRegion->toXMLList_func( a_writeInfo, indent2, false, false );
    m_fastRegion->toXMLList_func( a_writeInfo, indent2, false, false );

    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class ResonanceBackgroundRegion1d
 * Class for the GNDS <**resolvedRegion**> or GNDS <**fastRegion**> node.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

ResonanceBackgroundRegion1d::ResonanceBackgroundRegion1d( Construction::Settings const &a_construction, HAPI::Node const &a_node,
		        SetupInfo &a_setupInfo, Suite *a_parent ) :
        Function1dForm( a_construction, a_node, a_setupInfo, FormType::resonanceBackgroundRegion1d, a_parent ),
        m_function1d( data1dParse( a_construction, a_node.first_child( ), a_setupInfo, nullptr ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

ResonanceBackgroundRegion1d::~ResonanceBackgroundRegion1d( ) {

    delete m_function1d;
}

/* *********************************************************************************************************//**
 * Returns the domain minimum for the instance.
 *
 * @return          The domain minimum for the instance.
 ***********************************************************************************************************/

double ResonanceBackgroundRegion1d::domainMin( ) const {

    throw Exception( "ResonanceBackgroundRegion1d::domainMin: not implemented" );
}

/* *********************************************************************************************************//**
 * Returns the domain maximum for the instance.
 *
 * @return              The domain maximum for the instance.
 ***********************************************************************************************************/

double ResonanceBackgroundRegion1d::domainMax( ) const {

    throw Exception( "ResonanceBackgroundRegion1d::domainMax: not implemented" );
}

/* *********************************************************************************************************//**
 * The value of *y(x1)* at the point *a_x1*.
 * Currently not implemented.
 *
 * @param a_x1          [in]    The point for the *x1* axis.
 * @return                      The value of the function at the point *a_x1*.
 ***********************************************************************************************************/

double ResonanceBackgroundRegion1d::evaluate( LUPI_maybeUnused double a_x1 ) const {

    throw Exception( "ResonanceBackgroundRegion1d::evaluate: not implemented" );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 * @param       a_embedded          [in]        If *true*, *this* function is embedded in a higher dimensional function.
 * @param       a_inRegions         [in]        If *true*, *this* is in a Regions2d container.
 ***********************************************************************************************************/

void ResonanceBackgroundRegion1d::toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, LUPI_maybeUnused bool a_embedded, LUPI_maybeUnused bool a_inRegions ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    a_writeInfo.addNodeStarter( a_indent, moniker( ) );

    if( m_function1d != nullptr ) m_function1d->toXMLList_func( a_writeInfo, indent2, false, false );

    a_writeInfo.addNodeEnder( moniker( ) );
}

/* *********************************************************************************************************//**
 * Function that parses a node one-d function node. Called from a Suite::parse instance.
 *
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent          [in]    The parent GIDI::Suite that the returned Form will be added to.
 *
 * @return                          The parsed and constructed resonanceBackground region instance.
 ***********************************************************************************************************/

static ResonanceBackgroundRegion1d *resonanceBackground1dParse( Construction::Settings const &a_construction, HAPI::Node const &a_node,
		SetupInfo &a_setupInfo, Suite *a_parent ) {

    std::string name( a_node.name( ) );

    if( name == "" ) return( nullptr );

    return( new ResonanceBackgroundRegion1d( a_construction, a_node, a_setupInfo, a_parent ) );
}

}               // End namespace Functions.

}               // End namespace GIDI.
