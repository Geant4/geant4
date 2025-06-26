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

namespace Functions {

/*! \class URR_probabilityTables1d
 * Class for the **GNDS* **URR_probabilityTables1d** node.
 */

/* *********************************************************************************************************//**
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

URR_probabilityTables1d::URR_probabilityTables1d( Construction::Settings const &a_construction, HAPI::Node const &a_node,
		        SetupInfo &a_setupInfo, Suite *a_parent ) :
        Function1dForm( a_construction, a_node, a_setupInfo, FormType::URR_probabilityTables1d, a_parent ),
        m_function2d( data2dParse( a_construction, a_node.first_child( ), a_setupInfo, nullptr ) ) {

    if( m_function2d != nullptr ) m_function2d->setAncestor( this );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

URR_probabilityTables1d::~URR_probabilityTables1d( ) {

    delete m_function2d;
}

/* *********************************************************************************************************//**
 * Returns the minimum energy for the URR probability tables.
 ***********************************************************************************************************/

double URR_probabilityTables1d::domainMin( ) const {

    return( m_function2d->domainMin( ) );
}

/* *********************************************************************************************************//**
 * Returns the maximum energy for the URR probability tables.
 ***********************************************************************************************************/

double URR_probabilityTables1d::domainMax( ) const {

    return( m_function2d->domainMax( ) );
}

/* *********************************************************************************************************//**
 * Returns the average value for the URR probability tables at energy *a_x1*. Currently does not work; always returns 0.0.
 *
 * @param       a_x1        [in]     The projectile energy to evaluate the URR probability tables at.
 ***********************************************************************************************************/

double URR_probabilityTables1d::evaluate( LUPI_maybeUnused double a_x1 ) const {

    return( 0.0 );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 * @param       a_embedded          [in]        If *true*, *this* function is embedded in a higher dimensional function.
 * @param       a_inRegions         [in]        If *true*, *this* is in a Regions2d container.
 ***********************************************************************************************************/

void URR_probabilityTables1d::toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, LUPI_maybeUnused bool a_embedded, LUPI_maybeUnused bool a_inRegions ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string attributes;

    if( label( ) != "" ) attributes = a_writeInfo.addAttribute( GIDI_labelChars, label( ) );
    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );

    if( m_function2d != nullptr ) m_function2d->toXMLList_func( a_writeInfo, indent2, false, false );

    a_writeInfo.addNodeEnder( moniker( ) );
}

}               // End namespace Functions.

}               // End namespace GIDI.
