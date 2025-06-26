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

/*! \class Reference1d
 * Class for the GNDS <**reference**> node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the XYs2d.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

Reference1d::Reference1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Function1dForm( a_construction, a_node, a_setupInfo, FormType::reference1d, a_parent ),
        m_xlink( a_node.attribute_as_string( GIDI_hrefChars ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Reference1d::~Reference1d( ) {

}

/* *********************************************************************************************************//**
 * Returns the domain minimum for the instance.
 *
 * @return          The domain minimum for the instance.
 ***********************************************************************************************************/

double Reference1d::domainMin( ) const {

    throw Exception( "Reference1d::domainMin: not implemented" );
}

/* *********************************************************************************************************//**
 * Returns the domain maximum for the instance.
 *
 * @return              The domain maximum for the instance.
 ***********************************************************************************************************/

double Reference1d::domainMax( ) const {

    throw Exception( "Reference1d::domainMax: not implemented" );
}

/* *********************************************************************************************************//**
 * The **y** value of this at the domain value **a_x1**.
 * Currently not implemented.
 *
 * @param a_x1          [in]    Domain value to evaluate this at.
 * @return                      The value of this at the domain value **a_x1**.
 ***********************************************************************************************************/


double Reference1d::evaluate( LUPI_maybeUnused double a_x1 ) const {

    throw Exception( "Reference1d::evaluate: not implemented" );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Reference1d::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string attributes;

    attributes += a_writeInfo.addAttribute( GIDI_labelChars, label( ) );
    attributes += a_writeInfo.addAttribute( GIDI_hrefChars, m_xlink );
    a_writeInfo.addNodeStarterEnder( a_indent, moniker( ), attributes );
}

}               // End namespace Functions.

}               // End namespace GIDI.
