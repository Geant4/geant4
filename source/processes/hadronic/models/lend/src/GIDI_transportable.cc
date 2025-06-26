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

#define GIDI_conserveChars "conserve"

/*! \class Transportable
 * Class for the GNDS <**transportable**> node that resides under the <**transportables**> node.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed to construct a Transportable instance.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops            [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

Transportable::Transportable( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		        PoPI::Database const &a_pops, Suite *a_parent ) :
        Form( a_node, a_setupInfo, FormType::transportable, a_parent ),
        m_conserve( a_node.attribute_as_string( GIDI_conserveChars ) ),
        m_group( a_construction, a_node.child( GIDI_groupChars ), a_setupInfo, a_pops ) {
}

/* *********************************************************************************************************//**
 * Copy constructor.
 *
 * @param a_transportable   [in]    Transportable instance to copy.
 ***********************************************************************************************************/

Transportable::Transportable( Transportable const &a_transportable ) :
        Form( a_transportable ),
        m_conserve( a_transportable.conserve( ) ),
        m_group( a_transportable.group( ) ) {

}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Transportable::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string attributes = a_writeInfo.addAttribute( GIDI_labelChars, label( ) );

    attributes += a_writeInfo.addAttribute( GIDI_conserveChars, m_conserve );
    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );
    m_group.toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

}
