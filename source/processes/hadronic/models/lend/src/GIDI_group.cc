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

/*! \class Group
 * Class for the GNDS <**group**> node that resides under the <**transportable**> node.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed to construct a Group instance.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops            [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 ***********************************************************************************************************/

Group::Group( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, LUPI_maybeUnused PoPI::Database const &a_pops ) :
        Form( a_node, a_setupInfo, FormType::group ),
        m_grid( a_node.child( GIDI_gridChars ), a_setupInfo, a_construction.useSystem_strtod( ) ) {

}

/* *********************************************************************************************************//**
 * Copy constructor.
 *
 * @param a_group           [in]    Group instance to copy.
 ***********************************************************************************************************/

Group::Group( Group const &a_group ) :
        Form( a_group ),
        m_grid( a_group.grid( ) ) {

}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Group::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    a_writeInfo.addNodeStarter( a_indent, moniker( ), a_writeInfo.addAttribute( GIDI_labelChars, label( ) ) );
    m_grid.toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class Groups
 * Class for the GNDS <**groups**> node that contains a list of <**group**> nodes.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Groups::Groups( ) :
        Suite( GIDI_groupsChars ) {

}

/* *********************************************************************************************************//**
 * @param a_fileName            [in]    File containing a groups node to be parsed.
 ***********************************************************************************************************/

Groups::Groups( std::string const &a_fileName ) :
        Suite( GIDI_groupsChars ) {

    addFile( a_fileName );
}

/* *********************************************************************************************************//**
 * Adds the contents of the specified file to *this*.
 *
 * @param a_fileName            [in]    File containing a groups node to be parsed.
 ***********************************************************************************************************/

void Groups::addFile( std::string const &a_fileName ) {

    HAPI::File *doc = new HAPI::PugiXMLFile( a_fileName.c_str( ), "Groups::addFile" );

    HAPI::Node groups = doc->first_child( );

    std::string name( groups.name( ) );
    if( name != GIDI_groupsChars ) throw Exception( "Invalid groups node file: file node name is '" + name + "'." );

    Construction::Settings construction( Construction::ParseMode::all, GIDI::Construction::PhotoMode::atomicOnly );
    PoPI::Database pops;

    SetupInfo setupInfo( nullptr );

    std::string formatVersionString = groups.attribute_as_string( GIDI_formatChars );
    if( formatVersionString == "" ) formatVersionString = GNDS_formatVersion_1_10Chars;
    LUPI::FormatVersion formatVersion;
    formatVersion.setFormat( formatVersionString );
    if( !formatVersion.supported( ) ) throw Exception( "Unsupport GND format version" + formatVersionString );
    setupInfo.m_formatVersion = formatVersion;

    for( HAPI::Node child = groups.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
        Group *group = new Group( construction, child, setupInfo, pops );

        add( group );
    }
    delete doc;
}

}
