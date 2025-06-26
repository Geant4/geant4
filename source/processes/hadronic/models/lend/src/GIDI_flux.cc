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

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options for parsing.
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the Protare.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 ***********************************************************************************************************/

Flux::Flux( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo ) :
        Form( a_node, a_setupInfo, FormType::flux ),
        m_flux( data2dParse( a_construction, a_node.first_child( ), a_setupInfo, nullptr ) ) {

    if( m_flux != nullptr ) m_flux->setAncestor( this );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Flux::~Flux( ) {

    delete m_flux;
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Flux::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string attributes = a_writeInfo.addAttribute( GIDI_labelChars, label( ) );

    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );
    if( m_flux != nullptr ) m_flux->toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class Fluxes
 * Class for the GNDS <**fluxes**> node that contains a list of flux nodes each as a 3d function.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Fluxes::Fluxes( ) :
        Suite( GIDI_fluxesChars ) {

}

/* *********************************************************************************************************//**
 * @param a_fileName            [in]    File containing a fluxes node to be parsed.
 ***********************************************************************************************************/

Fluxes::Fluxes( std::string const &a_fileName ) :
        Suite( GIDI_fluxesChars ) {

    addFile( a_fileName );
}

/* *********************************************************************************************************//**
* Adds the contents of the specified file to *this*.
 *
 ***********************************************************************************************************/

void Fluxes::addFile( std::string const &a_fileName ) {

    HAPI::File *doc = new HAPI::PugiXMLFile( a_fileName.c_str( ), "Fluxes::addFile" );

    HAPI::Node fluxes = doc->first_child( );

    std::string name( fluxes.name( ) );
    Construction::Settings construction( Construction::ParseMode::all, GIDI::Construction::PhotoMode::atomicOnly );

    SetupInfo setupInfo( nullptr );

    std::string formatVersionString = fluxes.attribute_as_string( GIDI_formatChars );
    if( formatVersionString == "" ) formatVersionString = GNDS_formatVersion_1_10Chars;
    LUPI::FormatVersion formatVersion;
    formatVersion.setFormat( formatVersionString );
    if( !formatVersion.supported( ) ) throw Exception( "unsupport GND format version" );
    setupInfo.m_formatVersion = formatVersion;

    for( HAPI::Node child = fluxes.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
        Functions::Function3dForm *function3d = data3dParse( construction, child, setupInfo, nullptr );

        add( function3d );
    }
    delete doc;
}

}
