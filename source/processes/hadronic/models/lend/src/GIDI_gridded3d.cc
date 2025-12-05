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

/*! \class Gridded3d
 * Class for the GNDS <**gridded3d**> node.
 */

/* *********************************************************************************************************//**
 * Constructed from data in a <**product**> node.
 *
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the Gridded3d.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 ***********************************************************************************************************/

Gridded3d::Gridded3d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo ) :
        Function3dForm( a_construction, a_node, a_setupInfo, FormType::gridded3d ),
        m_data( a_node.child( GIDI_arrayChars ), a_setupInfo, a_construction.useSystem_strtod( ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Gridded3d::~Gridded3d( ) {

}

/* *********************************************************************************************************//**
 * Only for internal use. Called by ProtareTNSL instance to zero the lower energy multi-group data covered by the ProtareSingle that
 * contains the TNSL data covers the lower energy multi-group data.
 *
 * @param a_maxTNSL_index           [in]    All elements up to *a_maxTNSL_index* exclusive are zero-ed.
 ***********************************************************************************************************/

void Gridded3d::modifiedMultiGroupElasticForTNSL( LUPI_maybeUnused std::size_t a_maxTNSL_index ) {

}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Gridded3d::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    a_writeInfo.addNodeStarter( a_indent, moniker( ) );
    axes( ).toXMLList( a_writeInfo, indent2 );
    m_data.toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

}               // End namespace Functions.

}               // End namespace GIDI.
