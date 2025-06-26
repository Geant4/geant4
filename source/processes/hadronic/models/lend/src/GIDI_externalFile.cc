/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdlib.h>
#include "GIDI.hpp"

namespace GIDI {

/*! \class ExternalFile
 * This class represents the **GNDS** <**externalFile**> node.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

ExternalFile::ExternalFile( std::string const &a_label, std::string const &a_path ) :
        Form( GIDI_externalFileChars, FormType::externalFile, a_label ),
        m_path( a_path ) {

}

/* *********************************************************************************************************//**
 * @param a_node            [in]    The **HAPI::Node** to be parsed to construct a GeneralEvaporation2d instance.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

ExternalFile::ExternalFile( HAPI::Node const &a_node, SetupInfo &a_setupInfo, LUPI_maybeUnused GIDI::Suite *a_parent ) :
        Form( a_node, a_setupInfo, FormType::externalFile ),
        m_path( a_node.attribute_as_string( GIDI_pathChars ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

ExternalFile::~ExternalFile( ) {

}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML line that represent *this*.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void ExternalFile::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string attributes;

    attributes  = a_writeInfo.addAttribute( GIDI_labelChars, label( ) );
    attributes += a_writeInfo.addAttribute( GIDI_pathChars, path( ) );

    a_writeInfo.addNodeStarterEnder( a_indent, moniker( ), attributes );
}

}               // End of namespace GIDI.
