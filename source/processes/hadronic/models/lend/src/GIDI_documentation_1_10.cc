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

namespace Documentation_1_10 {

#define GIDI_nameChars "name"

/*! \class Documentation
 * This class supports storing **GNDS** 1.9 and 1.10 documentation.
 */

/*! \class Suite
 * This is the GIDI::Suite class but with a different parse function.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Suite::Suite( ) :
        GIDI::Suite( GIDI_documentations_1_10_Chars ) {

}

/* *********************************************************************************************************//**
 * Parse a GNDS 1.10 documentations node.
 *
 * @param a_node        [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo   [in]    Information create my the Protare constructor to help in parsing.
 ***********************************************************************************************************/

void Suite::parse( HAPI::Node const &a_node, SetupInfo &a_setupInfo ) {

    for( HAPI::Node child = a_node.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
        add( new Documentation( child, a_setupInfo, this ) );
    }
}

/*! \class Documentation
 */

/* *********************************************************************************************************//**
 *
 * @param a_node        [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo   [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent      [in]    The parent GIDI::Suite.
 * @return
 ***********************************************************************************************************/

Documentation::Documentation( HAPI::Node const &a_node, LUPI_maybeUnused SetupInfo &a_setupInfo, LUPI_maybeUnused GIDI::Suite *a_parent ) :
        Form( GIDI_documentationChars, FormType::generic, a_node.attribute_as_string( GIDI_nameChars ) ) {

    m_label = a_node.attribute_as_string( GIDI_nameChars );
    m_text = std::string( a_node.text().get() );
}

}

}
