/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "PoPI.hpp"

namespace PoPI {

#define PoPI_generationChars "generation"

/*! \class Lepton
 * This class represents **PoPs** lepton instance.
 */

/* *********************************************************************************************************//**
 * Constructor that parses an **HAPI** instance to create a **PoPs** lepton node.
 *
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_DB              [in]    The **PoPI::Database:: instance to add the constructed **Lepton** to.
 * @param a_parent          [in]    The parent suite that will contain *this*.
 ***********************************************************************************************************/

Lepton::Lepton( HAPI::Node const &a_node, Database *a_DB, LUPI_maybeUnused Database *a_parent ) :
        Particle( a_node, Particle_class::lepton, PoPI_leptonChars ),
        m_generation( a_node.attribute( PoPI_generationChars ).value( ) ) {

    if( ID( ).substr(0, 2) == IDs::electron ) setIntid( intidHelper( isAnti( ), Particle_class::lepton, 0 ) );

    addToDatabase( a_DB );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Lepton::~Lepton( ) {

}

/* *********************************************************************************************************//**
 * Returns the generation attribute.
 ***********************************************************************************************************/

std::string Lepton::toXMLListExtraAttributes( void ) const {

    return( std::string( " generation=\"" + m_generation + "\"" ) );
}

}
