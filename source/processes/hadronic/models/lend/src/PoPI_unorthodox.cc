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

/*! \class Unorthodox
 * This class represents **PoPs** unorthodox instance.
 */

/* *********************************************************************************************************//**
 * Constructor that parses an **HAPI** instance to create a **PoPs** unorthodox node.
 *
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_DB              [in]    The **PoPI::Database:: instance to add the constructed **Unorthodox** to.
 * @param a_parent          [in]    The parent suite that will contain *this*.
 ***********************************************************************************************************/

Unorthodox::Unorthodox( HAPI::Node const &a_node, Database *a_DB, LUPI_maybeUnused Database *a_parent ) :
        Particle( a_node, Particle_class::unorthodox, PoPI_unorthodoxChars ) {

    if( ID( ) == IDs::FissionProductENDL99120 ) {
        setHasNucleus( true );
        setIntid( intidHelper( isAnti( ), Particle_class::ENDL_fissionProduct, 99120 ) ); }
    if( ID( ) == IDs::FissionProductENDL99125 ) {
        setHasNucleus( true );
        setIntid( intidHelper( isAnti( ), Particle_class::ENDL_fissionProduct, 99125 ) );
    }

    addToDatabase( a_DB );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Unorthodox::~Unorthodox( ) {

}

}
