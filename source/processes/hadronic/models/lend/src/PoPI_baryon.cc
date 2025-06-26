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

/*! \class Baryon
 * This class represents **PoPs** baryon instance.
 */

/* *********************************************************************************************************//**
 * Constructor that parses an **HAPI** instance to create a **PoPs** baryon node.
 *
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_DB              [in]    The **PoPI::Database:: instance to add the constructed **Baryon** to.
 * @param a_parent          [in]    The parent suite that will contain *this*.
 ***********************************************************************************************************/

Baryon::Baryon( HAPI::Node const &a_node, Database *a_DB, LUPI_maybeUnused Database *a_parent ) :
        Particle( a_node, Particle_class::baryon, PoPI_baryonChars ) {

    int baryonIndex = -1;

    if( baseId( ) == IDs::neutron ) {
        baryonIndex = 0; }
    if( baseId( ) == IDs::proton ) {
        baryonIndex = 1;
    }
    if( baryonIndex != -1 ) setIntid( intidHelper( isAnti( ), Particle_class::baryon, baryonIndex ) );

    addToDatabase( a_DB );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Baryon::~Baryon( ) {

}

}
