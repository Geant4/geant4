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

/*! \class GaugeBoson
 * This class represents **PoPs** gaugeBoson instance.
 */

/* *********************************************************************************************************//**
 * Constructor that parses an **HAPI** instance to create a **PoPs** gaugeBoson node.
 *
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_DB              [in]    The **PoPI::Database:: instance to add the constructed **GaugeBoson** to.
 * @param a_parent          [in]    The parent suite that will contain *this*.
 ***********************************************************************************************************/

GaugeBoson::GaugeBoson( HAPI::Node const &a_node, Database *a_DB, LUPI_maybeUnused Database *a_parent ) :
        Particle( a_node, Particle_class::gaugeBoson, PoPI_gaugeBosonChars ) {

    if( ID( ) == IDs::photon ) setIntid( intidHelper( isAnti( ), Particle_class::gaugeBoson, 0 ) );

    addToDatabase( a_DB );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

GaugeBoson::~GaugeBoson( ) {

}

}
