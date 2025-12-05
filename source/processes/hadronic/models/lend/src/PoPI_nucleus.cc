/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdlib.h>
#include <stdexcept>

#include "PoPI.hpp"

namespace PoPI {

#define PoPI_energyChars "energy"

/*! \class Nucleus
 * This class represents **PoPs** nucleus instance.
 */

/* *********************************************************************************************************//**
 * Constructor that parses an **HAPI** instance to create a **PoPs** nucleus node.
 *
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_DB              [in]    The **PoPI::Database:: instance to add the constructed **Nucleus** to.
 * @param a_nuclide         [in]    This nuclide instance that will contain *this*.
 ***********************************************************************************************************/

Nucleus::Nucleus( HAPI::Node const &a_node, Database *a_DB, Nuclide *a_nuclide ) :
        Particle( a_node, Particle_class::nucleus, PoPI_nucleusChars, -1 ),
        m_nuclide( a_nuclide ),
        m_Z( a_nuclide->Z( ) ),
        m_A( a_nuclide->A( ) ),
        m_levelName( a_node.attribute( PoPI_indexChars ).value( ) ),                // The string version of m_levelIndex.
        m_levelIndex( a_node.attribute( PoPI_indexChars ).as_int( ) ),              // The int version of m_levelName.
        m_energy( a_node.child( PoPI_energyChars ) ) {

    if( a_node.empty( ) ) throw Exception( "nuclide is missing nucleus" );

    int sign = ( isAnti( ) ? -1 : 1 );
    setIntid( sign * ( 1000 * ( 1000 * (levelIndex( ) + 500) + Z( ) ) + A( ) ) );

    addToDatabase( a_DB );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Nucleus::~Nucleus( ) {

}

/* *********************************************************************************************************//**
 * Returns the atomic ID of the parent nuclide.
 *
 * @return                              The atomic ID of the parent nuclide.
 ***********************************************************************************************************/

std::string const &Nucleus::atomsID( void ) const {

    return( m_nuclide->atomsID( ) );
}

/* *********************************************************************************************************//**
 * Returns the mass of the nucleus in units of *a_unit*. Currently not fully implement  and does not support *a_unit*.
 *
 * @param a_unit                        [in]    The unit to return the mass in.
 *
 * @return                                      The mass in unit of *a_unit*.
 ***********************************************************************************************************/

double Nucleus::massValue( char const *a_unit ) const {

    if( mass( ).size( ) > 0 ) {
        PQ_double const *pq_mass = dynamic_cast<PQ_double const *>( mass( )[0] );

        if( pq_mass == nullptr ) throw Exception( "Particle does not have a PoPI::PQ_double mass." );
        return( pq_mass->value( a_unit ) );
    }

// FIXME: still need to correct for electron masses and binding energy. Currently, an approximation is done.
    double bindingEnergy = 0.0;
    if( m_Z == 1 ) {
        bindingEnergy = 13.5981e-6; }
    else if( m_Z == 2 ) {
        bindingEnergy = 79.005e-6;
    }

    return( m_nuclide->massValue( a_unit ) - ( m_Z * PoPI_electronMass_MeV_c2 - bindingEnergy ) / PoPI_AMU2MeV_c2 );
}

/* *********************************************************************************************************//**
 * Returns the excitation energy of the nucleus in units of *a_unit*. Currently not fully implement  and does not support *a_unit*.
 *
 * @param a_unit                        [in]    The unit to return the mass in.
 *
 * @return                                      The mass in unit of *a_unit*.
 ***********************************************************************************************************/

double Nucleus::energy( std::string const &a_unit ) const {

    if( m_energy.size( ) == 0 ) {
        if( m_levelIndex != 0 )
        std::cerr << std::endl << "Particle " << ID( ) << " missing energy node, please report to PoPs maintainer. Using 0.0 and continuing." << std::endl;
        return( 0.0 );
    }
    PQ_double *pq = static_cast<PQ_double *>( m_energy[0] );
    if( pq->unit( ) == "eV" ) return( pq->value( ) * 1e-6 );        // Kludge until units are functional.
    return( pq->value( a_unit ) );
}

/* *********************************************************************************************************//**
 * Returns the index attribute.
 ***********************************************************************************************************/

std::string Nucleus::toXMLListExtraAttributes( void ) const {

    return( std::string( " index=\"" + m_levelName + "\"" ) );
}

/* *********************************************************************************************************//**
 * Added the *m_energy* stuff to *a_XMLList*.
 *
 * @param a_XMLList                     [in]    The list to add an XML output representation of *this* to.
 * @param a_indent1                     [in]    The amount of indentation to added to each line added to *a_XMLList*.
 ***********************************************************************************************************/

void Nucleus::toXMLListExtraElements( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const {

    m_energy.toXMLList( a_XMLList, a_indent1 );
}

}
