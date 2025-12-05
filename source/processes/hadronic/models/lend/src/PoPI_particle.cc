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

/*! \class Particle
 * The base class representing a particle.
 */

/* *********************************************************************************************************//**
 * @param a_node                        [in]    The **HAPI::Node** node to be parsed.
 * @param a_class                       [in]    The class of the particle.
 * @param a_family                      [in]    The family of the particle.
 * @param a_hasNucleus                  [in]    Indicates if the particle is or contains a nucleus. 0 = no, -1 = yes and 1 = is nucleus.
 ***********************************************************************************************************/

Particle::Particle( HAPI::Node const &a_node, Particle_class a_class, std::string const &a_family, int a_hasNucleus ) :
        IDBase( a_node, a_class ),
        m_baseId( "" ),
        m_family( a_family ),
        m_anti( "" ),
        m_hasNucleus( a_hasNucleus ),
        m_mass( a_node.child( PoPI_massChars ) ),
        m_spin( a_node.child( PoPI_spinChars ) ),
        m_parity( a_node.child( PoPI_parityChars ) ),
        m_charge( a_node.child( PoPI_chargeChars ) ),
        m_halflife( a_node.child( PoPI_halflifeChars ) ),
        m_decayData( a_node.child( PoPI_decayDataChars ) ) {

    m_baseId = baseAntiQualifierFromID( ID( ), m_anti );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Particle::~Particle( ) {

}

/* *********************************************************************************************************//**
 * Returns the mass of the particle in units of *a_unit*. Currently not fully implement  and does not support *a_unit*.
 *
 * @param a_unit                        [in]    The unit to return the mass in.
 *
 * @return                                      The mass in unit of *a_unit*.
 ***********************************************************************************************************/

double Particle::massValue( char const *a_unit ) const {

    if( m_mass.size( ) == 0 ) throw Exception( "Particle '" + ID( ) + "' does not have any mass data." );

    PQ_double const *pq_mass = dynamic_cast<PQ_double const *>( mass( )[0] );

    if( pq_mass == nullptr ) throw Exception( "Particle '" + ID( ) + "' does not have a PoPI::PQ_double mass." );
    return( pq_mass->value( a_unit ) );
}

/* *********************************************************************************************************//**
 * Adds the contents of *this* to *a_XMLList* where each item in *a_XMLList* is one line (without linefeeds) to output as an XML representation of *this*.
 *
 * @param a_XMLList                     [in]    The list to add an XML output representation of *this* to.
 * @param a_indent1                     [in]    The amount of indentation to added to each line added to *a_XMLList*.
 ***********************************************************************************************************/

void Particle::toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const {

    std::string indent2 = a_indent1 + "  ";

    std::string header = a_indent1 + "<" + family( ) + " id=\"" + ID( ) + "\"" + toXMLListExtraAttributes( ) + ">";
    a_XMLList.push_back( std::move( header ) );

    m_mass.toXMLList( a_XMLList, indent2 );
    m_spin.toXMLList( a_XMLList, indent2 );
    m_parity.toXMLList( a_XMLList, indent2 );
    m_charge.toXMLList( a_XMLList, indent2 );
    m_halflife.toXMLList( a_XMLList, indent2 );
    toXMLListExtraElements( a_XMLList, indent2 );
    m_decayData.toXMLList( a_XMLList, indent2 );

    appendXMLEnd( a_XMLList, family( ) );
}

/* *********************************************************************************************************//**
 * Currently there are no extra attributes to add. Ergo, returns an empty string.
 ***********************************************************************************************************/

std::string Particle::toXMLListExtraAttributes( void ) const {

    return( "" );
}

/* *********************************************************************************************************//**
 * Currently there are no extra child nodes to add.
 *
 * @param a_XMLList                     [in]    The list to add an XML output representation of *this* to.
 * @param a_indent1                     [in]    The amount of indentation to added to each line added to *a_XMLList*.
 ***********************************************************************************************************/

void Particle::toXMLListExtraElements( LUPI_maybeUnused std::vector<std::string> &a_XMLList, LUPI_maybeUnused std::string const &a_indent1 ) const {

    return;
}

}
