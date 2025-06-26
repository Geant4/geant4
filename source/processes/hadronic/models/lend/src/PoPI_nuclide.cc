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

// FIXME - Must be removed once unit conversion is supported.
#define AMU2MeV 931.494028

/*! \class Nuclide
 * This class represents **PoPs** nuclide instance.
 */

/* *********************************************************************************************************//**
 * Constructor that parses an **HAPI** instance to create a **PoPs** nuclide node.
 *
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_DB              [in]    The **PoPI::Database:: instance to add the constructed **Nuclide** to.
 * @param a_isotope         [in]    This isotope instance that will contain *this*.
 ***********************************************************************************************************/

Nuclide::Nuclide( HAPI::Node const &a_node, Database *a_DB, Isotope *a_isotope ) :
        Particle( a_node, Particle_class::nuclide, PoPI_nuclideChars, -1 ),
        m_isotope( a_isotope ),
        m_nucleus( a_node.child( PoPI_nucleusChars ), a_DB, this ),
        m_gammaDecayData( a_node.child( PoPI_gammaDecayDataChars ) ) {

    int sign = ( isAnti( ) ? -1 : 1 );
    setIntid( sign * ( 1000 * ( 1000 * levelIndex( ) + Z( ) ) + A( ) ) );

    addToDatabase( a_DB );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Nuclide::~Nuclide( ) {

}

/* *********************************************************************************************************//**
 * Returns of atomic number of the parent isotope.
 *
 * @return                      The atomic number of *this*.
 ***********************************************************************************************************/

int Nuclide::Z( void ) const {

    return( m_isotope->Z( ) );
}

/* *********************************************************************************************************//**
 * Returns of atomic mass number of the parent isotope.
 *
 * @return                      The atomic mass number of *this*.
 ***********************************************************************************************************/

int Nuclide::A( void ) const {

    return( m_isotope->A( ) );
}

/* *********************************************************************************************************//**
 * Returns of atomic symbol of the parent isotope.
 *
 * @return                      The atomic symbol of *this*.
 ***********************************************************************************************************/

std::string const &Nuclide::atomsID( void ) const {

    return( m_isotope->symbol( ) );
}

/* *********************************************************************************************************//**
 * Returns the mass suite for the first nuclide in the isotope containing *this*.
 *
 * @return                      A const reference to a PQ_suite.
 ***********************************************************************************************************/

PQ_suite const &Nuclide::baseMass( void ) const {

    return( (*m_isotope).nuclides( )[0].mass( ) );
}

/* *********************************************************************************************************//**
 * Returns the mass of the nuclide in units of *a_unit* including the nucleus excitation energy. 
 * Currently not fully implement  and does not support *a_unit*.
 *
 * @param a_unit                        [in]    The unit to return the mass in.
 *
 * @return                                      The mass in unit of *a_unit*.
 ***********************************************************************************************************/

double Nuclide::massValue( char const *a_unit ) const {

    std::string unit_c2( a_unit );
    unit_c2 += " * c**2";
    PQ_double const *pq_mass;

    if( mass( ).size( ) > 0 ) {
        pq_mass = dynamic_cast<PQ_double const *>( mass( )[0] ); }
    else {
        if( baseMass( ).size( ) == 0 ) throw Exception( "nuclide::massValue: no mass in level 0 for particle '" + ID( ) + "'." );
        pq_mass = dynamic_cast<PQ_double const *>( baseMass( )[0] );
    }
    double _mass = pq_mass->value( a_unit );

    double v_levelEnergy = levelEnergy( unit_c2 ) / AMU2MeV;

    return( _mass + v_levelEnergy );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

void Nuclide::calculateNuclideGammaBranchStateInfos( PoPI::Database const &a_pops, NuclideGammaBranchStateInfos &a_nuclideGammaBranchStateInfos,
                    bool a_alwaysAdd ) const {

    if( a_nuclideGammaBranchStateInfos.find( ID( ) ) != nullptr )
        return;

    NuclideGammaBranchStateInfo *nuclideGammaBranchStateInfo = new NuclideGammaBranchStateInfo( ID( ), intid( ), kind( ), m_nucleus.energy( "MeV" ) );

    if( m_gammaDecayData.rows( ) > 0 ) {
        m_gammaDecayData.calculateNuclideGammaBranchStateInfo( a_pops, *nuclideGammaBranchStateInfo ); }
    else {
        decayData( ).calculateNuclideGammaBranchStateInfo( a_pops, *nuclideGammaBranchStateInfo );    
    }

    if( ( nuclideGammaBranchStateInfo->branches( ).size( ) > 0 ) || a_alwaysAdd ) {
        a_nuclideGammaBranchStateInfos.add( nuclideGammaBranchStateInfo ); }
    else {
        delete nuclideGammaBranchStateInfo;
    }
}

/* *********************************************************************************************************//**
 * Added *m_nucleus* stuff to *a_XMLList*.
 ***********************************************************************************************************/

void Nuclide::toXMLListExtraElements( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const {

    m_nucleus.toXMLList( a_XMLList, a_indent1 );
}

}
