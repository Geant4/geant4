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

/*! \class ParticleInfo
 * This class stores an abridged set of particle information from PoPI::Database as needed by GIDI.
 *
 * Note, the stored mass does not include the mass associated with nuclear excitation energy. This addition mass is
 * store in m_excitationEnergy.
 */

/* *********************************************************************************************************//**
 * @param a_ID                  [in]    The particle's PoPs ID.
 * @param a_pid                 [in]    The same as *a_id* unless particle is an alias, then the final particle's id.
 * @param a_mass                [in]    The particle's groud state mass. For nuclide and nucleus, this is the 
 * @param a_excitationEnergy    [in]    The particle's nuclear excitation energy.
 ***********************************************************************************************************/

ParticleInfo::ParticleInfo( std::string const &a_ID, std::string const &a_pid, double a_mass, double a_excitationEnergy ) :
        m_id( ParticleInfo::IDPortion( a_ID ) ),
        m_qualifier( ParticleInfo::qualifierPortion( a_ID ) ),
        m_pid( a_pid ),
        m_mass( a_mass, "amu" ),
        m_excitationEnergy( a_excitationEnergy, "MeV" ) {

}

/* *********************************************************************************************************//**
 * @param a_ID                      [in]    The particle's PoPs ID.
 * @param a_pops                    [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs            [in]    The internal PoPI::Database instance used to get particle indices and possibly other particle information.
 *                                          This is the <**PoPs**> node under the <**reactionSuite**> node.
 * @param a_requiredInGlobalPoPs    [in]    If *true*, the ID must be in *a_pops*.
 ***********************************************************************************************************/

ParticleInfo::ParticleInfo( std::string const &a_ID, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, bool a_requiredInGlobalPoPs ) :
        m_id( ParticleInfo::IDPortion( a_ID ) ),
        m_qualifier( ParticleInfo::qualifierPortion( a_ID ) ),
        m_pid( "" ),
        m_mass( -1, "amu" ),
        m_excitationEnergy( 0, "MeV" ) {

    PoPI::Base const *particleOrAlias = nullptr;       // Need to get the mass and nuclear excitation energy. Favor from internal PoPs if present.
    std::string energyUnit( "MeV" );

    if( a_pops.exists( m_id ) ) {
        particleOrAlias = &a_pops.get<PoPI::Base>( m_id );
        if( particleOrAlias->isAlias( ) ) {
            PoPI::Alias const *alias = static_cast<PoPI::Alias const *>( particleOrAlias );
            particleOrAlias = &a_pops.get<PoPI::Base>( alias->pid( ) );
        }
        m_pid = particleOrAlias->ID( ); }
    else {
        if( a_requiredInGlobalPoPs ) throw Exception( "ParticleInfo::ParticleInfo: required particle ID not in global PoPs: " + m_id );
    }

    if( a_internalPoPs.exists( m_id ) ) {
        particleOrAlias = &a_internalPoPs.get<PoPI::Base>( m_id );
        if( particleOrAlias->isAlias( ) ) {
            PoPI::Alias const *alias = static_cast<PoPI::Alias const *>( particleOrAlias );
            particleOrAlias = &a_pops.get<PoPI::Base>( alias->pidIndex( ) );
        }
    }

    if( particleOrAlias == nullptr ) throw Exception( "ParticleInfo::ParticleInfo: particle ID not in global PoPs: " + m_id );

    if( particleOrAlias->isParticle( ) ) {
        PoPI::Particle const &particle = static_cast<PoPI::Particle const &>( *particleOrAlias );

        try {
            m_mass = PhysicalQuantity( particle.massValue( "amu" ), "amu" ); }
        catch (...) {
            m_mass = PhysicalQuantity( -1, "amu" );
        }

        if( particle.isNuclide( ) ) {
            PoPI::Nuclide const &nuclide = static_cast<PoPI::Nuclide const &>( particle );

            m_excitationEnergy = PhysicalQuantity( nuclide.levelEnergy( energyUnit ), energyUnit );
        }
    }
}

/* *********************************************************************************************************//**
 * Copy constructor for ParticleInfo.
 *
 * @param a_particleInfo        [in]    ParticleInfo instance to copy.
 ***********************************************************************************************************/

ParticleInfo::ParticleInfo( ParticleInfo const &a_particleInfo ) :
        m_id( a_particleInfo.ID( ) ),
        m_qualifier( a_particleInfo.qualifier( ) ),
        m_pid( a_particleInfo.pid( ) ),
        m_mass( a_particleInfo.mass( ) ),
        m_excitationEnergy( a_particleInfo.excitationEnergy( ) ) {

}

/* *********************************************************************************************************//**
 * The assignment operator. This method sets the members of *this* to those of *a_rhs* except for those
 * not set by base classes.
 *
 * @param a_rhs                     [in]    Instance whose member are used to set the members of *this*.
 ***********************************************************************************************************/

ParticleInfo &ParticleInfo::operator=( ParticleInfo const &a_rhs ) {

    if( this != &a_rhs ) {
        m_id = a_rhs.ID( );
        m_qualifier = a_rhs.qualifier( );
        m_pid = a_rhs.pid( );
        m_mass = a_rhs.mass( );
        m_excitationEnergy = a_rhs.excitationEnergy( );
    }

    return( *this );
}

/* *********************************************************************************************************//**
 * Returns the particle's actual mass (i.e., its *m_mass* plus *m_excitationEnergy*) in unit of *a_unit*.
 *
 * @param a_unit            [in]    The requested unit for the returned mass.
 *
 * @return                          The mass in unit of *a_unit*.
 ***********************************************************************************************************/

double ParticleInfo::mass( LUPI_maybeUnused std::string const &a_unit ) const {

    return( PoPI_AMU2MeV_c2 * m_mass.value( ) );
}

/* *********************************************************************************************************//**
 * This static method returns the non-qualifier portion of a particle's id. For example, for the id "Th232{1s1/2}",
 * the string "Th232" is returned.
 *
 * @param a_ID              [in]    The id's whose non-qualifier portion is to be returned.
 *
 * @return                          The non-qualifier portion of the particle's id.
 ***********************************************************************************************************/

std::string const ParticleInfo::IDPortion( std::string const &a_ID ) {

    std::string::size_type index1 = a_ID.find( "{" );

    std::string ID( a_ID, 0, index1 );
    return( ID );
}
/* *********************************************************************************************************//**
 * This static method returns the qualifier portion of a particle's id. For example, for the id "Th232{1s1/2}",
 * the string "1s1/2" is returned.
 *
 * @param a_ID              [in]    The id's whose qualifier portion is to be returned.
  *
 * @return                          The non-qualifier portion of the particle's id.
 ***********************************************************************************************************/

std::string const ParticleInfo::qualifierPortion( std::string const &a_ID ) {

    std::string::size_type index1 = a_ID.find( "{" );

    std::string qualifier( "" );
    if( index1 == std::string::npos ) return( qualifier );

    std::string::size_type index2 = a_ID.find( "}" );
    qualifier = std::string( a_ID, index1 + 1, index2 - index1 - 1 );

    return( qualifier );
}

}                           // End of namespace GIDI.
