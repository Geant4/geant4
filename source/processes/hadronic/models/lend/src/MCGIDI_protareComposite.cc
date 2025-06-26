/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <map>

#include "MCGIDI.hpp"

namespace MCGIDI {

/*! \class ProtareComposite
 * Class representing a **GNDS** <**reactionSuite**> node with only data needed for Monte Carlo transport. The
 * data are also stored in a way that is better suited for Monte Carlo transport. For example, cross section data
 * for each reaction are not stored with its reaction, but within the HeatedCrossSections member of the Protare.
 */

/* *********************************************************************************************************//**
 * Default constructor used when broadcasting a Protare as needed by MPI or GPUs.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProtareComposite::ProtareComposite( ) :
        Protare( ProtareType::composite ),
        m_numberOfReactions( 0 ),
        m_numberOfOrphanProducts( 0 ),
        m_minimumEnergy( 0.0 ),
        m_maximumEnergy( 0.0 ) {

}

/* *********************************************************************************************************//**
 * @param a_smr                         [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_protare                     [in]    The GIDI::Protare whose data is to be used to construct *this*.
 * @param a_pops                        [in]    A PoPs Database instance used to get particle intids and possibly other particle information.
 * @param a_settings                    [in]    Used to pass user options to the *this* to instruct it which data are desired.
 * @param a_particles                   [in]    List of transporting particles and their information (e.g., multi-group boundaries and fluxes).
 * @param a_domainHash                  [in]    The hash data used when looking up a cross section.
 * @param a_temperatureInfos            [in]    The list of temperature data to extract from *a_protare*.
 * @param a_reactionsToExclude          [in]    A list of reaction to not include in the MCGIDI::Protare. This currently does not work for ProtareComposite.
 * @param a_reactionsToExcludeOffset    [in]    The starting index for the reactions in this ProtareSingle.
 * @param a_allowFixedGrid              [in]    For internal (i.e., MCGIDI) use only. Users must use the default value.
 ***********************************************************************************************************/

LUPI_HOST ProtareComposite::ProtareComposite( LUPI::StatusMessageReporting &a_smr, GIDI::ProtareComposite const &a_protare, PoPI::Database const &a_pops, 
                Transporting::MC &a_settings, GIDI::Transporting::Particles const &a_particles, DomainHash const &a_domainHash, 
                GIDI::Styles::TemperatureInfos const &a_temperatureInfos, std::set<int> const &a_reactionsToExclude, int a_reactionsToExcludeOffset, LUPI_maybeUnused bool a_allowFixedGrid ) :
        Protare( ProtareType::composite, a_protare, a_settings, a_pops ),
        m_numberOfReactions( 0 ),
        m_numberOfOrphanProducts( 0 ) {

    std::vector<GIDI::Protare *> &protares = static_cast<std::vector<GIDI::Protare *> &>( const_cast<GIDI::ProtareComposite &>( a_protare ).protares( ) );
    std::size_t length = static_cast<std::size_t>( protares.size( ) );

    std::set<int> product_intids;
    std::set<int> product_intids_transportable;
    std::set<int> product_indices;
    std::set<int> product_indices_transportable;

    m_protares.resize( length );
    for( std::size_t i1 = 0; i1 < length; ++i1 ) {
        GIDI::Protare const *protare = protares[i1];

        m_protares[i1] = static_cast<ProtareSingle *>( protareFromGIDIProtare( a_smr, *protare, a_pops, a_settings, a_particles, a_domainHash, a_temperatureInfos, a_reactionsToExclude, a_reactionsToExcludeOffset, false ) );

        m_numberOfReactions += m_protares[i1]->numberOfReactions( );
        m_numberOfOrphanProducts += m_protares[i1]->numberOfOrphanProducts( );

        if( i1 == 0 ) {
            m_minimumEnergy = m_protares[0]->minimumEnergy( );
            m_maximumEnergy = m_protares[0]->maximumEnergy( );
        }
        if( m_protares[i1]->minimumEnergy( ) < m_minimumEnergy ) m_minimumEnergy = m_protares[i1]->minimumEnergy( );
        if( m_protares[i1]->maximumEnergy( ) > m_maximumEnergy ) m_maximumEnergy = m_protares[i1]->maximumEnergy( );

        addVectorItemsToSet( m_protares[i1]->productIntids( false ), product_intids );
        addVectorItemsToSet( m_protares[i1]->productIntids( true  ), product_intids_transportable );
        addVectorItemsToSet( m_protares[i1]->productIndices( false ), product_indices );
        addVectorItemsToSet( m_protares[i1]->productIndices( true  ), product_indices_transportable );

        a_reactionsToExcludeOffset += protare->numberOfReactions( );
    }

    productIntidsAndIndices( product_intids, product_intids_transportable, product_indices, product_indices_transportable );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProtareComposite::~ProtareComposite( ) {

    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );

    for( std::size_t i1 = 0; i1 < length; ++i1 ) delete m_protares[i1];
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs index *a_particleIndex*.
 *
 * @param a_particleIndex       [in]    The PoPs index of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void ProtareComposite::setUserParticleIndex2( int a_particleIndex, int a_userParticleIndex ) {

    for( auto iter = m_protares.begin( ); iter != m_protares.end( ); ++iter ) (*iter)->setUserParticleIndex( a_particleIndex, a_userParticleIndex );
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs intid *a_particleIntid*.
 *
 * @param a_particleIntid       [in]    The PoPs intid of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/
    
LUPI_HOST void ProtareComposite::setUserParticleIndexViaIntid2( int a_particleIntid, int a_userParticleIndex ) {
    
    for( auto iter = m_protares.begin( ); iter != m_protares.end( ); ++iter ) (*iter)->setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
}

/* *********************************************************************************************************//**
 * Returns the pointer representing the (a_index - 1)th **ProtareSingle**.
 *
 * @param a_index               [in]    Index of the **ProtareSingle** to return.
 *
 * @return                              Pointer to the requested protare or nullptr if invalid *a_index*..
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProtareSingle const *ProtareComposite::protare( std::size_t a_index ) const {

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        std::size_t number = m_protares[i1]->numberOfProtares( );

        if( number > a_index ) return( m_protares[i1]->protare( a_index ) );
        a_index -= number;
    }

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Returns the pointer representing the (a_index - 1)th **ProtareSingle**.
 *
 * @param a_index               [in]    Index of the **ProtareSingle** to return.
 *
 * @return                              Pointer to the requested protare or nullptr if invalid *a_index*..
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProtareSingle *ProtareComposite::protare( std::size_t a_index ) {

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        std::size_t number = m_protares[i1]->numberOfProtares( );

        if( number > a_index ) return( m_protares[i1]->protare( a_index ) );
        a_index -= number;
    }

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Returns the pointer to the **ProtareSingle** that contains the (a_index - 1)th reaction.
 *
 * @param a_index               [in]    Index of the reaction.
 *
 * @return                              Pointer to the requested protare or nullptr if invalid *a_index*..
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProtareSingle const *ProtareComposite::protareWithReaction( int a_index ) const {

    if( a_index < 0 ) return( nullptr );

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        int numberOfReactions = m_protares[i1]->numberOfReactions( );

        if( a_index < numberOfReactions ) return( m_protares[i1] );
        a_index -= numberOfReactions;
    }

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Returns the list of temperatures for the requested ProtareSingle.
 *
 * @param a_index               [in]    Index of the reqested ProtareSingle.
 *
 * @return                              Vector of doubles.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Vector<double> ProtareComposite::temperatures( std::size_t a_index ) const {

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        std::size_t number = m_protares[i1]->numberOfProtares( );

        if( number > a_index ) return( m_protares[i1]->temperatures( a_index ) );
        a_index -= number;
    }

    LUPI_THROW( "ProtareSingle::temperatures: a_index not in range." );

    Vector<double> temps;                           // Only to stop compilers from complaining.
    return( temps );
}

/* *********************************************************************************************************//**
 * Returns the reaction at index *a_index*. If *a_index* is negative, the reaction of the TNSL protare at index -*a_index* is
 * returned; otherwise, the reaction from the regular protare at index *a_index* is returned.
 *
 * @param           a_index [in]    The index of the reaction to return.
 *
 * @return                          The reaction at index *a_index*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Reaction const *ProtareComposite::reaction( int a_index ) const {

    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );

    for( std::size_t i1 = 0; i1 < length; ++i1 ) {
        int numberOfReactions = m_protares[i1]->numberOfReactions( );

        if( a_index < numberOfReactions ) return( m_protares[i1]->reaction( a_index ) );
        a_index -= numberOfReactions;
    }

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Returns the reaction at index *a_index*. If *a_index* is negative, the reaction of the TNSL protare at index -*a_index* is
 * returned; otherwise, the reaction from the regular protare at index *a_index* is returned.
 *
 * @param           a_index [in]    The index of the reaction to return.
 *
 * @return                          The reaction at index *a_index*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Reaction const *ProtareComposite::orphanProduct( int a_index ) const {

    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );

    for( std::size_t i1 = 0; i1 < length; ++i1 ) {
        int numberOfReactions = m_protares[i1]->numberOfOrphanProducts( );

        if( a_index < numberOfReactions ) return( m_protares[i1]->orphanProduct( a_index ) );
        a_index -= numberOfReactions;
    }

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Returns *true* is one of the protares has a fission channel and *false* otherwise.
 *
 * @return                          *true* is one of the protares has a fission channel and *false* otherwise.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE bool ProtareComposite::hasFission( ) const {

    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );

    for( std::size_t i1 = 0; i1 < length; ++i1 ) {
        if( m_protares[i1]->hasFission( ) ) return( true );
    }

    return( false );
}

/* *********************************************************************************************************//**
 * Returns true if *this* has a photoatomic incoherent doppler broadened reaction and false otherwise.
 *
 * @return                          *true* is one of the protares has a fission channel and *false* otherwise.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE bool ProtareComposite::hasIncoherentDoppler( ) const {

    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );

    for( std::size_t i1 = 0; i1 < length; ++i1 ) {
        if( m_protares[i1]->hasIncoherentDoppler( ) ) return( true );
    }

    return( false );
}

/* *********************************************************************************************************//**
 * Returns true if one the the sub-protares of *this* has a unresolved resonance region (URR) data and false otherwise.
 *
 * @return                              true is if *this* has a URR data.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE bool ProtareComposite::hasURR_probabilityTables( ) const {
   
    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );

    for( std::size_t i1 = 0; i1 < length; ++i1 ) {
        if( m_protares[i1]->hasURR_probabilityTables( ) ) return( true );
    }

    return( false );
}

/* *********************************************************************************************************//**
 * Returns the minimum energy for the unresolved resonance region (URR) domain. If no URR data present, returns -1.
 *
 * @return                              The energy or -1 if not URR data present.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareComposite::URR_domainMin( ) const {

    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );
    double URR_domain_min = 1e32;

    for( std::size_t i1 = 0; i1 < length; ++i1 ) {
        if( m_protares[i1]->hasURR_probabilityTables( ) ) {
            if( URR_domain_min > m_protares[i1]->URR_domainMin( ) ) URR_domain_min = m_protares[i1]->URR_domainMin( );
        }
    }

    if( URR_domain_min == 1e32 ) URR_domain_min = -1.0;

    return( URR_domain_min );
}

/* *********************************************************************************************************//**
 * Returns the maximum energy for the unresolved resonance region (URR) domain. If no URR data present, returns -1.
 *
 * @return                              true is if *this* has a URR data.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareComposite::URR_domainMax( ) const {

    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );
    double URR_domain_max = -1.0;

    for( std::size_t i1 = 0; i1 < length; ++i1 ) {
        if( m_protares[i1]->hasURR_probabilityTables( ) ) {
            if( URR_domain_max < m_protares[i1]->URR_domainMax( ) ) URR_domain_max = m_protares[i1]->URR_domainMax( );
            
        }
    }

    return( URR_domain_max );
}

/* *********************************************************************************************************//**
 * Returns *true* if the reaction at index *a_index* has URR robability tables and false otherwise.
 *
 * @param           a_index [in]    The index of the reaction.
 *
 * @return                          *true* if the reaction has URR robability tables and false otherwise.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE bool ProtareComposite::reactionHasURR_probabilityTables( int a_index ) const {

    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );

    for( std::size_t i1 = 0; i1 < length; ++i1 ) {
        int numberOfReactions = m_protares[i1]->numberOfReactions( );

        if( a_index < numberOfReactions ) return( m_protares[i1]->reactionHasURR_probabilityTables( a_index ) );
        a_index -= numberOfReactions;
    }

    return( false );
}

/* *********************************************************************************************************//**
 * Returns the threshold for the reaction at index *a_index*. If *a_index* is negative, it is set to 0 before the
 * threshold in the regular protare is returned.
 *
 * @param           a_index [in]    The index of the reaction.
 *
 * @return                          The threshold for reaction at index *a_index*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareComposite::threshold( std::size_t a_index ) const {

    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );

    for( std::size_t i1 = 0; i1 < length; ++i1 ) {
        std::size_t numberOfReactions = m_protares[i1]->numberOfReactions( );

        if( a_index < numberOfReactions ) return( m_protares[i1]->threshold( a_index ) );
        a_index -= numberOfReactions;
    }

    return( 0.0 );
}

/* *********************************************************************************************************//**
 * Returns the total cross section.
 * 
 * @param   a_URR_protareInfos  [in]    URR information.
 * @param   a_hashIndex         [in]    The cross section hash index.
 * @param   a_temperature       [in]    The target temperature.
 * @param   a_energy            [in]    The projectile energy.
 * @param   a_sampling          [in]    Only used for multi-group cross sections. When sampling, the cross section in the group where threshold 
 *                                      is present the cross section is augmented.
 *
 * @return                              The total cross section.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareComposite::crossSection( URR_protareInfos const &a_URR_protareInfos, int a_hashIndex, double a_temperature, double a_energy, bool a_sampling ) const {

    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );
    double cross_section = 0.0;

    for( std::size_t i1 = 0; i1 < length; ++i1 ) cross_section += m_protares[i1]->crossSection( a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, a_sampling );

    return( cross_section );
}

/* *********************************************************************************************************//**
 * Adds the energy dependent, total cross section corresponding to the temperature *a_temperature* multiplied by *a_userFactor* to *a_crossSectionVector*.
 *
 * @param   a_temperature               [in]        Specifies the temperature of the material.
 * @param   a_userFactor                [in]        User factor which all cross sections are multiplied by.
 * @param   a_numberAllocated           [in]        The length of memory allocated for *a_crossSectionVector*.
 * @param   a_crossSectionVector        [in/out]    The energy dependent, total cross section to add cross section data to.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void ProtareComposite::crossSectionVector( double a_temperature, double a_userFactor, int a_numberAllocated, double *a_crossSectionVector ) const {

    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );

    for( std::size_t i1 = 0; i1 < length; ++i1 ) m_protares[i1]->crossSectionVector( a_temperature, a_userFactor, a_numberAllocated, a_crossSectionVector );
}

/* *********************************************************************************************************//**
 * Returns the cross section for reaction at index *a_reactionIndex*.
 *
 * @param   a_reactionIndex     [in]    The index of the reaction.
 * @param   a_URR_protareInfos  [in]    URR information.
 * @param   a_hashIndex         [in]    The cross section hash index.
 * @param   a_temperature       [in]    The target temperature.
 * @param   a_energy            [in]    The projectile energy.
 * @param   a_sampling          [in]    Only used for multi-group cross sections. When sampling, the cross section in the group where threshold 
 *                                      is present the cross section is augmented.
 *
 * @return                              The total cross section.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareComposite::reactionCrossSection( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, int a_hashIndex, double a_temperature, double a_energy, bool a_sampling ) const {

    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );
    double cross_section = 0.0;

    for( std::size_t i1 = 0; i1 < length; ++i1 ) {
        int numberOfReactions = m_protares[i1]->numberOfReactions( );

        if( a_reactionIndex < numberOfReactions ) {
            cross_section = m_protares[i1]->reactionCrossSection( a_reactionIndex, a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, a_sampling );
            break;
        }
        a_reactionIndex -= numberOfReactions;
    }

    return( cross_section );
}

/* *********************************************************************************************************//**
 * Returns the cross section for reaction at index *a_reactionIndex*.
 *
 * @param   a_reactionIndex     [in]    The index of the reaction.
 * @param   a_URR_protareInfos  [in]    URR information.
 * @param   a_temperature       [in]    The target temperature.
 * @param   a_energy            [in]    The projectile energy.
 *
 * @return                              The total cross section.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareComposite::reactionCrossSection( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, double a_temperature, double a_energy ) const {

    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );
    double cross_section = 0.0;

    for( std::size_t i1 = 0; i1 < length; ++i1 ) {
        int numberOfReactions = m_protares[i1]->numberOfReactions( );

        if( a_reactionIndex < numberOfReactions ) {
            cross_section = m_protares[i1]->reactionCrossSection( a_reactionIndex, a_URR_protareInfos, a_temperature, a_energy );
            break;
        }
        a_reactionIndex -= numberOfReactions;
    }

    return( cross_section );
}

/* *********************************************************************************************************//**
 * Returns the total deposition energy.
 *
 * @param   a_hashIndex     [in]    The cross section hash index.
 * @param   a_temperature   [in]    The target temperature.
 * @param   a_energy        [in]    The projectile energy.
 *
 * @return                          The total deposition energy.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareComposite::depositionEnergy( int a_hashIndex, double a_temperature, double a_energy ) const {

    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );
    double deposition_energy = 0.0;

    for( std::size_t i1 = 0; i1 < length; ++i1 ) deposition_energy += m_protares[i1]->depositionEnergy( a_hashIndex, a_temperature, a_energy );

    return( deposition_energy );
}

/* *********************************************************************************************************//**
 * Returns the total deposition momentum.
 *
 * @param   a_hashIndex     [in]    The cross section hash index.
 * @param   a_temperature   [in]    The target temperature.
 * @param   a_energy        [in]    The projectile energy.
 *
 * @return                          The total deposition momentum.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareComposite::depositionMomentum( int a_hashIndex, double a_temperature, double a_energy ) const {

    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );
    double deposition_momentum = 0.0;

    for( std::size_t i1 = 0; i1 < length; ++i1 ) deposition_momentum += m_protares[i1]->depositionMomentum( a_hashIndex, a_temperature, a_energy );

    return( deposition_momentum );
}

/* *********************************************************************************************************//**
 * Returns the total production energy.
 *
 * @param   a_hashIndex     [in]    The cross section hash index.
 * @param   a_temperature   [in]    The target temperature.
 * @param   a_energy        [in]    The projectile energy.
 *
 * @return                          The total production energy.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareComposite::productionEnergy( int a_hashIndex, double a_temperature, double a_energy ) const {

    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );
    double production_energy = 0.0;

    for( std::size_t i1 = 0; i1 < length; ++i1 ) production_energy += m_protares[i1]->productionEnergy( a_hashIndex, a_temperature, a_energy );

    return( production_energy );
}

/* *********************************************************************************************************//**
 * Returns the gain for particle with index *a_particleIndex*. 
 *
 * @param a_hashIndex           [in]    The continuous energy hash or multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The projectile energy.
 * @param a_particleIndex       [in]    The index of the particle whose gain is to be returned.
 *
 * @return                      [in]    A vector of the length of the number of multi-group groups.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareComposite::gain( int a_hashIndex, double a_temperature, double a_energy, int a_particleIndex ) const {

    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );
    double gain1 = m_protares[0]->gain( a_hashIndex, a_temperature, a_energy, a_particleIndex );

    for( std::size_t i1 = 1; i1 < length; ++i1 ) gain1 += m_protares[i1]->gain( a_hashIndex, a_temperature, a_energy, a_particleIndex );

    return( gain1 );
}

/* *********************************************************************************************************//**
 * Returns the gain for particle with intid *a_particleIntid*.
 *
 * @param a_hashIndex           [in]    The continuous energy hash or multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The projectile energy.
 * @param a_particleIntid       [in]    The intid of the particle whose gain is to be returned.
 *
 * @return                      [in]    A vector of the length of the number of multi-group groups.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareComposite::gainViaIntid( int a_hashIndex, double a_temperature, double a_energy, int a_particleIntid ) const {

    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );
    double gain1 = m_protares[0]->gainViaIntid( a_hashIndex, a_temperature, a_energy, a_particleIntid );

    for( std::size_t i1 = 1; i1 < length; ++i1 ) gain1 += m_protares[i1]->gainViaIntid( a_hashIndex, a_temperature, a_energy, a_particleIntid );

    return( gain1 );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void ProtareComposite::serialize2( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    std::size_t vectorSize = m_protares.size( );
    int vectorSizeInt = static_cast<int>( vectorSize );
    LUPI::DataBuffer *workingBuffer = &a_buffer;

    DATA_MEMBER_INT( m_numberOfReactions, a_buffer, a_mode );
    DATA_MEMBER_INT( m_numberOfOrphanProducts, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_minimumEnergy, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_maximumEnergy, a_buffer, a_mode );

    DATA_MEMBER_INT( vectorSizeInt, *workingBuffer, a_mode );
    vectorSize = static_cast<std::size_t>( vectorSizeInt );

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        m_protares.resize( vectorSize, &(workingBuffer->m_placement) );
        for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
            if( workingBuffer->m_placement != nullptr ) {
                m_protares[vectorIndex] = new(workingBuffer->m_placement) ProtareSingle;
                workingBuffer->incrementPlacement( sizeof( ProtareSingle ) ); }
            else {
                m_protares[vectorIndex] = new ProtareSingle;
            }
        }
    }
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        a_buffer.m_placement += m_protares.internalSize();
        a_buffer.incrementPlacement( sizeof( ProtareSingle ) * vectorSize );
    }

    for( std::size_t i1 = 0; i1 < vectorSize; ++i1 ) m_protares[i1]->serialize2( a_buffer, a_mode );
}

}
