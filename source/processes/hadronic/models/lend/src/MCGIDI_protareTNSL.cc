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

/* *********************************************************************************************************//**
 * For internal use only. Function to determine the GIDI::Styles::TemperatureInfos from a GIDI::ProtareSingle that contains TNSL data,
 * protare given the parameters specified in the ransporting::MC *a_settings*.
 *
 * @param a_protare             [in]    A GIDI::ProtareSingle that contains TNSL data.
 * @param a_settings            [in]    Used to pass user options to the *this* to instruct it which data are desired.
 *
 * @return              
 ***********************************************************************************************************/

static LUPI_HOST GIDI::Styles::TemperatureInfos TNSL_temperatureInfos( GIDI::ProtareSingle const &a_protare, LUPI_maybeUnused Transporting::MC &a_settings ) {

    return( a_protare.temperatures( ) );
}

/*! \class ProtareTNSL
 * Class representing a **GNDS** <**reactionSuite**> node with only data needed for Monte Carlo transport. The
 * data are also stored in a way that is better suited for Monte Carlo transport. For example, cross section data
 * for each reaction are not stored with its reaction, but within the HeatedCrossSections member of the Protare.
 */

/* *********************************************************************************************************//**
 * Default constructor used when broadcasting a Protare as needed by MPI or GPUs.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProtareTNSL::ProtareTNSL( ) :
        Protare( ProtareType::TNSL ),
        m_numberOfTNSLReactions( 0 ),
        m_TNSL_maximumEnergy( 0.0 ),
        m_TNSL_maximumTemperature( 0.0 ),
        m_protareWithElastic( nullptr ),
        m_TNSL( nullptr ),
        m_protareWithoutElastic( nullptr ) {

}

/* *********************************************************************************************************//**
 * @param a_smr                         [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_protare                     [in]    The GIDI::Protare whose data is to be used to construct *this*.
 * @param a_pops                        [in]    A PoPs Database instance used to get particle intids and possibly other particle information.
 * @param a_settings                    [in]    Used to pass user options to the *this* to instruct it which data are desired.
 * @param a_particles                   [in]    List of transporting particles and their information (e.g., multi-group boundaries and fluxes).
 * @param a_domainHash                  [in]    The hash data used when looking up a cross section.
 * @param a_temperatureInfos            [in]    The list of temperature data to extract from *a_protare*.
 * @param a_reactionsToExclude          [in]    A list of reaction to not include in the MCGIDI::Protare. This currently does not work for ProtareTNSL.
 * @param a_reactionsToExcludeOffset    [in]    The starting index for the reactions in this ProtareSingle.
 * @param a_allowFixedGrid              [in]    For internal (i.e., MCGIDI) use only. Users must use the default value.
 ***********************************************************************************************************/

LUPI_HOST ProtareTNSL::ProtareTNSL( LUPI::StatusMessageReporting &a_smr, GIDI::ProtareTNSL const &a_protare, PoPI::Database const &a_pops, Transporting::MC &a_settings, 
                GIDI::Transporting::Particles const &a_particles, DomainHash const &a_domainHash, GIDI::Styles::TemperatureInfos const &a_temperatureInfos,
                std::set<int> const &a_reactionsToExclude, int a_reactionsToExcludeOffset, LUPI_maybeUnused bool a_allowFixedGrid ) :
        Protare( ProtareType::TNSL, a_protare, a_settings, a_pops ),
        m_protareWithElastic( static_cast<ProtareSingle *>( protareFromGIDIProtare( a_smr, *a_protare.protare( ), a_pops, a_settings, a_particles, 
            a_domainHash, a_temperatureInfos, a_reactionsToExclude, a_reactionsToExcludeOffset, false ) ) ),
        m_TNSL( static_cast<ProtareSingle *>( protareFromGIDIProtare( a_smr, *a_protare.TNSL( ), a_pops, a_settings, a_particles, a_domainHash, 
            TNSL_temperatureInfos( *a_protare.TNSL( ), a_settings ), a_reactionsToExclude, a_reactionsToExcludeOffset + static_cast<int>( m_protareWithElastic->numberOfReactions( ) ), false ) ) ),
        m_protareWithoutElastic( nullptr )  {

    std::set<int> reactionsToExclude( a_reactionsToExclude );

    reactionsToExclude.insert( 0 );
    m_protareWithoutElastic = static_cast<ProtareSingle *>( protareFromGIDIProtare( a_smr, *a_protare.protare( ), a_pops, a_settings, a_particles, a_domainHash, a_temperatureInfos, reactionsToExclude ) );

    m_numberOfTNSLReactions = m_TNSL->numberOfReactions( );
    m_TNSL_maximumEnergy = m_TNSL->maximumEnergy( );
    m_TNSL_maximumTemperature = m_TNSL->temperatures( ).back( );

    std::set<int> product_intids;
    std::set<int> product_intids_transportable;
    std::set<int> product_indices;
    std::set<int> product_indices_transportable;

    addVectorItemsToSet( m_TNSL->productIntids( false ), product_intids );
    addVectorItemsToSet( m_protareWithElastic->productIntids( true  ), product_intids_transportable );
    addVectorItemsToSet( m_TNSL->productIndices( false ), product_indices );
    addVectorItemsToSet( m_protareWithElastic->productIndices( true  ), product_indices_transportable );
    productIntidsAndIndices( product_intids, product_intids_transportable, product_indices, product_indices_transportable );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProtareTNSL::~ProtareTNSL( ) {

    delete m_protareWithElastic;
    delete m_TNSL;
    delete m_protareWithoutElastic;
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs index *a_particleIndex*.
 *
 * @param a_particleIndex       [in]    The PoPs index of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/
 
LUPI_HOST void ProtareTNSL::setUserParticleIndex2( int a_particleIndex, int a_userParticleIndex ) {

    m_protareWithElastic->setUserParticleIndex( a_particleIndex, a_userParticleIndex );
    m_TNSL->setUserParticleIndex( a_particleIndex, a_userParticleIndex );
    m_protareWithoutElastic->setUserParticleIndex( a_particleIndex, a_userParticleIndex );
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs intid *a_particleIntid*.
 *
 * @param a_particleIntid       [in]    The PoPs intid of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void ProtareTNSL::setUserParticleIndexViaIntid2( int a_particleIntid, int a_userParticleIndex ) {

    m_protareWithElastic->setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
    m_TNSL->setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
    m_protareWithoutElastic->setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
}


/* *********************************************************************************************************//**
 * Returns the pointer representing the (a_index - 1)th **ProtareSingle**.
 *
 * @param a_index               [in]    Index of the **ProtareSingle** to return. Can only be 0 or 1.
 *
 * @return                              Pointer to the requested protare or nullptr if invalid *a_index*..
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProtareSingle const *ProtareTNSL::protare( std::size_t a_index ) const {

    if( a_index == 0 ) return( m_protareWithElastic );
    if( a_index == 1 ) return( m_TNSL );
    return( nullptr );
}

/* *********************************************************************************************************//**
 * Returns the pointer representing the (a_index - 1)th **ProtareSingle**.
 *
 * @param a_index               [in]    Index of the **ProtareSingle** to return. Can only be 0 or 1.
 *
 * @return                              Pointer to the requested protare or nullptr if invalid *a_index*..
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProtareSingle *ProtareTNSL::protare( std::size_t a_index ) {

    if( a_index == 0 ) return( m_protareWithElastic );
    if( a_index == 1 ) return( m_TNSL );
    return( nullptr );
}

/* *********************************************************************************************************//**
 * Returns the pointer to the **ProtareSingle** that contains the (a_index - 1)th reaction.
 *
 * @param a_index               [in]    Index of the reaction.
 *
 * @return                              Pointer to the requested protare or nullptr if invalid *a_index*..
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProtareSingle const *ProtareTNSL::protareWithReaction( int a_index ) const {

    int index = a_index - m_numberOfTNSLReactions;

    if( a_index < 0 ) return( nullptr );
    if( index < 0 ) return( m_TNSL );
    return( m_protareWithElastic->protareWithReaction( index ) );
}

/* *********************************************************************************************************//**
 * Returns the list of temperatures for the requested ProtareSingle.
 *
 * @param a_index               [in]    Index of the reqested ProtareSingle. 
 *
 * @return                              Vector of doubles.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Vector<double> ProtareTNSL::temperatures( std::size_t a_index ) const {

    if( a_index == 0 ) return( m_protareWithElastic->temperatures( 0 ) );
    if( a_index == 1 ) return( m_TNSL->temperatures( 0 ) );

    LUPI_THROW( "ProtareSingle::temperatures: a_index not 0 or 1." );

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

LUPI_HOST_DEVICE Reaction const *ProtareTNSL::reaction( int a_index ) const {

    int index = a_index - m_numberOfTNSLReactions;

    if( index < 0 ) return( m_TNSL->reaction( a_index ) );
    return( m_protareWithElastic->reaction( index ) );
}

/* *********************************************************************************************************//**
 * Returns *true* if the reaction at index *a_index* has URR robability tables and false otherwise.
 *
 * @param           a_index [in]    The index of the reaction.
 *
 * @return                          *true* if the reaction has URR robability tables and false otherwise.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE bool ProtareTNSL::reactionHasURR_probabilityTables( int a_index ) const {

    int index = a_index - m_numberOfTNSLReactions;

    if( index < 0 ) return( false );
    return( m_protareWithElastic->reactionHasURR_probabilityTables( index ) );
}

/* *********************************************************************************************************//**
 * Returns the threshold for the reaction at index *a_index*. If *a_index* is negative, it is set to 0 before the
 * threshold in the regular protare is returned.
 *
 * @param           a_index [in]    The index of the reaction.
 *
 * @return                          The threshold for reaction at index *a_index*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareTNSL::threshold( std::size_t a_index ) const {

    int index = a_index - m_numberOfTNSLReactions;

    if( index < 0 ) return( m_TNSL->threshold( a_index ) );
    return( m_protareWithElastic->threshold( index ) );
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

LUPI_HOST_DEVICE double ProtareTNSL::crossSection( URR_protareInfos const &a_URR_protareInfos, int a_hashIndex, double a_temperature, double a_energy, bool a_sampling ) const {

    double crossSection1 = 0.0;

    if( ( a_energy < m_TNSL_maximumEnergy ) && ( a_temperature <= m_TNSL_maximumTemperature ) ) {
        crossSection1 = m_TNSL->crossSection( a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, a_sampling ) +
                        m_protareWithoutElastic->crossSection( a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, a_sampling ); }
    else {
        crossSection1 = m_protareWithElastic->crossSection( a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, a_sampling );
    }

    return( crossSection1 );
}

/* *********************************************************************************************************//**
 * Adds the energy dependent, total cross section corresponding to the temperature *a_temperature* multiplied by *a_userFactor* to *a_crossSectionVector*.
 *
 * @param   a_temperature               [in]        Specifies the temperature of the material.
 * @param   a_userFactor                [in]        User factor which all cross sections are multiplied by.
 * @param   a_numberAllocated           [in]        The length of memory allocated for *a_crossSectionVector*.
 * @param   a_crossSectionVector        [in/out]    The energy dependent, total cross section to add cross section data to.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void ProtareTNSL::crossSectionVector( double a_temperature, double a_userFactor, int a_numberAllocated, double *a_crossSectionVector ) const {

    if( a_temperature <= m_TNSL_maximumTemperature ) {
        m_TNSL->crossSectionVector( a_temperature, a_userFactor, a_numberAllocated, a_crossSectionVector );
        m_protareWithoutElastic->crossSectionVector( a_temperature, a_userFactor, a_numberAllocated, a_crossSectionVector ); }
    else {
        m_protareWithElastic->crossSectionVector( a_temperature, a_userFactor, a_numberAllocated, a_crossSectionVector );
    }
}

/* *********************************************************************************************************//**
 * Returns the cross section for reaction at index *a_reactionIndex*, for target at temperature *a_temperature* and projectile of energy *a_energy*.
 *
 * @param   a_URR_protareInfos  [in]    URR information.
 * @param   a_reactionIndex     [in]    The index of the reaction.
 * @param   a_hashIndex         [in]    The cross section hash index.
 * @param   a_temperature       [in]    The target temperature.
 * @param   a_energy            [in]    The projectile energy.
 * @param   a_sampling          [in]    Only used for multi-group cross sections. When sampling, the cross section in the group where threshold 
 *                                      is present the cross section is augmented.
 *
 * @return                              The total cross section.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareTNSL::reactionCrossSection( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, int a_hashIndex, double a_temperature, double a_energy, bool a_sampling ) const {

    int index = a_reactionIndex - m_numberOfTNSLReactions;
    double crossSection1 = 0.0;

    if( ( a_energy < m_TNSL_maximumEnergy ) && ( a_temperature <= m_TNSL_maximumTemperature ) ) {
        if( index < 0 ) {
            crossSection1 = m_TNSL->reactionCrossSection( a_reactionIndex, a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, a_sampling ); }
        else {
            if( index > 0 ) crossSection1 = m_protareWithElastic->reactionCrossSection( index, a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, a_sampling );
        } }
    else {
        if( index >= 0 ) crossSection1 = m_protareWithElastic->reactionCrossSection( index, a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, a_sampling );
    }

    return( crossSection1 );
}

/* *********************************************************************************************************//**
 * Returns the cross section for reaction at index *a_reactionIndex*, for target at temperature *a_temperature* and projectile of energy *a_energy*.
 *
 * @param   a_URR_protareInfos  [in]    URR information.
 * @param   a_reactionIndex     [in]    The index of the reaction.
 * @param   a_temperature       [in]    The target temperature.
 * @param   a_energy            [in]    The projectile energy.
 *
 * @return                              The total cross section.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareTNSL::reactionCrossSection( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, double a_temperature, double a_energy ) const {

    int index = a_reactionIndex - m_numberOfTNSLReactions;
    double crossSection1 = 0.0;

    if( ( a_energy < m_TNSL_maximumEnergy ) && ( a_temperature <= m_TNSL_maximumTemperature ) ) {
        if( index < 0 ) {
            crossSection1 = m_TNSL->reactionCrossSection( a_reactionIndex, a_URR_protareInfos, a_temperature, a_energy ); }
        else {
            if( index > 0 ) crossSection1 = m_protareWithElastic->reactionCrossSection( index, a_URR_protareInfos, a_temperature, a_energy );
        } }
    else {
        if( index >= 0 ) crossSection1 = m_protareWithElastic->reactionCrossSection( index, a_URR_protareInfos, a_temperature, a_energy );
    }

    return( crossSection1 );
}

/* *********************************************************************************************************//**
 * Returns the total deposition energy.
 *
 * @param a_hashIndex     [in]    The cross section hash index.
 * @param a_temperature   [in]    The target temperature.
 * @param a_energy        [in]    The projectile energy.
 *
 * @return                          The total deposition energy.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareTNSL::depositionEnergy( int a_hashIndex, double a_temperature, double a_energy ) const {

    double deposition_energy = 0.0;

    if( ( a_energy < m_TNSL_maximumEnergy ) && ( a_temperature <= m_TNSL_maximumTemperature ) ) {
        deposition_energy = m_TNSL->depositionEnergy( a_hashIndex, a_temperature, a_energy ) +
                        m_protareWithoutElastic->depositionEnergy( a_hashIndex, a_temperature, a_energy ); }
    else {
        deposition_energy = m_protareWithElastic->depositionEnergy( a_hashIndex, a_temperature, a_energy );
    }
        
    return( deposition_energy );
}

/* *********************************************************************************************************//**
 * Returns the total deposition momentum.
 *
 * @param a_hashIndex     [in]    The cross section hash index.
 * @param a_temperature   [in]    The target temperature.
 * @param a_energy        [in]    The projectile energy.
 *
 * @return                          The total deposition momentum.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareTNSL::depositionMomentum( int a_hashIndex, double a_temperature, double a_energy ) const {

    double deposition_momentum = 0.0;

    if( ( a_energy < m_TNSL_maximumEnergy ) && ( a_temperature <= m_TNSL_maximumTemperature ) ) {
        deposition_momentum = m_TNSL->depositionMomentum( a_hashIndex, a_temperature, a_energy ) +
                        m_protareWithoutElastic->depositionMomentum( a_hashIndex, a_temperature, a_energy ); }
    else {
        deposition_momentum = m_protareWithElastic->depositionMomentum( a_hashIndex, a_temperature, a_energy );
    }

    return( deposition_momentum );
}

/* *********************************************************************************************************//**
 * Returns the total production energy.
 *
 * @param a_hashIndex     [in]    The cross section hash index.
 * @param a_temperature   [in]    The target temperature.
 * @param a_energy        [in]    The projectile energy.
 *
 * @return                          The total production energy.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareTNSL::productionEnergy( int a_hashIndex, double a_temperature, double a_energy ) const {

    double production_energy = 0.0;

    if( ( a_energy < m_TNSL_maximumEnergy ) && ( a_temperature <= m_TNSL_maximumTemperature ) ) {
        production_energy = m_TNSL->productionEnergy( a_hashIndex, a_temperature, a_energy ) +
                        m_protareWithoutElastic->productionEnergy( a_hashIndex, a_temperature, a_energy ); }
    else {
        production_energy = m_protareWithElastic->productionEnergy( a_hashIndex, a_temperature, a_energy );
    }

    return( production_energy );
}

/* *********************************************************************************************************//**
 * Returns the gain for particle with index *a_particleIndex*.
 *
 * @param a_hashIndex           [in]    The cross section hash index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The projectile energy.
 * @param a_particleIndex       [in]    The index of the particle whose gain is to be returned.
 *
 * @return                      [in]    A vector of the length of the number of multi-group groups.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareTNSL::gain( int a_hashIndex, double a_temperature, double a_energy, int a_particleIndex ) const {

    double gain1 = 0.0;

    if( ( a_energy < m_TNSL_maximumEnergy ) && ( a_temperature <= m_TNSL_maximumTemperature ) ) {
        gain1 = m_TNSL->gain( a_hashIndex, a_temperature, a_energy, a_particleIndex ) +
                        m_protareWithoutElastic->gain( a_hashIndex, a_temperature, a_energy, a_particleIndex ); }
    else {
        gain1 = m_protareWithElastic->gain( a_hashIndex, a_temperature, a_energy, a_particleIndex );
    }

    return( gain1 );
}

/* *********************************************************************************************************//**
 * Returns the gain for particle with intid *a_particleIntid*.
 *
 * @param a_hashIndex           [in]    The cross section hash index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The projectile energy.
 * @param a_particleIntid       [in]    The intid of the particle whose gain is to be returned.
 *
 * @return                      [in]    A vector of the length of the number of multi-group groups.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareTNSL::gainViaIntid( int a_hashIndex, double a_temperature, double a_energy, int a_particleIntid ) const {

    double gain1 = 0.0;

    if( ( a_energy < m_TNSL_maximumEnergy ) && ( a_temperature <= m_TNSL_maximumTemperature ) ) {
        gain1 = m_TNSL->gainViaIntid( a_hashIndex, a_temperature, a_energy, a_particleIntid ) +
                        m_protareWithoutElastic->gainViaIntid( a_hashIndex, a_temperature, a_energy, a_particleIntid ); }
    else {
        gain1 = m_protareWithElastic->gainViaIntid( a_hashIndex, a_temperature, a_energy, a_particleIntid );
    }

    return( gain1 );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void ProtareTNSL::serialize2( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    int numberOfTNSLReactions = static_cast<int>( m_numberOfTNSLReactions );
    DATA_MEMBER_INT( numberOfTNSLReactions, a_buffer, a_mode );
    m_numberOfTNSLReactions = static_cast<std::size_t>( numberOfTNSLReactions );

    DATA_MEMBER_DOUBLE( m_TNSL_maximumEnergy, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_TNSL_maximumTemperature, a_buffer, a_mode );

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        if( a_buffer.m_placement != nullptr ) {
            m_protareWithElastic = new(a_buffer.m_placement) ProtareSingle;
            a_buffer.incrementPlacement( sizeof( ProtareSingle ) );
            m_TNSL = new(a_buffer.m_placement) ProtareSingle;
            a_buffer.incrementPlacement( sizeof( ProtareSingle ) );
            m_protareWithoutElastic = new(a_buffer.m_placement) ProtareSingle;
            a_buffer.incrementPlacement( sizeof( ProtareSingle ) ); }
        else {
            m_protareWithElastic = new ProtareSingle( );
            m_TNSL = new ProtareSingle( );
            m_protareWithoutElastic = new ProtareSingle( );
        }
    }
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
            a_buffer.incrementPlacement( sizeof( ProtareSingle ) );
            a_buffer.incrementPlacement( sizeof( ProtareSingle ) );
            a_buffer.incrementPlacement( sizeof( ProtareSingle ) );
    }
    m_protareWithElastic->serialize2( a_buffer, a_mode );
    m_TNSL->serialize2( a_buffer, a_mode );
    m_protareWithoutElastic->serialize2( a_buffer, a_mode );
}

}
