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
 * Returns the proper **MCGIDI** protare base on the type of **GIDI** protare.
 *
 * @param a_smr                         [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_protare                     [in]    The GIDI::Protare whose data is to be used to construct *this*.
 * @param a_pops                        [in]    A PoPs Database instance used to get particle intids and possibly other particle information.
 * @param a_settings                    [in]    Used to pass user options to the *this* to instruct it which data are desired.
 * @param a_particles                   [in]    List of transporting particles and their information (e.g., multi-group boundaries and fluxes).
 * @param a_domainHash                  [in]    The hash data used when looking up a cross section.
 * @param a_temperatureInfos            [in]    The list of temperature data to extract from *a_protare*.
 * @param a_reactionsToExclude          [in]    A list of reaction to not include in the MCGIDI::Protare.
 * @param a_reactionsToExcludeOffset    [in]    The starting index for the reactions in this ProtareSingle.
 * @param a_allowFixedGrid              [in]    For internal (i.e., MCGIDI) use only. Users must use the default value.
 ***********************************************************************************************************/

LUPI_HOST Protare *protareFromGIDIProtare( LUPI::StatusMessageReporting &a_smr, GIDI::Protare const &a_protare, PoPI::Database const &a_pops, Transporting::MC &a_settings, 
                GIDI::Transporting::Particles const &a_particles, DomainHash const &a_domainHash, GIDI::Styles::TemperatureInfos const &a_temperatureInfos, 
                GIDI::ExcludeReactionsSet const &a_reactionsToExclude, std::size_t a_reactionsToExcludeOffset, bool a_allowFixedGrid ) {

    Protare *protare( nullptr );

    if( a_protare.protareType( ) == GIDI::ProtareType::single ) {
        protare = new ProtareSingle( a_smr, static_cast<GIDI::ProtareSingle const &>( a_protare ), a_pops, a_settings, a_particles, a_domainHash, 
                a_temperatureInfos, a_reactionsToExclude, a_reactionsToExcludeOffset, a_allowFixedGrid ); }
    else if( a_protare.protareType( ) == GIDI::ProtareType::composite ) {
        protare = new ProtareComposite( a_smr, static_cast<GIDI::ProtareComposite const &>( a_protare ), a_pops, a_settings, a_particles, a_domainHash,
                a_temperatureInfos, a_reactionsToExclude, a_reactionsToExcludeOffset, false ); }
    else if( a_protare.protareType( ) == GIDI::ProtareType::TNSL ) {
        protare = new ProtareTNSL( a_smr, static_cast<GIDI::ProtareTNSL const &>( a_protare ), a_pops, a_settings, a_particles, a_domainHash,
                a_temperatureInfos, a_reactionsToExclude, a_reactionsToExcludeOffset, false );
    }

    return( protare );
}

/*! \class Protare
 * Base class for the *MCGIDI* protare classes.
 */

/* *********************************************************************************************************//**
 * @param a_protareType         [in]    The enum for the type of Protare (i.e., single, composite or TNSL).
 *
 * Default constructor used when broadcasting a Protare as needed by MPI or GPUs.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Protare::Protare( ProtareType a_protareType ) :
        m_protareType( a_protareType ),
        m_projectileID( ),
        m_projectileIntid( -1 ),
        m_projectileIndex( -1 ),
        m_projectileUserIndex( -1 ),
        m_projectileMass( 0.0 ),
        m_projectileExcitationEnergy( 0.0 ),

        m_targetID( ),
        m_targetIntid( -1 ),
        m_targetIndex( -1 ),
        m_targetUserIndex( -1 ),
        m_targetMass( 0.0 ),
        m_targetExcitationEnergy( 0.0 ),

        m_neutronIndex( -1 ),
        m_userNeutronIndex( -1 ),
        m_photonIndex( -1 ),
        m_userPhotonIndex( -1 ),

        m_evaluation( ),
        m_projectileFrame( GIDI::Frame::lab ),

        m_isTNSL_ProtareSingle( false ) {

}

/* *********************************************************************************************************//**
 * Default base Protare constructor.
 *
 * @param a_protareType         [in]    The enum for the type of Protare (i.e., single, composite or TNSL).
 * @param a_protare             [in]    The GIDI::Protare whose data is to be used to construct *this*.
 * @param a_settings            [in]    Used to pass user options to the *this* to instruct it which data are desired.
 ***********************************************************************************************************/

LUPI_HOST Protare::Protare( ProtareType a_protareType, GIDI::Protare const &a_protare, LUPI_maybeUnused Transporting::MC const &a_settings, PoPI::Database const &a_pops ) :
        m_protareType( a_protareType ),
        m_projectileID( a_protare.projectile( ).ID( ).c_str( ) ),
        m_projectileIntid( -1 ),
        m_projectileIndex( MCGIDI_popsIndex( a_pops, m_projectileID.c_str( ) ) ),
        m_projectileUserIndex( -1 ),
        m_projectileMass( a_protare.projectile( ).mass( "MeV/c**2" ) ),          // Includes nuclear excitation energy.
        m_projectileExcitationEnergy( a_protare.projectile( ).excitationEnergy( ).value( ) ),

        m_targetID( a_protare.target( ).ID( ).c_str( ) ),
        m_targetIntid( -1 ),
        m_targetIndex( -1 ),
        m_targetUserIndex( -1 ),
        m_targetMass( a_protare.target( ).mass( "MeV/c**2" ) ),                  // Includes nuclear excitation energy.
        m_targetExcitationEnergy( a_protare.target( ).excitationEnergy( ).value( ) ),

        m_neutronIndex( MCGIDI_popsIndex( a_pops, PoPI::IDs::neutron ) ),
        m_userNeutronIndex( -1 ),
        m_photonIndex( MCGIDI_popsIndex( a_pops, PoPI::IDs::photon ) ),
        m_userPhotonIndex( -1 ),

        m_evaluation( a_protare.evaluation( ).c_str( ) ),
        m_projectileFrame( a_protare.projectileFrame( ) ),
        m_productIntids( 0 ),
        m_userProductIndices( 0 ),
        m_productIntidsTransportable( 0 ),
        m_userProductIndicesTransportable( 0 ),
        m_isTNSL_ProtareSingle( a_protare.isTNSL_ProtareSingle( ) ) {

    PoPI::Database const &pops = a_protare.protare( 0 )->internalPoPs( );
    m_projectileIntid = MCGIDI_popsIntid( pops, m_projectileID.c_str( ) );

    if( a_protare.protare( 0 )->interaction( ) != GIDI_MapInteractionTNSLChars ) {
        m_targetIntid = MCGIDI_popsIntid( pops, m_targetID.c_str( ) );
        m_targetIndex = MCGIDI_popsIndex( a_pops,  m_targetID.c_str( ) );
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Protare::~Protare( ) {

}

/* *********************************************************************************************************//**
 * Sets *this* members *m_productIntids* and *m_productIntidsTransportable* to *a_intids* and *a_transportableIntids* respectively.
 * And, sets *this* members *m_productIndices* and *m_productIndicesTransportable* to *a_indices* and *a_transportableIndices* respectively.
 *
 * @param a_intids                  [out]   The list of intids for the outgoing particles (i.e., products).
 * @param a_transportableIntids     [in]    The list of transportable intids for the outgoing particles (i.e., products).
 * @param a_indices                 [in]    The list of indices for the outgoing particles (i.e., products).
 * @param a_transportableIndices    [in]    The list of transportable indices for the outgoing particles (i.e., products).
 ***********************************************************************************************************/

LUPI_HOST void Protare::productIntidsAndIndices( std::set<int> const &a_intids, std::set<int> const &a_transportableIntids,
                std::set<int> const &a_indices, std::set<int> const &a_transportableIndices ) {

    m_productIntids.reserve( a_intids.size( ) );
    m_userProductIndices.reserve( a_intids.size( ) );
    for( std::set<int>::const_iterator iter = a_intids.begin( ); iter != a_intids.end( ); ++iter ) {
        m_productIntids.push_back( *iter );
        m_userProductIndices.push_back( -1 );
    }

    m_productIndices.reserve( a_indices.size( ) );
    for( auto iter = a_indices.begin( ); iter != a_indices.end( ); ++iter ) m_productIndices.push_back( *iter );

    m_productIntidsTransportable.reserve( a_transportableIntids.size( ) );
    m_userProductIndicesTransportable.reserve( a_transportableIntids.size( ) );
    for( std::set<int>::const_iterator iter = a_transportableIntids.begin( ); iter != a_transportableIntids.end( ); ++iter ) {
        m_productIntidsTransportable.push_back( *iter );
        m_userProductIndicesTransportable.push_back( -1 );
    }

    m_productIndicesTransportable.reserve( a_transportableIndices.size( ) );
    for( auto iter = a_transportableIndices.begin( ); iter != a_transportableIndices.end( ); ++iter ) m_productIndicesTransportable.push_back( *iter );
}

/* *********************************************************************************************************//**
 * Returns the list product intids. If *a_transportablesOnly* is true, the list only includes transportable particle.
 *
 * @param a_transportablesOnly  [in]    If **true**, a reference to *m_productIntidsTransportable* is returned; otherwise a reference to *m_productIntids* is returned.
 ***********************************************************************************************************/

LUPI_HOST Vector<int> const &Protare::productIntids( bool a_transportablesOnly ) const {

    if( a_transportablesOnly ) return( m_productIntidsTransportable );
    return( m_productIntids );
}

/* *********************************************************************************************************//**
 * Returns the list product indices. If *a_transportablesOnly* is true, the list only includes transportable particle.
 *
 * @param a_transportablesOnly  [in]    If **true**, a reference to *m_productIndicesTransportable* is returned; otherwise a reference to *m_productIndices* is returned.
 ***********************************************************************************************************/

LUPI_HOST Vector<int> const &Protare::productIndices( bool a_transportablesOnly ) const {

    if( a_transportablesOnly ) return( m_productIndicesTransportable );
    return( m_productIndices );
}

/* *********************************************************************************************************//**
 * Returns the list user product indices. If *a_transportablesOnly* is true, the list only includes transportable particle.
 *
 * @param a_transportablesOnly  [in]    If **true**, a reference to *m_userProductIndicesTransportable* is returned; otherwise a reference to *m_userProductIndices* is returned.
 ***********************************************************************************************************/

LUPI_HOST Vector<int> const &Protare::userProductIndices( bool a_transportablesOnly ) const {

    if( a_transportablesOnly ) return( m_userProductIndicesTransportable );
    return( m_userProductIndices );
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs index *a_particleIndex*.
 *
 * @param a_particleIndex       [in]    The PoPs index of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void Protare::setUserParticleIndex( int a_particleIndex, int a_userParticleIndex ) {

    if( m_projectileIndex == a_particleIndex ) m_projectileUserIndex = a_userParticleIndex;
    if( m_targetIndex == a_particleIndex ) m_targetUserIndex = a_userParticleIndex;

    if( m_photonIndex == a_particleIndex ) m_userPhotonIndex = a_userParticleIndex;

    for( std::size_t i1 = 0; i1 < m_productIndices.size( ); ++i1 ) {
        if( m_productIndices[i1] == a_particleIndex ) m_userProductIndices[i1] = a_userParticleIndex;
    }

    for( std::size_t i1 = 0; i1 < m_productIndicesTransportable.size( ); ++i1 ) {
        if( m_productIndicesTransportable[i1] == a_particleIndex ) m_userProductIndicesTransportable[i1] = a_userParticleIndex;
    }

    switch( m_protareType ) {
    case ProtareType::single:
        static_cast<ProtareSingle *>( this )->setUserParticleIndex2( a_particleIndex, a_userParticleIndex );
        break;
    case ProtareType::composite:
        static_cast<ProtareComposite *>( this )->setUserParticleIndex2( a_particleIndex, a_userParticleIndex );
        break;
    case ProtareType::TNSL:
        static_cast<ProtareTNSL *>( this )->setUserParticleIndex2( a_particleIndex, a_userParticleIndex );
        break;
    }
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs intid *a_particleIntid*.
 *
 * @param a_particleIntid       [in]    The PoPs intid of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void Protare::setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex ) {

    if( m_projectileIntid == a_particleIntid ) m_projectileUserIndex = a_userParticleIndex;
    if( m_targetIntid == a_particleIntid ) m_targetUserIndex = a_userParticleIndex;

    if( PoPI::Intids::photon == a_particleIntid ) m_userPhotonIndex = a_userParticleIndex;

    for( std::size_t i1 = 0; i1 < m_productIntids.size( ); ++i1 ) {
        if( m_productIntids[i1] == a_particleIntid ) m_userProductIndices[i1] = a_userParticleIndex;
    }

    for( std::size_t i1 = 0; i1 < m_productIntidsTransportable.size( ); ++i1 ) {
        if( m_productIntidsTransportable[i1] == a_particleIntid ) m_userProductIndicesTransportable[i1] = a_userParticleIndex;
    }

    switch( m_protareType ) {
    case ProtareType::single:
        static_cast<ProtareSingle *>( this )->setUserParticleIndexViaIntid2( a_particleIntid, a_userParticleIndex );
        break;
    case ProtareType::composite:
        static_cast<ProtareComposite *>( this )->setUserParticleIndexViaIntid2( a_particleIntid, a_userParticleIndex );
        break;
    case ProtareType::TNSL:
        static_cast<ProtareTNSL *>( this )->setUserParticleIndexViaIntid2( a_particleIntid, a_userParticleIndex );
        break;
    }
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Protare::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    serializeCommon( a_buffer, a_mode );

    switch( m_protareType ) {
    case ProtareType::single:
        static_cast<ProtareSingle *>( this )->serialize2( a_buffer, a_mode );
        break;
    case ProtareType::composite:
        static_cast<ProtareComposite *>( this )->serialize2( a_buffer, a_mode );
        break;
    case ProtareType::TNSL:
        static_cast<ProtareTNSL *>( this )->serialize2( a_buffer, a_mode );
        break;
    }
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Protare::serialize2( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    serializeCommon( a_buffer, a_mode );

    switch( m_protareType ) {
    case ProtareType::single:
        static_cast<ProtareSingle *>( this )->serialize2( a_buffer, a_mode );
        break;
    case ProtareType::composite:
        LUPI_THROW( "Protare::serialize2:: Oops1, this should not happend." );
        break;
    case ProtareType::TNSL:
        LUPI_THROW( "Protare::serialize2:: Oops2, this should not happend." );
        break;
    }
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Protare::serializeCommon( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    int protareType1 = 0;
    if( a_mode != LUPI::DataBuffer::Mode::Unpack ) {
        switch( m_protareType ) {
        case ProtareType::single :
            break;
        case ProtareType::composite :
            protareType1 = 1;
            break;
        case ProtareType::TNSL :
            protareType1 = 2;
            break;
        }
    }
    DATA_MEMBER_INT( protareType1, a_buffer, a_mode );
    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        switch( protareType1 ) {
        case 0 :
            m_protareType = ProtareType::single;
            break;
        case 1 :
            m_protareType = ProtareType::composite;
            break;
        case 2 :
            m_protareType = ProtareType::TNSL;
            break;
        }
    }

    DATA_MEMBER_STRING( m_projectileID, a_buffer, a_mode );
    DATA_MEMBER_INT( m_projectileIntid, a_buffer, a_mode );
    DATA_MEMBER_INT( m_projectileIndex, a_buffer, a_mode );
    DATA_MEMBER_INT( m_projectileUserIndex, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_projectileMass, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_projectileExcitationEnergy, a_buffer, a_mode );

    DATA_MEMBER_STRING( m_targetID, a_buffer, a_mode );
    DATA_MEMBER_INT( m_targetIntid, a_buffer, a_mode );
    DATA_MEMBER_INT( m_targetIndex, a_buffer, a_mode );
    DATA_MEMBER_INT( m_targetUserIndex, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_targetMass, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_targetExcitationEnergy, a_buffer, a_mode );

    DATA_MEMBER_INT( m_photonIndex, a_buffer, a_mode );
    DATA_MEMBER_INT( m_userPhotonIndex, a_buffer, a_mode );

    DATA_MEMBER_STRING( m_evaluation, a_buffer, a_mode );

    int frame = 0;
    if( m_projectileFrame == GIDI::Frame::centerOfMass ) frame = 1;
    DATA_MEMBER_INT( frame, a_buffer, a_mode );
    m_projectileFrame = GIDI::Frame::lab;
    if( frame == 1 ) m_projectileFrame = GIDI::Frame::centerOfMass;

    DATA_MEMBER_VECTOR_INT( m_productIntids, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_INT( m_productIndices, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_INT( m_userProductIndices, a_buffer, a_mode );

    DATA_MEMBER_VECTOR_INT( m_productIntidsTransportable, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_INT( m_productIndicesTransportable, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_INT( m_userProductIndicesTransportable, a_buffer, a_mode );

    DATA_MEMBER_CAST( m_isTNSL_ProtareSingle, a_buffer, a_mode, bool );
}

/* *********************************************************************************************************//**
 * Returns the number of memory bytes used by *this*.
 *
 * @return                      The number of bytes used by *this*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE long Protare::sizeOf( ) const {

    long sizeOf1 = 0;

    switch( m_protareType ) {
    case ProtareType::single:
        sizeOf1 = static_cast<ProtareSingle const *>( this )->sizeOf2( );
        break;
    case ProtareType::composite: 
        sizeOf1 = static_cast<ProtareComposite const *>( this )->sizeOf2( );
        break;
    case ProtareType::TNSL:     
        sizeOf1 = static_cast<ProtareTNSL const *>( this )->sizeOf2( );
        break;
    }

    return( sizeOf1 );
}

/* *********************************************************************************************************//**
 * This method counts the number of bytes of memory allocated by *this*. 
 * This is an improvement to the internalSize() method of getting memory size.
 ***********************************************************************************************************/
LUPI_HOST_DEVICE long Protare::memorySize( ) {

    LUPI::DataBuffer buf;
    // Written this way for debugger to modify buf.m_placementStart here for easier double checking.
    buf.m_placement = buf.m_placementStart + sizeOf();
    serialize(buf, LUPI::DataBuffer::Mode::Memory);
    return( ( buf.m_placement - buf.m_placementStart ) + ( buf.m_sharedPlacement - buf.m_sharedPlacementStart ) );
}

/* *********************************************************************************************************//**
 * This method counts the number of bytes of memory allocated by *this* and puts it into a_totalMemory.
 * If shared memory is used, the size of shared memory is a_sharedMemory. If using shared memory,
 * the host code only needs to allocate (a_totalMemory - a_sharedMemory) in main memory.
 ***********************************************************************************************************/
LUPI_HOST_DEVICE void Protare::incrementMemorySize( long &a_totalMemory, long &a_sharedMemory ) {

    LUPI::DataBuffer buf;         // Written this way for debugger to modify buf.m_placementStart here for easier double checking.

    buf.m_placement = buf.m_placementStart + sizeOf( );
    serialize( buf, LUPI::DataBuffer::Mode::Memory );
    a_totalMemory += buf.m_placement - buf.m_placementStart;
    a_sharedMemory += buf.m_sharedPlacement - buf.m_sharedPlacementStart;
}

/* *********************************************************************************************************//**
 * Returns the number of protares contained in *this*.
 *
 * @return                              Integer number of protares.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE std::size_t Protare::numberOfProtares( ) const {

    std::size_t numberOfProtares2 = 0;

    switch( protareType( ) ) {
    case ProtareType::single:
        numberOfProtares2 = static_cast<ProtareSingle const *>( this )->numberOfProtares( );
        break;
    case ProtareType::composite:
        numberOfProtares2 = static_cast<ProtareComposite const *>( this )->numberOfProtares( );
        break;
    case ProtareType::TNSL:
        numberOfProtares2 = static_cast<ProtareTNSL const *>( this )->numberOfProtares( );
        break;
    }

    return( numberOfProtares2 );
}

/* *********************************************************************************************************//**
 * Returns the const pointer representing the protare at index *a_index*.
 *
 * @param a_index               [in]    Index of protare in *this*.
 *
 * @return                              Returns the const pointer representing the protare.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProtareSingle const *Protare::protare( std::size_t a_index ) const {

    ProtareSingle const *protare1 = nullptr;

    switch( protareType( ) ) {
    case ProtareType::single:
        protare1 = static_cast<ProtareSingle const *>( this )->protare( a_index );
        break;
    case ProtareType::composite:
        protare1 = static_cast<ProtareComposite const *>( this )->protare( a_index );
        break;
    case ProtareType::TNSL:
        protare1 = static_cast<ProtareTNSL const *>( this )->protare( a_index );
        break;
    }

    return( protare1 );
}

/* *********************************************************************************************************//**
 * Returns the pointer representing the protare at index *a_index*.
 *
 * @param a_index               [in]    Index of protare in *this*.
 *
 * @return                              Returns the pointer representing the protare.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProtareSingle *Protare::protare( std::size_t a_index ) {

    ProtareSingle *protare1 = nullptr;

    switch( protareType( ) ) {
    case ProtareType::single:
        protare1 = static_cast<ProtareSingle *>( this )->protare( a_index );
        break;
    case ProtareType::composite:
        protare1 = static_cast<ProtareComposite *>( this )->protare( a_index );
        break;
    case ProtareType::TNSL:
        protare1 = static_cast<ProtareTNSL *>( this )->protare( a_index );
        break;
    }

    return( protare1 );
}

/* *********************************************************************************************************//**
 * Returns the pointer to the **ProtareSingle** that contains the (a_index - 1)th reaction.
 *
 * @param a_index               [in]    Index of the reaction.
 *
 * @return                              Pointer to the requested protare or nullptr if invalid *a_index*..
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProtareSingle const *Protare::protareWithReaction( std::size_t a_index ) const {

    ProtareSingle const *protare1 = nullptr;

    switch( protareType( ) ) {
    case ProtareType::single:
        protare1 = static_cast<ProtareSingle const *>( this )->protareWithReaction( a_index );
        break;
    case ProtareType::composite:
        protare1 = static_cast<ProtareComposite const *>( this )->protareWithReaction( a_index );
        break;
    case ProtareType::TNSL:
        protare1 = static_cast<ProtareTNSL const *>( this )->protareWithReaction( a_index );
        break;
    }

    return( protare1 );
}

/* *********************************************************************************************************//**
 * Returns the minimum cross section domain for all reaction..
 *
 * @return                              Returns the minimum cross section domain.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double Protare::minimumEnergy( ) const {

    double minimumEnergy1 = 0.0;

    switch( protareType( ) ) {
    case ProtareType::single: 
        minimumEnergy1 = static_cast<ProtareSingle const *>( this )->minimumEnergy( );
        break;
    case ProtareType::composite:
        minimumEnergy1 = static_cast<ProtareComposite const *>( this )->minimumEnergy( );
        break;
    case ProtareType::TNSL:
        minimumEnergy1 = static_cast<ProtareTNSL const *>( this )->minimumEnergy( );
        break;
    }

    return( minimumEnergy1 );
}

/* *********************************************************************************************************//**
 * Returns the maximum cross section domain for all reaction..
 *
 * @return                              Returns the maximum cross section domain.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double Protare::maximumEnergy( ) const {

    double maximumEnergy1 = 0.0;

    switch( protareType( ) ) {
    case ProtareType::single:
        maximumEnergy1 = static_cast<ProtareSingle const *>( this )->maximumEnergy( );
        break;
    case ProtareType::composite:
        maximumEnergy1 = static_cast<ProtareComposite const *>( this )->maximumEnergy( );
        break;
    case ProtareType::TNSL:
        maximumEnergy1 = static_cast<ProtareTNSL const *>( this )->maximumEnergy( );
        break;
    }

    return( maximumEnergy1 );
}

/* *********************************************************************************************************//**
 * Returns the list of temperatures for the requested ProtareSingle.
 *
 * @param a_index               [in]    Index of the reqested ProtareSingle.
 *
 * @return                              Vector of doubles.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Vector<double> Protare::temperatures( std::size_t a_index ) const {

    ProtareSingle const *protareSingle = protare( a_index );
    return( protareSingle->temperatures( ) );
}

/* *********************************************************************************************************//**
 * Returns the projectile's multi-group boundaries that were read from the file (i.e., pre-collapse).
 *
 * @param a_index               [in]    Index of the reqested ProtareSingle.
 *
 * @return                              Vector of doubles.
 ***********************************************************************************************************/

LUPI_HOST Vector<double> const &Protare::projectileMultiGroupBoundaries( ) const {

    ProtareSingle const *protareSingle = protare( 0 );
    return( protareSingle->projectileMultiGroupBoundaries( ) );
}


/* *********************************************************************************************************//**
 * Returns the projectile's collapsed multi-group boundaries (i.e., those used for transport).
 *
 * @param a_index               [in]    Index of the reqested ProtareSingle.
 *
 * @return                              Vector of doubles.
 ***********************************************************************************************************/

LUPI_HOST Vector<double> const &Protare::projectileMultiGroupBoundariesCollapsed( ) const {

    ProtareSingle const *protareSingle = protare( 0 );
    return( protareSingle->projectileMultiGroupBoundariesCollapsed( ) );
}

/* *********************************************************************************************************//**
 * Returns the number of reactions of *this*.  
 *
 * @return                              Number of reactions of *this*. 
 ***********************************************************************************************************/

LUPI_HOST_DEVICE std::size_t Protare::numberOfReactions( ) const {

    std::size_t numberOfReactions1 = 0;

    switch( protareType( ) ) {
    case ProtareType::single: 
        numberOfReactions1 = static_cast<ProtareSingle const *>( this )->numberOfReactions( );
        break;
    case ProtareType::composite:
        numberOfReactions1 = static_cast<ProtareComposite const *>( this )->numberOfReactions( );
        break;
    case ProtareType::TNSL:
        numberOfReactions1 = static_cast<ProtareTNSL const *>( this )->numberOfReactions( );
        break;
    }

    return( numberOfReactions1 );
}

/* *********************************************************************************************************//**
 * Returns the reaction at index *a_index*. 
 *
 * @param           a_index [in]    The index of the reaction to return.
 *
 * @return                          The reaction at index *a_index*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Reaction const *Protare::reaction( std::size_t a_index ) const {

    Reaction const *reaction1 = nullptr;

    switch( protareType( ) ) {
    case ProtareType::single: 
        reaction1 = static_cast<ProtareSingle const *>( this )->reaction( a_index );
        break;
    case ProtareType::composite:
        reaction1 = static_cast<ProtareComposite const *>( this )->reaction( a_index );
        break;
    case ProtareType::TNSL:
        reaction1 = static_cast<ProtareTNSL const *>( this )->reaction( a_index );
        break;
    }

    return( reaction1 );
}

/* *********************************************************************************************************//**
 * Returns the number of orphanProducts of *this*.  
 *
 * @return                              Number of orphanProducts of *this*. 
 ***********************************************************************************************************/

LUPI_HOST_DEVICE std::size_t Protare::numberOfOrphanProducts( ) const {

    std::size_t numberOfReactions1 = 0;

    switch( protareType( ) ) {
    case ProtareType::single: 
        numberOfReactions1 = static_cast<ProtareSingle const *>( this )->numberOfOrphanProducts( );
        break;
    case ProtareType::composite:
        numberOfReactions1 = static_cast<ProtareComposite const *>( this )->numberOfOrphanProducts( );
        break;
    case ProtareType::TNSL:
        numberOfReactions1 = static_cast<ProtareTNSL const *>( this )->numberOfOrphanProducts( );
        break;
    }

    return( numberOfReactions1 );
}

/* *********************************************************************************************************//**
 * Returns the orphanProduct at index *a_index*. 
 *
 * @param           a_index [in]    The index of the orphanProduct to return.
 *
 * @return                          The orphanProduct at index *a_index*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Reaction const *Protare::orphanProduct( std::size_t a_index ) const {

    Reaction const *orphanProduct1 = nullptr;

    switch( protareType( ) ) {
    case ProtareType::single: 
        orphanProduct1 = static_cast<ProtareSingle const *>( this )->orphanProduct( a_index );
        break;
    case ProtareType::composite:
        orphanProduct1 = static_cast<ProtareComposite const *>( this )->orphanProduct( a_index );
        break;
    case ProtareType::TNSL:
        orphanProduct1 = static_cast<ProtareTNSL const *>( this )->orphanProduct( a_index );
        break;
    }

    return( orphanProduct1 );
}

/* *********************************************************************************************************//**
 * Returns true if *this* has a fission reaction and false otherwise.
 *
 * @return                              true is if *this* has a fission reaction and false otherwise.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE bool Protare::hasFission( ) const {

    bool hasFission1 = false;

    switch( protareType( ) ) {
    case ProtareType::single: 
        hasFission1 = static_cast<ProtareSingle const *>( this )->hasFission( );
        break;
    case ProtareType::composite:
        hasFission1 = static_cast<ProtareComposite const *>( this )->hasFission( );
        break;
    case ProtareType::TNSL:
        hasFission1 = static_cast<ProtareTNSL const *>( this )->hasFission( );
        break;
    }

    return( hasFission1 );
}

/* *********************************************************************************************************//**
 * Returns true if *this* has a photoatomic incoherent doppler broadened reaction and false otherwise.
 *
 * @return                              true is if *this* has a specified reaction and false otherwise.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE bool Protare::hasIncoherentDoppler( ) const {

    bool hasIncoherentDoppler1 = false;

    switch( protareType( ) ) {
    case ProtareType::single:
        hasIncoherentDoppler1 = static_cast<ProtareSingle const *>( this )->hasIncoherentDoppler( );
        break;
    case ProtareType::composite:
        hasIncoherentDoppler1 = static_cast<ProtareComposite const *>( this )->hasIncoherentDoppler( );
        break;
    case ProtareType::TNSL:
        hasIncoherentDoppler1 = static_cast<ProtareTNSL const *>( this )->hasIncoherentDoppler( );
        break;
    }

    return( hasIncoherentDoppler1 );
}

/* *********************************************************************************************************//**
 * Returns URR index of *this*.
 *
 * @return                              Integer URR index of *this*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE int Protare::URR_index( ) const {

    int URR_index1 = 0;

    switch( protareType( ) ) {
    case ProtareType::single: 
        URR_index1 = static_cast<ProtareSingle const *>( this )->URR_index( );
        break;
    case ProtareType::composite:
        URR_index1 = static_cast<ProtareComposite const *>( this )->URR_index( );
        break;
    case ProtareType::TNSL:
        URR_index1 = static_cast<ProtareTNSL const *>( this )->URR_index( );
        break;
    }

    return( URR_index1 );
}

/* *********************************************************************************************************//**
 * Returns **true** if *this* has URR probability tables and **false** otherwise.
 *
 * @return                              boolean.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE bool Protare::hasURR_probabilityTables( ) const {

    bool hasURR_probabilityTables1 = false;

    switch( protareType( ) ) {
    case ProtareType::single: 
        hasURR_probabilityTables1 = static_cast<ProtareSingle const *>( this )->hasURR_probabilityTables( );
        break;
    case ProtareType::composite:
        hasURR_probabilityTables1 = static_cast<ProtareComposite const *>( this )->hasURR_probabilityTables( );
        break;
    case ProtareType::TNSL:
        hasURR_probabilityTables1 = static_cast<ProtareTNSL const *>( this )->hasURR_probabilityTables( );
        break;
    }

    return( hasURR_probabilityTables1 );
}

/* *********************************************************************************************************//**
 * Returns the minimum energy for the unresolved resonance region (URR) domain. If no URR data present, returns -1.
 *
 * @return                              Minimum energy for the unresolved resonance region (URR) domain.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double Protare::URR_domainMin( ) const {

    double  URR_domainMin1 = 0.0;

    switch( protareType( ) ) {
    case ProtareType::single: 
        URR_domainMin1 = static_cast<ProtareSingle const *>( this )->URR_domainMin( );
        break;
    case ProtareType::composite:
        URR_domainMin1 = static_cast<ProtareComposite const *>( this )->URR_domainMin( );
        break;
    case ProtareType::TNSL:
        URR_domainMin1 = static_cast<ProtareTNSL const *>( this )->URR_domainMin( );
        break;
    }

    return( URR_domainMin1 );
}

/* *********************************************************************************************************//**
 * Returns the maximum energy for the unresolved resonance region (URR) domain. If no URR data present, returns -1.
 *
 * @return                              Maximum energy for the unresolved resonance region (URR) domain.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double Protare::URR_domainMax( ) const {

    double  URR_domainMax1 = 0.0;

    switch( protareType( ) ) {
    case ProtareType::single:
        URR_domainMax1 = static_cast<ProtareSingle const *>( this )->URR_domainMax( );
        break;
    case ProtareType::composite:
        URR_domainMax1 = static_cast<ProtareComposite const *>( this )->URR_domainMax( );
        break;
    case ProtareType::TNSL:
        URR_domainMax1 = static_cast<ProtareTNSL const *>( this )->URR_domainMax( );
        break;
    }

    return( URR_domainMax1 );
}

/* *********************************************************************************************************//**
 * Returns *true* if the reaction at index *a_index* has URR robability tables and *false* otherwise.
 *
 * @param           a_index [in]    The index of the reaction.
 *
 * @return                          *true* if the reaction has URR robability tables and *false* otherwise.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE bool Protare::reactionHasURR_probabilityTables( std::size_t a_index ) const {

    bool reactionHasURR_probabilityTables1 = false;

    switch( protareType( ) ) {
    case ProtareType::single: 
        reactionHasURR_probabilityTables1 = static_cast<ProtareSingle const *>( this )->reactionHasURR_probabilityTables( a_index );
        break;
    case ProtareType::composite:
        reactionHasURR_probabilityTables1 = static_cast<ProtareComposite const *>( this )->reactionHasURR_probabilityTables( a_index );
        break;
    case ProtareType::TNSL:
        reactionHasURR_probabilityTables1 = static_cast<ProtareTNSL const *>( this )->reactionHasURR_probabilityTables( a_index );
        break;
    }

    return( reactionHasURR_probabilityTables1 );
}

/* *********************************************************************************************************//**
 * Returns the threshold for the reaction at index *a_index*. If *a_index* is negative, it is set to 0 before the
 * threshold in the regular protare is returned.
 *
 * @param           a_index [in]    The index of the reaction.
 *
 * @return                          The threshold for reaction at index *a_index*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double Protare::threshold( std::size_t a_index ) const {

    double threshold1 = 0.0;

    switch( protareType( ) ) {
    case ProtareType::single: 
        threshold1 = static_cast<ProtareSingle const *>( this )->threshold( a_index );
        break;
    case ProtareType::composite:
        threshold1 = static_cast<ProtareComposite const *>( this )->threshold( a_index );
        break;
    case ProtareType::TNSL:
        threshold1 = static_cast<ProtareTNSL const *>( this )->threshold( a_index );
        break;
    }

    return( threshold1 );
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

LUPI_HOST_DEVICE double Protare::crossSection( URR_protareInfos const &a_URR_protareInfos, std::size_t a_hashIndex, double a_temperature, double a_energy, bool a_sampling ) const {

    double crossSection1 = 0.0;

    switch( protareType( ) ) {
    case ProtareType::single: 
        crossSection1 = static_cast<ProtareSingle const *>( this )->crossSection( a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, a_sampling );
        break;
    case ProtareType::composite:
        crossSection1 = static_cast<ProtareComposite const *>( this )->crossSection( a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, a_sampling );
        break;
    case ProtareType::TNSL:
        crossSection1 = static_cast<ProtareTNSL const *>( this )->crossSection( a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, a_sampling );
        break;
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

LUPI_HOST_DEVICE void Protare::crossSectionVector( double a_temperature, double a_userFactor, std::size_t a_numberAllocated, 
                double *a_crossSectionVector ) const {

    switch( protareType( ) ) {
    case ProtareType::single: 
        static_cast<ProtareSingle const *>( this )->crossSectionVector( a_temperature, a_userFactor, a_numberAllocated, a_crossSectionVector );
        break;
    case ProtareType::composite:
        static_cast<ProtareComposite const *>( this )->crossSectionVector( a_temperature, a_userFactor, a_numberAllocated, a_crossSectionVector );
        break;
    case ProtareType::TNSL:
        static_cast<ProtareTNSL const *>( this )->crossSectionVector( a_temperature, a_userFactor, a_numberAllocated, a_crossSectionVector );
        break;
    }
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

LUPI_HOST_DEVICE double Protare::reactionCrossSection( std::size_t a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, std::size_t a_hashIndex,
                double a_temperature, double a_energy, bool a_sampling ) const {

    double reactionCrossSection1 = 0.0;

    switch( protareType( ) ) {
    case ProtareType::single: 
        reactionCrossSection1 = static_cast<ProtareSingle const *>( this )->reactionCrossSection( a_reactionIndex, a_URR_protareInfos, a_hashIndex,
                a_temperature, a_energy, a_sampling );
        break;
    case ProtareType::composite:
        reactionCrossSection1 = static_cast<ProtareComposite const *>( this )->reactionCrossSection( a_reactionIndex, a_URR_protareInfos, a_hashIndex,
                a_temperature, a_energy, a_sampling );
        break;
    case ProtareType::TNSL:
        reactionCrossSection1 = static_cast<ProtareTNSL const *>( this )->reactionCrossSection( a_reactionIndex, a_URR_protareInfos, a_hashIndex,
                a_temperature, a_energy, a_sampling );
        break;
    }

    return( reactionCrossSection1 );
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

LUPI_HOST_DEVICE double Protare::reactionCrossSection( std::size_t a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, double a_temperature, double a_energy ) const {

    double reactionCrossSection1 = 0.0;

    switch( protareType( ) ) {
    case ProtareType::single: 
        reactionCrossSection1 = static_cast<ProtareSingle const *>( this )->reactionCrossSection( a_reactionIndex, a_URR_protareInfos, a_temperature, a_energy );
        break;
    case ProtareType::composite:
        reactionCrossSection1 = static_cast<ProtareComposite const *>( this )->reactionCrossSection( a_reactionIndex, a_URR_protareInfos, a_temperature, a_energy );
        break;
    case ProtareType::TNSL:
        reactionCrossSection1 = static_cast<ProtareTNSL const *>( this )->reactionCrossSection( a_reactionIndex, a_URR_protareInfos, a_temperature, a_energy  );
        break;
    }

    return( reactionCrossSection1 );
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

LUPI_HOST_DEVICE double Protare::depositionEnergy( std::size_t a_hashIndex, double a_temperature, double a_energy ) const {

    double depositionEnergy1 = 0.0;

    switch( protareType( ) ) {
    case ProtareType::single: 
        depositionEnergy1 = static_cast<ProtareSingle const *>( this )->depositionEnergy( a_hashIndex, a_temperature, a_energy );
        break;
    case ProtareType::composite:
        depositionEnergy1 = static_cast<ProtareComposite const *>( this )->depositionEnergy( a_hashIndex, a_temperature, a_energy );
        break;
    case ProtareType::TNSL:
        depositionEnergy1 = static_cast<ProtareTNSL const *>( this )->depositionEnergy( a_hashIndex, a_temperature, a_energy );
        break;
    }

    return( depositionEnergy1 );
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

LUPI_HOST_DEVICE double Protare::depositionMomentum( std::size_t a_hashIndex, double a_temperature, double a_energy ) const {

    double depositionMomentum1 = 0.0;

    switch( protareType( ) ) {
    case ProtareType::single: 
        depositionMomentum1 = static_cast<ProtareSingle const *>( this )->depositionMomentum( a_hashIndex, a_temperature, a_energy );
        break;
    case ProtareType::composite:
        depositionMomentum1 = static_cast<ProtareComposite const *>( this )->depositionMomentum( a_hashIndex, a_temperature, a_energy );
        break;
    case ProtareType::TNSL:
        depositionMomentum1 = static_cast<ProtareTNSL const *>( this )->depositionMomentum( a_hashIndex, a_temperature, a_energy );
        break;
    }

    return( depositionMomentum1 );
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

LUPI_HOST_DEVICE double Protare::productionEnergy( std::size_t a_hashIndex, double a_temperature, double a_energy ) const {

    double productionEnergy1 = 0.0;

    switch( protareType( ) ) {
    case ProtareType::single: 
        productionEnergy1 = static_cast<ProtareSingle const *>( this )->productionEnergy( a_hashIndex, a_temperature, a_energy );
        break;
    case ProtareType::composite:
        productionEnergy1 = static_cast<ProtareComposite const *>( this )->productionEnergy( a_hashIndex, a_temperature, a_energy );
        break;
    case ProtareType::TNSL:
        productionEnergy1 = static_cast<ProtareTNSL const *>( this )->productionEnergy( a_hashIndex, a_temperature, a_energy );
        break;
    }

    return( productionEnergy1 );
}

/* *********************************************************************************************************//**
 * Returns the gain for particle with index *a_particleIndex*.
 *
 * @param a_hashIndex           [in]    The continuous energy hash or multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The projectile energy.
 * @param a_particleIndex       [in]    The index of the particle whose gain is to be returned.
 *
 * @return                              A vector of the length of the number of multi-group groups.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double Protare::gain( std::size_t a_hashIndex, double a_temperature, double a_energy, int a_particleIndex ) const {

    double gain1 = 0.0;

    switch( protareType( ) ) {
    case ProtareType::single: 
        gain1 = static_cast<ProtareSingle const *>( this )->gain( a_hashIndex, a_temperature, a_energy, a_particleIndex );
        break;
    case ProtareType::composite:
        gain1 = static_cast<ProtareComposite const *>( this )->gain( a_hashIndex, a_temperature, a_energy, a_particleIndex );
        break;
    case ProtareType::TNSL:
        gain1 = static_cast<ProtareTNSL const *>( this )->gain( a_hashIndex, a_temperature, a_energy, a_particleIndex );
        break;
    }

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
 * @return                              A vector of the length of the number of multi-group groups.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double Protare::gainViaIntid( std::size_t a_hashIndex, double a_temperature, double a_energy, int a_particleIntid ) const {

    double gain1 = 0.0;

    switch( protareType( ) ) {
    case ProtareType::single:
        gain1 = static_cast<ProtareSingle const *>( this )->gainViaIntid( a_hashIndex, a_temperature, a_energy, a_particleIntid );
        break;
    case ProtareType::composite:
        gain1 = static_cast<ProtareComposite const *>( this )->gainViaIntid( a_hashIndex, a_temperature, a_energy, a_particleIntid );
        break;
    case ProtareType::TNSL:
        gain1 = static_cast<ProtareTNSL const *>( this )->gainViaIntid( a_hashIndex, a_temperature, a_energy, a_particleIntid );
        break;
    }

    return( gain1 );
}

/* *********************************************************************************************************//**
 * Returns a reference to the **m_upscatterModelAGroupVelocities**.
 *
 * @return                          Reference to the upscatter model A group velocities.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Vector<double> const &Protare::upscatterModelAGroupVelocities( ) const {

    Vector<double> const *upscatterModelAGroupVelocities1 = nullptr;
    switch( protareType( ) ) {
    case ProtareType::single: 
        upscatterModelAGroupVelocities1 = &static_cast<ProtareSingle const *>( this )->upscatterModelAGroupVelocities( );
        break;
    case ProtareType::composite:
        upscatterModelAGroupVelocities1 = &static_cast<ProtareComposite const *>( this )->upscatterModelAGroupVelocities( );
        break;
    case ProtareType::TNSL:
        upscatterModelAGroupVelocities1 = &static_cast<ProtareTNSL const *>( this )->upscatterModelAGroupVelocities( );
        break;
    }

    return( *upscatterModelAGroupVelocities1 );
}

/*! \class ProtareSingle
 * Class representing a **GNDS** <**reactionSuite**> node with only data needed for Monte Carlo transport. The
 * data are also stored in a way that is better suited for Monte Carlo transport. For example, cross section data
 * for each reaction are not stored with its reaction, but within the HeatedCrossSections member of the Protare.
 */

/* *********************************************************************************************************//**
 * Default constructor used when broadcasting a Protare as needed by MPI or GPUs.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProtareSingle::ProtareSingle( ) :
        Protare( ProtareType::single ),
        m_URR_index( -1 ),
        m_hasURR_probabilityTables( false ),
        m_URR_domainMin( -1.0 ),
        m_URR_domainMax( -1.0 ),
        m_upscatterModelASupported( false ),
        m_projectileMultiGroupBoundaries( 0 ),
        m_projectileMultiGroupBoundariesCollapsed( 0 ),
        m_reactions( 0 ),
        m_orphanProducts( 0 ) {

}

/* *********************************************************************************************************//**
 * @param a_smr                         [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_protare                     [in]    The GIDI::Protare whose data is to be used to construct *this*.
 * @param a_pops                        [in]    A PoPs Database instance used to get particle intids and possibly other particle information.
 * @param a_settings                    [in]    Used to pass user options to the *this* to instruct it which data are desired.
 * @param a_particles                   [in]    List of transporting particles and their information (e.g., multi-group boundaries and fluxes).
 * @param a_domainHash                  [in]    The hash data used when looking up a cross section.
 * @param a_temperatureInfos            [in]    The list of temperature data to extract from *a_protare*.
 * @param a_reactionsToExclude          [in]    A list of reaction to not include in the MCGIDI::Protare.
 * @param a_reactionsToExcludeOffset    [in]    The starting index for the reactions in this ProtareSingle.
 * @param a_allowFixedGrid              [in]    For internal (i.e., MCGIDI) use only. Users must use the default value.
 ***********************************************************************************************************/

LUPI_HOST ProtareSingle::ProtareSingle( LUPI::StatusMessageReporting &a_smr, GIDI::ProtareSingle const &a_protare, PoPI::Database const &a_pops, 
                Transporting::MC &a_settings, GIDI::Transporting::Particles const &a_particles, DomainHash const &a_domainHash, 
                GIDI::Styles::TemperatureInfos const &a_temperatureInfos, GIDI::ExcludeReactionsSet const &a_reactionsToExclude, 
                std::size_t a_reactionsToExcludeOffset, bool a_allowFixedGrid ) :
        Protare( ProtareType::single, a_protare, a_settings, a_pops ),
        m_interaction( a_protare.interaction( ).c_str( ) ),
        m_URR_index( -1 ),
        m_hasURR_probabilityTables( false ),
        m_URR_domainMin( -1.0 ),
        m_URR_domainMax( -1.0 ),
        m_domainHash( a_domainHash ),
        m_upscatterModelASupported( ( projectileIntid( ) != PoPI::Intids::photon ) &&
                                    ( projectileIntid( ) != PoPI::Intids::electron ) &&
                                    !isTNSL_ProtareSingle( ) ),
        m_projectileMultiGroupBoundaries( 0 ),
        m_projectileMultiGroupBoundariesCollapsed( 0 ),
        m_reactions( 0 ),
        m_orphanProducts( 0 ),
        m_isPhotoAtomic( a_protare.isPhotoAtomic( ) ),
        m_heatedCrossSections( ),
        m_heatedMultigroupCrossSections( ) {

    a_protare.updateReactionIndices( 0 );           // This is not correct as the offset should be passed as an arguent.
    PoPI::Database const &pops = a_protare.protare( 0 )->internalPoPs( );

    if( !a_protare.isPhotoAtomic( ) ) {
        std::set<std::string> incompleteParticles;
        a_protare.incompleteParticles( a_settings, incompleteParticles );
        for( auto particle = a_particles.particles( ).begin( ); particle != a_particles.particles( ).end( ); ++particle ) {
            if( incompleteParticles.count( particle->first ) != 0 ) {
                std::string message = "Requested particle '" + particle->first + "' is incomplete in '" + a_protare.realFileName( ) + "'.";
                if( a_settings.throwOnError( ) ) {
                    throw std::runtime_error( message.c_str( ) ); }
                else {
                    smr_setReportError2p( a_smr.smr( ), 0, 0, message.c_str( ) );
                }
            }
        }
    }

    SetupInfo setupInfo( *this, a_protare, a_pops, pops );
    setupInfo.m_formatVersion = a_protare.formatVersion( );
    setupInfo.m_GRIN_continuumGammas = a_protare.GRIN_continuumGammas2( );

    GIDI::Transporting::Particles particles;
    for( std::map<std::string, GIDI::Transporting::Particle>::const_iterator particle = a_particles.particles( ).begin( ); particle != a_particles.particles( ).end( ); ++particle ) {
        setupInfo.m_particleIntids[particle->first] = MCGIDI_popsIntid( pops, particle->first );
        setupInfo.m_particleIndices[particle->first] = MCGIDI_popsIndex( a_pops, particle->first );

        if( ( m_interaction == GIDI_MapInteractionAtomicChars ) && 
                !( ( particle->first == PoPI::IDs::photon ) || ( particle->first == PoPI::IDs::electron ) ) ) continue;
        particles.add( particle->second );
    }

    GIDI::Transporting::MG multiGroupSettings( a_settings.projectileID( ), GIDI::Transporting::Mode::MonteCarloContinuousEnergy, a_settings.delayedNeutrons( ) );
    multiGroupSettings.setThrowOnError( a_settings.throwOnError( ) );

    setupInfo.m_distributionLabel = a_temperatureInfos[0].griddedCrossSection( );

    a_settings.styles( &a_protare.styles( ) );

    switch( a_settings.crossSectionLookupMode( ) ) {
    case Transporting::LookupMode::Data1d::continuousEnergy :
        m_continuousEnergy = true;
        break;
    case Transporting::LookupMode::Data1d::multiGroup :
        m_continuousEnergy = false;
        if( a_settings.upscatterModel( ) == Sampling::Upscatter::Model::B ) {
            multiGroupSettings.setMode( GIDI::Transporting::Mode::multiGroupWithSnElasticUpScatter ); }
        else {
            multiGroupSettings.setMode( GIDI::Transporting::Mode::multiGroup );
        }
        break;
    default :
        throw std::runtime_error( "ProtareSingle::ProtareSingle: invalid lookupMode" );
    }
    m_fixedGrid = a_allowFixedGrid && ( a_protare.projectile( ).ID( ) == PoPI::IDs::photon ) && ( a_settings.fixedGridPoints( ).size( ) > 0 );

    setupNuclideGammaBranchStateInfos( setupInfo, a_protare, a_settings.makePhotonEmissionProbabilitiesOne( ),
            a_settings.zeroNuclearLevelEnergyWidth( ) );
    convertACE_URR_probabilityTablesFromGIDI( a_protare, a_settings,  setupInfo );

    if( ( a_settings.crossSectionLookupMode( ) == Transporting::LookupMode::Data1d::multiGroup ) || 
        ( a_settings.other1dDataLookupMode( ) == Transporting::LookupMode::Data1d::multiGroup ) ) {

        GIDI::Suite const *transportables = nullptr;
        if( setupInfo.m_formatVersion.major( ) > 1 ) {
            GIDI::Styles::HeatedMultiGroup const &heatedMultiGroup = *a_protare.styles( ).get<GIDI::Styles::HeatedMultiGroup>( a_settings.label( ) );
           transportables = &heatedMultiGroup.transportables( ); }
        else {
            std::vector<GIDI::Suite::const_iterator> tags = a_protare.styles( ).findAllOfMoniker( GIDI_multiGroupStyleChars );

            if( tags.size( ) != 1 ) throw std::runtime_error( "MCGIDI::ProtareSingle::ProtareSingle: What is going on here?" );
            GIDI::Styles::MultiGroup const &multiGroup = static_cast<GIDI::Styles::MultiGroup const &>( **tags[0] );
            transportables = &multiGroup.transportables( );
        }

        GIDI::Transportable const transportable = *transportables->get<GIDI::Transportable>( a_protare.projectile( ).ID( ) );
        m_projectileMultiGroupBoundaries = transportable.groupBoundaries( );
        GIDI::Transporting::Particle const *particle = a_particles.particle( a_protare.projectile( ).ID( ) );
        m_projectileMultiGroupBoundariesCollapsed = particle->multiGroup( ).boundaries( );
    }

    std::vector<GIDI::Reaction const *> GIDI_reactions;
    std::set<std::string> product_ids;
    std::set<std::string> product_ids_transportable;
    GIDI::Reaction const *nuclearPlusCoulombInterferenceReaction = nullptr;
    if( a_settings.nuclearPlusCoulombInterferenceOnly( ) ) nuclearPlusCoulombInterferenceReaction = a_protare.nuclearPlusCoulombInterferenceOnlyReaction( );

    for( std::size_t reactionIndex = 0; reactionIndex < a_protare.reactions( ).size( ); ++reactionIndex ) {
        if( a_reactionsToExclude.find( reactionIndex + a_reactionsToExcludeOffset ) != a_reactionsToExclude.end( ) ) continue;

        GIDI::Reaction const *GIDI_reaction = a_protare.reaction( reactionIndex );

        if( !GIDI_reaction->active( ) ) continue;

        if( m_continuousEnergy ) {
            if( GIDI_reaction->crossSectionThreshold( ) >= a_settings.energyDomainMax( ) ) continue; }
        else {
            GIDI::Vector multi_group_cross_section = GIDI_reaction->multiGroupCrossSection( a_smr, multiGroupSettings, a_temperatureInfos[0] );
            GIDI::Vector vector = GIDI::collapse( multi_group_cross_section, a_settings, a_particles, 0.0 );

            std::size_t i1 = 0;
            for( ; i1 < vector.size( ); ++i1 ) if( vector[i1] != 0.0 ) break;
            if( i1 == vector.size( ) ) continue;
        }
        if( a_settings.ignoreENDF_MT5( ) && ( GIDI_reaction->ENDF_MT( ) == 5 ) && ( a_reactionsToExclude.size( ) == 0 ) ) continue;

        GIDI_reaction->productIDs( product_ids, particles, false );
        GIDI_reaction->productIDs( product_ids_transportable, particles, true );

        if( a_settings.nuclearPlusCoulombInterferenceOnly( ) && GIDI_reaction->RutherfordScatteringPresent( ) ) {
            if( nuclearPlusCoulombInterferenceReaction != nullptr ) GIDI_reactions.push_back( nuclearPlusCoulombInterferenceReaction ); }
        else {
            GIDI_reactions.push_back( GIDI_reaction );
        }
    }

    bool zeroReactions = GIDI_reactions.size( ) == 0;   // Happens when all reactions are skipped in the prior loop.
    if( zeroReactions ) GIDI_reactions.push_back( a_protare.reaction( 0 ) );    // Special case where no reaction in the protare is wanted so the first one is used but its cross section is set to 0.0 at all energies.

    setupInfo.m_reactionType = Transporting::Reaction::Type::Reactions;
    m_reactions.reserve( GIDI_reactions.size( ) );
    for( auto GIDI_reaction = GIDI_reactions.begin( ); GIDI_reaction != GIDI_reactions.end( ); ++GIDI_reaction ) {
        setupInfo.m_reaction = *GIDI_reaction;
        setupInfo.m_isPairProduction = (*GIDI_reaction)->isPairProduction( );
        setupInfo.m_isPhotoAtomicIncoherentScattering = (*GIDI_reaction)->isPhotoAtomicIncoherentScattering( );
        setupInfo.m_initialStateIndex = -1;
        Reaction *reaction = new Reaction( **GIDI_reaction, setupInfo, a_settings, particles, a_temperatureInfos );
        setupInfo.m_initialStateIndices[(*GIDI_reaction)->label( )] = setupInfo.m_initialStateIndex;
        reaction->updateProtareSingleInfo( this, m_reactions.size( ) );
        m_reactions.push_back( reaction );
    }

    std::set<int> product_intids;
    std::set<int> product_indices;
    for( std::set<std::string>::iterator iter = product_ids.begin( ); iter != product_ids.end( ); ++iter ) {
        product_intids.insert( MCGIDI_popsIntid( pops, *iter ) );
        product_indices.insert( MCGIDI_popsIndex( a_pops, *iter ) );
    }
    std::set<int> product_intids_transportable;
    std::set<int> product_indices_transportable;
    for( std::set<std::string>::iterator iter = product_ids_transportable.begin( ); iter != product_ids_transportable.end( ); ++iter ) {
        product_intids_transportable.insert( MCGIDI_popsIntid( pops, *iter ) );
        product_indices_transportable.insert( MCGIDI_popsIndex( a_pops, *iter ) );
    }
    productIntidsAndIndices( product_intids, product_intids_transportable, product_indices, product_indices_transportable );

    if( a_settings.sampleNonTransportingParticles( ) || particles.hasParticle( PoPI::IDs::photon ) ) {
        setupInfo.m_reactionType = Transporting::Reaction::Type::OrphanProducts;
        m_orphanProducts.reserve( a_protare.orphanProducts( ).size( ) );
        std::vector< std::vector<std::size_t> > associatedOrphanProductIndices( m_reactions.size( ) );

        for( std::size_t orphanProductIndex = 0; orphanProductIndex < a_protare.orphanProducts( ).size( ); ++orphanProductIndex ) {
            GIDI::Reaction const *GIDI_reaction = a_protare.orphanProduct( orphanProductIndex );

            if( GIDI_reaction->crossSectionThreshold( ) >= a_settings.energyDomainMax( ) ) continue;

            setupInfo.m_reaction = GIDI_reaction;
            Reaction *orphanProductReaction = new Reaction( *GIDI_reaction, setupInfo, a_settings, particles, a_temperatureInfos );
            orphanProductReaction->updateProtareSingleInfo( this, m_orphanProducts.size( ) );
            m_orphanProducts.push_back( orphanProductReaction );

            GIDI::Functions::Reference1d const *reference( GIDI_reaction->crossSection( ).get<GIDI::Functions::Reference1d>( 0 ) );
            std::string xlink = reference->xlink( );
            GUPI::Ancestry const *ancestry = a_protare.findInAncestry( xlink );
            if( ancestry == nullptr ) throw std::runtime_error( "Could not find xlink for orphan product - 1." );
            ancestry = ancestry->ancestor( );
            if( ancestry == nullptr ) throw std::runtime_error( "Could not find xlink for orphan product - 2." );
            if( ancestry->moniker( ) != GIDI_crossSectionSumChars ) {
                ancestry = ancestry->ancestor( );
                if( ancestry == nullptr ) throw std::runtime_error( "Could not find xlink for orphan product - 3." );
            }
            GIDI::Sums::CrossSectionSum const *crossSectionSum = static_cast<GIDI::Sums::CrossSectionSum const *>( ancestry );
            GIDI::Sums::Summands const &summands = crossSectionSum->summands( );
            for( std::size_t i1 = 0; i1 < summands.size( ); ++i1 ) {
                GIDI::Sums::Summand::Base const *summand = summands[i1];

                ancestry = a_protare.findInAncestry( summand->href( ) );
                if( ancestry == nullptr ) throw std::runtime_error( "Could not find href for summand - 1." );
                ancestry = ancestry->ancestor( );
                if( ancestry == nullptr ) throw std::runtime_error( "Could not find href for summand - 2." );

                GIDI::Reaction const *GIDI_reaction2 = static_cast<GIDI::Reaction const *>( ancestry );
                for( std::size_t reactionIndex = 0; reactionIndex < m_reactions.size( ); ++reactionIndex ) {
                    std::string label( m_reactions[reactionIndex]->label( ).c_str( ) );

                    if( label == GIDI_reaction2->label( ) ) {
                        associatedOrphanProductIndices[reactionIndex].push_back( m_orphanProducts.size( ) - 1 );
                        break;
                    }
                }
            }
        }

        for( std::size_t reactionIndex = 0; reactionIndex < m_reactions.size( ); ++reactionIndex ) {
            Reaction *reaction = m_reactions[reactionIndex];
            std::size_t size = associatedOrphanProductIndices[reactionIndex].size( );
            if( size > 0 ) {
                std::vector<Product *> associatedOrphanProducts;
                for( std::size_t index1 = 0; index1 < size; ++index1 ) {
                    std::size_t associatedOrphanProductIndex = associatedOrphanProductIndices[reactionIndex][index1];
                    m_orphanProducts[associatedOrphanProductIndex]->addOrphanProductToProductList( associatedOrphanProducts );
                }
                reaction->setOrphanProductData( associatedOrphanProductIndices[reactionIndex], associatedOrphanProducts );
            }
        }
    }

    std::vector<GIDI::Reaction const *> GIDI_orphanProducts;
    for( std::size_t reactionIndex = 0; reactionIndex < a_protare.orphanProducts( ).size( ); ++reactionIndex ) {
        GIDI::Reaction const *GIDI_reaction = a_protare.orphanProduct( reactionIndex );

        if( GIDI_reaction->crossSectionThreshold( ) >= a_settings.energyDomainMax( ) ) continue;
        GIDI_orphanProducts.push_back( GIDI_reaction );
    }

    bool removeContinuousEnergyData = false;
    if( m_continuousEnergy ) {
        m_heatedCrossSections.update( a_smr, setupInfo, a_settings, particles, a_domainHash, a_temperatureInfos, GIDI_reactions, GIDI_orphanProducts,
                m_fixedGrid, zeroReactions );
        m_hasURR_probabilityTables = m_heatedCrossSections.hasURR_probabilityTables( );
        m_URR_domainMin = m_heatedCrossSections.URR_domainMin( );
        m_URR_domainMax = m_heatedCrossSections.URR_domainMax( ); }
    else {
        m_heatedMultigroupCrossSections.update( a_smr, a_protare, setupInfo, a_settings, particles, a_temperatureInfos, GIDI_reactions, 
                GIDI_orphanProducts, zeroReactions, a_reactionsToExclude );

        if( ( a_settings.upscatterModelAGroupBoundaries().size( ) > 0 ) ||
                ( a_settings.upscatterModel( ) == Sampling::Upscatter::Model::DBRC ) ) { // Load pointwise data to recompute Model A cross sections on user-defined grid or for DBRC.
            removeContinuousEnergyData = true;
            m_heatedCrossSections.update( a_smr, setupInfo, a_settings, particles, a_domainHash, a_temperatureInfos, GIDI_reactions, GIDI_orphanProducts,
                    m_fixedGrid, zeroReactions );
        }
    }

    if( m_upscatterModelASupported && ( a_settings.upscatterModel( ) == Sampling::Upscatter::Model::A ) ) {
        std::vector<double> const &upscatterModelAGroupBoundaries = a_settings.upscatterModelAGroupBoundaries( );
        if( upscatterModelAGroupBoundaries.size( ) == 0 ) {
            GIDI::Styles::Base const *style = a_protare.styles( ).get<GIDI::Styles::Base>( a_temperatureInfos[0].heatedMultiGroup( ) );

            if( style->moniker( ) == GIDI_SnElasticUpScatterStyleChars ) style = a_protare.styles( ).get<GIDI::Styles::Base>( style->derivedStyle( ) );
            if( style->moniker( ) != GIDI_heatedMultiGroupStyleChars ) throw GIDI::Exception( "Label does not yield a heatedMultiGroup style." );

            GIDI::Styles::HeatedMultiGroup const &heatedMultiGroup = *static_cast<GIDI::Styles::HeatedMultiGroup const *>( style );
            std::vector<double> const &boundaries = heatedMultiGroup.groupBoundaries( a_protare.projectile( ).ID( ) );

            m_upscatterModelAGroupEnergies.resize( boundaries.size( ) );
            m_upscatterModelAGroupVelocities.resize( boundaries.size( ) );
            for( std::size_t i1 = 0; i1 < boundaries.size( ); ++i1 ) {
                m_upscatterModelAGroupEnergies[i1] = boundaries[i1];
                m_upscatterModelAGroupVelocities[i1] = MCGIDI_particleBeta( projectileMass( ), boundaries[i1] );
            }

            GIDI::ExcludeReactionsSet reactionsToExclude;
            auto upscatterModelACrossSectionForm = a_protare.multiGroupCrossSection( a_smr, multiGroupSettings, a_temperatureInfos[0], 
                    reactionsToExclude, a_temperatureInfos[0].heatedMultiGroup( ) );
            m_upscatterModelACrossSection.resize( upscatterModelACrossSectionForm.size( ) );
            for( std::size_t i1 = 0; i1 < upscatterModelACrossSectionForm.size( ); ++i1 ) 
                m_upscatterModelACrossSection[i1] = upscatterModelACrossSectionForm[i1]; }
        else {
            
            m_upscatterModelAGroupEnergies.reserve( upscatterModelAGroupBoundaries.size( ) );
            m_upscatterModelAGroupVelocities.reserve( upscatterModelAGroupBoundaries.size( ) );
            for( auto iter = upscatterModelAGroupBoundaries.begin( ); iter != upscatterModelAGroupBoundaries.end( ); ++iter ) {
                m_upscatterModelAGroupEnergies.push_back( *iter );
                m_upscatterModelAGroupVelocities.push_back( MCGIDI_particleBeta( projectileMass( ), *iter ) );
            }

            GIDI::Transporting::MultiGroup boundaries( "Model A", upscatterModelAGroupBoundaries );

            GIDI::Functions::XYs1d crossSectionXYs1d = m_heatedCrossSections.crossSectionAsGIDI_XYs1d( 0.0 );

            GIDI::Transporting::Flux flux( "Model A", 0.0 );
            std::vector<double> energies, fluxes;
            energies.push_back( upscatterModelAGroupBoundaries[0] );
            energies.push_back( upscatterModelAGroupBoundaries.back( ) );
            fluxes.push_back( 1.0 );
            fluxes.push_back( 1.0 );
            flux.addFluxOrder( GIDI::Transporting::Flux_order( 0, energies, fluxes ) );

            GIDI::Vector crossSectionVector = multiGroupXYs1d( boundaries, crossSectionXYs1d, flux );
            m_upscatterModelACrossSection.resize( crossSectionVector.size( ) );
            for( std::size_t index = 0; index < crossSectionVector.size( ); ++index )
                    m_upscatterModelACrossSection[index] = crossSectionVector[index];
        }
        if( !m_continuousEnergy ) m_multiGroupHash = MultiGroupHash( m_projectileMultiGroupBoundariesCollapsed );
    }

    if( ( PoPI::Intids::neutron  == projectileIntid( ) ) && ( a_settings.upscatterModel( ) == Sampling::Upscatter::Model::DBRC ) ) {
        std::size_t reactionIndex = 0;
        for( auto reactionIter = m_reactions.begin( ); reactionIter != m_reactions.end( ); ++reactionIter, ++reactionIndex ) {
            if( (*reactionIter)->ENDF_MT( ) == 2 ) {
                Reaction *reaction = *reactionIter;

                HeatedCrossSectionContinuousEnergy const *heatedCrossSectionContinuousEnergy = m_heatedCrossSections.heatedCrossSections( )[0];
                HeatedReactionCrossSectionContinuousEnergy const *heatedReactionCrossSectionContinuousEnergy 
                        = heatedCrossSectionContinuousEnergy->reactionCrossSection( reactionIndex );

                Vector<double> const &energies = heatedCrossSectionContinuousEnergy->energies( );
                Vector<MCGIDI_FLOAT> const &crossSectionsFloat = heatedReactionCrossSectionContinuousEnergy->crossSections( );
                Vector<double> crossSections( crossSectionsFloat.size( ) );
                std::size_t index = 0;
                for( auto iter = crossSectionsFloat.begin( ); iter != crossSectionsFloat.end( ); ++iter, ++index )
                    crossSections[index] = *iter;

                Sampling::Upscatter::ModelDBRC_data *modelDBRC_data = 
                        new Sampling::Upscatter::ModelDBRC_data( projectileMass( ), targetMass( ), energies, crossSections, a_domainHash );
                reaction->setModelDBRC_data( modelDBRC_data );
                break;
            }
        }
    }
    if( removeContinuousEnergyData ) {
        m_heatedCrossSections.clear( );
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProtareSingle::~ProtareSingle( ) {

    for( Vector<NuclideGammaBranchInfo *>::const_iterator iter = m_branches.begin( ); iter < m_branches.end( ); ++iter ) delete *iter;
    for( Vector<NuclideGammaBranchStateInfo *>::const_iterator iter = m_nuclideGammaBranchStateInfos.begin( ); iter < m_nuclideGammaBranchStateInfos.end( ); ++iter ) delete *iter;
    for( Vector<Reaction *>::const_iterator iter = m_reactions.begin( ); iter < m_reactions.end( ); ++iter ) delete *iter;
    for( Vector<Reaction *>::const_iterator iter = m_orphanProducts.begin( ); iter < m_orphanProducts.end( ); ++iter ) delete *iter;
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs index *a_particleIndex*.
 *
 * @param a_particleIndex       [in]    The PoPs index of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void ProtareSingle::setUserParticleIndex2( int a_particleIndex, int a_userParticleIndex ) {

    m_heatedCrossSections.setUserParticleIndex( a_particleIndex, a_userParticleIndex );
    m_heatedMultigroupCrossSections.setUserParticleIndex( a_particleIndex, a_userParticleIndex );
    for( auto iter = m_reactions.begin( ); iter < m_reactions.end( ); ++iter ) (*iter)->setUserParticleIndex( a_particleIndex, a_userParticleIndex );
    for( auto iter = m_orphanProducts.begin( ); iter < m_orphanProducts.end( ); ++iter ) (*iter)->setUserParticleIndex( a_particleIndex, a_userParticleIndex );
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs intid *a_particleIntid*.
 *  
 * @param a_particleIndex       [in]    The PoPs index of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/
 
LUPI_HOST void ProtareSingle::setUserParticleIndexViaIntid2( int a_particleIntid, int a_userParticleIndex ) {
 
    m_heatedCrossSections.setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
    m_heatedMultigroupCrossSections.setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
    for( auto iter = m_reactions.begin( ); iter < m_reactions.end( ); ++iter ) (*iter)->setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
    for( auto iter = m_orphanProducts.begin( ); iter < m_orphanProducts.end( ); ++iter ) (*iter)->setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
}

/* *********************************************************************************************************//**
 * Returns the pointer representing the protare (i.e., *this*) if *a_index* is 0 and nullptr otherwise.
 *
 * @param a_index               [in]    Must always be 0.
 *
 * @return                              Returns the pointer representing *this*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProtareSingle const *ProtareSingle::protare( std::size_t a_index ) const {

    if( a_index != 0 ) return( nullptr );
    return( this );
}

/* *********************************************************************************************************//**
 * Returns the pointer representing the protare (i.e., *this*) if *a_index* is 0 and nullptr otherwise.
 *
 * @param a_index               [in]    Must always be 0.
 *
 * @return                              Returns the pointer representing *this*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProtareSingle *ProtareSingle::protare( std::size_t a_index ) {

    if( a_index != 0 ) return( nullptr );
    return( this );
}

/* *********************************************************************************************************//**
 * Returns the pointer to the **this** if (*a_index* - 1)th is a value reaction index and nullptr otherwise.
 * 
 * @param a_index               [in]    Index of the reaction.
 * 
 * @return                              Pointer to the requested protare or nullptr if invalid *a_index*..
 ***********************************************************************************************************/
 
LUPI_HOST_DEVICE ProtareSingle const *ProtareSingle::protareWithReaction( std::size_t a_index ) const {
 
    if( static_cast<std::size_t>( a_index ) < numberOfReactions( ) ) return( this );
    return( nullptr );
}

/* *********************************************************************************************************//**
 * Returns the list of temperatures for *this*.
 *
 * @param a_index               [in]    Index of the reqested ProtareSingle. Must be 0.
 *
 * @return                              Vector of doubles.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Vector<double> ProtareSingle::temperatures( std::size_t a_index ) const {

    if( a_index != 0 ) LUPI_THROW( "ProtareSingle::temperatures: a_index not 0." );
    if( m_continuousEnergy ) return( m_heatedCrossSections.temperatures( ) );
    return( m_heatedMultigroupCrossSections.temperatures( ) );
}

/* *********************************************************************************************************//**
 * Sets up the nuclear gamma branching data needed to sample gamma decays.
 *
 * @param a_setupInfo                           [in]    Used internally when constructing a Protare to pass information to other constructors.
 * @param a_protare                             [in]    The **GIDI::Protare** whose data are to be used to construct gamma branching data.
 * @param a_makePhotonEmissionProbabilitiesOne  [in]    If true, all photon emission probabilities are set to 1.0 (i.e., all ICCs are set to 0.0).
 ***********************************************************************************************************/

LUPI_HOST void ProtareSingle::setupNuclideGammaBranchStateInfos( SetupInfo &a_setupInfo, GIDI::ProtareSingle const &a_protare, 
                bool a_makePhotonEmissionProbabilitiesOne, bool a_zeroNuclearLevelEnergyWidth ) {

    PoPI::NuclideGammaBranchStateInfos const &nuclideGammaBranchStateInfos = a_protare.nuclideGammaBranchStateInfos( );
    std::vector<NuclideGammaBranchInfo *> nuclideGammaBranchInfos;

    for( std::size_t i1 = 0; i1 < nuclideGammaBranchStateInfos.size( ); ++i1 ) {
        a_setupInfo.m_stateNamesToIndices[nuclideGammaBranchStateInfos[i1]->state( )] = (int) i1;
        PoPI::NuclideGammaBranchStateInfo const *nuclideGammaBranchStateInfo = nuclideGammaBranchStateInfos.find( nuclideGammaBranchStateInfos[i1]->state( ) );
        a_setupInfo.m_nuclearLevelEnergies[nuclideGammaBranchStateInfos[i1]->state( )] = nuclideGammaBranchStateInfo->nuclearLevelEnergy( );
    }

    m_nuclideGammaBranchStateInfos.reserve( nuclideGammaBranchStateInfos.size( ) );
    for( std::size_t i1 = 0; i1 < nuclideGammaBranchStateInfos.size( ); ++i1 ) {
        m_nuclideGammaBranchStateInfos.push_back( new NuclideGammaBranchStateInfo( *nuclideGammaBranchStateInfos[i1], nuclideGammaBranchInfos,
                a_setupInfo.m_stateNamesToIndices, a_makePhotonEmissionProbabilitiesOne, a_zeroNuclearLevelEnergyWidth ) );
    }

    m_branches.reserve( nuclideGammaBranchInfos.size( ) );
    for( std::size_t i1 = 0; i1 < nuclideGammaBranchInfos.size( ); ++i1 ) m_branches.push_back( nuclideGammaBranchInfos[i1] );
}

/* *********************************************************************************************************//**
 * Returns true if *this* has a fission reaction and false otherwise.
 *
 * @return                              true is if *this* has a fission reaction and false otherwise.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE bool ProtareSingle::hasFission( ) const {

    for( Vector<Reaction *>::const_iterator iter = m_reactions.begin( ); iter < m_reactions.end( ); ++iter ) {
        if( (*iter)->hasFission( ) ) return( true );
    }
    return( false );
}

/* *********************************************************************************************************//**
 * Returns true if *this* has an incoherent photoatomic doppler broadened reaction and false otherwise.
 *
 * @return                              true is if *this* has a specified reaction and false otherwise.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE bool ProtareSingle::hasIncoherentDoppler( ) const {

    for( Vector<Reaction *>::const_iterator iter = m_reactions.begin( ); iter < m_reactions.end( ); ++iter ) {
        if( (*iter)->ENDF_MT( ) == 1534 ) return( true );
    }
    return( false );
}

/* *********************************************************************************************************//**
 * Returns true if *a_energy* with unresolved resonance region (URR) of *this* and false otherwise.
 *
 * @return                              true is if *this* has a URR data.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE bool ProtareSingle::inURR( double a_energy ) const {

    if( a_energy < m_URR_domainMin ) return( false );
    if( a_energy > m_URR_domainMax ) return( false );

    return( true );
}

/* *********************************************************************************************************//**
 * Returns the total cross section for target temperature *a_temperature* and projectile energy *a_energy*. 
 * *a_sampling* is only used for multi-group cross section look up.
 *
 * @param a_URR_protareInfos    [in]    URR information.
 * @param a_hashIndex           [in]    Specifies the continuous energy or multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The energy of the projectile.
 * @param a_sampling            [in]    Used for multi-group look up. If *true*, use augmented cross sections.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareSingle::crossSection( URR_protareInfos const &a_URR_protareInfos, std::size_t a_hashIndex, double a_temperature, double a_energy, bool a_sampling ) const {

    if( m_continuousEnergy ) return( m_heatedCrossSections.crossSection( a_URR_protareInfos, m_URR_index, a_hashIndex, a_temperature, a_energy ) );

    return( m_heatedMultigroupCrossSections.crossSection( a_hashIndex, a_temperature, a_sampling ) );
}

/* *********************************************************************************************************//**
 * Adds the energy dependent, total cross section corresponding to the temperature *a_temperature* multiplied by *a_userFactor* to *a_crossSectionVector*.
 * 
 * @param   a_temperature               [in]        Specifies the temperature of the material.
 * @param   a_userFactor                [in]        User factor which all cross sections are multiplied by.
 * @param   a_numberAllocated           [in]        The length of memory allocated for *a_crossSectionVector*.
 * @param   a_crossSectionVector        [in/out]   The energy dependent, total cross section to add cross section data to.
 ***********************************************************************************************************/
 
LUPI_HOST_DEVICE void ProtareSingle::crossSectionVector( double a_temperature, double a_userFactor, std::size_t a_numberAllocated, 
                double *a_crossSectionVector ) const {

    if( m_continuousEnergy ) {
        if( !m_fixedGrid ) LUPI_THROW( "ProtareSingle::crossSectionVector: continuous energy cannot be supported." );
        m_heatedCrossSections.crossSectionVector( a_temperature, a_userFactor, a_numberAllocated, a_crossSectionVector ); }
    else {
        m_heatedMultigroupCrossSections.crossSectionVector( a_temperature, a_userFactor, a_numberAllocated, a_crossSectionVector );
    }
}

/* *********************************************************************************************************//**
 * Returns the reaction's cross section for the reaction at index *a_reactionIndex*, for target temperature *a_temperature* and projectile energy *a_energy*. 
 * *a_sampling* is only used for multi-group cross section look up.
 *
 * @param a_reactionIndex       [in]    The index of the reaction.
 * @param a_URR_protareInfos    [in]    URR information.
 * @param a_hashIndex           [in]    Specifies the continuous energy or multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The energy of the projectile.
 * @param a_sampling            [in]    Used for multi-group look up. If *true*, use augmented cross sections.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareSingle::reactionCrossSection( std::size_t a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, std::size_t a_hashIndex, 
                double a_temperature, double a_energy, bool a_sampling ) const {

    if( m_continuousEnergy ) return( m_heatedCrossSections.reactionCrossSection( a_reactionIndex, a_URR_protareInfos, m_URR_index, a_hashIndex, 
            a_temperature, a_energy ) );

    return( m_heatedMultigroupCrossSections.reactionCrossSection( a_reactionIndex, a_hashIndex, a_temperature, a_sampling ) );
}

/* *********************************************************************************************************//**
 * Returns the reaction's cross section for the reaction at index *a_reactionIndex*, for target temperature *a_temperature* and projectile energy *a_energy*.
 *
 * @param a_reactionIndex       [in]    The index of the reaction.
 * @param a_URR_protareInfos    [in]    URR information.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The energy of the projectile.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareSingle::reactionCrossSection( std::size_t a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, double a_temperature, double a_energy ) const {

    if( m_continuousEnergy ) return( m_heatedCrossSections.reactionCrossSection( a_reactionIndex, a_URR_protareInfos, m_URR_index, a_temperature, a_energy ) );

    return( m_heatedMultigroupCrossSections.reactionCrossSection( a_reactionIndex, a_temperature, a_energy ) );
}

/* *********************************************************************************************************//**
 * Returns the index of a sampled reaction for a target with termpature *a_temperature*, a projectile with energy *a_energy* and total cross section
 * *a_crossSection*. Random numbers are obtained via *a_userrng* and *a_rngState*.
 *
 * @param a_hashIndex           [in]    Specifies the continuous energy or multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The energy of the projectile.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareSingle::depositionEnergy( std::size_t a_hashIndex, double a_temperature, double a_energy ) const {

    if( m_continuousEnergy ) return( m_heatedCrossSections.depositionEnergy( a_hashIndex, a_temperature, a_energy ) );

    return( m_heatedMultigroupCrossSections.depositionEnergy( a_hashIndex, a_temperature ) );
}

/* *********************************************************************************************************//**
 * Returns the index of a sampled reaction for a target with termpature *a_temperature*, a projectile with energy *a_energy* and total cross section
 * *a_crossSection*. Random numbers are obtained via *a_userrng* and *a_rngState*.
 *
 * @param a_hashIndex           [in]    Specifies the continuous energy or multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The energy of the projectile.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareSingle::depositionMomentum( std::size_t a_hashIndex, double a_temperature, double a_energy ) const {

    if( m_continuousEnergy ) return( m_heatedCrossSections.depositionMomentum( a_hashIndex, a_temperature, a_energy ) );

    return( m_heatedMultigroupCrossSections.depositionMomentum( a_hashIndex, a_temperature ) );
}

/* *********************************************************************************************************//**
 * Returns the index of a sampled reaction for a target with termpature *a_temperature*, a projectile with energy *a_energy* and total cross section
 * *a_crossSection*. Random numbers are obtained via *a_userrng* and *a_rngState*.
 *
 * @param a_hashIndex           [in]    Specifies the continuous energy or multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The energy of the projectile.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareSingle::productionEnergy( std::size_t a_hashIndex, double a_temperature, double a_energy ) const {

    if( m_continuousEnergy ) return( m_heatedCrossSections.productionEnergy( a_hashIndex, a_temperature, a_energy ) );

    return( m_heatedMultigroupCrossSections.productionEnergy( a_hashIndex, a_temperature ) );
}

/* *********************************************************************************************************//**
 * Returns the index of a sampled reaction for a target with termpature *a_temperature*, a projectile with energy *a_energy* and total cross section
 * *a_crossSection*. Random numbers are obtained via *a_userrng* and *a_rngState*.
 *
 * @param a_hashIndex           [in]    Specifies the continuous energy or multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The energy of the projectile.
 * @param a_particleIndex       [in]    The index of the particle whose gain is to be returned.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareSingle::gain( std::size_t a_hashIndex, double a_temperature, double a_energy, int a_particleIndex ) const {

    if( m_continuousEnergy ) return( m_heatedCrossSections.gain( a_hashIndex, a_temperature, a_energy, a_particleIndex ) );

    return( m_heatedMultigroupCrossSections.gain( a_hashIndex, a_temperature, a_particleIndex ) );
}

/* *********************************************************************************************************//**
 * Returns the intid of a sampled reaction for a target with termpature *a_temperature*, a projectile with energy *a_energy* and total cross section
 * *a_crossSection*. Random numbers are obtained via *a_userrng* and *a_rngState*.
 *
 * @param a_hashIndex           [in]    Specifies the continuous energy or multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The energy of the projectile.
 * @param a_particleIntid       [in]    The intid of the particle whose gain is to be returned.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProtareSingle::gainViaIntid( std::size_t a_hashIndex, double a_temperature, double a_energy, int a_particleIntid ) const {

    if( m_continuousEnergy ) return( m_heatedCrossSections.gainViaIntid( a_hashIndex, a_temperature, a_energy, a_particleIntid ) );

    return( m_heatedMultigroupCrossSections.gainViaIntid( a_hashIndex, a_temperature, a_particleIntid ) );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void ProtareSingle::serialize2( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    std::size_t vectorSize;
    LUPI::DataBuffer *workingBuffer = &a_buffer;

    DATA_MEMBER_STRING( m_interaction, a_buffer, a_mode );
    DATA_MEMBER_INT( m_URR_index, a_buffer, a_mode );
    DATA_MEMBER_CAST( m_hasURR_probabilityTables, a_buffer, a_mode, bool );
    DATA_MEMBER_DOUBLE( m_URR_domainMin, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_URR_domainMax, a_buffer, a_mode );
    m_domainHash.serialize( a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_projectileMultiGroupBoundaries, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_projectileMultiGroupBoundariesCollapsed, a_buffer, a_mode );
    DATA_MEMBER_CAST( m_upscatterModelASupported, a_buffer, a_mode, bool );
    DATA_MEMBER_VECTOR_DOUBLE( m_upscatterModelAGroupEnergies, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_upscatterModelAGroupVelocities, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_upscatterModelACrossSection, a_buffer, a_mode );
    m_multiGroupHash.serialize( a_buffer, a_mode );

    vectorSize = m_nuclideGammaBranchStateInfos.size( );
    int vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, *workingBuffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        m_nuclideGammaBranchStateInfos.resize( vectorSize, &(workingBuffer->m_placement) );
        for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
            if (workingBuffer->m_placement != nullptr) {
                m_nuclideGammaBranchStateInfos[vectorIndex] = new(workingBuffer->m_placement) NuclideGammaBranchStateInfo;
                workingBuffer->incrementPlacement( sizeof( NuclideGammaBranchStateInfo ) );
            }
            else {
                m_nuclideGammaBranchStateInfos[vectorIndex] = new NuclideGammaBranchStateInfo;
            }
        }
    }
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        a_buffer.m_placement += m_nuclideGammaBranchStateInfos.internalSize();
        a_buffer.incrementPlacement( sizeof( NuclideGammaBranchStateInfo ) * vectorSize );
    }
    for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
        m_nuclideGammaBranchStateInfos[vectorIndex]->serialize( *workingBuffer, a_mode );
    }

    vectorSize = m_branches.size( );
    vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, *workingBuffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        m_branches.resize( vectorSize, &(workingBuffer->m_placement) );
        for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
            if (workingBuffer->m_placement != nullptr) {
                m_branches[vectorIndex] = new(workingBuffer->m_placement) NuclideGammaBranchInfo;
                workingBuffer->incrementPlacement( sizeof( NuclideGammaBranchInfo ) );
            }
            else {
                m_branches[vectorIndex] = new NuclideGammaBranchInfo;
            }
        }
    }
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        a_buffer.m_placement += m_branches.internalSize();
        workingBuffer->incrementPlacement( sizeof( NuclideGammaBranchInfo ) * vectorSize );
    }
    for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
        m_branches[vectorIndex]->serialize( *workingBuffer, a_mode );
    }

    vectorSize = m_reactions.size( );
    vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, *workingBuffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        m_reactions.resize( vectorSize, &(workingBuffer->m_placement) );
        for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
            if (workingBuffer->m_placement != nullptr) {
                m_reactions[vectorIndex] = new(workingBuffer->m_placement) Reaction;
                workingBuffer->incrementPlacement( sizeof(Reaction));
            }
            else {
                m_reactions[vectorIndex] = new Reaction;
            }
        }
    }
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        a_buffer.m_placement += m_reactions.internalSize();
        a_buffer.incrementPlacement( sizeof(Reaction) * vectorSize);
    }
    for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
        m_reactions[vectorIndex]->serialize( *workingBuffer, a_mode );
        m_reactions[vectorIndex]->updateProtareSingleInfo( this, vectorIndex );
    }

    vectorSize = m_orphanProducts.size( );
    vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, *workingBuffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        m_orphanProducts.resize( vectorSize, &(workingBuffer->m_placement) );
        for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
            if (workingBuffer->m_placement != nullptr) {
                m_orphanProducts[vectorIndex] = new(workingBuffer->m_placement) Reaction;
                workingBuffer->incrementPlacement( sizeof(Reaction));
            }
            else {
                m_orphanProducts[vectorIndex] = new Reaction;
            }
        }
    }

    if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        a_buffer.m_placement += m_orphanProducts.internalSize( );
        a_buffer.incrementPlacement( sizeof( Reaction ) * vectorSize );
    }

    for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
        m_orphanProducts[vectorIndex]->serialize( *workingBuffer, a_mode );
        m_orphanProducts[vectorIndex]->updateProtareSingleInfo( this, vectorIndex );
    }

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        for( auto reactionIter = m_reactions.begin( ); reactionIter != m_reactions.end( ); ++reactionIter ) {
            (*reactionIter)->addOrphanProductToProductList( m_orphanProducts );
        }
    }

    DATA_MEMBER_CAST( m_isPhotoAtomic, *workingBuffer, a_mode, bool );
    DATA_MEMBER_CAST( m_continuousEnergy, *workingBuffer, a_mode, bool );
    DATA_MEMBER_CAST( m_fixedGrid, *workingBuffer, a_mode, bool );
    m_heatedCrossSections.serialize( *workingBuffer, a_mode );
    m_heatedMultigroupCrossSections.serialize( *workingBuffer, a_mode );
}

}
