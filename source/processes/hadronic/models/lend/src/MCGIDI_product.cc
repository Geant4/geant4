/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "MCGIDI.hpp"

namespace MCGIDI {

/*! \class Product
 * This class represents a **GNDS** <**product**> node with only data needed for Monte Carlo transport.
 */

/* *********************************************************************************************************//**
 * Default constructor used when broadcasting a Protare as needed by MPI or GPUs.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Product::Product( ) :
        m_ID( ),
        m_intid( -1 ),
        m_index( -1 ),
        m_userParticleIndex( -1 ),
        m_mass( 0.0 ),
        m_excitationEnergy( 0.0 ),
        m_twoBodyOrder( TwoBodyOrder::notApplicable ),
        m_initialStateIndex( -1 ),
        m_multiplicity( nullptr ),
        m_distribution( nullptr ),
        m_outputChannel( nullptr ) {

}

/* *********************************************************************************************************//**
 * @param a_product             [in]    The GIDI::Product whose data is to be used to construct *this*.
 * @param a_setupInfo           [in]    Used internally when constructing a Protare to pass information to other constructors.
 * @param a_settings            [in]    Used to pass user options to the *this* to instruct it which data are desired.
 * @param a_particles           [in]    List of transporting particles and their information (e.g., multi-group boundaries and fluxes).
 * @param a_isFission           [in]    *true* if parent channel is a fission channel and *false* otherwise.
 ***********************************************************************************************************/

LUPI_HOST Product::Product( GIDI::Product const *a_product, SetupInfo &a_setupInfo, Transporting::MC const &a_settings, 
                GIDI::Transporting::Particles const &a_particles, bool a_isFission ) :
        m_ID( a_product->particle( ).ID( ).c_str( ) ),
        m_intid( MCGIDI_popsIntid( a_setupInfo.m_pops, a_product->particle( ).ID( ) ) ),
        m_index( MCGIDI_popsIndex( a_setupInfo.m_popsUser, a_product->particle( ).ID( ) ) ),
        m_userParticleIndex( -1 ),
        m_label( a_product->label( ).c_str( ) ),
        m_isCompleteParticle( a_product->isCompleteParticle( ) ),
        m_mass( a_product->particle( ).mass( "MeV/c**2" ) ),         // Includes nuclear excitation energy.
        m_excitationEnergy( a_product->particle( ).excitationEnergy( ).value( ) ),
        m_twoBodyOrder( a_setupInfo.m_twoBodyOrder ),
        m_initialStateIndex( -1 ),
        m_multiplicity( Functions::parseMultiplicityFunction1d( a_setupInfo, a_settings, a_product->multiplicity( ) ) ),
        m_distribution( nullptr ),
        m_outputChannel( nullptr ) {

    a_setupInfo.m_product1Mass = mass( );                           // Includes nuclear excitation energy.
    a_setupInfo.m_initialStateIndex = -1;
    m_distribution = Distributions::parseGIDI( a_product->distribution( ), a_setupInfo, a_settings );
    m_initialStateIndex = a_setupInfo.m_initialStateIndex;

    GIDI::OutputChannel const *output_channel = a_product->outputChannel( );
    if( output_channel != nullptr ) m_outputChannel = new OutputChannel( output_channel, a_setupInfo, a_settings, a_particles );

    if( a_isFission && ( m_intid == PoPI::Intids::neutron ) && a_settings.wantTerrellPromptNeutronDistribution( ) ) {
        Functions::Function1d_d1 *multiplicity1 = static_cast<Functions::Function1d_d1 *>( m_multiplicity );

        m_multiplicity = new Functions::TerrellFissionNeutronMultiplicityModel( -1.0, multiplicity1 );
    }
}

/* *********************************************************************************************************//**
 * @param a_pops                [in]    A PoPs Database instance used to get particle intids and possibly other particle information.
 * @param a_ID                  [in]    The PoPs id for the product.
 * @param a_label               [in]    The **GNDS** label for the product.
 ***********************************************************************************************************/

LUPI_HOST Product::Product( PoPI::Database const &a_pops, std::string const &a_ID, std::string const &a_label ) :
        m_ID( a_ID.c_str( ) ),
        m_intid( MCGIDI_popsIntid( a_pops, a_ID ) ),
        m_index( MCGIDI_popsIndex( a_pops, a_ID ) ),
        m_userParticleIndex( -1 ),
        m_label( a_label.c_str( ) ),
        m_mass( 0.0 ),                                  // FIXME, good for photon but nothing else. Still need to implement.
        m_excitationEnergy( 0.0 ),
        m_twoBodyOrder( TwoBodyOrder::notApplicable ),
        m_initialStateIndex( -1 ),
        m_multiplicity( nullptr ),
        m_distribution( nullptr ),
        m_outputChannel( nullptr ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Product::~Product( ) {

    delete m_multiplicity;

    Distributions::Type type = Distributions::Type::none;
    if( m_distribution != nullptr ) type = m_distribution->type( );
    switch( type ) {
    case Distributions::Type::none:
        break;
    case Distributions::Type::unspecified:
        delete static_cast<Distributions::Unspecified *>( m_distribution );
        break;
    case Distributions::Type::angularTwoBody:
        delete static_cast<Distributions::AngularTwoBody *>( m_distribution );
        break;
    case Distributions::Type::KalbachMann:
        delete static_cast<Distributions::KalbachMann *>( m_distribution );
        break;
    case Distributions::Type::uncorrelated:
        delete static_cast<Distributions::Uncorrelated *>( m_distribution );
        break;
    case Distributions::Type::branching3d:
        delete static_cast<Distributions::Branching3d *>( m_distribution );
        break;
    case Distributions::Type::energyAngularMC:
        delete static_cast<Distributions::EnergyAngularMC *>( m_distribution );
        break;
    case Distributions::Type::angularEnergyMC:
        delete static_cast<Distributions::AngularEnergyMC *>( m_distribution );
        break;
    case Distributions::Type::coherentPhotoAtomicScattering:
        delete static_cast<Distributions::CoherentPhotoAtomicScattering *>( m_distribution );
        break;
    case Distributions::Type::incoherentPhotoAtomicScattering:
        delete static_cast<Distributions::IncoherentPhotoAtomicScattering *>( m_distribution );
        break;
    case Distributions::Type::incoherentBoundToFreePhotoAtomicScattering:
        delete static_cast<Distributions::IncoherentBoundToFreePhotoAtomicScattering *>( m_distribution );
        break;
    case Distributions::Type::incoherentPhotoAtomicScatteringElectron:
        delete static_cast<Distributions::IncoherentPhotoAtomicScatteringElectron *>( m_distribution );
        break;
    case Distributions::Type::pairProductionGamma:
        delete static_cast<Distributions::PairProductionGamma *>( m_distribution );
        break;
    case Distributions::Type::coherentElasticTNSL:
        delete static_cast<Distributions::CoherentElasticTNSL *>( m_distribution );
        break;
    case Distributions::Type::incoherentElasticTNSL:
        delete static_cast<Distributions::IncoherentElasticTNSL *>( m_distribution );
        break;
    }

    delete m_outputChannel;
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs index *a_particleIndex*.
 *  
 * @param a_particleIndex       [in]    The PoPs index of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void Product::setUserParticleIndex( int a_particleIndex, int a_userParticleIndex ) {

    if( m_index == a_particleIndex ) m_userParticleIndex = a_userParticleIndex;
#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    if( m_outputChannel != nullptr ) m_outputChannel->setUserParticleIndex( a_particleIndex, a_userParticleIndex );
#endif
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs intid *a_particleIntid*.
 *
 * @param a_particleIntid       [in]    The PoPs intid of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void Product::setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex ) {

    if( m_intid == a_particleIntid ) m_userParticleIndex = a_userParticleIndex;
#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    if( m_outputChannel != nullptr ) m_outputChannel->setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
#endif
}

/* *********************************************************************************************************//**
 * This method calls the **setModelDBRC_data* method on the distribution of *this* with *a_modelDBRC_data*.
 *
 * @param a_modelDBRC_data      [in]    The instance storing data needed to treat the DRRC upscatter mode.
 ***********************************************************************************************************/

LUPI_HOST void Product::setModelDBRC_data( Sampling::Upscatter::ModelDBRC_data *a_modelDBRC_data ) {

    m_distribution->setModelDBRC_data( a_modelDBRC_data );
}

/* *********************************************************************************************************//**
 * This method returns the final Q for *this* by getting its output channel's finalQ.
 *
 * @param a_x1                  [in]    The energy of the projectile.
 *
 * @return                              The Q-value at product energy *a_x1*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double Product::finalQ( LUPI_maybeUnused double a_x1 ) const {

#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    if( m_outputChannel != nullptr ) return( m_outputChannel->finalQ( a_x1 ) );
#endif
    return( m_excitationEnergy );
}

/* *********************************************************************************************************//**
 * This method returns *true* if the output channel or any of its sub-output channels is a fission channel and *false* otherwise.
 *
 * @return                              *true* if any sub-output channel is a fission channel and *false* otherwise.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE bool Product::hasFission( ) const {

#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    if( m_outputChannel != nullptr ) return( m_outputChannel->hasFission( ) );
#endif
    return( false );
}

/* *********************************************************************************************************//**
 * Returns the energy dependent multiplicity for outgoing particle with pops index *a_index*. The returned value may not
 * be an integer. Energy dependent multiplicities mainly occurs for photons and fission neutrons.
 *
 * @param a_index                   [in]    The PoPs index of the requested particle.
 * @param a_projectileEnergy        [in]    The energy of the projectile.
 *
 * @return                                  The multiplicity value for the requested particle.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double Product::productAverageMultiplicity( int a_index, double a_projectileEnergy ) const {

    double multiplicity1 = 0.0;

    if( a_index == m_index ) {
        if( ( m_multiplicity->domainMin( ) <= a_projectileEnergy ) && ( m_multiplicity->domainMax( ) >= a_projectileEnergy ) )
            multiplicity1 += m_multiplicity->evaluate( a_projectileEnergy );
    }
#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    if( m_outputChannel != nullptr ) multiplicity1 += m_outputChannel->productAverageMultiplicity( a_index, a_projectileEnergy );
#endif

    return( multiplicity1 );
}

/* *********************************************************************************************************//**
 * Returns the energy dependent multiplicity for outgoing particle with pops intid *a_intid*. The returned value may not
 * be an integer. Energy dependent multiplicities mainly occurs for photons and fission neutrons.
 *
 * @param a_intid                   [in]    The PoPs intid of the requested particle.
 * @param a_projectileEnergy        [in]    The energy of the projectile.
 *
 * @return                                  The multiplicity value for the requested particle.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double Product::productAverageMultiplicityViaIntid( int a_intid, double a_projectileEnergy ) const {

    double multiplicity1 = 0.0;

    if( a_intid == m_intid ) {
        if( ( m_multiplicity->domainMin( ) <= a_projectileEnergy ) && ( m_multiplicity->domainMax( ) >= a_projectileEnergy ) )
            multiplicity1 += m_multiplicity->evaluate( a_projectileEnergy );
    }
#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    if( m_outputChannel != nullptr ) multiplicity1 += m_outputChannel->productAverageMultiplicityViaIntid( a_intid, a_projectileEnergy );
#endif

    return( multiplicity1 );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Product::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_STRING( m_ID, a_buffer, a_mode );
    DATA_MEMBER_INT( m_intid, a_buffer, a_mode );
    DATA_MEMBER_INT( m_index, a_buffer, a_mode );
    DATA_MEMBER_INT( m_userParticleIndex, a_buffer, a_mode );
    DATA_MEMBER_STRING( m_label, a_buffer, a_mode );
    DATA_MEMBER_CAST( m_isCompleteParticle, a_buffer, a_mode, bool );
    DATA_MEMBER_DOUBLE( m_mass, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_excitationEnergy, a_buffer, a_mode );

    int twoBodyOrder = 0;
    switch( m_twoBodyOrder ) {
    case TwoBodyOrder::notApplicable :
        break;
    case TwoBodyOrder::firstParticle :
        twoBodyOrder = 1;
        break;
    case TwoBodyOrder::secondParticle :
        twoBodyOrder = 2;
        break;
    }
    DATA_MEMBER_INT( twoBodyOrder , a_buffer, a_mode );
    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        switch( twoBodyOrder ) {
        case 0 :
            m_twoBodyOrder = TwoBodyOrder::notApplicable;
            break;
        case 1 :
            m_twoBodyOrder = TwoBodyOrder::firstParticle;
            break;
        case 2 :
            m_twoBodyOrder = TwoBodyOrder::secondParticle;
            break;
        }
    }

    DATA_MEMBER_INT( m_initialStateIndex, a_buffer, a_mode );

    m_multiplicity = serializeFunction1d( a_buffer, a_mode, m_multiplicity );
    m_distribution = serializeDistribution( a_buffer, a_mode, m_distribution );

#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    bool haveChannel = m_outputChannel != nullptr;
    DATA_MEMBER_CAST( haveChannel, a_buffer, a_mode, bool );
    if( haveChannel ) {
        if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
            if (a_buffer.m_placement != nullptr) {
                m_outputChannel = new(a_buffer.m_placement) OutputChannel();
                a_buffer.incrementPlacement( sizeof(OutputChannel));
            }
            else {
                m_outputChannel = new OutputChannel();
            }
        }
        if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
            a_buffer.incrementPlacement( sizeof(OutputChannel));
        }
        m_outputChannel->serialize( a_buffer, a_mode );
    }
#endif
}

}
