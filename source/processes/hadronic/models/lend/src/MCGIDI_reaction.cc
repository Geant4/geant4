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

/*! \class Reaction 
 * Class representing a **GNDS** <**reaction**> node with only data needed for Monte Carlo transport.
 */

/* *********************************************************************************************************//**
 * Default constructor used when broadcasting a Reaction as needed by MPI or GPUs.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Reaction::Reaction( ) :
        m_protareSingle( nullptr ),
        m_reactionIndex( MCGIDI_nullReaction ),
        m_GIDI_reactionIndex( MCGIDI_nullReaction ),
        m_label( ),
        m_ENDF_MT( 0 ),
        m_ENDL_C( 0 ),
        m_ENDL_S( 0 ),
        m_initialStateIndex( -1 ),
        m_neutronIndex( -1 ),
        m_hasFission( false ),
        m_projectileMass( 0.0 ),
        m_targetMass( 0.0 ),
        m_crossSectionThreshold( 0.0 ),
        m_twoBodyThreshold( 0.0 ),
        m_hasFinalStatePhotons( false ),
        m_fissionResiduaIntid( -1 ),
        m_fissionResiduaIndex( -1 ),
        m_fissionResiduaUserIndex( -1 ),
        m_fissionResiduals( GIDI::Construction::FissionResiduals::none ),
        m_fissionResidualMass( 0.0 ),
        m_totalDelayedNeutronMultiplicity( nullptr ),
#ifdef MCGIDI_USE_OUTPUT_CHANNEL
        m_outputChannel( nullptr ),
#endif

// GRIN extras:
        m_GRIN_specialSampleProducts( false ),
        m_GRIN_inelasticThreshold( 0.0 ),
        m_GRIN_maximumCaptureIncidentEnergy( 0.0 ),
        m_GRIN_inelastic( nullptr ),
        m_GRIN_capture( nullptr ) {
}

/* *********************************************************************************************************//**
 * @param a_reaction            [in]    The GIDI::Reaction whose data is to be used to construct *this*.
 * @param a_setupInfo           [in]    Used internally when constructing a Protare to pass information to other constructors.
 * @param a_settings            [in]    Used to pass user options to the *this* to instruct it which data are desired.
 * @param a_particles           [in]    List of transporting particles and their information (e.g., multi-group boundaries and fluxes).
 * @param a_temperatureInfos    [in]    The list of temperature data to extract from *a_protare*.
 ***********************************************************************************************************/

LUPI_HOST Reaction::Reaction( GIDI::Reaction const &a_reaction, SetupInfo &a_setupInfo, Transporting::MC const &a_settings, 
                GIDI::Transporting::Particles const &a_particles, LUPI_maybeUnused GIDI::Styles::TemperatureInfos const &a_temperatureInfos ) :
        m_protareSingle( nullptr ),
        m_reactionIndex( MCGIDI_nullReaction ),
        m_GIDI_reactionIndex( a_reaction.reactionIndex( ) ),
        m_label( a_reaction.label( ).c_str( ) ),
        m_ENDF_MT( a_reaction.ENDF_MT( ) ),
        m_ENDL_C( a_reaction.ENDL_C( ) ),
        m_ENDL_S( a_reaction.ENDL_S( ) ),
        m_initialStateIndex( -1 ),
        m_neutronIndex( a_setupInfo.m_neutronIndex ),
        m_hasFission( a_reaction.hasFission( ) ),
        m_projectileMass( a_setupInfo.m_protare.projectileMass( ) ),
        m_targetMass( a_setupInfo.m_protare.targetMass( ) ),
        m_crossSectionThreshold( a_reaction.crossSectionThreshold( ) ),
        m_twoBodyThreshold( a_reaction.twoBodyThreshold( ) ),
        m_fissionResiduaIntid( -1 ),
        m_fissionResiduaIndex( -1 ),
        m_fissionResiduaUserIndex( -1 ),
        m_fissionResiduals( GIDI::Construction::FissionResiduals::none ),
        m_fissionResidualMass( 0.0 ),

// GRIN extras:
        m_GRIN_specialSampleProducts( false ),
        m_GRIN_inelasticThreshold( 0.0 ),
        m_GRIN_maximumCaptureIncidentEnergy( 0.0 ),
        m_GRIN_inelastic( nullptr ),
        m_GRIN_capture( nullptr ) {

    a_setupInfo.m_hasFinalStatePhotons = false;
#ifndef MCGIDI_USE_OUTPUT_CHANNEL
    OutputChannel *m_outputChannel;
#endif
    m_outputChannel = new OutputChannel( a_reaction.outputChannel( ), a_setupInfo, a_settings, a_particles );

    std::set<std::string> product_ids;

    a_reaction.productIDs( product_ids, a_particles, false );
    m_productIntids.reserve( product_ids.size( ) );
    m_productIndices.reserve( product_ids.size( ) );
    m_userProductIndices.reserve( product_ids.size( ) );
    m_productMultiplicities.reserve( product_ids.size( ) );
    for( std::set<std::string>::iterator iter = product_ids.begin( ); iter != product_ids.end( ); ++iter ) {
        m_productIntids.push_back( MCGIDI_popsIntid( a_setupInfo.m_pops, *iter ) );
        m_productIndices.push_back( static_cast<int>( a_setupInfo.m_popsUser[*iter] ) );
        m_userProductIndices.push_back( -1 );
        m_productMultiplicities.push_back( a_reaction.productMultiplicity( *iter ) );
    }

    product_ids.clear( );
    a_reaction.productIDs( product_ids, a_particles, true );
    m_productIntidsTransportable.reserve( product_ids.size( ) );
    m_productIndicesTransportable.reserve( product_ids.size( ) );
    m_userProductIndicesTransportable.reserve( product_ids.size( ) );
    for( std::set<std::string>::iterator iter = product_ids.begin( ); iter != product_ids.end( ); ++iter ) {
        m_productIntidsTransportable.push_back( MCGIDI_popsIntid( a_setupInfo.m_pops, *iter ) );
        m_productIndicesTransportable.push_back( static_cast<int>( a_setupInfo.m_popsUser[*iter] ) );
        m_userProductIndicesTransportable.push_back( -1 );
    }

    m_hasFinalStatePhotons = a_setupInfo.m_hasFinalStatePhotons;
    m_fissionResiduals = a_reaction.outputChannel( )->fissionResiduals( );
    if(      m_fissionResiduals == GIDI::Construction::FissionResiduals::ENDL99120 ) {
        m_fissionResiduaIntid = PoPI::Intids::FissionProductENDL99120;
        m_fissionResiduaIndex = MCGIDI_popsIndex( a_setupInfo.m_popsUser, PoPI::IDs::FissionProductENDL99120 ); }
    else if( m_fissionResiduals == GIDI::Construction::FissionResiduals::ENDL99125 ) {
        m_fissionResiduaIntid = PoPI::Intids::FissionProductENDL99125;
        m_fissionResiduaIndex = MCGIDI_popsIndex( a_setupInfo.m_popsUser, PoPI::IDs::FissionProductENDL99125 );
    }
    m_fissionResidualMass = 117.5 * PoPI_AMU2MeV_c2;        // Hardwired for now as MeV. Only used if m_fissionResiduaIntid != -1.

#ifndef MCGIDI_USE_OUTPUT_CHANNEL
    std::vector<Product *> products;
    std::vector<DelayedNeutron *> delayedNeutrons;
    std::vector<Functions::Function1d_d1 *> Qs;

    m_totalDelayedNeutronMultiplicity = nullptr;
    m_outputChannel->moveProductsEtAlToReaction( products, &m_totalDelayedNeutronMultiplicity, delayedNeutrons, Qs );

    m_products.resize( products.size( ) );
    for( std::size_t index = 0; index < products.size( ); ++index ) m_products[index] = products[index];

    m_delayedNeutrons.resize( delayedNeutrons.size( ) );
    for( std::size_t index = 0; index < delayedNeutrons.size( ); ++index ) m_delayedNeutrons[index] = delayedNeutrons[index];

    m_Qs.resize( Qs.size( ) );
    for( std::size_t index = 0; index < Qs.size( ); ++index ) m_Qs[index] = Qs[index];

    delete m_outputChannel;
#endif

    GIDI::GRIN::GRIN_continuumGammas const *GRIN_continuumGammas = a_setupInfo.m_GRIN_continuumGammas;
    if( GRIN_continuumGammas != nullptr ) {
        if( m_ENDF_MT == 102 ) {
            if( GRIN_continuumGammas->captureLevelProbabilities( ).size( ) > 0 ) {
                m_GRIN_specialSampleProducts = true;
                m_GRIN_maximumCaptureIncidentEnergy = GRIN_continuumGammas->maximumCaptureIncidentEnergy( ).value( );
                m_GRIN_capture = new GRIN_capture( a_setupInfo, *GRIN_continuumGammas );
            } }
        else if( m_ENDF_MT == 91 ) {
            GIDI::Suite const &inelasticIncidentEnergies = GRIN_continuumGammas->inelasticIncidentEnergies( );
            if( inelasticIncidentEnergies.size( ) > 0 ) {
                m_GRIN_specialSampleProducts = true;
                GIDI::GRIN::InelasticIncidentEnergy const *inelasticIncidentEnergy = inelasticIncidentEnergies.get<GIDI::GRIN::InelasticIncidentEnergy const>( 0 );
                m_GRIN_inelasticThreshold = inelasticIncidentEnergy->energy( );
                m_GRIN_inelastic = new GRIN_inelastic( a_setupInfo, *GRIN_continuumGammas );
            }
        }
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Reaction::~Reaction( ) {

#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    delete m_outputChannel;
#else
    delete m_totalDelayedNeutronMultiplicity;
    for( auto iter = m_products.begin( ); iter != m_products.end( ); ++iter ) delete *iter;
    for( auto iter = m_delayedNeutrons.begin( ); iter != m_delayedNeutrons.end( ); ++iter ) delete *iter;
    for( auto iter = m_Qs.begin( ); iter != m_Qs.end( ); ++iter ) delete *iter;
#endif
}
/* *********************************************************************************************************//**
 * Returns the Q-value for projectile energy *a_energy*. 
 *
 * @param a_URR_protareInfos    [in]    URR information.
 * @param a_hashIndex           [in]    Specifies the continuous energy or multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy_in           [in]    The energy of the projectile.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double Reaction::finalQ( double a_energy ) const { 

#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    return( m_outputChannel->finalQ( a_energy ) );
#else
    double Q = 0.0;
    for( auto Q_iter = m_Qs.begin( ); Q_iter != m_Qs.end( ); ++Q_iter ) Q += (*Q_iter)->evaluate( a_energy );

    return( Q );
#endif
}

/* *********************************************************************************************************//**
 * Returns the reaction's cross section for target temperature *a_temperature* and projectile energy *a_energy_in*.
 *
 * @param a_URR_protareInfos    [in]    URR information.
 * @param a_hashIndex           [in]    Specifies the continuous energy or multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy_in           [in]    The energy of the projectile.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double Reaction::crossSection( URR_protareInfos const &a_URR_protareInfos, std::size_t a_hashIndex, double a_temperature, double a_energy_in ) const {

    return( m_protareSingle->reactionCrossSection( m_reactionIndex, a_URR_protareInfos, a_hashIndex, a_temperature, a_energy_in, false ) );
}

/* *********************************************************************************************************//**
 * Returns the reaction's cross section for target temperature *a_temperature* and projectile energy *a_energy_in*.
 *
 * @param a_URR_protareInfos    [in]    URR information.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy_in           [in]    The energy of the projectile.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double Reaction::crossSection( URR_protareInfos const &a_URR_protareInfos, double a_temperature, double a_energy_in ) const {

    return( m_protareSingle->reactionCrossSection( m_reactionIndex, a_URR_protareInfos, a_temperature, a_energy_in ) );
}

/* *********************************************************************************************************//**
 * Returns the reaction's cross section as a pointer to a GIDI::Functions::XYs1d instance.
 *
 * @returns                 A GIDI::Functions::XYs1d instance.
 ***********************************************************************************************************/

LUPI_HOST GIDI::Functions::XYs1d Reaction::crossSectionAsGIDI_XYs1d( double a_temperature ) const {

    return( m_protareSingle->heatedCrossSections( ).reactionCrossSectionAsGIDI_XYs1d( m_reactionIndex, a_temperature ) );
}

/* *********************************************************************************************************//**
 * Returns the multiplicity for outgoing particle with pops index *a_index*. If the multiplicity is energy dependent,
 * the returned value is -1. For energy dependent multiplicities it is better to use the method **productAverageMultiplicity**.
 *
 * @param a_index                   [in]    The PoPs index of the requested particle.
 *
 * @return                                  The multiplicity value for the requested particle.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE int Reaction::productMultiplicity( int a_index ) const {

    std::size_t i1 = 0;

    for( Vector<int>::iterator iter = m_productIndices.begin( ); iter != m_productIndices.end( ); ++iter, ++i1 ) {
        if( *iter == a_index ) return( m_productMultiplicities[i1] );
    }

    return( 0 );
}
/* *********************************************************************************************************//**
 * Returns the multiplicity for outgoing particle with pops intid *a_intid*. If the multiplicity is energy dependent,
 * the returned value is -1. For energy dependent multiplicities it is better to use the method **productAverageMultiplicity**.
 *
 * @param a_intid                   [in]    The PoPs intid of the requested particle.
 *
 * @return                                  The multiplicity value for the requested particle.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE int Reaction::productMultiplicityViaIntid( int a_intid ) const {

    std::size_t i1 = 0;

    for( Vector<int>::iterator iter = m_productIntids.begin( ); iter != m_productIntids.end( ); ++iter, ++i1 ) {
        if( *iter == a_intid ) return( m_productMultiplicities[i1] );
    }

    return( 0 );
}

/* *********************************************************************************************************//**
 * Returns the energy dependent multiplicity for outgoing particle with pops index *a_index*. The returned value may not
 * be an integer. Energy dependent multiplicity mainly occurs for photons and fission neutrons.
 *
 * @param a_index                   [in]    The PoPs index of the requested particle.
 * @param a_projectileEnergy        [in]    The energy of the projectile.
 *
 * @return                                  The multiplicity value for the requested particle.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double Reaction::productAverageMultiplicity( int a_index, double a_projectileEnergy ) const {

    double multiplicity = 0.0;

    if( m_crossSectionThreshold > a_projectileEnergy ) return( multiplicity );

    std::size_t i1 = 0;
    for( Vector<int>::iterator iter = m_productIndices.begin( ); iter != m_productIndices.end( ); ++iter, ++i1 ) {
        if( *iter == a_index ) {
            multiplicity = m_productMultiplicities[i1];
            break;
        }
    }

    if( multiplicity < 0 ) {
#ifdef MCGIDI_USE_OUTPUT_CHANNEL
        multiplicity = m_outputChannel->productAverageMultiplicity( a_index, a_projectileEnergy );
#else
        multiplicity = 0.0;
        for( auto productIter = m_products.begin( ); productIter != m_products.end( ); ++productIter ) {
            multiplicity += (*productIter)->productAverageMultiplicity( a_index, a_projectileEnergy );
        }

        if( ( m_totalDelayedNeutronMultiplicity != nullptr ) && ( a_index == m_neutronIndex ) ) {
            multiplicity += m_totalDelayedNeutronMultiplicity->evaluate( a_projectileEnergy );
        }
#endif
    }

    return( multiplicity );
}

/* *********************************************************************************************************//**
 * Returns the energy dependent multiplicity for outgoing particle with pops intid *a_intid*. The returned value may not
 * be an integer. Energy dependent multiplicity mainly occurs for photons and fission neutrons.
 *
 * @param a_intid                   [in]    The PoPs intid of the requested particle.
 * @param a_projectileEnergy        [in]    The energy of the projectile.
 *
 * @return                                  The multiplicity value for the requested particle.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double Reaction::productAverageMultiplicityViaIntid( int a_intid, double a_projectileEnergy ) const {

    double multiplicity = 0.0;

    if( m_crossSectionThreshold > a_projectileEnergy ) return( multiplicity );

    std::size_t i1 = 0;
    for( Vector<int>::iterator iter = m_productIntids.begin( ); iter != m_productIntids.end( ); ++iter, ++i1 ) {
        if( *iter == a_intid ) {
            multiplicity = m_productMultiplicities[i1];
            break;
        }
    }

    if( multiplicity < 0 ) {
#ifdef MCGIDI_USE_OUTPUT_CHANNEL
        multiplicity = m_outputChannel->productAverageMultiplicityViaIntid( a_intid, a_projectileEnergy );
#else
        multiplicity = 0.0;
        for( auto productIter = m_products.begin( ); productIter != m_products.end( ); ++productIter ) {
            multiplicity += (*productIter)->productAverageMultiplicityViaIntid( a_intid, a_projectileEnergy );
        }

        if( ( m_totalDelayedNeutronMultiplicity != nullptr ) && ( a_intid == PoPI::Intids::neutron ) ) {
            multiplicity += m_totalDelayedNeutronMultiplicity->evaluate( a_projectileEnergy );
        }
#endif
    }

    return( multiplicity );
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs index *a_particleIndex*.
 *
 * @param a_particleIndex       [in]    The PoPs index of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void Reaction::setUserParticleIndex( int a_particleIndex, int a_userParticleIndex ) {

#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    m_outputChannel->setUserParticleIndex( a_particleIndex, a_userParticleIndex );
#else
    for( auto productIter = m_products.begin( ); productIter != m_products.end( ); ++productIter ) {
        (*productIter)->setUserParticleIndex( a_particleIndex, a_userParticleIndex );
    }

    for( auto iter = m_delayedNeutrons.begin( ); iter != m_delayedNeutrons.end( ); ++iter ) 
        (*iter)->setUserParticleIndex( a_particleIndex, a_userParticleIndex );
#endif

    for( std::size_t i1 = 0; i1 < m_productIndices.size( ); ++i1 ) {
        if( m_productIndices[i1] == a_particleIndex ) m_userProductIndices[i1] = a_userParticleIndex;
    }

    for( std::size_t i1 = 0; i1 < m_productIndicesTransportable.size( ); ++i1 ) {
        if( m_productIndicesTransportable[i1] == a_particleIndex ) m_userProductIndicesTransportable[i1] = a_userParticleIndex;
    }

    if( a_particleIndex == m_fissionResiduaIndex ) m_fissionResiduaUserIndex = a_userParticleIndex;

    if( m_GRIN_capture != nullptr ) m_GRIN_capture->setUserParticleIndex( a_particleIndex, a_userParticleIndex );
    if( m_GRIN_inelastic != nullptr ) m_GRIN_inelastic->setUserParticleIndex( a_particleIndex, a_userParticleIndex );
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs index *a_particleIntid*.
 *
 * @param a_particleIndex       [in]    The PoPs intid of the particle whose user intid is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void Reaction::setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex ) {

#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    m_outputChannel->setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
#else
    for( auto productIter = m_products.begin( ); productIter != m_products.end( ); ++productIter ) {
        (*productIter)->setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
    }

    for( auto iter = m_delayedNeutrons.begin( ); iter != m_delayedNeutrons.end( ); ++iter )
        (*iter)->setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
#endif

    for( std::size_t i1 = 0; i1 < m_productIntids.size( ); ++i1 ) {
        if( m_productIntids[i1] == a_particleIntid ) m_userProductIndices[i1] = a_userParticleIndex;
    }

    for( std::size_t i1 = 0; i1 < m_productIntidsTransportable.size( ); ++i1 ) {
        if( m_productIntidsTransportable[i1] == a_particleIntid ) m_userProductIndicesTransportable[i1] = a_userParticleIndex;
    }

    if( a_particleIntid == m_fissionResiduaIntid ) m_fissionResiduaUserIndex = a_userParticleIndex;

    if( m_GRIN_capture != nullptr ) m_GRIN_capture->setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
    if( m_GRIN_inelastic != nullptr ) m_GRIN_inelastic->setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
}

/* *********************************************************************************************************//**
 * This method calls the **setModelDBRC_data* method on the first product of *this* with *a_modelDBRC_data*.
 *
 * @param a_modelDBRC_data      [in]    The instance storing data needed to treat the DRRC upscatter mode.
 ***********************************************************************************************************/

LUPI_HOST void Reaction::setModelDBRC_data( Sampling::Upscatter::ModelDBRC_data *a_modelDBRC_data ) {

#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    m_outputChannel->setModelDBRC_data( a_modelDBRC_data );
#else
    m_products[0]->setModelDBRC_data( a_modelDBRC_data );
#endif
}

/* *********************************************************************************************************//**
 * Adds the associated orphan products of an orphan product reaction to *a_associatedOrphanProducts*.
 *
 * @param a_associatedOrphanProducts    [in]    A list where the associated orphan products are added to.
 ***********************************************************************************************************/

LUPI_HOST void Reaction::addOrphanProductToProductList( std::vector<Product *> &a_associatedOrphanProducts ) const {

#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    m_outputChannel->addOrphanProductToProductList( a_associatedOrphanProducts );
#else
    for( auto productIter = m_products.begin( ); productIter != m_products.end( ); ++productIter ) {
        a_associatedOrphanProducts.push_back( *productIter );
    }
#endif
}

/* *********************************************************************************************************//**
 * Adds the associated orphan products of an orphan product reaction to *a_associatedOrphanProducts*.
 *
 * @param a_associatedOrphanProducts    [in]    A list where the associated orphan products are added to.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Reaction::addOrphanProductToProductList( Vector<Product *> &a_associatedOrphanProducts ) const {

#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    m_outputChannel->addOrphanProductToProductList( a_associatedOrphanProducts );
#else
    for( auto productIter = m_products.begin( ); productIter != m_products.end( ); ++productIter ) {
        a_associatedOrphanProducts.push_back( *productIter );
    }
#endif
}

/* *********************************************************************************************************//**
 * Adds the associated orphan products of an orphan product to *a_associatedOrphanProducts*.
 *
 * @param a_associatedOrphanProducts    [in]    A list where the associated orphan products are added to.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Reaction::addOrphanProductToProductList( Vector<Reaction *> &a_orphanProducts ) {

    for( auto associatedOrphanProductIndex = m_associatedOrphanProductIndices.begin( );
            associatedOrphanProductIndex != m_associatedOrphanProductIndices.end( ); ++associatedOrphanProductIndex ) {
        Reaction *orphanProduct = a_orphanProducts[*associatedOrphanProductIndex];
        orphanProduct->addOrphanProductToProductList( m_associatedOrphanProducts );
    }
}

/* *********************************************************************************************************//**
 * Adds the contents of **a_associatedOrphanProductIndcies** and **a_associatedOrphanProducts** to this instance.
 *
 * @param a_associatedOrphanProductIndcies  [in]    The list of indices of the orphanProduct reaction that make up the product in **a_associatedOrphanProducts**.
 * @param a_associatedOrphanProducts        [in]    The list of pointers to the associated orphan products.
 ***********************************************************************************************************/

LUPI_HOST void Reaction::setOrphanProductData( std::vector<std::size_t> const &a_associatedOrphanProductIndcies,
                std::vector<Product *> const &a_associatedOrphanProducts ) {

    m_associatedOrphanProductIndices.reserve( a_associatedOrphanProductIndcies.size( ) );
    for( auto iter = a_associatedOrphanProductIndcies.begin( ); iter != a_associatedOrphanProductIndcies.end( ); ++iter )
        m_associatedOrphanProductIndices.push_back( *iter );

    m_associatedOrphanProducts.reserve( a_associatedOrphanProducts.size( ) );
    for( auto iter = a_associatedOrphanProducts.begin( ); iter != a_associatedOrphanProducts.end( ); ++iter )
        m_associatedOrphanProducts.push_back( *iter );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Reaction::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {
    
    DATA_MEMBER_SIZE_T( m_GIDI_reactionIndex, a_buffer, a_mode );
    DATA_MEMBER_STRING( m_label, a_buffer, a_mode );
    DATA_MEMBER_INT( m_ENDF_MT, a_buffer, a_mode );
    DATA_MEMBER_INT( m_ENDL_C, a_buffer, a_mode );
    DATA_MEMBER_INT( m_ENDL_S, a_buffer, a_mode );
    DATA_MEMBER_INT( m_initialStateIndex, a_buffer, a_mode );
    DATA_MEMBER_INT( m_neutronIndex, a_buffer, a_mode );
    DATA_MEMBER_CAST( m_hasFission, a_buffer, a_mode, bool );
    DATA_MEMBER_DOUBLE( m_projectileMass, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_targetMass, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_crossSectionThreshold, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_twoBodyThreshold, a_buffer, a_mode );
    DATA_MEMBER_CAST( m_hasFinalStatePhotons, a_buffer, a_mode, bool );
    DATA_MEMBER_INT( m_fissionResiduaIntid, a_buffer, a_mode );
    DATA_MEMBER_INT( m_fissionResiduaIndex, a_buffer, a_mode );
    DATA_MEMBER_INT( m_fissionResiduaUserIndex, a_buffer, a_mode );
    serializeFissionResiduals( m_fissionResiduals, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_fissionResidualMass, a_buffer, a_mode );

    DATA_MEMBER_VECTOR_INT( m_productIntids, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_INT( m_productIndices, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_INT( m_userProductIndices, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_INT( m_productMultiplicities, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_INT( m_productIntidsTransportable, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_INT( m_productIndicesTransportable, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_INT( m_userProductIndicesTransportable, a_buffer, a_mode );

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
            } }
        else if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
            a_buffer.incrementPlacement( sizeof( OutputChannel ) );
        }
        m_outputChannel->serialize( a_buffer, a_mode );
    }
#else
    serializeQs( a_buffer, a_mode, m_Qs );
    serializeProducts( a_buffer, a_mode, m_products );
    m_totalDelayedNeutronMultiplicity = serializeFunction1d( a_buffer, a_mode, m_totalDelayedNeutronMultiplicity );
    serializeDelayedNeutrons( a_buffer, a_mode, m_delayedNeutrons );
#endif

    DATA_MEMBER_VECTOR_SIZE_T( m_associatedOrphanProductIndices, a_buffer, a_mode );

    std::size_t vectorSize = m_associatedOrphanProducts.size( );
    int vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;
    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        m_associatedOrphanProducts.reserve( vectorSize, &a_buffer.m_placement ); }
    else if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        a_buffer.m_placement += m_associatedOrphanProducts.internalSize( );
    }

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        m_protareSingle = nullptr;
        m_reactionIndex = MCGIDI_nullReaction;
    }

    DATA_MEMBER_CAST( m_GRIN_specialSampleProducts, a_buffer, a_mode, bool );
    DATA_MEMBER_DOUBLE( m_GRIN_inelasticThreshold, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_GRIN_maximumCaptureIncidentEnergy, a_buffer, a_mode );

    if( m_GRIN_specialSampleProducts ) {
        if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
            if( a_buffer.m_placement != nullptr ) {
                if( m_ENDF_MT == 91 ) {
                    m_GRIN_inelastic = new(a_buffer.m_placement) GRIN_inelastic;
                    a_buffer.incrementPlacement( sizeof( GRIN_inelastic ) ); }
                else {
                    m_GRIN_capture = new(a_buffer.m_placement) GRIN_capture;
                    a_buffer.incrementPlacement( sizeof( GRIN_capture ) );
                } }
            else {
                if( m_ENDF_MT == 91 ) {
                    m_GRIN_inelastic = new GRIN_inelastic( ); }
                else {
                    m_GRIN_capture = new GRIN_capture( );
                }
            }
        }

        if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
            if( m_ENDF_MT == 91 ) {
                a_buffer.incrementPlacement( sizeof( GRIN_inelastic ) ); }
            else {
                a_buffer.incrementPlacement( sizeof( GRIN_capture ) );
            }
        }

        if( m_ENDF_MT == 91 ) {
            m_GRIN_inelastic->serialize( a_buffer, a_mode ); }
        else {
            m_GRIN_capture->serialize( a_buffer, a_mode );
        }
    }
}

}
