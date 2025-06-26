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
 * This class represents a **GNDS** <**outputChannel**> node with only data needed for Monte Carlo transport.
 */


/* *********************************************************************************************************//**
 * Default constructor used when broadcasting a Protare as needed by MPI or GPUs.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE OutputChannel::OutputChannel( ) :
        m_channelType( ChannelType::none ),
        m_neutronIndex( -1 ),
        m_isFission( false ),
        m_hasFinalStatePhotons( false ),
        m_Q( nullptr ),
        m_products( ),
        m_totalDelayedNeutronMultiplicity( nullptr ) {

}

/* *********************************************************************************************************//**
 * @param a_outputChannel       [in]    The GIDI::OutputChannel whose data is to be used to construct *this*.
 * @param a_setupInfo           [in]    Used internally when constructing a Protare to pass information to other constructors.
 * @param a_settings            [in]    Used to pass user options to the *this* to instruct it which data are desired.
 * @param a_particles           [in]    List of transporting particles and their information (e.g., multi-group boundaries and fluxes).
 ***********************************************************************************************************/

LUPI_HOST OutputChannel::OutputChannel( GIDI::OutputChannel const *a_outputChannel, SetupInfo &a_setupInfo, Transporting::MC const &a_settings, 
                GIDI::Transporting::Particles const &a_particles ) :
        m_channelType( ChannelType::none ),
        m_neutronIndex( a_setupInfo.m_neutronIndex ),
        m_isFission( false ),
        m_hasFinalStatePhotons( false ),
        m_Q( nullptr ),
        m_products( ),
        m_totalDelayedNeutronMultiplicity( nullptr ) {

    if( a_outputChannel != nullptr ) {
        m_channelType = a_outputChannel->twoBody( ) ? ChannelType::twoBody : ChannelType::uncorrelatedBodies;
        m_isFission = a_outputChannel->isFission( );

        m_Q = Functions::parseFunction1d_d1( a_outputChannel->Q( ).get<GIDI::Functions::Function1dForm>( 0 ) );
        if( a_setupInfo.m_isPairProduction ) {
            double domainMin = m_Q->domainMin( ), domainMax = m_Q->domainMax( );

            delete m_Q;
            m_Q = new Functions::Constant1d( domainMin, domainMax, 0.0 );
        }
        a_setupInfo.m_Q = m_Q->evaluate( 0 );                                  // Needed for NBodyPhaseSpace.

        GIDI::Suite const &products = a_outputChannel->products( );
        if( m_channelType == ChannelType::twoBody ) {
            if( !a_setupInfo.m_protare.isTNSL_ProtareSingle( ) ) {
                GIDI::Product const *product = products.get<GIDI::Product>( 1 );
                a_setupInfo.m_product2Mass = product->particle( ).mass( "MeV/c**2" );         // Includes nuclear excitation energy.
            }
        }

        bool electronPresent = false;
        std::size_t size = 0;
        std::set<std::size_t> productsToDo;
        for( std::size_t i1 = 0; i1 < a_outputChannel->products( ).size( ); ++i1 ) {
            GIDI::Product const *product = products.get<GIDI::Product>( i1 );

            if( product->particle( ).ID( ) == PoPI::IDs::electron ) electronPresent = true;
            if( electronPresent && !a_setupInfo.m_isPhotoAtomicIncoherentScattering ) continue;
            if( !electronPresent && !product->isCompleteParticle( ) && ( product->outputChannel( ) == nullptr ) ) continue;

            if( ( product->outputChannel( ) != nullptr ) || a_settings.sampleNonTransportingParticles( ) || a_particles.hasParticle( product->particle( ).ID( ) ) )
                productsToDo.insert( i1 );
        }
        size = productsToDo.size( );
        if( a_setupInfo.m_isPairProduction ) {
            size += 2;
            size = 2;                               // This is a kludge until the ENDL to GNDS translator is fixed.
        }

        bool addIncoherentPhotoAtomicScatteringElectron = false;
        if( a_setupInfo.m_isPhotoAtomicIncoherentScattering && a_particles.hasParticle( PoPI::IDs::electron ) ) {   // May need to add electron for legacy GNDS files.
            if( !electronPresent ) {
// FIXME: BRB 7/Nov/2024, Why is an electron added, this is incoherent atomic scattering which does not emit an electron?
                addIncoherentPhotoAtomicScatteringElectron = true;
                ++size;
            }
        }
        m_products.reserve( size );

        if( a_setupInfo.m_isPairProduction ) {
            std::string ID( PoPI::IDs::photon );
            std::string label = ID;

            Product *product = new Product( a_setupInfo.m_popsUser, ID, label );
            product->setMultiplicity( new Functions::Constant1d( a_setupInfo.m_domainMin, a_setupInfo.m_domainMax, 1.0, 0.0 ) );
            product->distribution( new Distributions::PairProductionGamma( a_setupInfo, true ) );
            m_products.push_back( product );

            label += "__a";
            product = new Product( a_setupInfo.m_popsUser, ID, label );
            product->setMultiplicity( new Functions::Constant1d( a_setupInfo.m_domainMin, a_setupInfo.m_domainMax, 1.0, 0.0 ) );
            product->distribution( new Distributions::PairProductionGamma( a_setupInfo, false ) );
            m_products.push_back( product );
        }

        for( std::size_t i1 = 0; i1 < a_outputChannel->products( ).size( ); ++i1 ) {
            if( productsToDo.find( i1 ) == productsToDo.end( ) ) continue;

            GIDI::Product const *product = products.get<GIDI::Product>( i1 );

            if( a_setupInfo.m_isPairProduction ) {
                if( !a_settings.sampleNonTransportingParticles( ) ) continue;
                if( a_setupInfo.m_protare.targetIntid( ) != MCGIDI_popsIntid( a_setupInfo.m_pops, product->particle( ).ID( ) ) ) continue;
            }
            a_setupInfo.m_twoBodyOrder = TwoBodyOrder::notApplicable;
            if( m_channelType == ChannelType::twoBody ) a_setupInfo.m_twoBodyOrder = ( ( i1 == 0 ? TwoBodyOrder::firstParticle : TwoBodyOrder::secondParticle ) );
            m_products.push_back( new Product( product, a_setupInfo, a_settings, a_particles, m_isFission ) );

            if( addIncoherentPhotoAtomicScatteringElectron && ( product->particle( ).ID( ) == PoPI::IDs::photon ) ) {
                addIncoherentPhotoAtomicScatteringElectron = false;

                Product *product2 = new Product( a_setupInfo.m_pops, PoPI::IDs::electron, PoPI::IDs::electron );
                product2->setMultiplicity( new Functions::Constant1d( a_setupInfo.m_domainMin, a_setupInfo.m_domainMax, 1.0, 0.0 ) );
                product2->distribution( new Distributions::IncoherentPhotoAtomicScatteringElectron( a_setupInfo ) );
                m_products.push_back( product2 );
            }
        }

        if( ( a_settings.delayedNeutrons( ) == GIDI::Transporting::DelayedNeutrons::on ) && a_particles.hasParticle( PoPI::IDs::neutron ) ) {
            GIDI::FissionFragmentData const &fissionFragmentData = a_outputChannel->fissionFragmentData( );
            GIDI::Suite const &delayedNeutrons = fissionFragmentData.delayedNeutrons( );

            if( delayedNeutrons.size( ) > 0 ) {
                bool missingData = false;
                GIDI::Axes axes;
                GIDI::Functions::XYs1d totalDelayedNeutronMultiplicity( axes, ptwXY_interpolationLinLin );

                m_delayedNeutrons.reserve( delayedNeutrons.size( ) );
                for( std::size_t i1 = 0; i1 < delayedNeutrons.size( ); ++i1 ) {
                    GIDI::DelayedNeutron const *delayedNeutron = delayedNeutrons.get<GIDI::DelayedNeutron>( i1 );
                    GIDI::Product const &product = delayedNeutron->product( );
                    GIDI::Suite const &multiplicity = product.multiplicity( );

                    GIDI::Functions::Function1dForm const *form1d = multiplicity.get<GIDI::Functions::Function1dForm>( 0 );

                    if( form1d->type( ) == GIDI::FormType::unspecified1d ) {
                        missingData = true;
                        break;
                    }

                    if( form1d->type( ) != GIDI::FormType::XYs1d ) {
                        std::cerr << "OutputChannel::OutputChannel: GIDI::DelayedNeutron multiplicity type != GIDI::FormType::XYs1d" << std::endl;
                        missingData = true;
                        break;
                    }

                    GIDI::Functions::XYs1d const *multiplicityXYs1d = static_cast<GIDI::Functions::XYs1d const *>( form1d );
                    totalDelayedNeutronMultiplicity += *multiplicityXYs1d;

                    m_delayedNeutrons.push_back( new DelayedNeutron( static_cast<int>( i1 ), delayedNeutron, a_setupInfo, a_settings, a_particles ) );
                }
                if( !missingData ) m_totalDelayedNeutronMultiplicity = new Functions::XYs1d( totalDelayedNeutronMultiplicity );
            }
        }

        m_hasFinalStatePhotons = a_setupInfo.m_hasFinalStatePhotons;
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE OutputChannel::~OutputChannel( ) {

    delete m_Q;
    for( std::size_t i1 = 0; i1 < m_products.size( ); ++i1 ) delete m_products[i1];

    delete m_totalDelayedNeutronMultiplicity;
    for( std::size_t i1 = 0; i1 < m_delayedNeutrons.size( ); ++i1 ) delete m_delayedNeutrons[i1];
}

/* *********************************************************************************************************//**
 * This method returns the final Q for *this* by getting its final Q plus any sub-output channel's finalQ.
 *
 * @param a_x1                  [in]    The energy of the projectile.
 *
 * @return                              The Q-value at product energy *a_x1*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double OutputChannel::finalQ( double a_x1 ) const {

    double final_Q = m_Q->evaluate( a_x1 );

    for( std::size_t i1 = 0; i1 < m_products.size( ); ++i1 ) final_Q += m_products[i1]->finalQ( a_x1 );
    return( final_Q );
}

/* *********************************************************************************************************//**
 * This method returns *true* if the output channel or any of its sub-output channels is a fission channel and *false* otherwise.
 *
 * @return                              *true* if *this* or any sub-output channel is a fission channel and *false* otherwise.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE bool OutputChannel::hasFission( ) const {

    if( m_isFission ) return( true );
    for( std::size_t i1 = 0; i1 < m_products.size( ); ++i1 ) {
        if( m_products[i1]->hasFission( ) ) return( true );
    }
    return( false );
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs index *a_particleIndex*.
 *  
 * @param a_particleIndex       [in]    The PoPs index of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void OutputChannel::setUserParticleIndex( int a_particleIndex, int a_userParticleIndex ) {

    for( auto iter = m_products.begin( ); iter != m_products.end( ); ++iter ) (*iter)->setUserParticleIndex( a_particleIndex, a_userParticleIndex );
    for( auto iter = m_delayedNeutrons.begin( ); iter != m_delayedNeutrons.end( ); ++iter ) (*iter)->setUserParticleIndex( a_particleIndex, a_userParticleIndex );
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs intid *a_particleIntid*.
 *  
 * @param a_particleIntid       [in]    The PoPs intid of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/
 
LUPI_HOST void OutputChannel::setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex ) {

    for( auto iter = m_products.begin( ); iter != m_products.end( ); ++iter ) (*iter)->setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
    for( auto iter = m_delayedNeutrons.begin( ); iter != m_delayedNeutrons.end( ); ++iter ) (*iter)->setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
} 

/* *********************************************************************************************************//**
 * This method calls the **setModelDBRC_data* method on the first product of *this* with *a_modelDBRC_data*.
 *
 * @param a_modelDBRC_data      [in]    The instance storing data needed to treat the DRRC upscatter mode.
 ***********************************************************************************************************/

LUPI_HOST void OutputChannel::setModelDBRC_data( Sampling::Upscatter::ModelDBRC_data *a_modelDBRC_data ) {

    m_products[0]->setModelDBRC_data( a_modelDBRC_data );
}

/* *********************************************************************************************************//**
 * Returns the energy dependent multiplicity for outgoing particle with pops id *a_id*. The returned value may not
 * be an integer. Energy dependent multiplicity mainly occurs for photons and fission neutrons.
 *
 * @param a_products                            [in]    The std::vector instance to add products to.
 * @param a_totalDelayedNeutronMultiplicity     [in]    The std::vector instance to add delayed neutron multiplicities to.
 * @param a_delayedNeutrons                     [in]    The std::vector instance to add delayed neutrons to.
 * @param a_Qs                                  [in]    The std::vector instance to add Q functions to.
 ***********************************************************************************************************/

LUPI_HOST void OutputChannel::moveProductsEtAlToReaction( std::vector<Product *> &a_products, Functions::Function1d **a_totalDelayedNeutronMultiplicity, 
                std::vector<DelayedNeutron *> &a_delayedNeutrons, std::vector<Functions::Function1d_d1 *> &a_Qs  ) {

    if( a_totalDelayedNeutronMultiplicity != nullptr ) {    /* This will not work if fission is a nested channel. Ergo "n + (R -> fission)". */
        *a_totalDelayedNeutronMultiplicity = m_totalDelayedNeutronMultiplicity;
        m_totalDelayedNeutronMultiplicity = nullptr;
        for( std::size_t index = 0; index < m_delayedNeutrons.size( ); ++index ) {
            a_delayedNeutrons.push_back( m_delayedNeutrons[index] );
            m_delayedNeutrons[index] = nullptr;
        }
    }

    a_Qs.push_back( m_Q );
    m_Q = nullptr;
    for( std::size_t productIndex = 0; productIndex < m_products.size( ); ++productIndex ) {
        Product *product = m_products[productIndex];

        if( product->outputChannel( ) != nullptr ) {
            product->outputChannel( )->moveProductsEtAlToReaction( a_products, nullptr, a_delayedNeutrons, a_Qs );
            delete product; }
        else {
            a_products.push_back( product );
        }
        m_products[productIndex] = nullptr;
    }
}

#ifdef MCGIDI_USE_OUTPUT_CHANNEL
/* *********************************************************************************************************//**
 * Adds the associated orphan products of an orphan product to *a_products*.
 *
 * @param a_products                [in]    The std::vector instance to add products to.
 *
 * @return                                  This does not seem to be working. Needs work.
 ***********************************************************************************************************/

LUPI_HOST void OutputChannel::addOrphanProductToProductList( std::vector<Product *> &a_products ) const {

    for( int productIndex = 0; productIndex < m_products.size( ); ++productIndex ) {
        Product *product = m_products[productIndex];

        if( product->outputChannel( ) != nullptr ) {
            product->outputChannel( )->addOrphanProductToProductList( a_products ); }
        else {
            a_products.push_back( product );
        }
    }

}

/* *********************************************************************************************************//**
 * Adds the associated orphan products of an orphan product to *a_products*.
 *
 * @param a_products                [in]    The std::vector instance to add products to.
 *
 * @return                                  This does not seem to be working. Needs work.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void OutputChannel::addOrphanProductToProductList( Vector<Product *> &a_products ) const {

    for( int productIndex = 0; productIndex < m_products.size( ); ++productIndex ) {
        Product *product = m_products[productIndex];

        if( product->outputChannel( ) != nullptr ) {
            product->outputChannel( )->addOrphanProductToProductList( a_products ); }
        else {
            a_products.push_back( product );
        }
    }

}

/* *********************************************************************************************************//**
 * Adds the associated orphan products of an orphan product to *a_products*.
 *
 * @param a_products                [in]    The std::vector instance to add products to.
 *
 * @return                                  This does not seem to be working. Needs work.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void OutputChannel::addOrphanProductToProductList( Vector<Product *> &a_products ) const {

    for( int productIndex = 0; productIndex < m_products.size( ); ++productIndex ) {
        Product *product = m_products[productIndex];

        if( product->outputChannel( ) != nullptr ) {
            product->outputChannel( )->addOrphanProductToProductList( a_products ); }
        else {
            a_products.push_back( product );
        }
    }

}
#endif

/* *********************************************************************************************************//**
 * Returns the energy dependent multiplicity for outgoing particle with pops index *a_index*. The returned value may not
 * be an integer. Energy dependent multiplicity mainly occurs for photons and fission neutrons.
 *
 * @param a_index                   [in]    The index of the requested particle.
 * @param a_projectileEnergy        [in]    The energy of the projectile.
 *
 * @return                                  The multiplicity value for the requested particle.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double OutputChannel::productAverageMultiplicity( int a_index, double a_projectileEnergy ) const {

    double multiplicity = 0.0;

    for( Vector<Product *>::const_iterator iter = m_products.begin( ); iter != m_products.end( ); ++iter ) {
        multiplicity += (*iter)->productAverageMultiplicity( a_index, a_projectileEnergy );
    }

    if( m_totalDelayedNeutronMultiplicity != nullptr ) {
        if( a_index == m_neutronIndex ) multiplicity += m_totalDelayedNeutronMultiplicity->evaluate( a_projectileEnergy );
    }

    return( multiplicity );
}

/* *********************************************************************************************************//**
 * Returns the energy dependent multiplicity for outgoing particle with pops intid *a_intid*. The returned value may not
 * be an integer. Energy dependent multiplicity mainly occurs for photons and fission neutrons.
 *
 * @param a_intid                   [in]    The intid of the requested particle.
 * @param a_projectileEnergy        [in]    The energy of the projectile.
 *
 * @return                                  The multiplicity value for the requested particle.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double OutputChannel::productAverageMultiplicityViaIntid( int a_intid, double a_projectileEnergy ) const {

    double multiplicity = 0.0;

    for( Vector<Product *>::const_iterator iter = m_products.begin( ); iter != m_products.end( ); ++iter ) {
        multiplicity += (*iter)->productAverageMultiplicityViaIntid( a_intid, a_projectileEnergy );
    }

    if( m_totalDelayedNeutronMultiplicity != nullptr ) {
        if( a_intid == PoPI::Intids::neutron ) multiplicity += m_totalDelayedNeutronMultiplicity->evaluate( a_projectileEnergy );
    }

    return( multiplicity );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void OutputChannel::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    int channelType = 0;
    switch( m_channelType ) {
    case ChannelType::none :
        break;
    case ChannelType::twoBody :
        channelType = 1;
        break;
    case ChannelType::uncorrelatedBodies :
        channelType = 2;
        break;
    }
    DATA_MEMBER_INT( channelType, a_buffer, a_mode );
    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        switch( channelType ) {
        case 0 :
            m_channelType = ChannelType::none;
            break;
        case 1 :
            m_channelType = ChannelType::twoBody;
            break;
        case 2 :
            m_channelType = ChannelType::uncorrelatedBodies;
            break;
        }
    }

    DATA_MEMBER_INT( m_neutronIndex, a_buffer, a_mode );
    DATA_MEMBER_CAST( m_isFission, a_buffer, a_mode, bool );
    DATA_MEMBER_CAST( m_hasFinalStatePhotons, a_buffer, a_mode, bool );

    m_Q = serializeFunction1d_d1( a_buffer, a_mode, m_Q );
    serializeProducts( a_buffer, a_mode, m_products );
    m_totalDelayedNeutronMultiplicity = serializeFunction1d( a_buffer, a_mode, m_totalDelayedNeutronMultiplicity );
    serializeDelayedNeutrons( a_buffer, a_mode, m_delayedNeutrons );
}

}
