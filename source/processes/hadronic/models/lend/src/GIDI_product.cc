/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdlib.h>
#include "GIDI.hpp"
#include <HAPI.hpp>

namespace GIDI {

/*! \class Product
 * Class to store a GNDS <**product**> node.
 */

/* *********************************************************************************************************//**
 * Constructed 
 ***********************************************************************************************************/

Product::Product( PoPI::Database const &a_pops, std::string const &a_productID, std::string const &a_label ) :
        Form( FormType::product ),
        m_particle( ParticleInfo( a_productID, a_pops, a_pops, true ) ),
        m_GNDS_particle( ParticleInfo( a_productID, a_pops, a_pops, true ) ),
        m_productMultiplicity( 0 ),
        m_treatProductAsIfInfinityMass( false ),
        m_multiplicity( GIDI_multiplicityChars, GIDI_labelChars ),
        m_distribution( GIDI_distributionChars, GIDI_labelChars ),
        m_averageEnergy( GIDI_averageEnergyChars, GIDI_labelChars ),
        m_averageMomentum( GIDI_averageMomentumChars, GIDI_labelChars ),
        m_outputChannel( nullptr ) {

    setMoniker( GIDI_productChars );
    setLabel( a_label );

    m_multiplicity.setAncestor( this );
    m_distribution.setAncestor( this );
    m_averageEnergy.setAncestor( this );
    m_averageMomentum.setAncestor( this );
}

/* *********************************************************************************************************//**
 * Constructed from data in a <**product**> node.
 *
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The HAPI::Node to be parsed and used to construct the Product.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops            [in]    The *external* PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs    [in]    The *internal* PoPI::Database instance used to get particle indices and possibly other particle information.
 *                                  This is the <**PoPs**> node under the <**reactionSuite**> node.
 * @param a_parent          [in]    The **m_products** member of GIDI::OutputChannel this product belongs to.
 * @param a_styles          [in]    The <**styles**> node under the <**reactionSuite**> node.
 ***********************************************************************************************************/

Product::Product( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, 
                PoPI::Database const &a_internalPoPs, Suite *a_parent, Styles::Suite const *a_styles ) :
        Form( a_node, a_setupInfo, FormType::product, a_parent ),
        m_particle( a_node.attribute_as_string( GIDI_pidChars ), a_pops, a_internalPoPs, false ),
        m_GNDS_particle( a_node.attribute_as_string( GIDI_pidChars ), a_pops, a_internalPoPs, false ),
        m_productMultiplicity( 0 ),
        m_treatProductAsIfInfinityMass( false ),
        m_multiplicity( a_construction, GIDI_multiplicityChars, GIDI_labelChars, a_node, a_setupInfo, a_pops, a_internalPoPs, parseMultiplicitySuite, a_styles ),
        m_distribution( a_construction, GIDI_distributionChars, GIDI_labelChars, a_node, a_setupInfo, a_pops, a_internalPoPs, parseDistributionSuite, a_styles ),
        m_averageEnergy( a_construction, GIDI_averageEnergyChars, GIDI_labelChars, a_node, a_setupInfo, a_pops, a_internalPoPs, parseAverageEnergySuite, a_styles ),
        m_averageMomentum( a_construction, GIDI_averageMomentumChars, GIDI_labelChars, a_node, a_setupInfo, a_pops, a_internalPoPs, parseAverageMomentumSuite, a_styles ),
        m_outputChannel( nullptr ) {

    if( a_setupInfo.m_outputChannelLevel == 0 ) a_setupInfo.m_initialState = m_particle.ID( );
    m_multiplicity.setAncestor( this );
    m_distribution.setAncestor( this );
    m_averageEnergy.setAncestor( this );
    m_averageMomentum.setAncestor( this );

    auto iter = a_setupInfo.m_particleSubstitution->find( m_GNDS_particle.ID( ) );
    if( iter != a_setupInfo.m_particleSubstitution->end( ) ) m_particle = iter->second;

    if( a_setupInfo.m_protare->projectile( ).ID( ) == PoPI::IDs::photon ) {
        if( m_particle.ID( ) == a_setupInfo.m_protare->GNDS_target( ).ID( ) ) m_particle = a_setupInfo.m_protare->target( );
    }

    if( a_setupInfo.m_protare->isPhotoAtomic( ) || a_setupInfo.m_protare->isTNSL_ProtareSingle( ) ) {
        m_treatProductAsIfInfinityMass = m_GNDS_particle.ID( ) == a_setupInfo.m_protare->GNDS_target( ).ID( );
    }

    HAPI::Node const _outputChannel = a_node.child( GIDI_outputChannelChars );
    a_setupInfo.m_outputChannelLevel += 1;
    if( ! _outputChannel.empty( ) ) 
        m_outputChannel = new OutputChannel( a_construction, _outputChannel, a_setupInfo, a_pops, a_internalPoPs, a_styles, false, false );
    a_setupInfo.m_outputChannelLevel -= 1;

    if( m_outputChannel == nullptr ) {
        if( m_multiplicity.size( ) > 0 ) {
            GIDI::Functions::Function1dForm const *function1d = m_multiplicity.get<GIDI::Functions::Function1dForm>( 0 );

            if( function1d->type( ) == FormType::constant1d ) {
                m_productMultiplicity = static_cast<int>( function1d->evaluate( 0.0 ) ); }
            else if( function1d->type( ) != FormType::unspecified1d ) {
                m_productMultiplicity = -1;
            }
        } }
    else {
        m_outputChannel->setAncestor( this );
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Product::~Product( ) {

    if( m_outputChannel != nullptr ) delete m_outputChannel;
}

/* *********************************************************************************************************//**
 * Returns the maximum product depth for this product.
 *
 * @return The maximum product depth.
 ***********************************************************************************************************/

int Product::depth( ) const {

    if( m_outputChannel == nullptr ) return( 0 );
    return( m_outputChannel->depth( ) );
}

/* *********************************************************************************************************//**
 * Only for internal use. Called by ProtareTNSL instance to zero the lower energy multi-group data covered by the ProtareSingle that
 * contains the TNSL data covers the lower energy multi-group data.
 *
 * @param a_maximumTNSL_MultiGroupIndex     [in]    A map that contains labels for heated multi-group data and the last valid group boundary
 *                                                  for the TNSL data for that boundary.
 ***********************************************************************************************************/

void Product::modifiedMultiGroupElasticForTNSL( std::map<std::string,std::size_t> const &a_maximumTNSL_MultiGroupIndex ) {

    m_multiplicity.modifiedMultiGroupElasticForTNSL( a_maximumTNSL_MultiGroupIndex );
    m_distribution.modifiedMultiGroupElasticForTNSL( a_maximumTNSL_MultiGroupIndex );
    m_averageEnergy.modifiedMultiGroupElasticForTNSL( a_maximumTNSL_MultiGroupIndex );
    m_averageMomentum.modifiedMultiGroupElasticForTNSL( a_maximumTNSL_MultiGroupIndex );
}

/* *********************************************************************************************************//**
 * Returns true if the product has an output channel and its output channel hasFission returns true, and false otherwise.
 *
 * @return  true if at least one output channel is a fission channel.
 ***********************************************************************************************************/

bool Product::hasFission( ) const {

    if( m_outputChannel != nullptr ) return( m_outputChannel->hasFission( ) );
    return( false );
}

/* *********************************************************************************************************//**
 * Returns **false* if *this* has delayed fission neutrons and they are not complete; otherwise, returns **true**.
 *
 * @param       a_isDelayedNeutron  [in]    **true** is called from FissionFragmentData::isDelayedFissionNeutronComplete and **false** otherwise.
 *
 * @return      bool
 ***********************************************************************************************************/

bool Product::isDelayedFissionNeutronComplete( bool a_isDelayedNeutron ) const {

    if( a_isDelayedNeutron ) {
        return( isCompleteParticle( ) ); }
    else {
        if( m_outputChannel != nullptr ) return( m_outputChannel->isDelayedFissionNeutronComplete( ) );
    }

    return( true );
}

/* *********************************************************************************************************//**
 * Returns **true** if all outgoing particles (i.e., products) are specifed in *a_particles*. That is, the user
 * will be tracking all products of *this* reaction.
 *
 * @param a_particles           [in]    The list of particles to be transported.
 *
 * @return                              bool.
 ***********************************************************************************************************/

bool Product::areAllProductsTracked( Transporting::Particles const &a_particles ) const {

    if( m_outputChannel != nullptr ) return( m_outputChannel->areAllProductsTracked( a_particles ) );

    return( a_particles.hasParticle( m_particle.ID( ) ) );
}

/* *********************************************************************************************************//**
 * Used by GUPI::Ancestry to tranverse GNDS nodes. This method returns a pointer to a derived class' a_item member or nullptr if none exists.
 *
 * @param a_item    [in]    The name of the class member whose pointer is to be return.
 * @return                  The pointer to the class member or nullptr if class does not have a member named a_item.
 ***********************************************************************************************************/

GUPI::Ancestry *Product::findInAncestry3( std::string const &a_item ) {

    if( a_item == GIDI_multiplicityChars ) return( &m_multiplicity );
    if( a_item == GIDI_distributionChars ) return( &m_distribution );
    if( a_item == GIDI_averageEnergyChars ) return( &m_averageEnergy );
    if( a_item == GIDI_averageMomentumChars ) return( &m_averageMomentum );
    if( a_item == GIDI_outputChannelChars ) return( m_outputChannel );

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Used by GUPI::Ancestry to tranverse GNDS nodes. This method returns a pointer to a derived class' a_item member or nullptr if none exists.
 *
 * @param a_item    [in]    The name of the class member whose pointer is to be return.
 * @return                  The pointer to the class member or nullptr if class does not have a member named a_item.
 ***********************************************************************************************************/

GUPI::Ancestry const *Product::findInAncestry3( std::string const &a_item ) const {

    if( a_item == GIDI_multiplicityChars ) return( &m_multiplicity );
    if( a_item == GIDI_distributionChars ) return( &m_distribution );
    if( a_item == GIDI_averageEnergyChars ) return( &m_averageEnergy );
    if( a_item == GIDI_averageMomentumChars ) return( &m_averageMomentum );
    if( a_item == GIDI_outputChannelChars ) return( m_outputChannel );

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Insert a std::set with the products index and any product in in its output channel.
 * If a_transportablesOnly is true, only transportable product indices are return.
 *
 * @param a_indices                 [out]   The unique list of product indices.
 * @param a_particles               [in]    The list of particles to be transported.
 * @param a_transportablesOnly      [in]    If true, only transportable product indices are added in the list.
 ***********************************************************************************************************/

void Product::productIDs( std::set<std::string> &a_indices, Transporting::Particles const &a_particles, bool a_transportablesOnly ) const {

    if( m_outputChannel == nullptr ) {
        if( a_transportablesOnly && !a_particles.hasParticle( m_particle.ID( ) ) ) return;
        if( m_particle.ID( ) != "" ) a_indices.insert( m_particle.ID( ) ); }
    else {
        m_outputChannel->productIDs( a_indices, a_particles, a_transportablesOnly );
    }
}

/* *********************************************************************************************************//**
 * Returns the product multiplicity (e.g., 0, 1, 2, ...) or -1 if energy dependent or not an integer for particle with id *a_productID*.
 *
 * @param a_productID;              [in]    The id of the requested particle.
 *
 * @return                                  The multiplicity for the requested particle.
 ***********************************************************************************************************/

int Product::productMultiplicity( std::string const &a_productID ) const {

    if( m_outputChannel != nullptr ) return( m_outputChannel->productMultiplicity( a_productID ) );

    if( PoPI::compareSpecialParticleIDs( a_productID, m_particle.ID( ) ) ) return( m_productMultiplicity );
    return( 0 );
}

/* *********************************************************************************************************//**
 * Determines the maximum Legredre order present in the multi-group transfer matrix for a this product and any sub-products for a give label.
 *
 * @param a_smr             [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings        [in]    Specifies the requested label.
 * @param a_temperatureInfo [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID       [in]    Particle id for the requested product.
 *
 * @return                          The maximum Legredre order. If no transfer matrix data are present for the requested product, -1 is returned.
 ***********************************************************************************************************/

int Product::maximumLegendreOrder( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    int _maximumLegendreOrder = -1;

    if( m_outputChannel == nullptr ) {
        if( PoPI::compareSpecialParticleIDs( m_particle.ID( ), a_productID ) ) {
            Distributions::MultiGroup3d const *form = dynamic_cast<Distributions::MultiGroup3d const*>( a_settings.form( 
                    a_smr, m_distribution, a_temperatureInfo, "distribution for maximumLegendreOrder" ) );
            if( form != nullptr ) {
                Functions::Gridded3d const &gdata = form->data( );
                Array3d const &data = gdata.data( );
                _maximumLegendreOrder = static_cast<int>( data.size( ) ) - 1;
            }
        } }
    else {
        int __maximumLegendreOrder = m_outputChannel->maximumLegendreOrder( a_smr, a_settings, a_temperatureInfo, a_productID );
        if( __maximumLegendreOrder > _maximumLegendreOrder ) _maximumLegendreOrder = __maximumLegendreOrder;
    }

    return( _maximumLegendreOrder );
}

/* *********************************************************************************************************//**
 * Returns the sum of the multi-group, multiplicity for the requested label for the this product and any sub-product. 
 * This is a cross section weighted multiplicity.
 *
 * @param a_smr             [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings        [in]    Specifies the requested label.
 * @param a_temperatureInfo [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID       [in]    Particle id for the requested product.
 *
 * @return                          The requested multi-group multiplicity as a GIDI::Vector.
 ***********************************************************************************************************/

Vector Product::multiGroupMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    Vector vector( 0 );

    if( m_outputChannel == nullptr ) {
        if( ( PoPI::compareSpecialParticleIDs( m_particle.ID( ), a_productID ) ) && ( !m_treatProductAsIfInfinityMass ) ) {
            Functions::Gridded1d const *form = dynamic_cast<Functions::Gridded1d const*>( a_settings.form( 
                    a_smr, m_multiplicity, a_temperatureInfo, "multiplicity" ) );
            if( form != nullptr ) vector += form->data( );
        } }
    else {
        vector += m_outputChannel->multiGroupMultiplicity( a_smr, a_settings, a_temperatureInfo, a_productID );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the sum of the multi-group, Q for the requested label for the this product and any sub-product . This is a cross section weighted Q.
 * If a_final is false, only the Q for the products output channel is returned, otherwise, the Q for all output channels
 * summed, including output channels for each products.
 *
 * @param a_smr             [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings        [in]    Specifies the requested label.
 * @param a_temperatureInfo [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_final           [in]    If true, the Q is calculated for all output channels, including those for products.
 *
 * @return                          The requested multi-group Q as a GIDI::Vector.
 ***********************************************************************************************************/

Vector Product::multiGroupQ( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, bool a_final ) const {
    
    Vector _vector( 0 );

    if( m_outputChannel != nullptr ) _vector += m_outputChannel->multiGroupQ( a_smr, a_settings, a_temperatureInfo, a_final );

    return( _vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, product matrix for the requested label for the requested product index for the requested Legendre order.
 * If no data are found, an empty GIDI::Matrix is returned.
 *
 * @param a_smr             [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings        [in]    Specifies the requested label and if delayed neutrons should be included.
 * @param a_temperatureInfo [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_particles       [in]    The list of particles to be transported.
 * @param a_productID       [in]    Particle id for the requested product.
 * @param a_order           [in]    Requested product matrix, Legendre order.
 *
 * @return                          The requested multi-group product matrix as a GIDI::Matrix.
 ***********************************************************************************************************/

Matrix Product::multiGroupProductMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, std::string const &a_productID, 
                std::size_t a_order ) const {

    Matrix matrix( 0, 0 );

    if( m_outputChannel == nullptr ) {
        if( ( PoPI::compareSpecialParticleIDs( m_particle.ID( ), a_productID ) ) && ( !m_treatProductAsIfInfinityMass ) ) {
            Distributions::MultiGroup3d const *form = dynamic_cast<Distributions::MultiGroup3d const*>( a_settings.form( 
                    a_smr, m_distribution, a_temperatureInfo, "distribution for product matrix" ) );
            if( form != nullptr ) {
                Functions::Gridded3d const &gdata = form->data( );
                Array3d const &data = gdata.data( );
                matrix = data.matrix( a_order );
            }
        } }
    else {
        matrix += m_outputChannel->multiGroupProductMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, a_productID, a_order );
    }

    return( matrix );
}

/* *********************************************************************************************************//**
 * Returns the sum of the multi-group, average energy for the requested label for the requested product. This is a cross section weighted average energy
 * summed over this and all sub-products.
 *
 * @param a_smr             [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings        [in]    Specifies the requested label.
 * @param a_temperatureInfo [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID       [in]    Particle id for the requested product.
 *
 * @return                          The requested multi-group average energy as a GIDI::Vector.
 ***********************************************************************************************************/

Vector Product::multiGroupAverageEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    Vector vector( 0 );

    if( m_outputChannel == nullptr ) {
        if( ( PoPI::compareSpecialParticleIDs( m_particle.ID( ), a_productID ) ) && ( !m_treatProductAsIfInfinityMass ) ) {
            Functions::Gridded1d const *form = dynamic_cast<Functions::Gridded1d const*>( a_settings.form( 
                    a_smr, m_averageEnergy, a_temperatureInfo, "average product energy" ) );
            if( form != nullptr ) vector += form->data( );
        } }
    else {
        vector += m_outputChannel->multiGroupAverageEnergy( a_smr, a_settings, a_temperatureInfo, a_productID );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the sum of the multi-group, average momentum for the requested label for the requested product. This is a cross section weighted average momentum 
 * summed over this and all sub-products.
 *
 * @param a_smr             [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings        [in]    Specifies the requested label.
 * @param a_temperatureInfo [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID       [in]    Particle id for the requested product.
 *
 * @return                          The requested multi-group average momentum as a GIDI::Vector.
 ***********************************************************************************************************/

Vector Product::multiGroupAverageMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    Vector vector( 0 );

    if( m_outputChannel == nullptr ) {
        if( ( PoPI::compareSpecialParticleIDs( m_particle.ID( ), a_productID ) ) && ( !m_treatProductAsIfInfinityMass ) ) {
            Functions::Gridded1d const *form = dynamic_cast<Functions::Gridded1d const*>( a_settings.form( 
                    a_smr, m_averageMomentum, a_temperatureInfo, "average product momentum" ) );
            if( form != nullptr ) vector += form->data( );
        } }
    else {
        vector += m_outputChannel->multiGroupAverageMomentum( a_smr, a_settings, a_temperatureInfo, a_productID );
    }
    return( vector );
}

/* *********************************************************************************************************//**
 * Returns, via arguments, the average energy and momentum, and gain for product with particle id *a_particleID*.
 *
 * @param a_settings                    [in]    Specifies the requested label.
 * @param a_particleID                  [in]    The particle id of the product.
 * @param a_energy                      [in]    The energy of the projectile.
 * @param a_productEnergy               [in]    The average energy of the product.
 * @param a_productMomentum             [in]    The average momentum of the product.
 * @param a_productGain                 [in]    The gain of the product.
 * @param a_ignoreIncompleteParticles   [in]    If *true*, incomplete particles are ignore, otherwise a *throw* is executed.
 ***********************************************************************************************************/

void Product::continuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID, double a_energy, 
                double &a_productEnergy, double &a_productMomentum, double &a_productGain, bool a_ignoreIncompleteParticles ) const {

    if( m_outputChannel == nullptr ) {
        if( a_particleID == m_particle.ID( ) ) {
            if( !isCompleteParticle( ) ) {
                if( a_ignoreIncompleteParticles ) return;
                throw Exception( "GIDI::Product::continuousEnergyProductData: particle is incomplete: " + toXLink( ) + "." );
            }
            a_productEnergy += averageEnergy( ).get<GIDI::Functions::Function1dForm>( 0 )->evaluate( a_energy );
            a_productMomentum += averageMomentum( ).get<GIDI::Functions::Function1dForm>( 0 )->evaluate( a_energy );
            a_productGain += multiplicity( ).get<GIDI::Functions::Function1dForm>( 0 )->evaluate( a_energy ); } }
    else {
        m_outputChannel->continuousEnergyProductData( a_settings, a_particleID, a_energy, a_productEnergy, a_productMomentum, a_productGain, a_ignoreIncompleteParticles );
    }
}

/* *********************************************************************************************************//**
 * Modifies the average product energies, momenta and gains for product with particle id *a_particleID*.
 * 
 * @param a_settings                    [in]    Specifies user options.
 * @param a_particleID                  [in]    The particle id of the product.
 * @param a_energies                    [in]    The vector of energies to map the data to.
 * @param a_offset                      [in]    The index of the first energy whose data are to be added to the vectors.
 * @param a_productEnergies             [out]   The vector of average energies of the product.
 * @param a_productMomenta              [out]   The vector of average momenta of the product.
 * @param a_productGains                [out]   The vector of gain of the product.
 * @param a_ignoreIncompleteParticles   [in]    If *true*, incomplete particles are ignore, otherwise a *throw* is executed.
 ***********************************************************************************************************/
 
void Product::mapContinuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID, 
                std::vector<double> const &a_energies, std::size_t a_offset, std::vector<double> &a_productEnergies, std::vector<double> &a_productMomenta, 
                std::vector<double> &a_productGains, bool a_ignoreIncompleteParticles ) const {

    if( m_outputChannel == nullptr ) {
        if( a_particleID == m_particle.ID( ) ) {
            if( !isCompleteParticle( ) ) {
                if( a_ignoreIncompleteParticles ) return;
                throw Exception( "GIDI::Product::mapContinuousEnergyProductData: particle is incomplete: " + toXLink( ) + "." );
            }
            GIDI::Functions::Function1dForm const *multiplicity1 = multiplicity( ).get<GIDI::Functions::Function1dForm>( 0 );
            if( multiplicity1->type( ) == FormType::branching1d ) return;

            averageEnergy( ).get<GIDI::Functions::Function1dForm>( 0 )->mapToXsAndAdd( a_offset, a_energies, a_productEnergies, -1.0 );
            averageMomentum( ).get<GIDI::Functions::Function1dForm>( 0 )->mapToXsAndAdd( a_offset, a_energies, a_productMomenta, -1.0 );
            multiplicity1->mapToXsAndAdd( a_offset, a_energies, a_productGains, 1.0 );
        } }
    else {
        m_outputChannel->mapContinuousEnergyProductData( a_settings, a_particleID, a_energies, a_offset, a_productEnergies, a_productMomenta, 
                a_productGains, a_ignoreIncompleteParticles );
    }
}

/* *********************************************************************************************************//**
 * Returns *true* if the product is complete and *false* otherwise. A product is complete if its multiplicity and/or distribution 
 * are other than *unspecified*.
 *
 * @return      *true* if the product is complete and *false* otherwise.
 ***********************************************************************************************************/

bool Product::isCompleteParticle( ) const {

    if( m_multiplicity.size( ) == 0 ) return( false );
    if( m_multiplicity.get<Form>( 0 )->type( ) == FormType::unspecified ) return( false );
    if( m_distribution.size( ) == 0 ) return( false );
    if( m_distribution.get<Form>( 0 )->type( ) == FormType::unspecified ) return( false );

    return( true );
}

/* *********************************************************************************************************//**
 * If the product has a distribution and its first distribution form is not **unspecified**, its PoPs id is added to *a_incompleteParticles*.
 *
 * @return      Returns *true* if the product has distribution data and *false* otherwise.
 ***********************************************************************************************************/

void Product::incompleteParticles( Transporting::Settings const &a_settings, std::set<std::string> &a_incompleteParticles ) const {

    if( m_outputChannel == nullptr ) {
        if( !isCompleteParticle( ) ) a_incompleteParticles.insert( particle( ).ID( ) ); }
    else {
        m_outputChannel->incompleteParticles( a_settings, a_incompleteParticles );
    }
}

/* *********************************************************************************************************//**
 * This methods calculates multi-group data for all needed components and adds each component's multi-group with label *a_heatedMultiGroupLabel*.
 *
 * @param   a_temperatureInfo                   [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param   a_heatedMultiGroupLabel             [in]    The label of the style for the multi-group data being added.
 * @param   a_multiGroupCalulationInformation   [in]    Store multi-group boundary and flux data used for multi-grouping.
 * @param   a_crossSectionXYs1d                 [in[    The cross section weight.
 ***********************************************************************************************************/

void Product::calculateMultiGroupData( ProtareSingle const *a_protare, Styles::TemperatureInfo const &a_temperatureInfo, 
                std::string const &a_heatedMultiGroupLabel, MultiGroupCalulationInformation const &a_multiGroupCalulationInformation, 
                Functions::XYs1d const &a_crossSectionXYs1d ) {

// FIXME, need to calculateMultiGroupData for the distribution.

    if( isCompleteParticle( ) ) {
        if( m_multiplicity.find( a_heatedMultiGroupLabel ) != m_multiplicity.end( ) ) {
            calculate1dMultiGroupDataInComponent( a_protare, a_heatedMultiGroupLabel, a_multiGroupCalulationInformation, m_multiplicity, a_crossSectionXYs1d );
            calculate1dMultiGroupDataInComponent( a_protare, a_heatedMultiGroupLabel, a_multiGroupCalulationInformation, m_averageEnergy, a_crossSectionXYs1d );
            calculate1dMultiGroupDataInComponent( a_protare, a_heatedMultiGroupLabel, a_multiGroupCalulationInformation, m_averageMomentum, a_crossSectionXYs1d );
        }
    }

    if( m_outputChannel != nullptr ) 
        m_outputChannel->calculateMultiGroupData( a_protare, a_temperatureInfo, a_heatedMultiGroupLabel, a_multiGroupCalulationInformation, a_crossSectionXYs1d );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Product::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string attributes;

    attributes += a_writeInfo.addAttribute( GIDI_pidChars, m_GNDS_particle.ID( ) );
    if( label( ) != "" ) attributes += a_writeInfo.addAttribute( GIDI_labelChars, label( ) );
    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );

    m_multiplicity.toXMLList( a_writeInfo, indent2 );
    m_distribution.toXMLList( a_writeInfo, indent2 );
    m_averageEnergy.toXMLList( a_writeInfo, indent2 );
    m_averageMomentum.toXMLList( a_writeInfo, indent2 );
    if( m_outputChannel != nullptr ) m_outputChannel->toXMLList( a_writeInfo, indent2 );

    a_writeInfo.addNodeEnder(  moniker( ) );
}

}
