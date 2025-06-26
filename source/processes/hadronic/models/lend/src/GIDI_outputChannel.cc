/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "GIDI.hpp"
#include <HAPI.hpp>

namespace GIDI {

/*! \class OutputChannel
 * This class represents a **GNDS** outputChannel.
*/

OutputChannel::OutputChannel( bool a_twoBody, bool a_fissions, std::string a_process ) :
        GUPI::Ancestry( GIDI_outputChannelChars ),
        m_twoBody( a_twoBody ),
        m_fissions( a_fissions ),
        m_process( a_process ),
        m_Q( GIDI_QChars, GIDI_labelChars ),
        m_products( GIDI_productsChars, GIDI_labelChars ),
        m_fissionFragmentData( ),
        m_fissionResiduals( Construction::FissionResiduals::none ) {

    m_Q.setAncestor( this );
    m_products.setAncestor( this );
    m_fissionFragmentData.setAncestor( this );
}

/* *********************************************************************************************************//**
 * Constructed from data in a <**outputChannel**> node.
 *
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The reaction HAPI::Node to be parsed and used to construct the reaction.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops            [in]    The *external* PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs    [in]    The *internal* PoPI::Database instance used to get particle indices and possibly other particle information.
 *                                  This is the <**PoPs**> node under the <**reactionSuite**> node.
 * @param a_styles          [in]    The <**styles**> node under the <**reactionSuite**> node.
 * @param a_isFission       [in]    Boolean indicating if output channel is a fission channel (true) or not (false).
 ***********************************************************************************************************/

OutputChannel::OutputChannel( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, 
                PoPI::Database const &a_internalPoPs, Styles::Suite const *a_styles, bool a_isFission, LUPI_maybeUnused bool a_addFissionResiduals ) :
        GUPI::Ancestry( a_node.name( ) ),
        m_twoBody( std::string( a_node.attribute_as_string( GIDI_genreChars ) ) == GIDI_twoBodyChars ),
        m_fissions( a_isFission ),
        m_process( std::string( a_node.attribute_as_string( GIDI_processChars ) ) ),
        m_Q( a_construction, GIDI_QChars, GIDI_labelChars, a_node, a_setupInfo, a_pops, a_internalPoPs, parseQSuite, a_styles ),
        m_products( a_construction, GIDI_productsChars, GIDI_labelChars, a_node, a_setupInfo, a_pops, a_internalPoPs, parseProductSuite, a_styles ),
        m_fissionFragmentData( a_construction, a_node.child( GIDI_fissionFragmentDataChars ), a_setupInfo, a_pops, a_internalPoPs, a_styles ),
        m_fissionResiduals( Construction::FissionResiduals::none ) {

    if( a_isFission ) m_fissionResiduals = a_construction.fissionResiduals( );

    m_Q.setAncestor( this );
    m_products.setAncestor( this );
    m_fissionFragmentData.setAncestor( this );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

OutputChannel::~OutputChannel( ) {

}

/* *********************************************************************************************************//**
 * Returns the maximum product depth for this output channel.
 *
 * @return The maximum product depth.
 ***********************************************************************************************************/

int OutputChannel::depth( ) const {

    int _depth = 0;
    std::size_t size = m_products.size( );

    for( std::size_t index = 0; index < size; ++index ) {
        Product const &product = *m_products.get<Product>( index );

        int productDepth = product.depth( );
        if( productDepth > _depth ) _depth = productDepth;
    }
    return( _depth + 1 );
}

/* *********************************************************************************************************//**
 * Returns **true** if all outgoing particles (i.e., products) are specifed in *a_particles*. That is, the user
 * will be tracking all products of *this* reaction.
 *
 * @param a_particles           [in]    The list of particles to be transported.
 *
 * @return                              bool.
 ***********************************************************************************************************/

bool OutputChannel::areAllProductsTracked( Transporting::Particles const &a_particles ) const {
// Does not check m_fissionFragmentData as its will only have neutrons which should already be in m_products at least for now.

    if( isFission( ) ) return( false );

    for( auto iter = m_products.begin( ); iter != m_products.end( ); ++iter ) {
        Product *product = static_cast<Product *>( *iter );

        if( !product->areAllProductsTracked( a_particles ) ) return( false );
    }

    return( true );
}

/* *********************************************************************************************************//**
 * Only for internal use. Called by a ProtareTNSL instance to zero the lower energy multi-group data covered by the TNSL ProtareSingle.
 *
 * @param a_maximumTNSL_MultiGroupIndex     [in]    A map that contains labels for heated multi-group data and the last valid group boundary
 *                                                  for the TNSL data for that boundary.
 ***********************************************************************************************************/

void OutputChannel::modifiedMultiGroupElasticForTNSL( std::map<std::string,std::size_t> a_maximumTNSL_MultiGroupIndex ) {

        // No need to fix m_Q as it is all 0.0's for elastic scattering.
    for( auto iter = m_products.begin( ); iter != m_products.end( ); ++iter ) {
        Product *product = static_cast<Product *>( *iter );

        product->modifiedMultiGroupElasticForTNSL( a_maximumTNSL_MultiGroupIndex );
    }
}

/* *********************************************************************************************************//**
 * Used by GUPI::Ancestry to tranverse GNDS nodes. This method returns a pointer to a derived class' a_item member or nullptr if none exists.
 *
 * @param a_item    [in]    The name of the class member whose pointer is to be return.
 * @return                  The pointer to the class member or nullptr if class does not have a member named a_item.
 ***********************************************************************************************************/

GUPI::Ancestry *OutputChannel::findInAncestry3( std::string const &a_item ) {

    if( a_item == GIDI_QChars ) return( &m_Q );
    if( a_item == GIDI_productsChars ) return( &m_products );
    if( a_item == GIDI_fissionFragmentDataChars ) return( &m_fissionFragmentData );

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Used by GUPI::Ancestry to tranverse GNDS nodes. This method returns a pointer to a derived class' a_item member or nullptr if none exists.
 *
 * @param a_item    [in]    The name of the class member whose pointer is to be return.
 * @return                  The pointer to the class member or nullptr if class does not have a member named a_item.
 ***********************************************************************************************************/

GUPI::Ancestry const *OutputChannel::findInAncestry3( std::string const &a_item ) const {

    if( a_item == GIDI_QChars ) return( &m_Q );
    if( a_item == GIDI_productsChars ) return( &m_products );
    if( a_item == GIDI_fissionFragmentDataChars ) return( &m_fissionFragmentData );

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Returns true if the product has an output channel and its output channel hasFission returns true, and false otherwise.
 *
 * @return  true if at least one output channel is a fission channel.
 ***********************************************************************************************************/

bool OutputChannel::hasFission( ) const {

    if( m_fissions ) return( true );

    std::size_t size = m_products.size( );
    for( std::size_t index = 0; index < size; ++index ) {
        Product const &product = *m_products.get<Product>( index );

        if( product.hasFission( ) ) return( true );
    }
    return( false );
}

/* *********************************************************************************************************//**
 * Returns **false* if outputChannel has delayed fission neutrons and they are not complete; otherwise, returns **true**.
 *
 * @return      bool
 ***********************************************************************************************************/

bool OutputChannel::isDelayedFissionNeutronComplete( ) const {

    if( !m_fissionFragmentData.isDelayedFissionNeutronComplete( ) ) return( false );

    std::size_t size = m_products.size( );
    for( std::size_t index = 0; index < size; ++index ) {
        Product const &product = *m_products.get<Product>( index );

        if( !product.isDelayedFissionNeutronComplete( false ) ) return( false );
    }

    return( true );
}

/* *********************************************************************************************************//**
 * Insert a std::set with the products id and any product in in its output channel.
 * If a_transportablesOnly is true, only transportable product indices are return.
 *
 * @param a_indices                 [out]   The unique list of product indices.
 * @param a_particles               [in]    The list of particles to be transported.
 * @param a_transportablesOnly      [in]    If true, only transportable product indices are added in the list.
 ***********************************************************************************************************/

void OutputChannel::productIDs( std::set<std::string> &a_indices, Transporting::Particles const &a_particles, bool a_transportablesOnly ) const {

    std::size_t size = m_products.size( );

    for( std::size_t index = 0; index < size; ++index ) {
        Product const &product = *m_products.get<Product>( index );

        product.productIDs( a_indices, a_particles, a_transportablesOnly );
    }

    m_fissionFragmentData.productIDs( a_indices, a_particles, a_transportablesOnly );

    if( !a_transportablesOnly && isFission( ) ) {
        if(      m_fissionResiduals == Construction::FissionResiduals::ENDL99120 ) {
            a_indices.insert( PoPI::IDs::FissionProductENDL99120 ); }
        else if( m_fissionResiduals == Construction::FissionResiduals::ENDL99125 ) {
            a_indices.insert( PoPI::IDs::FissionProductENDL99125 );
        }
    }
}

/* *********************************************************************************************************//**
 * Returns the product multiplicity (e.g., 0, 1, 2, ...) or -1 if energy dependent or not an integer for particle with id *a_productID*.
 *
 * @param a_productID;              [in]    The id of the requested particle.
 *
 * @return                                  The multiplicity for the requested particle.
 ***********************************************************************************************************/

int OutputChannel::productMultiplicity( std::string const &a_productID ) const {

    int total_multiplicity = 0;
    std::size_t size = m_products.size( );

    if( isFission( ) ) {
        if(      ( a_productID == PoPI::IDs::FissionProductENDL99120 ) && ( m_fissionResiduals == Construction::FissionResiduals::ENDL99120 ) ) {
            return( 2 ); }
        else if( ( a_productID == PoPI::IDs::FissionProductENDL99125 ) && ( m_fissionResiduals == Construction::FissionResiduals::ENDL99125 ) ) {
            return( 2 );
        }
    }

    for( std::size_t index = 0; index < size; ++index ) {
        Product const &product = *m_products.get<Product>( index );
        int multiplicity = product.productMultiplicity( a_productID );

        if( multiplicity < 0 ) return( -1 );
        total_multiplicity += multiplicity;
    }

    int multiplicity = m_fissionFragmentData.productMultiplicity( a_productID );
    if( multiplicity < 0 ) return( -1 );

    return( total_multiplicity + multiplicity );
}

/* *********************************************************************************************************//**
 * Determines the maximum Legredre order present in the multi-group transfer matrix for the specified products of this output channel.
 *
 * @param a_smr             [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings        [in]    Specifies the requested label.
 * @param a_temperatureInfo [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID       [in]    Particle id of the requested product.
 *
 * @return                          The maximum Legredre order. If no transfer matrix data are present for the requested product, -1 is returned.
 ***********************************************************************************************************/

int OutputChannel::maximumLegendreOrder( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    std::size_t size = m_products.size( );
    int _maximumLegendreOrder = -1;

    for( std::size_t index = 0; index < size; ++index ) {
        Product const &product = *m_products.get<Product>( index );
        int r_maximumLegendreOrder = product.maximumLegendreOrder( a_smr, a_settings, a_temperatureInfo, a_productID );

        if( r_maximumLegendreOrder > _maximumLegendreOrder ) _maximumLegendreOrder = r_maximumLegendreOrder;
    }

    int r_maximumLegendreOrder = m_fissionFragmentData.maximumLegendreOrder( a_smr, a_settings, a_temperatureInfo, a_productID );
    if( r_maximumLegendreOrder > _maximumLegendreOrder ) _maximumLegendreOrder = r_maximumLegendreOrder;

    return( _maximumLegendreOrder );
}

/* *********************************************************************************************************//**
 * Returns the sum of the multi-group multiplicity for the requested label for the request product of this output channel. 
 * This is a cross section weighted multiplicity.
 *
 * @param a_smr             [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings        [in]    Specifies the requested label.
 * @param a_temperatureInfo [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID       [in]    Particle id for the requested product.
 *
 * @return                          The requested multi-group multiplicity as a GIDI::Vector.
 ***********************************************************************************************************/

Vector OutputChannel::multiGroupMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    Vector vector( 0 );

    for( std::size_t index = 0; index < m_products.size( ); ++index ) {
        Product const &product = *m_products.get<Product>( index );

        vector += product.multiGroupMultiplicity( a_smr, a_settings, a_temperatureInfo, a_productID );
    }

    vector += m_fissionFragmentData.multiGroupMultiplicity( a_smr, a_settings, a_temperatureInfo, a_productID );

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the sum of the multi-group, Q for the requested label for the this output channel. This is a cross section weighted Q.
 * If a_final is false, only the Q for the output channels directly under each reaction is summed. Otherwise, the Q for all output channels
 * summed, including output channels for each products.
 *
 * @param a_smr             [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings        [in]    Specifies the requested label.
 * @param a_temperatureInfo [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_final           [in]    If true, the Q is calculated for all output channels, including those for products.
 *
 * @return                          The requested multi-group Q as a GIDI::Vector.
 ***********************************************************************************************************/

Vector OutputChannel::multiGroupQ( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, bool a_final ) const {

    Vector vector( 0 );

    Functions::Gridded1d const *form = dynamic_cast<Functions::Gridded1d const*>( a_settings.form( a_smr, m_Q, a_temperatureInfo, "Q-value" ) );

    if( form != nullptr ) vector += form->data( );

    if( a_final ) {
        for( std::size_t index = 0; index < m_products.size( ); ++index ) {
            Product const &product1 = *m_products.get<Product>( index );

            vector += product1.multiGroupQ( a_smr, a_settings, a_temperatureInfo, a_final );
        }
    }

    vector += m_fissionFragmentData.multiGroupQ( a_smr, a_settings, a_temperatureInfo, a_final );

    return( vector );
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

Matrix OutputChannel::multiGroupProductMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, std::string const &a_productID, int a_order ) const {

    Matrix matrix( 0, 0 );

    for( std::size_t index = 0; index < m_products.size( ); ++index ) {
        Product const &product = *m_products.get<Product>( index );

        matrix += product.multiGroupProductMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, a_productID, a_order );
    }

    matrix += m_fissionFragmentData.multiGroupProductMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, a_productID, a_order ); 

    return( matrix );
}

/* *********************************************************************************************************//**
 * Returns the sum of the multi-group, average energy for the requested label for the requested product. This is a cross section weighted average energy.
 *
 * @param a_smr             [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings        [in]    Specifies the requested label.
 * @param a_temperatureInfo [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID       [in]    Particle id for the requested product.
 *
 * @return                          The requested multi-group average energy as a GIDI::Vector.
 ***********************************************************************************************************/

Vector OutputChannel::multiGroupAverageEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    Vector vector( 0 );

    for( std::size_t index = 0; index < m_products.size( ); ++index ) {
        Product const &product = *m_products.get<Product>( index );

        vector += product.multiGroupAverageEnergy( a_smr, a_settings, a_temperatureInfo, a_productID );
    }

    vector += m_fissionFragmentData.multiGroupAverageEnergy( a_smr, a_settings, a_temperatureInfo, a_productID ); 

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the sum of the multi-group, average momentum for the requested label for the requested product. This is a cross section weighted average momentum.
 *
 * @param a_smr             [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings        [in]    Specifies the requested label.
 * @param a_temperatureInfo [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID       [in]    Particle id for the requested product.
 *
 * @return                          The requested multi-group average momentum as a GIDI::Vector.
 ***********************************************************************************************************/

Vector OutputChannel::multiGroupAverageMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    Vector vector( 0 );

    for( std::size_t index = 0; index < m_products.size( ); ++index ) {
        Product const &product = *m_products.get<Product>( index );

        vector += product.multiGroupAverageMomentum( a_smr, a_settings, a_temperatureInfo, a_productID );
    }

    vector += m_fissionFragmentData.multiGroupAverageMomentum( a_smr, a_settings, a_temperatureInfo, a_productID ); 

    return( vector );
}

/* *********************************************************************************************************//**
 * Loops over the instances in **m_products** calling their **incompleteParticles** methods and calls the :**incompleteParticles** method
 * for the **m_fissionFragmentData** member.
 *
 * @param a_settings                [in]    Specifies the requested label.
 * @param a_incompleteParticles     [out]   The list of particles whose **completeParticle** method returns *false*.
 ***********************************************************************************************************/

void OutputChannel::incompleteParticles( Transporting::Settings const &a_settings, std::set<std::string> &a_incompleteParticles ) const {

    for( std::size_t index = 0; index < m_products.size( ); ++index ) {
        Product const &product = *m_products.get<Product>( index );

        product.incompleteParticles( a_settings, a_incompleteParticles );
    }

    m_fissionFragmentData.incompleteParticles( a_settings, a_incompleteParticles );
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

void OutputChannel::continuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID, double a_energy, 
                double &a_productEnergy, double &a_productMomentum, double &a_productGain, bool a_ignoreIncompleteParticles ) const {

    for( std::size_t index = 0; index < m_products.size( ); ++index ) {
        Product const &product = *m_products.get<Product>( index );

        product.continuousEnergyProductData( a_settings, a_particleID, a_energy, a_productEnergy, a_productMomentum, a_productGain, a_ignoreIncompleteParticles );
    }

    m_fissionFragmentData.continuousEnergyProductData( a_settings, a_particleID, a_energy, a_productEnergy, a_productMomentum, a_productGain, 
            a_ignoreIncompleteParticles );
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
 
void OutputChannel::mapContinuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID, 
                std::vector<double> const &a_energies, int a_offset, std::vector<double> &a_productEnergies, std::vector<double> &a_productMomenta, 
                std::vector<double> &a_productGains, bool a_ignoreIncompleteParticles ) const {

    for( std::size_t index = 0; index < m_products.size( ); ++index ) {
        Product const &product = *m_products.get<Product>( index );

        product.mapContinuousEnergyProductData( a_settings, a_particleID, a_energies, a_offset, a_productEnergies, a_productMomenta, 
                a_productGains, a_ignoreIncompleteParticles );
    }

    m_fissionFragmentData.mapContinuousEnergyProductData( a_settings, a_particleID, a_energies, a_offset, a_productEnergies, a_productMomenta,
                a_productGains, a_ignoreIncompleteParticles );
}

/* *********************************************************************************************************//**
 * This methods calculates multi-group data for all needed components and adds each component's multi-group with label *a_heatedMultiGroupLabel*.
 *
 * @param   a_temperatureInfo                   [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param   a_heatedMultiGroupLabel             [in]    The label of the style for the multi-group data being added.
 * @param   a_multiGroupCalulationInformation   [in]    Store multi-group boundary and flux data used for multi-grouping.
 * @param   a_crossSectionXYs1d                 [in[    The cross section weight.
 ***********************************************************************************************************/

void OutputChannel::calculateMultiGroupData( ProtareSingle const *a_protare, Styles::TemperatureInfo const &a_temperatureInfo, 
                std::string const &a_heatedMultiGroupLabel, MultiGroupCalulationInformation const &a_multiGroupCalulationInformation, 
                Functions::XYs1d const &a_crossSectionXYs1d ) {

    calculate1dMultiGroupDataInComponent( a_protare, a_heatedMultiGroupLabel, a_multiGroupCalulationInformation, m_Q, a_crossSectionXYs1d );

    for( std::size_t index = 0; index < m_products.size( ); ++index ) {
        Product &product = *m_products.get<Product>( index );

        product.calculateMultiGroupData( a_protare, a_temperatureInfo, a_heatedMultiGroupLabel, a_multiGroupCalulationInformation, a_crossSectionXYs1d );
    }

    m_fissionFragmentData.calculateMultiGroupData( a_protare, a_temperatureInfo, a_heatedMultiGroupLabel, a_multiGroupCalulationInformation, a_crossSectionXYs1d );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *  
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/
 
void OutputChannel::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {
    
    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string attributes;

    if( m_twoBody ) {
        attributes = a_writeInfo.addAttribute( GIDI_genreChars, GIDI_twoBodyChars ); }
    else {
        attributes = a_writeInfo.addAttribute( GIDI_genreChars, GIDI_NBodyChars );
    }

    if( m_process != "" ) attributes += a_writeInfo.addAttribute( GIDI_processChars, m_process );

    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );

    m_Q.toXMLList( a_writeInfo, indent2 ); 
    m_products.toXMLList( a_writeInfo, indent2 );
    m_fissionFragmentData.toXMLList( a_writeInfo, indent2 );

    a_writeInfo.addNodeEnder( moniker( ) );
}

}
