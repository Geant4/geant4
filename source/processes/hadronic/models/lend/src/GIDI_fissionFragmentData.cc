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

/*! \class FissionFragmentData
 * This class represents a **GNDS** fissionFragmentData.
*/

/* *********************************************************************************************************//**
 * Default constructor for FissionFragmentData.
 ***********************************************************************************************************/

FissionFragmentData::FissionFragmentData( ) :
        GUPI::Ancestry( GIDI_fissionFragmentDataChars ),
        m_delayedNeutrons( GIDI_delayedNeutronsChars, GIDI_labelChars ),
        m_fissionEnergyReleases( GIDI_fissionEnergyReleasesChars, GIDI_labelChars ) {

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
 ***********************************************************************************************************/

FissionFragmentData::FissionFragmentData( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		        PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, Styles::Suite const *a_styles ) :
        GUPI::Ancestry( a_node.name( ) ),
        m_delayedNeutrons( a_construction, GIDI_delayedNeutronsChars, GIDI_labelChars, a_node, a_setupInfo, a_pops, a_internalPoPs, 
                parseDelayedNeutronsSuite, a_styles ),
        m_fissionEnergyReleases( a_construction, GIDI_fissionEnergyReleasesChars, GIDI_labelChars, a_node, a_setupInfo, a_pops, a_internalPoPs, 
                parseFissionEnergyReleasesSuite, a_styles ) {

    m_delayedNeutrons.setAncestor( this );
    m_fissionEnergyReleases.setAncestor( this );

    for( std::size_t i1 = 0; i1 < m_delayedNeutrons.size( ); ++i1 ) {
        DelayedNeutron *delayedNeutron = m_delayedNeutrons.get<DelayedNeutron>( i1 );

        delayedNeutron->setDelayedNeutronIndex( i1 );
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

FissionFragmentData::~FissionFragmentData( ) {

}

/* *********************************************************************************************************//**
 * Returns **false* if *this* has delayed fission neutrons and they are not complete; otherwise, returns **true**.
 *
 * @return      bool
 ***********************************************************************************************************/

bool FissionFragmentData::isDelayedFissionNeutronComplete( ) const {

    for( std::size_t index = 0; index < m_delayedNeutrons.size( ); ++index ) {
        DelayedNeutron const &delayedNeutrons1 = *m_delayedNeutrons.get<DelayedNeutron>( index );

        Product const &product = delayedNeutrons1.product( );
        if( !product.isDelayedFissionNeutronComplete( true ) ) return( false );
    }

    return( true );
}

/* *********************************************************************************************************//**
 * Used by GUPI::Ancestry to tranverse GNDS nodes. This method returns a pointer to a derived class' a_item member or nullptr if none exists.
 *
 * @param a_item    [in]    The name of the class member whose pointer is to be return.
 * @return                  The pointer to the class member or nullptr if class does not have a member named a_item.
 ***********************************************************************************************************/

GUPI::Ancestry *FissionFragmentData::findInAncestry3( std::string const &a_item ) {

    if( a_item == GIDI_delayedNeutronsChars ) return( &m_delayedNeutrons );
    if( a_item == GIDI_fissionEnergyReleasesChars ) return( &m_fissionEnergyReleases );

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Used by GUPI::Ancestry to tranverse GNDS nodes. This method returns a pointer to a derived class' a_item member or nullptr if none exists.
 *
 * @param a_item    [in]    The name of the class member whose pointer is to be return.
 * @return                  The pointer to the class member or nullptr if class does not have a member named a_item.
 ***********************************************************************************************************/

GUPI::Ancestry const *FissionFragmentData::findInAncestry3( std::string const &a_item ) const {

    if( a_item == GIDI_delayedNeutronsChars ) return( &m_delayedNeutrons );
    if( a_item == GIDI_fissionEnergyReleasesChars ) return( &m_fissionEnergyReleases );

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Insert a std::set with the products id and any product in in its output channel.
 * If a_transportablesOnly is true, only transportable product indices are return.
 *
 * @param a_indices                 [out]   The unique list of product indices.
 * @param a_particles               [in]    The list of particles to be transported.
 * @param a_transportablesOnly      [in]    If true, only transportable product indices are added in the list.
 ***********************************************************************************************************/

void FissionFragmentData::productIDs( std::set<std::string> &a_indices, Transporting::Particles const &a_particles, bool a_transportablesOnly ) const {

    for( std::size_t index = 0; index < m_delayedNeutrons.size( ); ++index ) {
        DelayedNeutron const &delayedNeutrons1 = *m_delayedNeutrons.get<DelayedNeutron>( index );

        delayedNeutrons1.productIDs( a_indices, a_particles, a_transportablesOnly );
    }
}

/* *********************************************************************************************************//**
 * Returns the product multiplicity (e.g., 0, 1, 2, ...) or -1 if energy dependent or not an integer for particle with id *a_productID*.
 * 
 * @param a_productID;              [in]    The id of the requested particle.
 *
 * @return                                  The multiplicity for the requested particle.
 ***********************************************************************************************************/
 
int FissionFragmentData::productMultiplicity( std::string const &a_productID ) const {

    int total_multiplicity = 0;

    for( std::size_t index = 0; index < m_delayedNeutrons.size( ); ++index ) {
        DelayedNeutron const &delayedNeutrons1 = *m_delayedNeutrons.get<DelayedNeutron>( index );

        int multiplicity = delayedNeutrons1.productMultiplicity( a_productID );
        if( multiplicity < 0 ) return( -1 );
        total_multiplicity += multiplicity;
    }

    return( total_multiplicity );
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

int FissionFragmentData::maximumLegendreOrder( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    int _maximumLegendreOrder = -1;

    if( a_settings.delayedNeutrons( ) == Transporting::DelayedNeutrons::on ) {
        for( std::size_t index = 0; index < m_delayedNeutrons.size( ); ++index ) {
            DelayedNeutron const &delayedNeutrons1 = *m_delayedNeutrons.get<DelayedNeutron>( index );
            int r_maximumLegendreOrder = delayedNeutrons1.maximumLegendreOrder( a_smr, a_settings, a_temperatureInfo, a_productID );

            if( r_maximumLegendreOrder > _maximumLegendreOrder ) _maximumLegendreOrder = r_maximumLegendreOrder;
        }
    }

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

Vector FissionFragmentData::multiGroupMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    Vector vector( 0 );

    if( a_settings.delayedNeutrons( ) != Transporting::DelayedNeutrons::on ) return( vector );

    for( std::size_t index = 0; index < m_delayedNeutrons.size( ); ++index ) {
        DelayedNeutron const &delayedNeutrons1 = *m_delayedNeutrons.get<DelayedNeutron>( index );

        vector += delayedNeutrons1.multiGroupMultiplicity( a_smr, a_settings, a_temperatureInfo, a_productID );
    }

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

Vector FissionFragmentData::multiGroupQ( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, LUPI_maybeUnused bool a_final ) const {

    Vector vector( 0 );

    if( a_settings.delayedNeutrons( ) != Transporting::DelayedNeutrons::on ) return( vector );

    if( m_fissionEnergyReleases.size( ) == 0 ) return( vector );

    Functions::FissionEnergyRelease const *form = dynamic_cast<Functions::FissionEnergyRelease const *>( a_settings.form( 
            a_smr, m_fissionEnergyReleases, a_temperatureInfo, "Q-value" ) );

    if( form != nullptr ) vector += form->multiGroupQ( a_smr, a_settings, a_temperatureInfo );

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

Matrix FissionFragmentData::multiGroupProductMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, std::string const &a_productID, int a_order ) const {

    Matrix matrix( 0, 0 );

    if( a_settings.delayedNeutrons( ) != Transporting::DelayedNeutrons::on ) return( matrix );

    for( std::size_t index = 0; index < m_delayedNeutrons.size( ); ++index ) {
        DelayedNeutron const &delayedNeutrons1 = *m_delayedNeutrons.get<DelayedNeutron>( index );

        matrix += delayedNeutrons1.multiGroupProductMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, a_productID, a_order );
    }

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

Vector FissionFragmentData::multiGroupAverageEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    Vector vector( 0 );

    if( a_settings.delayedNeutrons( ) != Transporting::DelayedNeutrons::on ) return( vector );

    for( std::size_t index = 0; index < m_delayedNeutrons.size( ); ++index ) {
        DelayedNeutron const &delayedNeutrons1 = *m_delayedNeutrons.get<DelayedNeutron>( index );

        vector += delayedNeutrons1.multiGroupAverageEnergy( a_smr, a_settings, a_temperatureInfo, a_productID );
    }

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

Vector FissionFragmentData::multiGroupAverageMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    Vector vector( 0 );

    if( a_settings.delayedNeutrons( ) != Transporting::DelayedNeutrons::on ) return( vector );

    for( std::size_t index = 0; index < m_delayedNeutrons.size( ); ++index ) {
        DelayedNeutron const &delayedNeutrons1 = *m_delayedNeutrons.get<DelayedNeutron>( index );

        vector += delayedNeutrons1.multiGroupAverageMomentum( a_smr, a_settings, a_temperatureInfo, a_productID );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Appends a DelayedNeutronProduct instance for each delayed neutron in *m_delayedNeutrons*.
 *
 * @param       a_delayedNeutronProducts    [in]        The list to append the delayed neutrons to.
 ***********************************************************************************************************/

void FissionFragmentData::delayedNeutronProducts( DelayedNeutronProducts &a_delayedNeutronProducts ) const {

    for( std::size_t index = 0; index < m_delayedNeutrons.size( ); ++index ) {
        DelayedNeutron const &delayedNeutron = *m_delayedNeutrons.get<DelayedNeutron>( index );

        PhysicalQuantity rate = *delayedNeutron.rate( ).get<PhysicalQuantity>( 0 );
        a_delayedNeutronProducts.push_back( DelayedNeutronProduct( delayedNeutron.delayedNeutronIndex( ), rate, &delayedNeutron.product( ) ) );
    }
}

/* *********************************************************************************************************//**
 * Loops over the instances in ** m_delayedNeutrons** calling their incompleteParticles method.
 *
 * @param a_settings                [in]    Specifies the requested label.
 * @param a_incompleteParticles     [out]   The list of particles whose **completeParticle** method returns *false*.
 ***********************************************************************************************************/

void FissionFragmentData::incompleteParticles( Transporting::Settings const &a_settings, std::set<std::string> &a_incompleteParticles ) const {

    if( a_settings.delayedNeutrons( ) == Transporting::DelayedNeutrons::on ) {
        for( std::size_t index = 0; index < m_delayedNeutrons.size( ); ++index ) {
            DelayedNeutron const &delayedNeutron = *m_delayedNeutrons.get<DelayedNeutron>( index );

            delayedNeutron.incompleteParticles( a_settings, a_incompleteParticles );
        }
    }
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

void FissionFragmentData::continuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID, double a_energy, 
                double &a_productEnergy, double &a_productMomentum, double &a_productGain, bool a_ignoreIncompleteParticles ) const {

    if( a_settings.delayedNeutrons( ) != Transporting::DelayedNeutrons::on ) return;

    for( std::size_t index = 0; index < m_delayedNeutrons.size( ); ++index ) {
        DelayedNeutron const &delayedNeutrons1 = *m_delayedNeutrons.get<DelayedNeutron>( index );

        delayedNeutrons1.continuousEnergyProductData( a_settings, a_particleID, a_energy, a_productEnergy, a_productMomentum, a_productGain, 
                a_ignoreIncompleteParticles );
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

void FissionFragmentData::mapContinuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID,

                std::vector<double> const &a_energies, int a_offset, std::vector<double> &a_productEnergies, std::vector<double> &a_productMomenta,
                std::vector<double> &a_productGains, bool a_ignoreIncompleteParticles ) const {

    if( a_settings.delayedNeutrons( ) != Transporting::DelayedNeutrons::on ) return;

    for( std::size_t index = 0; index < m_delayedNeutrons.size( ); ++index ) {
        DelayedNeutron const &delayedNeutrons1 = *m_delayedNeutrons.get<DelayedNeutron>( index );

        delayedNeutrons1.mapContinuousEnergyProductData( a_settings, a_particleID, a_energies, a_offset, a_productEnergies, a_productMomenta,
                a_productGains, a_ignoreIncompleteParticles );
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

void FissionFragmentData::calculateMultiGroupData( ProtareSingle const *a_protare, Styles::TemperatureInfo const &a_temperatureInfo, 
                std::string const &a_heatedMultiGroupLabel, MultiGroupCalulationInformation const &a_multiGroupCalulationInformation, 
                Functions::XYs1d const &a_crossSectionXYs1d ) {

    for( std::size_t index = 0; index < m_delayedNeutrons.size( ); ++index ) {
        DelayedNeutron &delayedNeutrons1 = *m_delayedNeutrons.get<DelayedNeutron>( index );

        delayedNeutrons1.calculateMultiGroupData( a_protare, a_temperatureInfo, a_heatedMultiGroupLabel, a_multiGroupCalulationInformation, a_crossSectionXYs1d );
    }

    if( ( m_fissionEnergyReleases.size( ) > 0 ) && ( m_fissionEnergyReleases.find( a_heatedMultiGroupLabel ) !=  m_fissionEnergyReleases.end( ) ) ) {
        Functions::FissionEnergyRelease const *fissionEnergyRelease = m_fissionEnergyReleases.get<Functions::FissionEnergyRelease>( 0 );
        Functions::FissionEnergyRelease *multiGroupFissionEnergyRelease = m_fissionEnergyReleases.get<Functions::FissionEnergyRelease>( a_heatedMultiGroupLabel );

        calculate1dMultiGroupFissionEnergyRelease( a_multiGroupCalulationInformation, a_crossSectionXYs1d, 
                fissionEnergyRelease->promptProductKE( ), multiGroupFissionEnergyRelease->promptProductKE( ) );
        calculate1dMultiGroupFissionEnergyRelease( a_multiGroupCalulationInformation, a_crossSectionXYs1d, 
                fissionEnergyRelease->promptNeutronKE( ), multiGroupFissionEnergyRelease->promptNeutronKE( ) );
        calculate1dMultiGroupFissionEnergyRelease( a_multiGroupCalulationInformation, a_crossSectionXYs1d, 
                fissionEnergyRelease->delayedNeutronKE( ), multiGroupFissionEnergyRelease->delayedNeutronKE( ) );
        calculate1dMultiGroupFissionEnergyRelease( a_multiGroupCalulationInformation, a_crossSectionXYs1d, 
                fissionEnergyRelease->promptGammaEnergy( ), multiGroupFissionEnergyRelease->promptGammaEnergy( ) );
        calculate1dMultiGroupFissionEnergyRelease( a_multiGroupCalulationInformation, a_crossSectionXYs1d, 
                fissionEnergyRelease->delayedGammaEnergy( ), multiGroupFissionEnergyRelease->delayedGammaEnergy( ) );
        calculate1dMultiGroupFissionEnergyRelease( a_multiGroupCalulationInformation, a_crossSectionXYs1d, 
                fissionEnergyRelease->delayedBetaEnergy( ), multiGroupFissionEnergyRelease->delayedBetaEnergy( ) );
        calculate1dMultiGroupFissionEnergyRelease( a_multiGroupCalulationInformation, a_crossSectionXYs1d, 
                fissionEnergyRelease->neutrinoEnergy( ), multiGroupFissionEnergyRelease->neutrinoEnergy( ) );
        calculate1dMultiGroupFissionEnergyRelease( a_multiGroupCalulationInformation, a_crossSectionXYs1d, 
                fissionEnergyRelease->nonNeutrinoEnergy( ), multiGroupFissionEnergyRelease->nonNeutrinoEnergy( ) );
        calculate1dMultiGroupFissionEnergyRelease( a_multiGroupCalulationInformation, a_crossSectionXYs1d, 
                fissionEnergyRelease->totalEnergy( ), multiGroupFissionEnergyRelease->totalEnergy( ) );
    }
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void FissionFragmentData::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    if( moniker( ) == "" ) return;
    if( ( m_delayedNeutrons.size( ) == 0 ) && ( m_fissionEnergyReleases.size( ) == 0 ) ) return;
    a_writeInfo.addNodeStarter( a_indent, moniker( ), "" );
    m_delayedNeutrons.toXMLList( a_writeInfo, indent2 );
    m_fissionEnergyReleases.toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

}
