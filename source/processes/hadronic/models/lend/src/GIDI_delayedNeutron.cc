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

/*! \class DelayedNeutron
 * This class represents a **GNDS** delayedNeutron.
*/

/* *********************************************************************************************************//**
 * Constructed from data in a <**outputChannel**> node.
 *
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The reaction HAPI::Node to be parsed and used to construct the reaction.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops            [in]    The *external* PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs    [in]    The *internal* PoPI::Database instance used to get particle indices and possibly other particle information.
 *                                  This is the <**PoPs**> node under the <**reactionSuite**> node.
 * @param a_parent          [in]    The parent GIDI::Suite.
 * @param a_styles          [in]    The <**styles**> node under the <**reactionSuite**> node.
 ***********************************************************************************************************/

DelayedNeutron::DelayedNeutron( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, 
                PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, Suite *a_parent, Styles::Suite const *a_styles ) :
        Form( a_node, a_setupInfo, FormType::delayedNeutron, a_parent ),
        m_delayedNeutronIndex( 0 ),
        m_rate( a_construction, GIDI_rateChars, GIDI_labelChars, a_node, a_setupInfo, a_pops, a_internalPoPs, parsePhysicalQuantitySuite, a_styles ),
        m_product( a_construction, a_node.child( GIDI_productChars ), a_setupInfo, a_pops, a_internalPoPs, nullptr, a_styles ) {

    m_rate.setAncestor( this );
    m_product.setAncestor( this );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

DelayedNeutron::~DelayedNeutron( ) {

}

/* *********************************************************************************************************//**
 * Used by GUPI::Ancestry to tranverse GNDS nodes. This method returns a pointer to a derived class' a_item member or nullptr if none exists.
 *
 * @param a_item    [in]    The name of the class member whose pointer is to be return.
 * @return                  The pointer to the class member or nullptr if class does not have a member named a_item.
 ***********************************************************************************************************/

GUPI::Ancestry *DelayedNeutron::findInAncestry3( std::string const &a_item ) {

    if( a_item == GIDI_rateChars ) return( &m_rate );
    if( a_item == GIDI_productChars ) return( &m_product );

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Used by GUPI::Ancestry to tranverse GNDS nodes. This method returns a pointer to a derived class' a_item member or nullptr if none exists.
 * 
 * @param a_item    [in]    The name of the class member whose pointer is to be return.
 * @return                  The pointer to the class member or nullptr if class does not have a member named a_item.
 ***********************************************************************************************************/

GUPI::Ancestry const *DelayedNeutron::findInAncestry3( std::string const &a_item ) const {

    if( a_item == GIDI_rateChars ) return( &m_rate );
    if( a_item == GIDI_productChars ) return( &m_product );

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

void DelayedNeutron::productIDs( std::set<std::string> &a_indices, Transporting::Particles const &a_particles, bool a_transportablesOnly ) const {

    m_product.productIDs( a_indices, a_particles, a_transportablesOnly );
}

/* *********************************************************************************************************//**
 * Returns the product multiplicity (e.g., 0, 1, 2, ...) or -1 if energy dependent or not an integer for particle with id *a_productID*.
 *
 * @param a_productID;              [in]    The id of the requested particle.
 *
 * @return                                  The multiplicity for the requested particle.
 ***********************************************************************************************************/

int DelayedNeutron::productMultiplicity( std::string const &a_productID ) const {

    return( m_product.productMultiplicity( a_productID ) );
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

int DelayedNeutron::maximumLegendreOrder( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    return( m_product.maximumLegendreOrder( a_smr, a_settings, a_temperatureInfo, a_productID ) );
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

Vector DelayedNeutron::multiGroupMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    return( m_product.multiGroupMultiplicity( a_smr, a_settings, a_temperatureInfo, a_productID ) );
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

Matrix DelayedNeutron::multiGroupProductMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, std::string const &a_productID, 
                std::size_t a_order ) const {

    return( m_product.multiGroupProductMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, a_productID, a_order ) );
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

Vector DelayedNeutron::multiGroupAverageEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    return( m_product.multiGroupAverageEnergy( a_smr, a_settings, a_temperatureInfo, a_productID ) );
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

Vector DelayedNeutron::multiGroupAverageMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    return( m_product.multiGroupAverageMomentum( a_smr, a_settings, a_temperatureInfo, a_productID ) );
}

/* *********************************************************************************************************//**
 * Added the product to *a_incompleteParticles* if the product's completeParticle returns *false*.
 *
 * @param a_settings                [in]    Specifies the requested label.
 * @param a_incompleteParticles     [out]   The list of particles whose **completeParticle** method returns *false*.
 ***********************************************************************************************************/

void DelayedNeutron::incompleteParticles( Transporting::Settings const &a_settings, std::set<std::string> &a_incompleteParticles ) const {

    m_product.incompleteParticles( a_settings, a_incompleteParticles );
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

void DelayedNeutron::continuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID, double a_energy,
                double &a_productEnergy, double &a_productMomentum, double &a_productGain, bool a_ignoreIncompleteParticles ) const {

    m_product.continuousEnergyProductData( a_settings, a_particleID, a_energy, a_productEnergy, a_productMomentum, a_productGain, 
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

void DelayedNeutron::mapContinuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID,
                std::vector<double> const &a_energies, std::size_t a_offset, std::vector<double> &a_productEnergies, std::vector<double> &a_productMomenta,
                std::vector<double> &a_productGains, bool a_ignoreIncompleteParticles ) const {
    
    m_product.mapContinuousEnergyProductData( a_settings, a_particleID, a_energies, a_offset, a_productEnergies, a_productMomenta,
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

void DelayedNeutron::calculateMultiGroupData( ProtareSingle const *a_protare, Styles::TemperatureInfo const &a_temperatureInfo, 
                std::string const &a_heatedMultiGroupLabel, MultiGroupCalulationInformation const &a_multiGroupCalulationInformation, 
                Functions::XYs1d const &a_crossSectionXYs1d ) {

    m_product.calculateMultiGroupData( a_protare, a_temperatureInfo, a_heatedMultiGroupLabel, a_multiGroupCalulationInformation, a_crossSectionXYs1d );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void DelayedNeutron::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {
    
    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    
    a_writeInfo.addNodeStarter( a_indent, moniker( ), a_writeInfo.addAttribute( GIDI_labelChars, label( ) ) );
    m_rate.toXMLList( a_writeInfo, indent2 );
    m_product.toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

}
