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

namespace Transporting {

/*! \class MC
 * Class to store user defined preferences for creating an MCGIDI::Protare instance.
 */

/* *********************************************************************************************************//**
 * Class to store user defined preferences for creating an MCGIDI::Protare instance.
 *
 * @param a_pops                        [in]    A PoPs Database instance used to get particle intids and possibly other particle information.
 * @param a_projectileID                [in]    The PoPs id for the projectile.
 * @param a_styles                      [in]    The styles child node of the GIDI::Protare.
 * @param a_label                       [in]    
 * @param a_delayedNeutrons             [in]    Sets whether delayed neutron data will be load or not if available.
 * @param a_energyDomainMax             [in]    The maximum projectile energy for which data should be loaded.
 ***********************************************************************************************************/

LUPI_HOST MC::MC( LUPI_maybeUnused PoPI::Database const &a_pops, std::string const &a_projectileID, GIDI::Styles::Suite const *a_styles, std::string const &a_label, 
                GIDI::Transporting::DelayedNeutrons a_delayedNeutrons, double a_energyDomainMax ) :
        GIDI::Transporting::Settings( a_projectileID, a_delayedNeutrons ),
        m_styles( a_styles ),
        m_label( a_label ),
        m_energyDomainMax( a_energyDomainMax ),
        m_ignoreENDF_MT5( false ),
        m_sampleNonTransportingParticles( false ),
        m_useSlowerContinuousEnergyConversion( false ),
        m_addExpectedValueData( true ),
        m_crossSectionLookupMode( LookupMode::Data1d::continuousEnergy ),
        m_other1dDataLookupMode( LookupMode::Data1d::continuousEnergy ),
        m_distributionLookupMode( LookupMode::Distribution::pdf_cdf ),
        m_upscatterModel( Sampling::Upscatter::Model::none ),
        m_URR_mode( URR_mode::none ),
        m_wantTerrellPromptNeutronDistribution( false ),
        m_wantRawTNSL_distributionSampling( true ),
        m_makePhotonEmissionProbabilitiesOne( false ),
        m_zeroNuclearLevelEnergyWidth( false )  {

}

/* *********************************************************************************************************//**
 * Class to store user defined preferences for creating an MCGIDI::Protare instance.
 *
 * @param a_pops                        [in]    A PoPs Database instance used to get particle intids and possibly other particle information.
 * @param a_protare                     [in]    GIDI::Protare whose information is used to fill *this*.
 * @param a_label                       [in]    
 * @param a_delayedNeutrons             [in]    Sets whether delayed neutron data will be load or not if available.
 * @param a_energyDomainMax             [in]    The maximum projectile energy for which data should be loaded.
 ***********************************************************************************************************/

LUPI_HOST MC::MC( LUPI_maybeUnused PoPI::Database const &a_pops, GIDI::Protare const &a_protare, std::string const &a_label, 
                GIDI::Transporting::DelayedNeutrons a_delayedNeutrons, double a_energyDomainMax ) :
        GIDI::Transporting::Settings( a_protare.projectile( ).ID( ), a_delayedNeutrons ),
        m_styles( &a_protare.styles( ) ),
        m_label( a_label ),
        m_energyDomainMax( a_energyDomainMax ),
        m_ignoreENDF_MT5( false ),
        m_sampleNonTransportingParticles( false ),
        m_useSlowerContinuousEnergyConversion( false ),
        m_addExpectedValueData( true ),
        m_crossSectionLookupMode( LookupMode::Data1d::continuousEnergy ),
        m_other1dDataLookupMode( LookupMode::Data1d::continuousEnergy ),
        m_distributionLookupMode( LookupMode::Distribution::pdf_cdf ),
        m_upscatterModel( Sampling::Upscatter::Model::none ),
        m_URR_mode( URR_mode::none ),
        m_wantTerrellPromptNeutronDistribution( false ),
        m_wantRawTNSL_distributionSampling( true ),
        m_makePhotonEmissionProbabilitiesOne( false ),
        m_zeroNuclearLevelEnergyWidth( false )  {

}

/* *********************************************************************************************************//**
 * Sets the *m_crossSectionLookupMode* member of *this* to *a_crossSectionLookupMode*.
 *
 * @param a_crossSectionLookupMode      [in]    The *LookupMode::Data1d* data mode.
 ***********************************************************************************************************/

LUPI_HOST void MC::setCrossSectionLookupMode( LookupMode::Data1d a_crossSectionLookupMode ) {

    if( ( a_crossSectionLookupMode != LookupMode::Data1d::continuousEnergy ) && 
        ( a_crossSectionLookupMode != LookupMode::Data1d::multiGroup ) ) {
        throw( "Invalided cross section mode request." );
    }
    m_crossSectionLookupMode = a_crossSectionLookupMode;
}

/* *********************************************************************************************************//**
 * Sets the *m_other1dDataLookupMode* member of *this* to *a_other1dDataLookupMode*.
 *
 * @param a_other1dDataLookupMode       [in]    The *LookupMode::Data1d* data mode.
 ***********************************************************************************************************/

LUPI_HOST void MC::setOther1dDataLookupMode( LookupMode::Data1d a_other1dDataLookupMode ) {

    if( a_other1dDataLookupMode != LookupMode::Data1d::continuousEnergy ) throw( "Invalided other mode request." );
    m_other1dDataLookupMode = a_other1dDataLookupMode;
}

/* *********************************************************************************************************//**
 * Sets the *m_distributionLookupMode* member of *this* to *a_distributionLookupMode*.
 *
 * @param a_distributionLookupMode      [in]    The *LookupMode::Data1d* data mode.
 ***********************************************************************************************************/

LUPI_HOST void MC::setDistributionLookupMode( LookupMode::Distribution a_distributionLookupMode ) {

    if( a_distributionLookupMode != LookupMode::Distribution::pdf_cdf ) throw( "Invalided distribution mode request." );
    m_distributionLookupMode = a_distributionLookupMode;
}

/* *********************************************************************************************************//**
 * This method sets the member *m_upscatterModelAGroupBoundaries* to *a_groupBoundaries*. It also checks that
 * the groups are in ascending order and executes a throw if they are not.
 *      
 * @param a_groupBoundaries     [in]    List of multi-group boundaries.
 ***********************************************************************************************************/
        
LUPI_HOST void MC::setUpscatterModelAGroupBoundaries( std::vector<double> const &a_groupBoundaries ) {

    double priorValue = 0.0;

    for( std::size_t index = 0; index < a_groupBoundaries.size( ); ++index ) {
        double value = a_groupBoundaries[index];

        if( index != 0 ) {
            if( value <= priorValue ) throw( "MC::setUpscatterModelAGroupBoundaries: group boundaries not in ascending order/" );
        }
        priorValue = value;
    }

    m_upscatterModelAGroupBoundaries = a_groupBoundaries;
}

}

}
