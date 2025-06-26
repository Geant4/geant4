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

namespace Functions {

/*! \class FissionEnergyRelease
 * Class to store GNDS <**fissionEnergyRelease**> node.
 */

/* *********************************************************************************************************//**
 * @param a_construction        [in]    Used to pass user options for parsing.
 * @param a_node                [in]    The HAPI::Node to be parsed to construct the instance.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

FissionEnergyRelease::FissionEnergyRelease( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Function1dForm( a_construction, a_node, a_setupInfo, FormType::fissionEnergyRelease1d, a_parent ) {

    m_promptProductKE = data1dParse( a_construction, a_node.child( GIDI_promptProductKEChars ).first_child( ), a_setupInfo, nullptr );
    m_promptNeutronKE = data1dParse( a_construction, a_node.child( GIDI_promptNeutronKEChars ).first_child( ), a_setupInfo, nullptr );
    m_delayedNeutronKE = data1dParse( a_construction, a_node.child( GIDI_delayedNeutronKEChars ).first_child( ), a_setupInfo, nullptr );
    m_promptGammaEnergy = data1dParse( a_construction, a_node.child( GIDI_promptGammaEnergyChars ).first_child( ), a_setupInfo, nullptr );
    m_delayedGammaEnergy = data1dParse( a_construction, a_node.child( GIDI_delayedGammaEnergyChars ).first_child( ), a_setupInfo, nullptr );
    m_delayedBetaEnergy = data1dParse( a_construction, a_node.child( GIDI_delayedBetaEnergyChars ).first_child( ), a_setupInfo, nullptr );
    m_neutrinoEnergy = data1dParse( a_construction, a_node.child( GIDI_neutrinoEnergyChars ).first_child( ), a_setupInfo, nullptr );
    m_nonNeutrinoEnergy = data1dParse( a_construction, a_node.child( GIDI_nonNeutrinoEnergyChars ).first_child( ), a_setupInfo, nullptr );
    m_totalEnergy = data1dParse( a_construction, a_node.child( GIDI_totalEnergyChars ).first_child( ), a_setupInfo, nullptr );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

FissionEnergyRelease::~FissionEnergyRelease( ) {

    if( m_promptProductKE != nullptr ) delete m_promptProductKE;
    if( m_promptNeutronKE != nullptr ) delete m_promptNeutronKE;
    if( m_delayedNeutronKE != nullptr ) delete m_delayedNeutronKE;
    if( m_promptGammaEnergy != nullptr ) delete m_promptGammaEnergy;
    if( m_delayedGammaEnergy != nullptr ) delete m_delayedGammaEnergy;
    if( m_delayedBetaEnergy != nullptr ) delete m_delayedBetaEnergy;
    if( m_neutrinoEnergy != nullptr ) delete m_neutrinoEnergy;
    if( m_nonNeutrinoEnergy != nullptr ) delete m_nonNeutrinoEnergy;
    if( m_totalEnergy != nullptr ) delete m_totalEnergy;
}

/* *********************************************************************************************************//**
 * Returns the multi-group Q-value.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 *
 * @return                              Multi-group Q-value.
 ***********************************************************************************************************/

Vector FissionEnergyRelease::multiGroupQ( LUPI_maybeUnused LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, LUPI_maybeUnused Styles::TemperatureInfo const &a_temperatureInfo ) const {

    Vector vector( 0 );

    if( a_settings.delayedNeutrons( ) == Transporting::DelayedNeutrons::on ) {
        Gridded1d const *gridded1d = dynamic_cast<Gridded1d const *>( m_delayedNeutronKE );

        vector += gridded1d->data( );
        gridded1d = dynamic_cast<Gridded1d const *>( m_delayedGammaEnergy );
        vector += gridded1d->data( );
        gridded1d = dynamic_cast<Gridded1d const *>( m_delayedBetaEnergy );
        vector += gridded1d->data( );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void FissionEnergyRelease::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    a_writeInfo.addNodeStarter( a_indent, moniker( ), a_writeInfo.addAttribute( GIDI_labelChars, label( ) ) );

    energyReleaseToXMLList( a_writeInfo, GIDI_promptProductKEChars, indent2, m_promptProductKE );
    energyReleaseToXMLList( a_writeInfo, GIDI_promptNeutronKEChars, indent2, m_promptNeutronKE );
    energyReleaseToXMLList( a_writeInfo, GIDI_delayedNeutronKEChars, indent2, m_delayedNeutronKE );
    energyReleaseToXMLList( a_writeInfo, GIDI_promptGammaEnergyChars, indent2, m_promptGammaEnergy );
    energyReleaseToXMLList( a_writeInfo, GIDI_delayedGammaEnergyChars, indent2, m_delayedGammaEnergy );
    energyReleaseToXMLList( a_writeInfo, GIDI_delayedBetaEnergyChars, indent2, m_delayedBetaEnergy );
    energyReleaseToXMLList( a_writeInfo, GIDI_neutrinoEnergyChars, indent2, m_neutrinoEnergy );
    energyReleaseToXMLList( a_writeInfo, GIDI_nonNeutrinoEnergyChars, indent2, m_nonNeutrinoEnergy );
    energyReleaseToXMLList( a_writeInfo, GIDI_totalEnergyChars, indent2, m_totalEnergy );

    a_writeInfo.addNodeEnder( moniker( ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_moniker           [in]        The moniker for the energy type.
 * @param       a_indent            [in]        The amount to indent *this* node.
 * @param       a_function          [in]        The component of the energy released in fission.
 ***********************************************************************************************************/

void FissionEnergyRelease::energyReleaseToXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_moniker, std::string const &a_indent, Function1dForm *a_function ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    if( a_function == nullptr ) return;

    a_writeInfo.addNodeStarter( a_indent, a_moniker, "" );
    a_function->toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( a_moniker );
}

}               // End namespace Functions.

}               // End namespace GIDI.
