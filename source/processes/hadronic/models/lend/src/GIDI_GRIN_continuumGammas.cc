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

namespace GRIN {

static Form *parseInelasticIncidentEnergySuite( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, 
        SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, std::string const &a_name, 
        Styles::Suite const *a_styles );
static Form *parseCaptureLevelProbabilitySuite( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, 
        SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, std::string const &a_name, 
        Styles::Suite const *a_styles );

/*! \class GRIN_continuumGammas
 * Base class for the protare sub-classes.
 */

/* *********************************************************************************************************//**
 * Base Protare constructor.
 ***********************************************************************************************************/

GRIN_continuumGammas::GRIN_continuumGammas( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, 
                PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, LUPI_maybeUnused ProtareSingle const &a_protare, Styles::Suite const *a_styles ) :
        GUPI::Ancestry( "" ),
        m_captureNeutronSeparationEnergy( a_node.child( GIDI_captureNeutronSeparationEnergyChars ), a_setupInfo ),
        m_maximumCaptureIncidentEnergy( a_node.child( GIDI_maximumIncidentEnergyChars ), a_setupInfo ),
        m_pops( ),
        m_inelasticIncidentEnergies( a_construction, GIDI_inelasticIncidentEnergiesChars, GIDI_labelChars, a_node, a_setupInfo, a_pops, 
                a_internalPoPs, parseInelasticIncidentEnergySuite, a_styles ),
        m_captureLevelProbabilities( a_construction, GIDI_captureLevelProbabilitiesChars, GIDI_labelChars, a_node, a_setupInfo, a_pops, 
                a_internalPoPs, parseCaptureLevelProbabilitySuite, a_styles ),
        m_captureResidualIntid( -1 ),
        m_captureResidualIndex( -1 ),
        m_captureResidualMass( 0.0 ) {

    m_captureNeutronSeparationEnergy.setAncestor( this );
    m_maximumCaptureIncidentEnergy.setAncestor( this );
    m_inelasticIncidentEnergies.setAncestor( this );
    m_captureLevelProbabilities.setAncestor( this );

    m_pops.addDatabase( a_node.child( GIDI_PoPsChars ), true );

    PoPI::Nuclide const &target = a_pops.get<PoPI::Nuclide>( a_setupInfo.m_protare->target( ).pid( ) );
    std::string captureResidualId = target.isotope( )->chemicalElement( )->symbol( ) + std::to_string( target.A( ) + 1 );
    PoPI::Nuclide const &captureResidual = a_pops.get<PoPI::Nuclide>( captureResidualId );
    m_captureResidualId = captureResidualId;
    m_captureResidualIntid = captureResidual.intid( );
    m_captureResidualIndex = captureResidual.index( );
    m_captureResidualMass = captureResidual.massValue( "MeV/c**2" );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

GRIN_continuumGammas::~GRIN_continuumGammas( ) {
}

/* *********************************************************************************************************//**
 * Returns a pointer to the member whose moniker is *a_item*.
 *
 * @param a_item            [in]    The moniker of the member to return.
 * @return                          Returns the pointer to the member of nullptr if it does not exists.
 ***********************************************************************************************************/

GUPI::Ancestry *GRIN_continuumGammas::findInAncestry3( std::string const &a_item ) {

    if( a_item == GIDI_captureNeutronSeparationEnergyChars ) return( &m_captureNeutronSeparationEnergy );
    if( a_item == GIDI_maximumIncidentEnergyChars ) return( &m_maximumCaptureIncidentEnergy );
    if( a_item == GIDI_inelasticIncidentEnergiesChars ) return( &m_inelasticIncidentEnergies );
    if( a_item == GIDI_captureLevelProbabilitiesChars ) return( &m_captureLevelProbabilities );
// The following does not work as PoPI::Database does not yet inherent from GUPI::Ancestry. This needs to be fixed.
//    if( a_item == PoPI_PoPsChars ) return( &m_pops );

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Returns a pointer to the member whose moniker is *a_item*.
 *
 * @param a_item            [in]    The moniker of the member to return.
 * @return                          Returns the pointer to the member of nullptr if it does not exists.
 ***********************************************************************************************************/

GUPI::Ancestry const *GRIN_continuumGammas::findInAncestry3( std::string const &a_item ) const {
    
    if( a_item == GIDI_captureNeutronSeparationEnergyChars ) return( &m_captureNeutronSeparationEnergy );
    if( a_item == GIDI_maximumIncidentEnergyChars ) return( &m_maximumCaptureIncidentEnergy );
    if( a_item == GIDI_inelasticIncidentEnergiesChars ) return( &m_inelasticIncidentEnergies );
    if( a_item == GIDI_captureLevelProbabilitiesChars ) return( &m_captureLevelProbabilities );
// The following does not work as PoPI::Database does not yet inherent from GUPI::Ancestry. This needs to be fixed.
//    if( a_item == PoPI_PoPsChars ) return( &m_pops );

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Function that parses a <**inelasticIncidentEnergy**> node. Called from a Suite::parse instance.
 *
 * @param a_construction            [in]    Used to pass user options for parsing.
 * @param a_parent                  [in]    The parent GIDI::Suite that the returned Form will be added to.
 * @param a_node                    [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo               [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops                    [in]    A PoPs Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs            [in]    The *internal* PoPI::Database instance used to get particle indices and possibly other particle information.
 *                                          This is the <**PoPs**> node under the <**reactionSuite**> node.
 * @param a_name                    [in]    The moniker for the node to be parsed.
 * @param a_styles                  [in]    A pointer to the <**styles**> node.
 *
 * @return                                  The parsed and constructed GIDI::Form or nullptr if the node is not supported.
 ***********************************************************************************************************/

static Form *parseInelasticIncidentEnergySuite( Construction::Settings const &a_construction, LUPI_maybeUnused Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
        LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, std::string const &a_name, LUPI_maybeUnused Styles::Suite const *a_styles ) {

    Form *form = nullptr;

    if( a_name == GIDI_inelasticIncidentEnergyChars ) {
        form = new InelasticIncidentEnergy( a_construction, a_node, a_setupInfo ); }
    else {
        std::cout << "parseInelasticIncidentEnergySuite: Ignoring unsupported Form '" << a_name << "'." << std::endl;
    }

    return( form );
}

/* *********************************************************************************************************//**
 * Function that parses a <**captureLevelProbability**> node. Called from a Suite::parse instance.
 *
 * @param a_construction            [in]    Used to pass user options for parsing.
 * @param a_parent                  [in]    The parent GIDI::Suite that the returned Form will be added to.
 * @param a_node                    [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo               [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops                    [in]    A PoPs Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs            [in]    The *internal* PoPI::Database instance used to get particle indices and possibly other particle information.
 *                                          This is the <**PoPs**> node under the <**reactionSuite**> node.
 * @param a_name                    [in]    The moniker for the node to be parsed.
 * @param a_styles                  [in]    A pointer to the <**styles**> node.
 *
 * @return                                  The parsed and constructed GIDI::Form or nullptr if the node is not supported.
 ***********************************************************************************************************/

static Form *parseCaptureLevelProbabilitySuite( Construction::Settings const &a_construction, LUPI_maybeUnused Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
        LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, std::string const &a_name, LUPI_maybeUnused Styles::Suite const *a_styles ) {

    Form *form = nullptr;

    if( a_name == GIDI_captureLevelProbabilityChars ) {
        form = new CaptureLevelProbability( a_construction, a_node, a_setupInfo ); }
    else {
        std::cout << "parseCaptureLevelProbabilitySuite: Ignoring unsupported Form '" << a_name << "'." << std::endl;
    }

    return( form );
}
/* *********************************************************************************************************//**
 * Function that parses a <**inelasticIncidentEnergy**> node. Called from a Suite::parse instance.
 *
 * @param a_construction            [in]    Used to pass user options for parsing.
 * @param a_node                    [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo               [in]    Information create my the Protare constructor to help in parsing.
 ***********************************************************************************************************/

InelasticIncidentEnergy::InelasticIncidentEnergy( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo ) :
        Form( a_node, a_setupInfo, FormType::GRIN_inelasticIncidentEnergy ),
        m_energy( a_node.attribute_as_double( GIDI_energyChars ) ),
        m_unit( a_node.attribute_as_string( GIDI_unitChars ) ),
        m_table( a_construction, a_node.child( GIDI_tableChars ), a_setupInfo ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

InelasticIncidentEnergy::~InelasticIncidentEnergy( ) {

}

/* *********************************************************************************************************//**
 * Function that parses a <**captureLevelProbability**> node. Called from a Suite::parse instance.
 *
 * @param a_construction            [in]    Used to pass user options for parsing.
 * @param a_node                    [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo               [in]    Information create my the Protare constructor to help in parsing.
 ***********************************************************************************************************/

CaptureLevelProbability::CaptureLevelProbability( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo ) :
        Form( a_node, a_setupInfo, FormType::GRIN_captureLevelProbability ),
        m_probabilty( a_node.attribute_as_double( GIDI_probabilityChars ) ), 
        m_spin( a_node.attribute_as_double( PoPI_spinChars ) ),
        m_spinUnit( a_node.attribute_as_string( GIDI_spinUnitChars ) ),
        m_parity( a_node.attribute_as_int( PoPI_parityChars ) ),
        m_capturePrimaryToContinua( a_node.attribute_as_string( GIDI_capturePrimaryToContinuaChars ) ),
        m_table( a_construction, a_node.child( GIDI_tableChars ), a_setupInfo ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

CaptureLevelProbability::~CaptureLevelProbability( ) {

}

}               // End namespace GRIN.

}               // End namespace GIDI.
