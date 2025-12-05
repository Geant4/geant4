/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <math.h>

#include "GIDI.hpp"
#include <HAPI.hpp>

namespace GIDI {

/* *********************************************************************************************************//**
 * Function that parses a <**style**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::Form or nullptr if the node is not supported.
 ***********************************************************************************************************/

Form *parseExternalFilesSuite( LUPI_maybeUnused Construction::Settings const &a_construction, GIDI::Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, std::string const &a_name, LUPI_maybeUnused Styles::Suite const *a_styles ) {

    Form *form = nullptr;

    if( a_name == GIDI_externalFileChars ) {
        form = new ExternalFile( a_node, a_setupInfo, a_parent ); }
    else {
        std::cout << "parseExternalFilesSuite: Ignoring unsupported externalFile = '" << a_name << "'." << std::endl;
    }

    return( form );
}

/* *********************************************************************************************************//**
 * Function that parses a <**style**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::Form or nullptr if the node is not supported.
 ***********************************************************************************************************/

Form *parseStylesSuite( Construction::Settings const &a_construction, GIDI::Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, std::string const &a_name, LUPI_maybeUnused Styles::Suite const *a_styles ) {

    Form *form = nullptr;

//  Styles not parsed are angularDistributionReconstructed.

    if(      a_name == GIDI_evaluatedStyleChars ) {
        form = new Styles::Evaluated( a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_crossSectionReconstructedStyleChars ) {
        form = new Styles::CrossSectionReconstructed( a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_CoulombPlusNuclearElasticMuCutoffStyleChars ) {
        form = new Styles::CoulombPlusNuclearElasticMuCutoff( a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_realizationChars ) {
        form = new Styles::Realization( a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_averageProductDataStyleChars ) {
        form = new Styles::AverageProductData( a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_MonteCarlo_cdfStyleChars ) {
        form = new Styles::MonteCarlo_cdf( a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_multiGroupStyleChars ) {
        form = new Styles::MultiGroup( a_construction, a_node, a_setupInfo, a_pops, a_internalPoPs, a_parent );
        if( a_setupInfo.m_multiGroup == nullptr ) {
            a_setupInfo.m_multiGroup = static_cast<Styles::MultiGroup *>( form ); }
        else {
            std::cout << "Multiple multiGroup style instances found which is not supported. Ignoring all but first instance." << std::endl;
        } }
    else if( a_name == GIDI_heatedStyleChars ) {
        form = new Styles::Heated( a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_heatedMultiGroupStyleChars ) {
        form = new Styles::HeatedMultiGroup( a_construction, a_node, a_setupInfo, a_pops, a_parent ); }
    else if( a_name == GIDI_SnElasticUpScatterStyleChars ) {
        form = new Styles::SnElasticUpScatter( a_node, a_setupInfo, a_pops, a_parent ); }
    else if( a_name == GIDI_griddedCrossSectionStyleChars ) {
        form = new Styles::GriddedCrossSection( a_construction, a_node, a_setupInfo, a_pops, a_parent ); }
    else if( a_name == GIDI_URR_probabilityTablesStyleChars ) {
        form = new Styles::URR_probabilityTables( a_construction, a_node, a_setupInfo, a_pops, a_parent ); }
    else {
        std::cout << "parseStylesSuite: Ignoring unsupported style = '" << a_name << "'." << std::endl;
    }

    return( form );
}

/* *********************************************************************************************************//**
 * Function that parses a <**transportable**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::Form or nullptr if the node is not supported.
 ***********************************************************************************************************/

Form *parseTransportablesSuite( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, std::string const &a_name, LUPI_maybeUnused Styles::Suite const *a_styles ) {

    Form *form = nullptr;

    if( a_name == GIDI_transportableChars ) {
        form = new Transportable( a_construction, a_node, a_setupInfo, a_pops, a_parent ); }
    else {
        std::cout << "parseTransportablesSuite: Ignoring unsupported Form '" << a_name << "'." << std::endl;
    }

    return( form );
}

/* *********************************************************************************************************//**
 * Function that parses a <**reaction**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::Reaction instance.
 ***********************************************************************************************************/

Form *parseReaction( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, std::string const &a_name, Styles::Suite const *a_styles ) {

    return( parseReactionType( GIDI_reactionChars, a_construction, a_parent, a_node, a_setupInfo, a_pops, a_internalPoPs, a_name, a_styles ) );
}

/* *********************************************************************************************************//**
 * Function that parses an <**orphanProduct**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::Reaction instance.
 ***********************************************************************************************************/

Form *parseOrphanProduct( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, std::string const &a_name, Styles::Suite const *a_styles ) {

    return( parseReactionType( GIDI_orphanProductChars, a_construction, a_parent, a_node, a_setupInfo, a_pops, a_internalPoPs, a_name, a_styles ) );
}

/* *********************************************************************************************************//**
 * Function that parses an <**orphanProduct**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::Reaction instance.
 ***********************************************************************************************************/

Form *parseFissionComponent( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, std::string const &a_name, Styles::Suite const *a_styles ) {

    return( parseReactionType( GIDI_fissionComponentChars, a_construction, a_parent, a_node, a_setupInfo, a_pops, a_internalPoPs, a_name, a_styles ) );
}

/* *********************************************************************************************************//**
 * Function that parses a <**reaction**> or an <**orphanProduct**> node. Called from a Suite::parse instance.
 *
 * @param a_moniker                 [in]    The moniker for the form to parse.
 * @param a_construction            [in]    Used to pass user options for parsing.
 * @param a_parent                  [in]    The parent GIDI::Suite that the returned Form will be added to.
 * @param a_node                    [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo               [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops                    [in]    A PoPs Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs            [in]    The *internal* PoPI::Database instance used to get particle indices and possibly other particle information.
 *                                          This is the <**PoPs**> node under the <**reactionSuite**> node.
 * @param a_name                    [in]    The moniker for the node to be parsed.
 * @param a_styles                  [in]    A pointer to the <**styles**> node.
 * @return                                  The parsed and constructed GIDI::Reaction instance.
 ***********************************************************************************************************/

Form *parseReactionType( std::string const &a_moniker, Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node,
		SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, std::string const &a_name,
		Styles::Suite const *a_styles ) {

    Form *form = nullptr;

    if( a_name == a_moniker ) {
        Protare const &protare( *static_cast<Protare const *>( a_parent->root( ) ) );
        form = new Reaction( a_construction, a_node, a_setupInfo, a_pops, a_internalPoPs, protare, a_styles ); }
    else {                                  // This should never happend.
        std::cout << "parseReactionType: Ignoring '" << a_moniker << "' unsupported form '" << a_name << "'." << std::endl;
    }

    return( form );
}

/* *********************************************************************************************************//**
 * Function that parses a <**crossSectionSum**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::CrossSectionSum instance.
 ***********************************************************************************************************/

Form *parseSumsCrossSectionsSuite( Construction::Settings const &a_construction, LUPI_maybeUnused Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, std::string const &a_name, LUPI_maybeUnused Styles::Suite const *a_styles ) {

    Form *form = nullptr;

    if( a_name == GIDI_crossSectionSumChars ) {
        form = new Sums::CrossSectionSum( a_construction, a_node, a_setupInfo, a_pops, a_internalPoPs ); }
    else {                                  // This should never happend.
        std::cout << "parseSumsCrossSectionsSuite: Ignoring unsupported Form '" << a_name << "'." << std::endl;
    }

    return( form );
}

/* *********************************************************************************************************//**
 * Function that parses a <**multiplicitySum**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::MultiplicitySum instance.
 ***********************************************************************************************************/

Form *parseSumsMultiplicitiesSuite( Construction::Settings const &a_construction, LUPI_maybeUnused Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, std::string const &a_name, LUPI_maybeUnused Styles::Suite const *a_styles ) {

    Form *form = nullptr;

    if( a_construction.parseMode( ) == Construction::ParseMode::outline ) return( nullptr );

    if( a_name == GIDI_multiplicitySumChars ) {
        form = new Sums::MultiplicitySum( a_construction, a_node, a_setupInfo, a_pops, a_internalPoPs ); }
    else {                                  // This should never happend.
        std::cout << "parseSumsMultiplicitiesSuite: Ignoring unsupported Form '" << a_name << "'." << std::endl;
    }

    return( form );
}

/* *********************************************************************************************************//**
 * Function that parses a node under the <**doubleDifferentialCrossSection**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::Form or nullptr if the node is not supported.
 ***********************************************************************************************************/

Form *parseDoubleDifferentialCrossSectionSuite( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node,
		SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, std::string const &a_name,
		LUPI_maybeUnused Styles::Suite const *a_styles ) {

    if( a_construction.parseMode( ) == Construction::ParseMode::outline ) return( nullptr );
    if( a_construction.parseMode( ) == Construction::ParseMode::multiGroupOnly ) return( nullptr );

    Form *form = nullptr;

    if( a_name == GIDI_coherentPhotonScatteringChars ) {
        form = new DoubleDifferentialCrossSection::CoherentPhotoAtomicScattering( a_construction, a_node, a_setupInfo, a_pops, a_internalPoPs, a_parent ); }
    else if( a_name == GIDI_incoherentPhotonScatteringChars ) {
        form = new DoubleDifferentialCrossSection::IncoherentPhotoAtomicScattering( a_construction, a_node, a_setupInfo, a_pops, a_internalPoPs, a_parent ); }
    else if( a_name == GIDI_incoherentBoundToFreePhotonScatteringChars ) {
        form = new DoubleDifferentialCrossSection::IncoherentBoundToFreePhotoAtomicScattering( a_construction, a_node, a_setupInfo, a_pops, a_internalPoPs, a_parent ); }
    else if( a_name == GIDI_TNSL_coherentElasticChars ) {
        form = new DoubleDifferentialCrossSection::n_ThermalNeutronScatteringLaw::CoherentElastic( a_construction, a_node, a_setupInfo, a_pops, a_internalPoPs, a_parent ); }
    else if( a_name == GIDI_TNSL_incoherentElasticChars ) {
        form = new DoubleDifferentialCrossSection::n_ThermalNeutronScatteringLaw::IncoherentElastic( a_construction, a_node, a_setupInfo, a_pops, a_internalPoPs, a_parent ); }
    else if( a_name == GIDI_TNSL_incoherentInelasticChars ) {
        form = new DoubleDifferentialCrossSection::n_ThermalNeutronScatteringLaw::IncoherentInelastic( a_construction, a_node, a_setupInfo, a_pops, a_internalPoPs, a_parent ); }
    else if( a_name == GIDI_CoulombPlusNuclearElasticChars ) { 
        }
    else {
        std::cout << "parseDoubleDifferentialCrossSectionSuite: Ignoring unsupported Form '" << a_name << "'." << std::endl;
    }

    return( form );
}

/* *********************************************************************************************************//**
 * Function that parses a node under the <**doubleDifferentialCrossSection**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::Form or nullptr if the node is not supported.
 ***********************************************************************************************************/

Form *parseScatteringAtom( Construction::Settings const &a_construction, LUPI_maybeUnused Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, LUPI_maybeUnused std::string const &a_name, LUPI_maybeUnused Styles::Suite const *a_styles ) {

    return( new DoubleDifferentialCrossSection::n_ThermalNeutronScatteringLaw::ScatteringAtom( a_construction, a_node, a_setupInfo ) );
}

/* *********************************************************************************************************//**
 * Function that parses a node under the <**crossSection**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::Form or nullptr if the node is not supported.
 ***********************************************************************************************************/

Form *parseCrossSectionSuite( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, std::string const &a_name, LUPI_maybeUnused Styles::Suite const *a_styles ) {

    if( a_construction.parseMode( ) == Construction::ParseMode::outline ) return( nullptr );
    if( ( a_construction.parseMode( ) == Construction::ParseMode::multiGroupOnly ) && ( a_name != GIDI_gridded1dChars ) ) return( nullptr );
    if( ( a_construction.parseMode( ) == Construction::ParseMode::MonteCarloContinuousEnergy ) && ( a_name != GIDI_Ys1dChars ) ) return( nullptr );

// Form not parsed is CoulombPlusNuclearElastic.
    Form *form = nullptr;

    if( a_name == GIDI_resonancesWithBackgroundChars ) {
        return( new Functions::ResonancesWithBackground1d( a_construction, a_node, a_setupInfo, a_parent ) ); }
    else if( a_name == GIDI_TNSL1dChars ) {
        form = new Functions::ThermalNeutronScatteringLaw1d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_URR_probabilityTables1dChars ) {
        form = new Functions::URR_probabilityTables1d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_CoulombPlusNuclearElasticChars ) {
        }
    else {
        form = data1dParse( a_construction, a_node, a_setupInfo, a_parent );
    }

    return( form );
}

/* *********************************************************************************************************//**
 * Function that parses a node under the <**DelayedNeutrons**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::Form or nullptr if the node is not supported.
 ***********************************************************************************************************/

Form *parseDelayedNeutronsSuite( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, std::string const &a_name, Styles::Suite const *a_styles ) {

    if( a_name != GIDI_delayedNeutronChars ) throw Exception( std::string( "Invalid " ) + GIDI_delayedNeutronsChars + " child node of moniker " + a_name );
    return( new DelayedNeutron( a_construction, a_node, a_setupInfo, a_pops, a_internalPoPs, a_parent, a_styles ) );
}

/* *********************************************************************************************************//**
 * Function that parses a node under the <**FissionEnergyReleases**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::Form or nullptr if the node is not supported.
 ***********************************************************************************************************/

Form *parseFissionEnergyReleasesSuite( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, std::string const &a_name, LUPI_maybeUnused Styles::Suite const *a_styles ) {

    if( a_name != GIDI_fissionEnergyReleaseChars ) throw Exception( std::string( "Invalid " ) + GIDI_fissionEnergyReleasesChars " child node of moniker " + a_name );
    return( new Functions::FissionEnergyRelease( a_construction, a_node, a_setupInfo, a_parent ) );
}

/* *********************************************************************************************************//**
 * Function that parses a node under the <**rate**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::Form or nullptr if the node is not supported.
 ***********************************************************************************************************/

Form *parsePhysicalQuantitySuite( LUPI_maybeUnused Construction::Settings const &a_construction, LUPI_maybeUnused Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, LUPI_maybeUnused std::string const &a_name, LUPI_maybeUnused Styles::Suite const *a_styles ) {

    return( new PhysicalQuantity( a_node, a_setupInfo ) );
}

/* *********************************************************************************************************//**
 * Function that parses a node under an <**availableEnergy**> or <**availableMomentum**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::Function1d instance.
 ***********************************************************************************************************/

Form *parseAvailableSuite( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, std::string const &a_name, LUPI_maybeUnused Styles::Suite const *a_styles ) {

    if( a_construction.parseMode( ) == Construction::ParseMode::outline ) return( nullptr );
    if( ( a_construction.parseMode( ) == Construction::ParseMode::multiGroupOnly ) && ( a_name != GIDI_gridded1dChars ) ) return( nullptr );

    return( data1dParse( a_construction, a_node, a_setupInfo, a_parent ) );
}

/* *********************************************************************************************************//**
 * Function that parses a node under a <**Q**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::Function1d instance.
 ***********************************************************************************************************/

Form *parseQSuite( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, LUPI_maybeUnused std::string const &a_name, LUPI_maybeUnused Styles::Suite const *a_styles ) {

    Form *form = nullptr;

    if( a_construction.parseMode( ) == Construction::ParseMode::outline ) return( nullptr );

    form = data1dParse( a_construction, a_node, a_setupInfo, a_parent );

    return( form );
}

/* *********************************************************************************************************//**
 * Function that parses a node under a <**products**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::Product instance.
 ***********************************************************************************************************/

Form *parseProductSuite( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, std::string const &a_name, Styles::Suite const *a_styles ) {

    Form *form = nullptr;

    if( a_name == GIDI_productChars ) {
        form = new Product( a_construction, a_node, a_setupInfo, a_pops, a_internalPoPs, a_parent, a_styles ); }
    else {
        std::cout << "parseProductSuite: Ignoring unsupported element in products " << a_node.name( ) << std::endl;
    }

    return( form );
}

/* *********************************************************************************************************//**
 * Function that parses a node under a <**multiplicity**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::Function1d instance.
 ***********************************************************************************************************/

Form *parseMultiplicitySuite( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, std::string const &a_name, LUPI_maybeUnused Styles::Suite const *a_styles ) {

    if( a_construction.parseMode( ) == Construction::ParseMode::outline ) return( nullptr );

    if( a_name == GIDI_branching1dChars ) return( new Functions::Branching1d( a_construction, a_node, a_setupInfo, a_parent ) );

    return( data1dParse( a_construction, a_node, a_setupInfo, a_parent ) );
}

/* *********************************************************************************************************//**
 * Function that parses a node under a <**distribution**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::Form or nullptr if the node is not supported.
 ***********************************************************************************************************/

Form *parseDistributionSuite( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, std::string const &a_name, LUPI_maybeUnused Styles::Suite const *a_styles ) {

    if( a_construction.parseMode( ) == Construction::ParseMode::outline ) return( nullptr );
    if( a_name == GIDI_multiGroup3dChars ) {
        if( ( a_construction.parseMode( ) == Construction::ParseMode::MonteCarloContinuousEnergy ) || 
            ( a_construction.parseMode( ) == Construction::ParseMode::excludeProductMatrices ) ) return( nullptr ); }
    else {
        if( a_construction.parseMode( ) == Construction::ParseMode::multiGroupOnly ) return( nullptr );
    }

//  Distributions not parsed are GIDI_LLNLLegendreChars.
    Form *form = nullptr;

    if( a_name == GIDI_multiGroup3dChars ) {
        form = new Distributions::MultiGroup3d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_angularTwoBodyChars ) {
        form = new Distributions::AngularTwoBody( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_uncorrelatedChars ) {
        form = new Distributions::Uncorrelated( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_KalbachMannChars ) {
        form = new Distributions::KalbachMann( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_energyAngularChars ) {
        form = new Distributions::EnergyAngular( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_energyAngularMCChars ) {
        form = new Distributions::EnergyAngularMC( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_angularEnergyChars ) {
        form = new Distributions::AngularEnergy( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_angularEnergyMCChars ) {
        form = new Distributions::AngularEnergyMC( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_LLNLAngularEnergyChars ) {
        form = new Distributions::LLNLAngularEnergy( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_coherentPhotonScatteringChars ) {
        form = new Distributions::CoherentPhotoAtomicScattering( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_incoherentPhotonScatteringChars ) {
        form = new Distributions::IncoherentPhotoAtomicScattering( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_incoherentBoundToFreePhotonScatteringChars ) {
        form = new Distributions::IncoherentBoundToFreePhotoAtomicScattering( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_thermalNeutronScatteringLawChars ) {
        form = new Distributions::ThermalNeutronScatteringLaw( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_branching3dChars ) {
        form = new Distributions::Branching3d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_referenceChars ) {
        form = new Distributions::Reference3d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_unspecifiedChars ) {
        form = new Distributions::Unspecified( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_CoulombPlusNuclearElasticChars ) {
        form = new Distributions::CoulombPlusNuclearElastic( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( a_name == GIDI_LLNLLegendreChars ) {
        form = new Distributions::LLNLLegendre( a_construction, a_node, a_setupInfo, a_parent ); }
    else {
        std::cout << "parseDistributionSuite: Ignoring unsupported distribution " << a_node.name( ) << std::endl;
    }

    return( form );
}

/* *********************************************************************************************************//**
 * Function that parses a node under an <**averageEnergy**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::Function1d instance.
 ***********************************************************************************************************/

Form *parseAverageEnergySuite( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, std::string const &a_name, LUPI_maybeUnused Styles::Suite const *a_styles ) {

    if( a_construction.parseMode( ) == Construction::ParseMode::outline ) return( nullptr );
    if( ( a_construction.parseMode( ) == Construction::ParseMode::multiGroupOnly ) && ( a_name != GIDI_gridded1dChars ) ) return( nullptr );

    return( data1dParse( a_construction, a_node, a_setupInfo, a_parent ) );
}

/* *********************************************************************************************************//**
 * Function that parses a node under an <**averageMomentum**> node. Called from a Suite::parse instance.
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
 * @return                                  The parsed and constructed GIDI::Function1d instance.
 ***********************************************************************************************************/

Form *parseAverageMomentumSuite( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		        LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, std::string const &a_name, LUPI_maybeUnused Styles::Suite const *a_styles ) {

    if( a_construction.parseMode( ) == Construction::ParseMode::outline ) return( nullptr );
    if( ( a_construction.parseMode( ) == Construction::ParseMode::multiGroupOnly ) && ( a_name != GIDI_gridded1dChars ) ) return( nullptr );

    return( data1dParse( a_construction, a_node, a_setupInfo, a_parent ) );
}

/* *********************************************************************************************************//**
 * This function parses a **probabilityTable** child node of a **probabilityTables** node.
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
 * @return                                  The parsed and constructed GIDI::Function1d instance.
 ***********************************************************************************************************/

Form *parseACE_URR_probabilityTables( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
                LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, std::string const &a_name, LUPI_maybeUnused Styles::Suite const *a_styles ) {

    if( a_name != GIDI_ACE_URR_probabilityTableChars ) throw Exception( std::string( "Invalid " ) + GIDI_ACE_URR_probabilityTablesChars " child node of moniker " + a_name );

    return( new ACE_URR::ProbabilityTable( a_construction, a_node, a_setupInfo, a_parent ) );
}

/* *********************************************************************************************************//**
 * This function parses a **column** child node of a **columnHeaders** node.
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
 * @return                                  The parsed and constructed GIDI::Function1d instance.
 ***********************************************************************************************************/

Form *parseColumnHeaders( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
                LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, std::string const &a_name, LUPI_maybeUnused Styles::Suite const *a_styles ) {

    if( a_name != GIDI_columnChars ) throw Exception( std::string( "Invalid " ) + GIDI_columnHeadersChars + " child node of moniker " + a_name );

    Table::Column *column = new Table::Column( a_construction, a_node, a_setupInfo, a_parent );
    column->setKeyName( GIDI_indexChars );

    return column;
}

/* *********************************************************************************************************//**
 * Function that parses a node one-d function node. Called from a Suite::parse instance.
 *
 * @param a_construction            [in]    Used to pass user options for parsing.
 * @param a_node                    [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo               [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent                  [in]    The parent GIDI::Suite that the returned Form will be added to.
 *
 * @return                          The parsed and constructed GIDI::Function1d instance.
 ***********************************************************************************************************/

Functions::Function1dForm *data1dParse( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) {

    Functions::Function1dForm *form = nullptr;
    std::string name( a_node.name( ) );

    if( name == GIDI_constant1dChars ) {
        form = new Functions::Constant1d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_XYs1dChars ) {
        form = new Functions::XYs1d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_Ys1dChars ) {
        form = new Functions::Ys1d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_polynomial1dChars ) {
        form = new Functions::Polynomial1d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_LegendreChars ) {
        form = new Functions::Legendre1d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_regions1dChars ) {
        form = new Functions::Regions1d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_gridded1dChars ) {
        form = new Functions::Gridded1d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_referenceChars ) {
        form = new Functions::Reference1d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_xs_pdf_cdf1dChars ) {
        form = new Functions::Xs_pdf_cdf1d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_unspecifiedChars ) {
        form = new Functions::Unspecified1d( a_construction, a_node, a_setupInfo, a_parent ); }
    else {
        std::cout << "data1dParse: Ignoring unsupported 1d function = '" << name << "'" << std::endl;
    }

    return( form );
}

/* *********************************************************************************************************//**
 * Function that parses a node one-d function node. Called from a Suite::parse instance. If no node exists, returns nullptr.
 *
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent          [in]    The parent GIDI::Suite that the returned Form will be added to.
 *
 * @return                          The parsed and constructed GIDI::Function1d instance.
 ***********************************************************************************************************/

Functions::Function1dForm *data1dParseAllowEmpty( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		Suite *a_parent ) {

    std::string name( a_node.name( ) );

    if( name == "" ) return( nullptr );
    return( data1dParse( a_construction, a_node, a_setupInfo, a_parent ) );
}

/* *********************************************************************************************************//**
 * Function that parses the list of 1d function nodes contained in *a_node*.
 *
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_function1ds     [in]    The object to fill with the list of parsed 1d functions.
 ***********************************************************************************************************/

void data1dListParse( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
                std::vector<Functions::Function1dForm *> &a_function1ds ) {

    for( HAPI::Node child = a_node.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
        Functions::Function1dForm *form = data1dParse( a_construction, child, a_setupInfo, nullptr );

        if( form == nullptr ) throw Exception( "data1dListParse data1dParse returned nullptr." );
        a_function1ds.push_back( form );
    }
}

/* *********************************************************************************************************//**
 * Function that parses a node two-d function node. Called from a Suite::parse instance.
 *
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent          [in]    The parent GIDI::Suite that the returned Form will be added to.
 * @return                          The parsed and constructed GIDI::Function2d instance.
 ***********************************************************************************************************/
 
Functions::Function2dForm *data2dParse( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) {

    Functions::Function2dForm *form = nullptr;
    std::string name( a_node.name( ) );

    if( name == GIDI_XYs2dChars ) {
        form = new Functions::XYs2d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_recoilChars ) {
        form = new Functions::Recoil2d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_isotropic2dChars ) {
        form = new Functions::Isotropic2d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_discreteGammaChars ) {
        form = new Functions::DiscreteGamma2d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_primaryGammaChars ) {
        form = new Functions::PrimaryGamma2d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_generalEvaporationChars ) {
        form = new Functions::GeneralEvaporation2d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_simpleMaxwellianFissionChars ) {
        form = new Functions::SimpleMaxwellianFission2d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_evaporationChars ) {
        form = new Functions::Evaporation2d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_WattChars ) {
        form = new Functions::Watt2d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_MadlandNixChars ) {
        form = new Functions::MadlandNix2d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_weightedFunctionalsChars ) {
        form = new Functions::WeightedFunctionals2d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_NBodyPhaseSpaceChars ) {
        form = new Functions::NBodyPhaseSpace2d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_regions2dChars ) {
        form = new Functions::Regions2d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_gridded2dChars ) {
        form = new Functions::Gridded2d( a_construction, a_node, a_setupInfo, a_parent ); }
    else {
        std::cout << "data2dParse: Ignoring unsupported 2d function = '" << name << "'" << std::endl;
    }

    return( form );
}

/* *********************************************************************************************************//**
 * Function that parses the list of 2d function nodes contained in *a_node*.
 *
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_function2ds     [in]    The object to fill with the list of parsed 2d functions.
 ***********************************************************************************************************/

void data2dListParse( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		std::vector<Functions::Function2dForm *> &a_function2ds ) {

    for( HAPI::Node child = a_node.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
        Functions::Function2dForm *form = data2dParse( a_construction, child, a_setupInfo, nullptr );

        if( form == nullptr ) throw Exception( "data2dListParse data2dParse returned nullptr." );
        a_function2ds.push_back( form );
    }
}

/* *********************************************************************************************************//**
 * Function that parses a node three-d function node. Called from a Suite::parse instance.
 *
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent          [in]    The parent GIDI::Suite that the returned Form will be added to.
 *
 * @return                          The parsed and constructed GIDI::Function3d instance.
 ***********************************************************************************************************/
 
Functions::Function3dForm *data3dParse( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) {

    Functions::Function3dForm *form = nullptr;
    std::string name( a_node.name( ) );

    if( name == GIDI_XYs3dChars ) {
        form = new Functions::XYs3d( a_construction, a_node, a_setupInfo, a_parent ); }
    else if( name == GIDI_gridded3dChars ) {
        form = new Functions::Gridded3d( a_construction, a_node, a_setupInfo ); }
    else {
        std::cout << "data3dParse: Ignoring unsupported 3d function = '" << name << "'" << std::endl;
    }
    
    return( form );
}

/* *********************************************************************************************************//**
 * Function that checks that the outerDomainValue values of the elements of *a_functions* are increasing and fills *a_Xs* with the outerDomainValue values.
 *
 * @param a_functions       [in]    List of functions to check.
 * @param a_Xs              [in]    A list of doubles that is filled with the outerDomainValues from the list of functions in *a_functions*.
 ***********************************************************************************************************/

void checkOuterDomainValues1d( std::vector<Functions::Function1dForm *> &a_functions, std::vector<double> &a_Xs ) {

    for( auto iter = a_functions.begin( ); iter != a_functions.end( ); ++iter ) {
        if( a_Xs.size( ) > 0 ) {
            if( (*iter)->outerDomainValue( ) <= a_Xs.back( ) ) throw Exception( "checkOuterDomainValues1d: next outerDomainValue <= current outerDomainValue." );
        }
        a_Xs.push_back( (*iter)->outerDomainValue( ) );
    }
}

/* *********************************************************************************************************//**
 * Function that checks that the outerDomainValue values of the elements of *a_functions* are increasing and fills *a_Xs* with the outerDomainValue values.
 *
 * @param a_functions       [in]    List of functions to check.
 * @param a_Xs              [in]    A list of doubles that is filled with the outerDomainValues from the list of functions in *a_functions*.
 ***********************************************************************************************************/

void checkOuterDomainValues2d( std::vector<Functions::Function2dForm *> &a_functions, std::vector<double> &a_Xs ) {

    for( auto iter = a_functions.begin( ); iter != a_functions.end( ); ++iter ) {
        if( a_Xs.size( ) > 0 ) {
            if( (*iter)->outerDomainValue( ) <= a_Xs.back( ) ) throw Exception( "checkOuterDomainValues2d: next outerDomainValue <= current outerDomainValue." );
        }
        a_Xs.push_back( (*iter)->outerDomainValue( ) );
    }
}

/* *********************************************************************************************************//**
 * Function that checks that the domain overlap of the elements of *a_functions* are .
 * The domain minimum values from each function and the domain maximum value are filled into *a_Xs*.
 *
 * @param a_functions       [in]    The list of functions whose domain limits are checked.
 * @param a_Xs              [in]    A std::vector<double> that is with the domain minimum values from each function and the domain maximum value.
 ***********************************************************************************************************/

void checkSequentialDomainLimits1d( std::vector<Functions::Function1dForm *> &a_functions, std::vector<double> &a_Xs ) {

    double domainMax = -1;

    for( auto iter = a_functions.begin( ); iter != a_functions.end( ); ++iter ) {
        double domainMin = (*iter)->domainMin( );

        if( a_Xs.size( ) > 0 ) {
            if( fabs( domainMax - domainMin ) > 1e-8 * ( fabs( domainMin ) + fabs( domainMax ) ) ) 
                    throw Exception( "checkSequentialDomainLimits1d: domains not abutting." );
        }
        a_Xs.push_back( domainMin );
        domainMax = (*iter)->domainMax( );
    }
    if( a_Xs.size( ) > 0 ) a_Xs.push_back( domainMax );
}

/* *********************************************************************************************************//**
 * Function that checks that the domain overlap of the elements of *a_functions* are .
 * The domain minimum values from each function and the domain maximum value are filled into *a_Xs*.
 *
 * @param a_functions       [in]    The list of functions whose domain limits are checked.
 * @param a_Xs              [in]    A std::vector<double> that is with the domain minimum values from each function and the domain maximum value.
 ***********************************************************************************************************/

void checkSequentialDomainLimits2d( std::vector<Functions::Function2dForm *> &a_functions, std::vector<double> &a_Xs ) {

    double domainMax = -1;

    for( auto iter = a_functions.begin( ); iter != a_functions.end( ); ++iter ) {
        double domainMin = (*iter)->domainMin( );

        if( a_Xs.size( ) > 0 ) {
            if( fabs( domainMax - domainMin ) > 1e-8 * ( fabs( domainMin ) + fabs( domainMax ) ) ) 
                    throw Exception( "checkSequentialDomainLimits2d: domains not abutting." );
        }
        a_Xs.push_back( domainMin );
        domainMax = (*iter)->domainMax( );
    }
    if( a_Xs.size( ) > 0 ) a_Xs.push_back( a_functions.back( )->domainMax( ) );
}

}
