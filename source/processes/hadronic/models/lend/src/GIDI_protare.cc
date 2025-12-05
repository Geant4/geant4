/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdlib.h>
#include <algorithm>
#include <cmath>

#include "GIDI.hpp"
#include <HAPI.hpp>

namespace GIDI {

static bool sortTemperatures( Styles::TemperatureInfo const &lhs, Styles::TemperatureInfo const &rhs );

/*! \class Protare
 * Base class for the protare sub-classes.
 */

/* *********************************************************************************************************//**
 * Base Protare constructor.
 ***********************************************************************************************************/

Protare::Protare( ) :
        GUPI::Ancestry( "" ),
        m_projectile( "", "", -1.0 ),
        m_target( "", "", -1.0 ),
        m_GNDS_target( "", "", -1.0 ) {
    
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Protare::~Protare( ) {

}

/* *********************************************************************************************************//**
 * Called by the constructs. This method does most of the parsing.
 *
 * @param a_node                        [in]    The protare (i.e., reactionSuite) node to be parsed and used to construct a Protare.
 * @param a_setupInfo                   [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops                        [in]    A PoPs Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs                [in]    The internal PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_targetRequiredInGlobalPoPs  [in]    If *true*, the target is required to be in **a_pops**.
 * @param a_requiredInPoPs              [in]    If *true*, particle is required to be in **a_pops**.
 ***********************************************************************************************************/

void Protare::initialize( HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs,
		bool a_targetRequiredInGlobalPoPs, bool a_requiredInPoPs ) {

    setMoniker( a_node.name( ) );    

    std::string projectileID = a_node.attribute_as_string( GIDI_projectileChars );
    m_projectile = ParticleInfo( projectileID, a_pops, a_internalPoPs, a_requiredInPoPs );

    std::string targetID = a_node.attribute_as_string( GIDI_targetChars );
    m_GNDS_target = ParticleInfo( targetID, a_pops, a_internalPoPs, a_targetRequiredInGlobalPoPs && a_requiredInPoPs );

    auto iter = a_setupInfo.m_particleSubstitution->find( m_GNDS_target.ID( ) );
    if( iter != a_setupInfo.m_particleSubstitution->end( ) ) {
        m_target = iter->second; }
    else {
        m_target = m_GNDS_target;
    }
}

/* *********************************************************************************************************//**
 * If the protare is a ProtareTNSL then summing over all reactions will include the standard protare's elastic cross section
 * in the domain of the TNSL data. The standard elastic cross section should not be added in this domain.
 * If needed, this function corrects the cross section for this over counting of the elastic cross section.
 * This method does nothing unless overwritten by the ProtareTNSL class.
 *
 * @param       a_label                     [in]    The label of the elastic cross section data to use if over counting needs to be corrected.
 * @param       a_crossSectionSum           [in]    The cross section to correct.
 ***********************************************************************************************************/

void Protare::TNSL_crossSectionSumCorrection( LUPI_maybeUnused std::string const &a_label, LUPI_maybeUnused Functions::XYs1d &a_crossSectionSum ) {

}

/* *********************************************************************************************************//**
 * If the protare is a ProtareTNSL then summing over all reactions will include the standard protare's elastic cross section
 * in the domain of the TNSL data. The standard elastic cross section should not be added in this domain.
 * If needed, this function corrects the cross section for this over counting of the elastic cross section.
 * This method does nothing unless overwritten by the ProtareTNSL class.
 *
 * @param       a_label                     [in]    The label of the elastic cross section data to use if over counting needs to be corrected.
 * @param       a_crossSectionSum           [in]    The cross section to correct.
 ***********************************************************************************************************/

void Protare::TNSL_crossSectionSumCorrection( LUPI_maybeUnused std::string const &a_label, LUPI_maybeUnused Functions::Ys1d &a_crossSectionSum ) {

}

/* *********************************************************************************************************//**
 * If the protare is a ProtareTNSL then summing over all reactions will include the standard protare's elastic cross section
 * in the domain of the TNSL data. The standard elastic cross section should not be added in this domain.
 * If needed, this function corrects the cross section for this over counting of the elastic cross section.
 * This method does nothing as the multi-group cross section data for the standard protare's elastic cross section are
 * zeroed when the data are read in. However, this method is added so codes do not have to check the type of data they are accessing.
 *
 * @param       a_label                     [in]    The label of the elastic cross section data to use if over counting needs to be corrected.
 * @param       a_crossSectionSum           [in]    The cross section to correct.
 ***********************************************************************************************************/

void Protare::TNSL_crossSectionSumCorrection( LUPI_maybeUnused std::string const &a_label, LUPI_maybeUnused Vector &a_crossSectionSum ) {

}

/* *********************************************************************************************************//**
 * Returns a list of all reaction indices whose ENDL C value is in the set *a_CValues*.
 *
 * @param       a_CValues                       [in]    A list of ENDL C values.
 * @param       a_checkActiveState              [in]    If true, all reactions whose active state is false are not included in the returned set even if their CValue match on in the list.
 *
 * @return                                      The list of reaction indices.
 ***********************************************************************************************************/

ExcludeReactionsSet Protare::reactionIndicesMatchingENDLCValues( std::set<int> const &a_CValues, bool a_checkActiveState ) {

    ExcludeReactionsSet indices;

    for( std::size_t i1 = 0; i1 < numberOfReactions( ); ++i1 ) {
        Reaction *reaction1 = reaction( i1 );

        if( a_checkActiveState && !reaction1->active( ) ) continue;
        if( a_CValues.find( reaction1->ENDL_C( ) ) != a_CValues.end( ) ) indices.insert( i1 );
    }

    return( indices );
}

/*! \class ProtareSingle
 * Class to store a GNDS <**reactionSuite**> node.
 */

/* *********************************************************************************************************//**
 * Parses a GNDS file to construct the Protare instance. Calls the initialize method which does most of the work.
 *
 * @param a_pops            [in]    A PoPs Database instance used to get particle indices and possibly other particle information.
 * @param a_projectileID    [in]    The PoPs id of the projectile.
 * @param a_targetID        [in]    The PoPs id of the target.
 * @param a_evaluation      [in]    The evaluation string.
 * @param a_interaction     [in]    The interaction flag for the protare.
 * @param a_formatVersion   [in]    The GNDS format to use.
 ***********************************************************************************************************/

ProtareSingle::ProtareSingle( PoPI::Database const &a_pops, std::string const &a_projectileID, std::string const &a_targetID, std::string const &a_evaluation,
                std::string const &a_interaction, std::string const &a_formatVersion ) :
        m_doc( nullptr ),
        m_dataManager( nullptr ),
        m_numberOfLazyParsingHelperForms( 0 ),
        m_numberOfLazyParsingHelperFormsReplaced( 0 ),
        m_formatVersion( a_formatVersion ),
        m_evaluation( a_evaluation ),
        m_interaction( a_interaction ),
        m_projectileFrame( Frame::lab ),
        m_decayPositronium( false ),
        m_thresholdFactor( 0.0 ),
        m_nuclearPlusCoulombInterferenceOnlyReaction( nullptr ),
        m_pointwiseAverageProductEnergy( GIDI_averageEnergyChars, GIDI_labelChars ) {

    setMoniker( GIDI_topLevelChars );
    initialize( );

    setProjectile( ParticleInfo( a_projectileID, a_pops, a_pops, true ) );
    setTarget( ParticleInfo( a_targetID, a_pops, a_pops, true ) );
}

/* *********************************************************************************************************//**
 * Parses a GNDS file to construct the Protare instance. Calls the initialize method which does most of the work.
 *
 * @param a_construction                [in]    Used to pass user options to the constructor.
 * @param a_fileName                    [in]    File containing a protare (i.e., reactionSuite) node that is parsed and used to construct the Protare.
 * @param a_fileType                    [in]    File type of a_fileType. Currently, only GIDI::FileType::XML and GIDI::FileType::HDF are supported.
 * @param a_pops                        [in]    A PoPs Database instance used to get particle indices and possibly other particle information.
 * @param a_particleSubstitution        [in]    Map of particles to substitute with another particles.
 * @param a_libraries                   [in]    The list of libraries that were searched to find *this*.
 * @param a_interaction                 [in]    The interaction flag for the protare.
 * @param a_targetRequiredInGlobalPoPs  [in]    If *true*, the target is required to be in **a_pops**.
 * @param a_requiredInPoPs              [in]    If *true*, no particle is required to be in **a_pops**.
 ***********************************************************************************************************/

ProtareSingle::ProtareSingle( Construction::Settings const &a_construction, std::string const &a_fileName, FileType a_fileType, 
                PoPI::Database const &a_pops, ParticleSubstitution const &a_particleSubstitution, std::vector<std::string> const &a_libraries, 
                std::string const &a_interaction, bool a_targetRequiredInGlobalPoPs, bool a_requiredInPoPs ) :
        Protare( ),
        m_doc( nullptr ),
        m_dataManager( nullptr ),
        m_numberOfLazyParsingHelperForms( 0 ),
        m_numberOfLazyParsingHelperFormsReplaced( 0 ),
        m_libraries( a_libraries ),
        m_interaction( a_interaction ),
        m_fileName( a_fileName ),
        m_realFileName( LUPI::FileInfo::realPath( a_fileName ) ),
        m_decayPositronium( a_construction.decayPositronium( ) ),
        m_pointwiseAverageProductEnergy( GIDI_averageEnergyChars, GIDI_labelChars ) {

#ifdef HAPI_USE_PUGIXML
    if( a_fileType == GIDI::FileType::XML ) {
        m_doc = new HAPI::PugiXMLFile( a_fileName.c_str( ), "ProtareSingle::ProtareSingle" );
    }
#endif
#ifdef HAPI_USE_HDF5
    if( a_fileType == GIDI::FileType::HDF ) {
        m_doc = new HAPI::HDFFile( a_fileName.c_str( ) );
    }
#endif
    if( m_doc == nullptr ) {
        throw std::runtime_error( "Only XML/HDF file types supported." );
    }

#ifdef GIDIP_TEST_PARSING
    if( a_construction.parseMode( ) != Construction::ParseMode::noParsing ) {
#endif
        HAPI::Node protare = m_doc->first_child( );

        SetupInfo setupInfo( this );
        ParticleSubstitution particleSubstitution( a_particleSubstitution );
        setupInfo.m_particleSubstitution = &particleSubstitution;

        initialize( a_construction, protare, setupInfo, a_pops, a_targetRequiredInGlobalPoPs, a_requiredInPoPs );
#ifdef GIDIP_TEST_PARSING
    }
#endif
}

/* *********************************************************************************************************//**
 * Parses a GNDS HAPI::Node instance to construct the Protare instance. Calls the ProtareSingle::initialize method which does most of the work.
 *
 * @param a_construction                [in]    Used to pass user options to the constructor.
 * @param a_node                        [in]    The **HAPI::Node** to be parsed and used to construct the Protare.
 * @param a_pops                        [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_particleSubstitution        [in]    Map of particles to substitute with another particles.
 * @param a_libraries                   [in]    The list of libraries that were searched to find *this*.
 * @param a_targetRequiredInGlobalPoPs  [in]    If *true*, the target is required to be in **a_pops**.
 * @param a_interaction                 [in]    The interaction between the projectile and target.
 * @param a_requiredInPoPs              [in]    If *true*, no particle is required to be in **a_pops**.
 ***********************************************************************************************************/

ProtareSingle::ProtareSingle( Construction::Settings const &a_construction, HAPI::Node const &a_node, PoPI::Database const &a_pops,
                ParticleSubstitution const &a_particleSubstitution, std::vector<std::string> const &a_libraries, 
                LUPI_maybeUnused std::string const &a_interaction, bool a_targetRequiredInGlobalPoPs, bool a_requiredInPoPs ) :
        Protare( ),
        m_doc( nullptr ),
        m_dataManager( nullptr ),
        m_numberOfLazyParsingHelperForms( 0 ),
        m_numberOfLazyParsingHelperFormsReplaced( 0 ),
        m_libraries( a_libraries ),
        m_pointwiseAverageProductEnergy( GIDI_averageEnergyChars, GIDI_labelChars ) {

    SetupInfo setupInfo( this );
    ParticleSubstitution particleSubstitution( a_particleSubstitution );
    setupInfo.m_particleSubstitution = &particleSubstitution;

    initialize( a_construction, a_node, setupInfo, a_pops, a_targetRequiredInGlobalPoPs, a_requiredInPoPs );
}

/* *********************************************************************************************************//**
 * Base initializer to be called by all constructors (directly or indirectly).
 ***********************************************************************************************************/

void ProtareSingle::initialize( ) {

    m_externalFiles.setMoniker( GIDI_externalFilesChars );
    m_externalFiles.setAncestor( this );

    m_styles.setMoniker( GIDI_stylesChars );
    m_styles.setAncestor( this );

    m_documentations.setAncestor( this );
    m_documentations.setMoniker( GIDI_documentations_1_10_Chars );

    m_reactions.setAncestor( this );
    m_reactions.setMoniker( GIDI_reactionsChars );

    m_orphanProducts.setAncestor( this );
    m_orphanProducts.setMoniker( GIDI_orphanProductsChars );

    m_incompleteReactions.setAncestor( this );
    m_incompleteReactions.setMoniker( GIDI_reactionsChars );

    m_sums.setAncestor( this );

    m_fissionComponents.setAncestor( this );
    m_fissionComponents.setMoniker( GIDI_fissionComponentsChars );

    m_RutherfordScatteringPresent = false;
    m_onlyRutherfordScatteringPresent = false;
    m_nuclearPlusCoulombInterferenceOnlyReaction = nullptr;
    m_multiGroupSummedReaction = nullptr;
    m_multiGroupSummedDelayedNeutrons = nullptr;

    m_ACE_URR_probabilityTables.setAncestor( this );
    m_ACE_URR_probabilityTables.setMoniker( GIDI_ACE_URR_probabilityTablesChars );

    m_photoAtomicIncoherentDoppler.setAncestor( this );
    m_photoAtomicIncoherentDoppler.setMoniker( GIDI_LLNL_photoAtomicIncoherentDoppler_Chars );

    m_pointwiseAverageProductEnergy.setAncestor( this );
    m_GRIN_continuumGammas = nullptr;
}

/* *********************************************************************************************************//**
 * Called by the constructs. This method does most of the parsing.
 *
 * @param a_construction                [in]    Used to pass user options to the constructor.
 * @param a_node                        [in]    The protare (i.e., reactionSuite) node to be parsed and used to construct a Protare.
 * @param a_setupInfo                   [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops                        [in]    A PoPs Database instance used to get particle indices and possibly other particle information.
 * @param a_targetRequiredInGlobalPoPs  [in]    If *true*, the target is required to be in **a_pops**.
 * @param a_requiredInPoPs              [in]    If *true*, no particle is required to be in **a_pops**.
 ***********************************************************************************************************/

void ProtareSingle::initialize( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		PoPI::Database const &a_pops, bool a_targetRequiredInGlobalPoPs, bool a_requiredInPoPs ) {

    if( a_node.name( ) != GIDI_topLevelChars ) throw Exception( "Node '" + a_node.name( ) + "' is not a 'reactionSuite' node." );

    m_externalFiles.parse( a_construction, a_node.child( GIDI_externalFilesChars ), a_setupInfo, a_pops, m_internalPoPs, parseExternalFilesSuite, nullptr );
    std::string parentDir = m_fileName.substr( 0, m_fileName.find_last_of( "/" ) );
    m_externalFiles.registerBinaryFiles( parentDir, a_setupInfo );

    HAPI::Node internalPoPs = a_node.child( GIDI_PoPsChars );
    m_internalPoPs.addDatabase( internalPoPs, true );
    std::vector<PoPI::Alias *> const aliases = m_internalPoPs.aliases( );
    for( auto alias = aliases.begin( ); alias != aliases.end( ); ++alias ) {
        a_setupInfo.m_particleSubstitution->insert( { (*alias)->pid( ), ParticleInfo( (*alias)->ID( ), a_pops, a_pops, true ) } );
    }

    m_interaction = a_node.attribute_as_string( GIDI_interactionChars );
    if( m_interaction == GIDI_MapInteractionTNSLChars ) a_targetRequiredInGlobalPoPs = false;

    Protare::initialize( a_node, a_setupInfo, a_pops, m_internalPoPs, a_targetRequiredInGlobalPoPs, a_requiredInPoPs );
    initialize( );

    m_formatVersion.setFormat( a_node.attribute_as_string( GIDI_formatChars ) );
    if( !m_formatVersion.supported( ) ) throw Exception( "unsupported GNDS format version" );
    a_setupInfo.m_formatVersion = m_formatVersion;

    if( m_formatVersion.major( ) > 1 ) {
        m_interaction = a_node.attribute_as_string( GIDI_interactionChars ); }
    else {
        HAPI::Node const firstReaction = a_node.child( GIDI_reactionsChars ).first_child( );
        if( firstReaction.attribute_as_int( GIDI_ENDF_MT_Chars ) == 502 ) m_interaction = GIDI_MapInteractionAtomicChars;
    }
    m_isPhotoAtomic = ( m_interaction == GIDI_MapInteractionAtomicChars ) && ( PoPI::IDs::photon == projectile( ).ID( ) );

    PoPI::Database const *GRIN_pops = nullptr;
    HAPI::Node const &applicationData = a_node.child( GIDI_applicationDataChars );
    std::vector<std::string> extraGammaBranchStates;
    for( HAPI::Node child1 = applicationData.first_child( ); !child1.empty( ); child1.to_next_sibling( ) ) {
        std::string nodeName( child1.name( ) );
        if( nodeName == GIDI_institutionChars ) {
            std::string label = child1.attribute_as_string( GIDI_labelChars );
            if( label == GIDI_LLNL_GRIN_continuumGammas ) {
                HAPI::Node child2 = child1.child( GIDI_GRIN_continuumGammasChars );
                if( a_construction.GRIN_continuumGammas( ) ) {
                    m_GRIN_continuumGammas = new GRIN::GRIN_continuumGammas( a_construction, child2, a_setupInfo, a_pops, m_internalPoPs, *this, &m_styles );
                    GRIN_pops = &(m_GRIN_continuumGammas->pops( ));
                    extraGammaBranchStates.push_back( m_GRIN_continuumGammas->captureResidualId( ) );
                }
            }
        }
    }
    m_internalPoPs.calculateNuclideGammaBranchStateInfos( m_nuclideGammaBranchStateInfos, GRIN_pops, extraGammaBranchStates );

    m_isTNSL_ProtareSingle = false;
    if( m_interaction == GIDI_TNSLChars ) m_interaction = GIDI_MapInteractionTNSLChars;
    if( m_interaction == GIDI_MapInteractionTNSLChars ) {
        m_isTNSL_ProtareSingle = true; }
    else {                                  // For some legacy GNDS 1.10 files.
        std::string name( a_node.child( GIDI_reactionsChars ).first_child( ).child( GIDI_doubleDifferentialCrossSectionChars ).first_child( ).name( ) );
        m_isTNSL_ProtareSingle = name.find( "thermalNeutronScatteringLaw" ) != std::string::npos;
        if( m_isTNSL_ProtareSingle ) m_interaction = GIDI_MapInteractionTNSLChars;
    }

    m_thresholdFactor = 1.0;
    if( a_pops.exists( target( ).pid( ) ) ) {
        m_thresholdFactor = 1.0 + projectile( ).mass( "amu" ) / target( ).mass( "amu" );            // BRB FIXME, I think only this statement needs to be in this if section.
    }

    m_evaluation = a_node.attribute_as_string( GIDI_evaluationChars );

    m_projectileFrame = parseFrame( a_node, a_setupInfo, GIDI_projectileFrameChars );

    m_styles.parse( a_construction, a_node.child( GIDI_stylesChars ), a_setupInfo, a_pops, m_internalPoPs, parseStylesSuite, nullptr );
    m_styles.updateChainEnds( );

    Styles::Evaluated *evaluated = m_styles.get<Styles::Evaluated>( 0 );

    m_projectileEnergyMin = evaluated->projectileEnergyDomain( ).minimum( );
    m_projectileEnergyMax = evaluated->projectileEnergyDomain( ).maximum( );

    m_reactions.parse( a_construction, a_node.child( GIDI_reactionsChars ), a_setupInfo, a_pops, m_internalPoPs, parseReaction, &m_styles );
    m_orphanProducts.parse( a_construction, a_node.child( GIDI_orphanProductsChars ), a_setupInfo, a_pops, m_internalPoPs, parseOrphanProduct, &m_styles );
    m_incompleteReactions.parse( a_construction, a_node.child( GIDI_incompleteReactionsChars ), a_setupInfo, a_pops, m_internalPoPs, parseReaction, &m_styles );

    m_sums.parse( a_construction, a_node.child( GIDI_sumsChars ), a_setupInfo, a_pops, m_internalPoPs );
    m_fissionComponents.parse( a_construction, a_node.child( GIDI_fissionComponentsChars ), a_setupInfo, a_pops, m_internalPoPs, parseFissionComponent, &m_styles );

    m_RutherfordScatteringPresent = false;
    m_onlyRutherfordScatteringPresent = false;
    for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
        Reaction const *reaction1 = m_reactions.get<Reaction>( i1 );
        if( reaction1->RutherfordScatteringPresent( ) ) {
            m_RutherfordScatteringPresent = true;
            m_onlyRutherfordScatteringPresent = reaction1->onlyRutherfordScatteringPresent( );
            break;
        }
    }

    for( HAPI::Node child1 = applicationData.first_child( ); !child1.empty( ); child1.to_next_sibling( ) ) {
        std::string nodeName( child1.name( ) );

        if( nodeName == GIDI_institutionChars ) {
            std::string label = child1.attribute_as_string( GIDI_labelChars );
            if( label == GIDI_LLNL_Chars ) {
                for( HAPI::Node child2 = child1.first_child( ); !child2.empty( ); child2.to_next_sibling( ) ) {
                    if( child2.name( ) == GIDI_nuclearPlusCoulombInterferenceChars ) {
                        HAPI::Node const reactionNode = child2.child( GIDI_reactionChars );
                        a_setupInfo.m_isENDL_C_9 = true;
                        m_nuclearPlusCoulombInterferenceOnlyReaction = new Reaction( a_construction, reactionNode, a_setupInfo, a_pops, m_internalPoPs, *this, &m_styles );
                        a_setupInfo.m_isENDL_C_9 = false;
                        bool dropC_9 = false;
                        auto &crossSectionSuite = m_nuclearPlusCoulombInterferenceOnlyReaction->crossSection( );
                        for( auto crossSectionIter = crossSectionSuite.begin( ); crossSectionIter != crossSectionSuite.end( ); ++crossSectionIter ) {
                            auto iter = crossSectionSuite.checkLazyParsingHelperFormIterator( crossSectionIter );
                            if( (*iter)->type( ) == FormType::XYs1d ) {
                                Functions::XYs1d *xys1d = static_cast<Functions::XYs1d *>( *iter );
                                auto ys = xys1d->ys( );
                                for( auto yIter = ys.begin( ); yIter != ys.end( ); ++yIter ) {
                                    if( *yIter < 0.0 ) {
                                        dropC_9 = true;
                                        break;
                                    }
                                }
                                break;
                            }
                        }
                        if( dropC_9 ) {
                            delete m_nuclearPlusCoulombInterferenceOnlyReaction;
                            m_nuclearPlusCoulombInterferenceOnlyReaction = nullptr;
                        }
                    }
                } }
            else if( label == GIDI_LLNL_multiGroupReactions_Chars ) {
                HAPI::Node child2 = child1.child( GIDI_reactionChars );
                m_multiGroupSummedReaction = new Reaction( a_construction, child2, a_setupInfo, a_pops, m_internalPoPs, *this, &m_styles ); }
            else if( label == GIDI_LLNL_multiGroupDelayedNeutrons_Chars ) {
                HAPI::Node child2 = child1.child( GIDI_outputChannelChars );
                m_multiGroupSummedDelayedNeutrons = new OutputChannel( a_construction, child2, a_setupInfo, a_pops, m_internalPoPs, &m_styles, true, false ); }
            else if( label == GIDI_LLNL_URR_probability_tables_Chars ) {
                m_ACE_URR_probabilityTables.parse( a_construction, child1.child( GIDI_ACE_URR_probabilityTablesChars ), a_setupInfo, a_pops, 
                        m_internalPoPs, parseACE_URR_probabilityTables, &m_styles ); }
            else if( label == GIDI_LLNL_photoAtomicIncoherentDoppler_Chars ) {
                HAPI::Node child2 = child1.child( GIDI_reactionsChars );
                m_photoAtomicIncoherentDoppler.parse( a_construction, child2, a_setupInfo, a_pops, m_internalPoPs, parseReaction, &m_styles );
                if( a_construction.usePhotoAtomicIncoherentDoppler( ) && m_photoAtomicIncoherentDoppler.size( ) > 0 ) {
                                                                                        // Add the consistent doppler broadened incoherent reactions to the list.
                    while( m_photoAtomicIncoherentDoppler.size( ) > 0 ) {
                        Reaction *photoAtomicIncoherentDopplerReaction = m_photoAtomicIncoherentDoppler.pop<Reaction>( 0 );
                        m_reactions.add( photoAtomicIncoherentDopplerReaction );
                    }
                    for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {         // Set standard incoherent reaction to inactive
                        Reaction *reaction1 = m_reactions.get<Reaction>( i1 );
                        if( reaction1->ENDF_MT( ) == 504 ) {
                            reaction1->setActive(false);
                        }
                    }
                }
            }
            else if( label == GIDI_LLNL_pointwiseAverageProductEnergies ) {
                m_pointwiseAverageProductEnergy.parse( a_construction, child1.child( GIDI_averageEnergyChars ), a_setupInfo, a_pops, 
                        m_internalPoPs, parseAverageEnergySuite, &m_styles ); }
            else if( label == GIDI_LLNL_GRIN_continuumGammas ) {        // Already parsed above.
                continue;
            } }
        else {
            std::cout << "parseStylesSuite: Ignoring unsupported style = '" << nodeName << "'." << std::endl;
        }
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

ProtareSingle::~ProtareSingle( ) {

    delete m_doc;
    delete m_dataManager;
    delete m_nuclearPlusCoulombInterferenceOnlyReaction;
    delete m_multiGroupSummedReaction;
    delete m_multiGroupSummedDelayedNeutrons;
}

/* *********************************************************************************************************//**
 * Checks various things to determine if it is okay to use summed multi-group data or not.
 *
 * @param a_settings        [in]    Specifies user options.
 *
 * @return                          Returns **true** if summed data are to be returned and **false** otherwise.
 ***********************************************************************************************************/

bool ProtareSingle::useMultiGroupSummedData( Transporting::MG const &a_settings, ExcludeReactionsSet const &a_reactionsToExclude ) const {

    if( ( numberOfInactiveReactions( ) > 0 ) || ( a_reactionsToExclude.size( ) > 0 ) ) return( false );
    if( m_multiGroupSummedReaction == nullptr ) return( false );
    if( a_settings.nuclearPlusCoulombInterferenceOnly( ) && m_RutherfordScatteringPresent ) return( false );
    if( m_decayPositronium ) return( false );                   // Need to handle decay of positronium by reaction.

    return( a_settings.useMultiGroupSummedData( ) );
}

/* *********************************************************************************************************//**
 * Returns **true** if delayed neutrons are to be included and **false** otherwise.
 *
 * @param a_settings        [in]    Specifies user options.
 *
 * @return                          Returns **true** if delayed neutrons are to be included and **false** otherwise.
 ***********************************************************************************************************/

bool ProtareSingle::useMultiGroupSummedDelayedNeutronsData( Transporting::MG const &a_settings ) const {

    return( ( a_settings.delayedNeutrons( ) == Transporting::DelayedNeutrons::on ) && ( m_multiGroupSummedDelayedNeutrons != nullptr ) );
}

/* *********************************************************************************************************//**
 * Returns *a_reaction* except when Rutherford scattering present and only nuclear + Coulomb interference is wanted. Otherwise
 * returns *m_nuclearPlusCoulombInterferenceOnlyReaction* which may be a **nullptr**.
 *
 * @param a_settings        [in]    Specifies user options.
 * @param a_reaction        [in]    Reaction pointer to return if check passes.
 *
 * @return                          Const reaction pointer or *m_nuclearPlusCoulombInterferenceOnlyReaction*.
 ***********************************************************************************************************/

Reaction const *ProtareSingle::checkIf_nuclearPlusCoulombInterferenceWanted( Transporting::MG const &a_settings, Reaction const *a_reaction ) const {

    if( a_reaction->RutherfordScatteringPresent( ) && a_settings.nuclearPlusCoulombInterferenceOnly( ) ) {
        return( m_nuclearPlusCoulombInterferenceOnlyReaction );
    }

    return( a_reaction );
}

/* *********************************************************************************************************//**
 * Returns *a_reaction* except when Rutherford scattering present and only nuclear + Coulomb interference is wanted. Otherwise
 * returns *m_nuclearPlusCoulombInterferenceOnlyReaction* which may be a **nullptr**.
 *
 * @param a_settings        [in]    Specifies user options.
 * @param a_reaction        [in]    Reaction pointer to return if check passes.
 *
 * @return                          Const reaction pointer or *m_nuclearPlusCoulombInterferenceOnlyReaction*.
 ***********************************************************************************************************/

Reaction const *ProtareSingle::reactionToMultiGroup( Transporting::MG const &a_settings, std::size_t a_index, 
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Reaction const *reaction1 = m_reactions.get<Reaction>( a_index );

    if( !reaction1->active( ) ) return( nullptr );
    if( a_reactionsToExclude.find( a_index ) != a_reactionsToExclude.end( ) ) return( nullptr );

    return( checkIf_nuclearPlusCoulombInterferenceWanted( a_settings, reaction1 ) );
}

/* *********************************************************************************************************//**
 * Returns the pointer representing the protare (i.e., *this*) if *a_index* is 0 and nullptr otherwise.
 *
 * @param a_index               [in]    Must always be 0.
 *
 * @return                              Returns the pointer representing *this*.
 ***********************************************************************************************************/

ProtareSingle *ProtareSingle::protare( std::size_t a_index ) {

    if( a_index != 0 ) return( nullptr );
    return( this );
}

/* *********************************************************************************************************//**
 * Returns the pointer representing the protare (i.e., *this*) if *a_index* is 0 and nullptr otherwise.
 *
 * @param a_index               [in]    Must always be 0.
 *
 * @return                              Returns the pointer representing *this*.
 ***********************************************************************************************************/

ProtareSingle const *ProtareSingle::protare( std::size_t a_index ) const {

    if( a_index != 0 ) return( nullptr );
    return( this );
}

/* *********************************************************************************************************//**
 * Returns the intid for the requested particle or -1 if the particle is not in *this* PoPs database.
 *
 * @param a_id                 [in]    The GNDS PoPs id for particle whose intd is requested.
 *
 * @return                             C++ int for the requested particle or -1 if particle is not in PoPs.
 ******************************************************************/

int ProtareSingle::intid( std::string const &a_id ) const {

    return( m_internalPoPs.intid( a_id ) );
}

/* *********************************************************************************************************//**
 * Fills in a std::set with a unique list of all product indices produced by reactions and orphanProducts. 
 * If a_transportablesOnly is true, only transportable product incides are return.
 *
 * @param a_ids                 [out]   The unique list of product indices.
 * @param a_particles           [in]    The list of particles to be transported.
 * @param a_transportablesOnly  [in]    If true, only transportable product indices are added in the list.
 ***********************************************************************************************************/

void ProtareSingle::productIDs( std::set<std::string> &a_ids, Transporting::Particles const &a_particles, bool a_transportablesOnly ) const {

    for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
        Reaction const *reaction1 = m_reactions.get<Reaction>( i1 );

        if( !reaction1->active( ) ) continue;
        reaction1->productIDs( a_ids, a_particles, a_transportablesOnly );
    }

    for( std::size_t i1 = 0; i1 < m_orphanProducts.size( ); ++i1 ) {
        Reaction const *reaction1 = m_orphanProducts.get<Reaction>( i1 );

        if( !reaction1->active( ) ) continue;
        reaction1->productIDs( a_ids, a_particles, a_transportablesOnly );
    }
}

/* *********************************************************************************************************//**
 * Determines the maximum Legredre order present in the multi-group transfer matrix for a give product for a give label.
 *
 * @param a_smr             [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings        [in]    Specifies the requested label.
 * @param a_temperatureInfo [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID       [in]    The id of the requested product.
 *
 * @return                          The maximum Legredre order. If no transfer matrix data are present for the requested product, -1 is returned.
 ***********************************************************************************************************/

int ProtareSingle::maximumLegendreOrder( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    int _maximumLegendreOrder = -1;
    ExcludeReactionsSet excludeReactionsSet;

    if( useMultiGroupSummedData( a_settings, excludeReactionsSet ) ) {
        _maximumLegendreOrder = m_multiGroupSummedReaction->maximumLegendreOrder( a_smr, a_settings, a_temperatureInfo, a_productID );
        if( useMultiGroupSummedDelayedNeutronsData( a_settings ) ) {
            _maximumLegendreOrder = std::max( _maximumLegendreOrder, m_multiGroupSummedDelayedNeutrons->maximumLegendreOrder( 
                    a_smr, a_settings, a_temperatureInfo, a_productID ) );
        } }
    else {
        for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
            Reaction const *reaction1 = m_reactions.get<Reaction>( i1 );

            if( !reaction1->active( ) ) continue;
            int r_maximumLegendreOrder = reaction1->maximumLegendreOrder( a_smr, a_settings, a_temperatureInfo, a_productID );
            if( r_maximumLegendreOrder > _maximumLegendreOrder ) _maximumLegendreOrder = r_maximumLegendreOrder;
        }

        for( std::size_t i1 = 0; i1 < m_orphanProducts.size( ); ++i1 ) {
            Reaction const *reaction1 = m_orphanProducts.get<Reaction>( i1 );

            if( !reaction1->active( ) ) continue;
            int r_maximumLegendreOrder = reaction1->maximumLegendreOrder( a_smr, a_settings, a_temperatureInfo, a_productID );
            if( r_maximumLegendreOrder > _maximumLegendreOrder ) _maximumLegendreOrder = r_maximumLegendreOrder;
        }
    }

    return( _maximumLegendreOrder );
}

/* *********************************************************************************************************//**
 * Returns a list of all process temperature data. For each temeprature, the labels for its 
 *
 *   + heated cross section data,
 *   + gridded cross section data,
 *   + multi-group data, and
 *   + multi-group upscatter data.
 *
 * are returned. If no data are present for a give data type (e.g., gridded cross section, multi-group upscatter), its label is an empty std::string.
 *
 * @return  The list of temperatures and their labels via an Styles::TemperatureInfos instance. The Styles::TemperatureInfos class
 *          has many (if not all) the method of a std::vector.
 ***********************************************************************************************************/

Styles::TemperatureInfos ProtareSingle::temperatures( ) const {

    std::size_t size( m_styles.size( ) );
    Styles::TemperatureInfos temperature_infos;

    for( std::size_t i1 = 0; i1 < size; ++i1 ) {
        Styles::Base const *style1 = m_styles.get<Styles::Base>( i1 );

        if( style1->moniker( ) == GIDI_heatedStyleChars ) {
            PhysicalQuantity const &temperature = style1->temperature( );
            std::string heated_cross_section( style1->label( ) );
            std::string gridded_cross_section( "" );
            std::string URR_probability_tables( "" );
            std::string heated_multi_group( "" );
            std::string Sn_elastic_upscatter( "" );

            for( std::size_t i2 = 0; i2 < size; ++i2 ) {
                Styles::Base const *style2 = m_styles.get<Styles::Base>( i2 );

                if( style2->moniker( ) == GIDI_multiGroupStyleChars ) continue;
                if( style2->temperature( ).value( ) != temperature.value( ) ) continue;

                if( style2->moniker( ) == GIDI_griddedCrossSectionStyleChars ) {
                    gridded_cross_section = style2->label( ); }
                else if( style2->moniker( ) == GIDI_URR_probabilityTablesStyleChars ) {
                    URR_probability_tables = style2->label( ); }
                else if( style2->moniker( ) == GIDI_SnElasticUpScatterStyleChars ) {
                    Sn_elastic_upscatter = style2->label( ); }
                else if( style2->moniker( ) == GIDI_heatedMultiGroupStyleChars ) {
                    heated_multi_group = style2->label( );
                }
            }
            temperature_infos.push_back( Styles::TemperatureInfo( temperature, heated_cross_section, gridded_cross_section, URR_probability_tables,
                    heated_multi_group, Sn_elastic_upscatter ) );
        }
    }

    std::sort( temperature_infos.begin( ), temperature_infos.end( ), sortTemperatures );

    return( temperature_infos );
}

/* *********************************************************************************************************//**
 * FOR INTERNAL USE ONLY.
 *
 * Determines if the temperature of lhs is less than that of rhs, or not.
 *
 * @param lhs   [in]
 * @param rhs   [in]
 * @return      true if temperature of lhs < rhs and false otherwise.
 ***********************************************************************************************************/

bool sortTemperatures( Styles::TemperatureInfo const &lhs, Styles::TemperatureInfo const &rhs ) {

    return( lhs.temperature( ).value( ) < rhs.temperature( ).value( ) );
}

/* *********************************************************************************************************//**
 * Returns the (*a_index*+1)th reaction or **nullptr**. If the indexed reaction is deactivated or exlucded, 
 * a **nullptr** is returned.
 *
 * @param a_index               [in]    The index of the requested reaction.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the cross section.
 *
 * @return                          The (*a_index*+1)th reaction or **nullptr**.
 ***********************************************************************************************************/

Reaction const *ProtareSingle::reaction( std::size_t a_index, Transporting::MG const &a_settings,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    return( reactionToMultiGroup( a_settings, a_index, a_reactionsToExclude ) );
}

/* *********************************************************************************************************//**
 * The method returns the number of reactions of *this* that have been deactivated.
 *
 * @return              The number of deactivated reaction of *this*.
 ***********************************************************************************************************/

std::size_t ProtareSingle::numberOfInactiveReactions( ) const {

    std::size_t count = 0;

    for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
        Reaction const *reaction1 = m_reactions.get<Reaction>( i1 );

        if( !reaction1->active( ) ) ++count;
    }

    return( count );
}

/* *********************************************************************************************************//**
 * Re-indexs the reactions in the reactions, orphanProducts and fissionComponents suites.
 *
 ***********************************************************************************************************/

void ProtareSingle::updateReactionIndices( std::size_t a_offset ) const {

    for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
        Reaction const *reaction1 = m_reactions.get<Reaction>( i1 + a_offset );

        reaction1->setReactionIndex( i1 );
    }
    for( std::size_t i1 = 0; i1 < m_orphanProducts.size( ); ++i1 ) {
        Reaction const *reaction1 = m_orphanProducts.get<Reaction>( i1 );
    
        reaction1->setReactionIndex( i1 );
    }
    for( std::size_t i1 = 0; i1 < m_orphanProducts.size( ); ++i1 ) {
        Reaction const *reaction1 = m_orphanProducts.get<Reaction>( i1 );
    
        reaction1->setReactionIndex( i1 );
    }
}

/* *********************************************************************************************************//**
 * Returns the multi-group boundaries for the requested label and product.
 *
 * @param a_settings        [in]    Specifies the requested label.
 * @param a_temperatureInfo [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID       [in]    ID for the requested product.
 *
 * @return                          List of multi-group boundaries.
 ***********************************************************************************************************/

std::vector<double> ProtareSingle::groupBoundaries( LUPI_maybeUnused Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    Styles::HeatedMultiGroup const *heatedMultiGroupStyle1 = m_styles.get<Styles::HeatedMultiGroup>( a_temperatureInfo.heatedMultiGroup( ) );

    return( heatedMultiGroupStyle1->groupBoundaries( a_productID ) );
}

/* *********************************************************************************************************//**
 * Returns the inverse speeds for the requested label. The label must be for a heated multi-group style.
 *
 * @param a_smr             [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings        [in]    Specifies the requested label.
 * @param a_temperatureInfo [in]    Specifies the temperature and labels use to lookup the requested data.
 *
 * @return                          List of inverse speeds.
 ***********************************************************************************************************/

Vector ProtareSingle::multiGroupInverseSpeed( LUPI_maybeUnused LUPI::StatusMessageReporting &a_smr, LUPI_maybeUnused Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo ) const {

    Styles::HeatedMultiGroup const *heatedMultiGroupStyle1 = m_styles.get<Styles::HeatedMultiGroup>( a_temperatureInfo.heatedMultiGroup( ) );

    return( heatedMultiGroupStyle1->inverseSpeedData( ) );
}

/* *********************************************************************************************************//**
 * Returns true if at least one reaction contains a fission channel.
 *
 * @return  true if at least one reaction contains a fission channel.
 ***********************************************************************************************************/

bool ProtareSingle::hasFission( ) const {

    for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
        Reaction const *reaction1 = m_reactions.get<Reaction>( i1 );

        if( !reaction1->active( ) ) continue;
        if( reaction1->hasFission( ) ) return( true );
    }
    return( false );
}

/* *********************************************************************************************************//**
 * Returns **false* if protare has delayed fission neutrons for an active reaction and they are not complete; otherwise, returns **true**.
 *
 * @return      bool
 ***********************************************************************************************************/

bool ProtareSingle::isDelayedFissionNeutronComplete( ) const {

    for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
        Reaction const *reaction1 = m_reactions.get<Reaction>( i1 );

        if( !reaction1->active( ) ) continue;
        if( reaction1->hasFission( ) ) {
            OutputChannel const *outputChannel = reaction1->outputChannel( );
            if( !outputChannel->isDelayedFissionNeutronComplete( ) ) return( false );
        }
    }

    return( true );
}

/* *********************************************************************************************************//**
 * Used by GUPI::Ancestry to tranverse GNDS nodes. This method returns a pointer to a derived class' a_item member or nullptr if none exists.
 *
 * @param a_item    [in]    The name of the class member whose pointer is to be return.
 * @return                  The pointer to the class member or nullptr if class does not have a member named a_item.
 ***********************************************************************************************************/

GUPI::Ancestry *ProtareSingle::findInAncestry3( std::string const &a_item ) {

    if( a_item == GIDI_stylesChars ) return( &m_styles );
    if( a_item == GIDI_reactionsChars ) return( &m_reactions );
    if( a_item == GIDI_orphanProductsChars ) return( &m_orphanProducts );
    if( a_item == GIDI_sumsChars ) return( &m_sums );
    if( a_item == GIDI_fissionComponentsChars ) return( &m_fissionComponents );
    if( a_item == GIDI_ACE_URR_probabilityTablesChars ) return( &m_ACE_URR_probabilityTables );

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Used by GUPI::Ancestry to tranverse GNDS nodes. This method returns a pointer to a derived class' a_item member or nullptr if none exists.
 *
 * @param a_item    [in]    The name of the class member whose pointer is to be return.
 * @return                  The pointer to the class member or nullptr if class does not have a member named a_item.
 ***********************************************************************************************************/

GUPI::Ancestry const *ProtareSingle::findInAncestry3( std::string const &a_item ) const {

    if( a_item == GIDI_stylesChars ) return( &m_styles );
    if( a_item == GIDI_reactionsChars ) return( &m_reactions );
    if( a_item == GIDI_orphanProductsChars ) return( &m_orphanProducts );
    if( a_item == GIDI_sumsChars ) return( &m_sums );
    if( a_item == GIDI_fissionComponentsChars ) return( &m_fissionComponents );
    if( a_item == GIDI_ACE_URR_probabilityTablesChars ) return( &m_ACE_URR_probabilityTables );

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total cross section for the requested label. This is summed over all reactions.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the cross section.
 * @param a_label               [in]    If not an empty string, this is used as the label for the form to return and the *a_temperatureInfo* labels are ignored.
 *
 * @return                          The requested multi-group cross section as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareSingle::multiGroupCrossSection( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude, std::string const &a_label ) const {

    Vector vector;

    if( useMultiGroupSummedData( a_settings, a_reactionsToExclude ) ) {
        vector = m_multiGroupSummedReaction->multiGroupCrossSection( a_smr, a_settings, a_temperatureInfo, a_label ); }
    else {
        for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
            Reaction const *reaction1 = reactionToMultiGroup( a_settings, i1, a_reactionsToExclude );

            if( reaction1 != nullptr ) vector += reaction1->multiGroupCrossSection( a_smr, a_settings, a_temperatureInfo, a_label );
        }
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total multiplicity for the requested label for the requested product. This is a cross section weighted multiplicity.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID           [in]    Id for the requested product.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the multiplicity.
 *
 * @return                          The requested multi-group multiplicity as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareSingle::multiGroupMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector;

    if( useMultiGroupSummedData( a_settings, a_reactionsToExclude ) ) {
        vector = m_multiGroupSummedReaction->multiGroupMultiplicity( a_smr, a_settings, a_temperatureInfo, a_productID );
        if( useMultiGroupSummedDelayedNeutronsData( a_settings ) ) {
            vector += m_multiGroupSummedDelayedNeutrons->multiGroupMultiplicity( a_smr, a_settings, a_temperatureInfo, a_productID );
        } }
    else {
        for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
            Reaction const *reaction1 = reactionToMultiGroup( a_settings, i1, a_reactionsToExclude );

            if( reaction1 != nullptr ) vector += reaction1->multiGroupMultiplicity( a_smr, a_settings, a_temperatureInfo, a_productID );
        }

        for( std::size_t i1 = 0; i1 < m_orphanProducts.size( ); ++i1 ) {
            Reaction const *reaction1 = m_orphanProducts.get<Reaction>( i1 );

            if( !reaction1->active( ) ) continue;
            vector += reaction1->multiGroupMultiplicity( a_smr, a_settings, a_temperatureInfo, a_productID );
        }
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total fission neutron multiplicity for the requested label. This is a cross section weighted multiplicity.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the cross multiplicity.
 *
 * @return                          The requested multi-group fission neutron multiplicity as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareSingle::multiGroupFissionNeutronMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, LUPI_maybeUnused ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector( 0 );

    for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
        Reaction const *reaction1 = m_reactions.get<Reaction>( i1 );

        if( !reaction1->active( ) ) continue;
        if( reaction1->hasFission( ) ) vector += reaction1->multiGroupMultiplicity( a_smr, a_settings, a_temperatureInfo, PoPI::IDs::neutron );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total fission gamma multiplicity for the requested label. This is a cross section weighted multiplicity.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the multiplicity.
 *
 * @return                          The requested multi-group fission neutron multiplicity as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareSingle::multiGroupFissionGammaMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, LUPI_maybeUnused ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector( 0 );

    for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
        Reaction const *reaction1 = m_reactions.get<Reaction>( i1 );

        if( !reaction1->active( ) ) continue;
        if( reaction1->hasFission( ) ) vector += reaction1->multiGroupMultiplicity( a_smr, a_settings, a_temperatureInfo, PoPI::IDs::photon );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total Q for the requested label. This is a cross section weighted Q
 * summed over all reactions
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_final               [in]    If false, only the Q for the primary reactions are return, otherwise, the Q for the final reactions.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the Q.
 *
 * @return                          The requested multi-group Q as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareSingle::multiGroupQ( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, bool a_final, LUPI_maybeUnused bool a_effectivePhotoAtomic,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector( 0 );

    if( useMultiGroupSummedData( a_settings, a_reactionsToExclude ) ) {
        vector = m_multiGroupSummedReaction->multiGroupQ( a_smr, a_settings, a_temperatureInfo, a_final );
        if( useMultiGroupSummedDelayedNeutronsData( a_settings ) ) {
            vector += m_multiGroupSummedDelayedNeutrons->multiGroupQ( a_smr, a_settings, a_temperatureInfo, a_final );
        } }
    else {
        for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
            Reaction const *reaction1 = reactionToMultiGroup( a_settings, i1, a_reactionsToExclude );

            if( reaction1 != nullptr ) vector += reaction1->multiGroupQ( a_smr, a_settings, a_temperatureInfo, a_final );
        }
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total product matrix for the requested label for the requested product id for the requested Legendre order.
 * If no data are found, an empty GIDI::Matrix is returned.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label and if delayed neutrons should be included.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_particles           [in]    The list of particles to be transported.
 * @param a_productID           [in]    PoPs id for the requested product.
 * @param a_order               [in]    Requested product matrix, Legendre order.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the product matrix.
 *
 * @return                          The requested multi-group product matrix as a GIDI::Matrix.
 ***********************************************************************************************************/

Matrix ProtareSingle::multiGroupProductMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, 
                std::string const &a_productID, std::size_t a_order, ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Matrix matrix( 0, 0 );

    if( useMultiGroupSummedData( a_settings, a_reactionsToExclude ) ) {
        matrix = m_multiGroupSummedReaction->multiGroupProductMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, a_productID, a_order );
        if( useMultiGroupSummedDelayedNeutronsData( a_settings ) ) {
            matrix += m_multiGroupSummedDelayedNeutrons->multiGroupProductMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, a_productID, a_order );
        } }
    else {
        for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
            Reaction const *reaction1 = reactionToMultiGroup( a_settings, i1, a_reactionsToExclude );

            if( reaction1 != nullptr ) matrix += reaction1->multiGroupProductMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, a_productID, a_order );
        }

        for( std::size_t i1 = 0; i1 < m_orphanProducts.size( ); ++i1 ) {
            Reaction const *reaction1 = m_orphanProducts.get<Reaction>( i1 );

            if( !reaction1->active( ) ) continue;
            matrix += reaction1->multiGroupProductMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, a_productID, a_order );
        }
    }

    return( matrix );
}

/* *********************************************************************************************************//**
 * Like ProtareSingle::multiGroupProductMatrix, but only returns the fission neutron, transfer matrix.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label and if delayed neutrons should be included.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_particles           [in]    The list of particles to be transported.
 * @param a_order               [in]    Requested product matrix, Legendre order.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the fission matrix.
 *
 * @return                          The requested multi-group neutron fission matrix as a GIDI::Matrix.
 ***********************************************************************************************************/

Matrix ProtareSingle::multiGroupFissionMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, std::size_t a_order,
                LUPI_maybeUnused ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Matrix matrix( 0, 0 );

    for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
        Reaction const *reaction1 = m_reactions.get<Reaction>( i1 );

        if( !reaction1->active( ) ) continue;
        if( reaction1->hasFission( ) ) matrix += reaction1->multiGroupFissionMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, a_order );
    }

    return( matrix );
}

/* *********************************************************************************************************//**
 * Returns the multi-group transport correction for the requested label. The transport correction is calculated from the transfer matrix
 * for the projectile id for the Legendre order of *a_order + 1*.
 *
 * @param a_smr                     [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings                [in]    Specifies the requested label.
 * @param a_temperatureInfo         [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_particles               [in]    The list of particles to be transported.
 * @param a_order                   [in]    Maximum Legendre order for transport. The returned transport correction is for the next higher Legender order.
 * @param a_transportCorrectionType [in]    Requested transport correction type.
 * @param a_temperature             [in]    The temperature of the flux to use when collapsing. Pass to the GIDI::collapse method.
 * @param a_reactionsToExclude      [in]    A list of reaction indices that are to be ignored when calculating the transport correction.
 *
 * @return                                  The requested multi-group transport correction as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareSingle::multiGroupTransportCorrection( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, std::size_t a_order, 
                TransportCorrectionType a_transportCorrectionType, double a_temperature, LUPI_maybeUnused ExcludeReactionsSet const &a_reactionsToExclude ) const {

    if( a_transportCorrectionType == TransportCorrectionType::None ) return( Vector( 0 ) );

    Matrix matrix( multiGroupProductMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, projectile( ).ID( ), a_order + 1 ) );
    Matrix matrixCollapsed = collapse( matrix, a_settings, a_particles, a_temperature, projectile( ).ID( ) );
    std::size_t size = matrixCollapsed.size( );
    std::vector<double> transportCorrection1( size, 0 );

    if( a_transportCorrectionType == TransportCorrectionType::None ) {
        }
    else if( a_transportCorrectionType == TransportCorrectionType::Pendlebury ) {
        for( std::size_t index = 0; index < size; ++index ) transportCorrection1[index] = matrixCollapsed[index][index]; }
    else {
        throw Exception( "Unsupported transport correction: only None and Pendlebury (i.e., Pendlebury/Underhill) are currently supported." );
    }
    return( Vector( transportCorrection1 ) );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total available energy for the requested label. This is a cross section weighted available energy
 * summed over all reactions.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the available energy.
 *
 * @return                          The requested multi-group available energy as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareSingle::multiGroupAvailableEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector( 0 );

    if( useMultiGroupSummedData( a_settings, a_reactionsToExclude ) ) {
        vector = m_multiGroupSummedReaction->multiGroupAvailableEnergy( a_smr, a_settings, a_temperatureInfo ); }
    else {
        for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
            Reaction const *reaction1 = reactionToMultiGroup( a_settings, i1, a_reactionsToExclude );

            if( reaction1 != nullptr ) vector += reaction1->multiGroupAvailableEnergy( a_smr, a_settings, a_temperatureInfo );
        }
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total average energy for the requested label for the requested product. This is a cross section weighted average energy
 * summed over all reactions.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID           [in]    Particle id for the requested product.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the average energy.
 *
 * @return                          The requested multi-group average energy as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareSingle::multiGroupAverageEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector( 0 );

    if( useMultiGroupSummedData( a_settings, a_reactionsToExclude ) ) {
        vector = m_multiGroupSummedReaction->multiGroupAverageEnergy( a_smr, a_settings, a_temperatureInfo, a_productID );
        if( useMultiGroupSummedDelayedNeutronsData( a_settings ) ) {
            vector += m_multiGroupSummedDelayedNeutrons->multiGroupAverageEnergy( a_smr, a_settings, a_temperatureInfo, a_productID );
        } }
    else {
        for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
            Reaction const *reaction1 = reactionToMultiGroup( a_settings, i1, a_reactionsToExclude );

            if( reaction1 != nullptr ) vector += reaction1->multiGroupAverageEnergy( a_smr, a_settings, a_temperatureInfo, a_productID );
        }

        for( std::size_t i1 = 0; i1 < m_orphanProducts.size( ); ++i1 ) {
            Reaction const *reaction1 = m_orphanProducts.get<Reaction>( i1 );

            if( !reaction1->active( ) ) continue;
            vector += reaction1->multiGroupAverageEnergy( a_smr, a_settings, a_temperatureInfo, a_productID );
        }
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total deposition energy for the requested label. This is a cross section weighted deposition energy
 * summed over all reactions. The deposition energy is calculated by subtracting the average energy from each transportable particle
 * from the available energy. The list of transportable particles is specified via the list of particle specified in the *a_settings* argument.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label and the products that are transported.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_particles           [in]    The list of particles to be transported.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the deposition energy.
 *
 * @return                      The requested multi-group deposition energy as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareSingle::multiGroupDepositionEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    bool atLeastOneReactionHasAllParticlesTracked = false;
    Vector vector;

    if( a_settings.zeroDepositionIfAllProductsTracked( ) ) {
        for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
            Reaction const *reaction1 = reactionToMultiGroup( a_settings, i1, a_reactionsToExclude );

            if( reaction1 != nullptr ) {
                if( reaction1->areAllProductsTracked( a_particles ) ) {
                    atLeastOneReactionHasAllParticlesTracked = true;
                    break;
                }
            }
        }
    }

    if( atLeastOneReactionHasAllParticlesTracked ) {
        for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
            Reaction const *reaction1 = reactionToMultiGroup( a_settings, i1, a_reactionsToExclude );

            if( reaction1 != nullptr ) vector += reaction1->multiGroupDepositionEnergy( a_smr, a_settings, a_temperatureInfo, a_particles );
        }
        for( std::size_t i1 = 0; i1 < m_orphanProducts.size( ); ++i1 ) {
            Reaction const *reaction1 = m_orphanProducts.get<Reaction>( i1 );

            if( !reaction1->active( ) ) continue;
            vector += reaction1->multiGroupDepositionEnergy( a_smr, a_settings, a_temperatureInfo, a_particles );
        } }
    else {
        std::map<std::string, Transporting::Particle> const &products( a_particles.particles( ) );
        vector = multiGroupAvailableEnergy( a_smr, a_settings, a_temperatureInfo, a_reactionsToExclude );
        Vector availableEnergy( vector );

        for( std::map<std::string, Transporting::Particle>::const_iterator iter = products.begin( ); iter != products.end( ); ++iter ) {
            vector -= multiGroupAverageEnergy( a_smr, a_settings, a_temperatureInfo, iter->first, a_reactionsToExclude );
        }

        for( std::size_t index = 0; index < availableEnergy.size( ); ++index ) {        // Check for values that should probably be 0.0.
            if( std::fabs( vector[index] ) < std::fabs( 1e-14 * availableEnergy[index] ) ) vector[index] = 0.0;
        }
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total available momentum for the requested label. This is a cross section weighted available momentum
 * summed over all reactions.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the available momentum.
 *
 * @return                          The requested multi-group available momentum as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareSingle::multiGroupAvailableMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector( 0 );

    if( useMultiGroupSummedData( a_settings, a_reactionsToExclude ) ) {
        vector = m_multiGroupSummedReaction->multiGroupAvailableMomentum( a_smr, a_settings, a_temperatureInfo ); }
    else {
        for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
            Reaction const *reaction1 = reactionToMultiGroup( a_settings, i1, a_reactionsToExclude );

            if( reaction1 != nullptr ) vector += reaction1->multiGroupAvailableMomentum( a_smr, a_settings, a_temperatureInfo );
        }
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total average momentum for the requested label for the requested product. This is a cross section weighted average momentum
 * summed over all reactions.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID           [in]    Particle id for the requested product.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the average momentum.
 *
 * @return                          The requested multi-group average momentum as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareSingle::multiGroupAverageMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector( 0 );

    if( useMultiGroupSummedData( a_settings, a_reactionsToExclude ) ) {
        vector = m_multiGroupSummedReaction->multiGroupAverageMomentum( a_smr, a_settings, a_temperatureInfo, a_productID );
        if( useMultiGroupSummedDelayedNeutronsData( a_settings ) ) {
            vector += m_multiGroupSummedDelayedNeutrons->multiGroupAverageMomentum( a_smr, a_settings, a_temperatureInfo, a_productID );
        } }
    else {
        for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
            Reaction const *reaction1 = reactionToMultiGroup( a_settings, i1, a_reactionsToExclude );

            if( reaction1 != nullptr ) vector += reaction1->multiGroupAverageMomentum( a_smr, a_settings, a_temperatureInfo, a_productID );
        }

        for( std::size_t i1 = 0; i1 < m_orphanProducts.size( ); ++i1 ) {
            Reaction const *reaction1 = m_orphanProducts.get<Reaction>( i1 );

            if( !reaction1->active( ) ) continue;
            vector += reaction1->multiGroupAverageMomentum( a_smr, a_settings, a_temperatureInfo, a_productID );
        }
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total deposition momentum for the requested label. This is a cross section weighted deposition momentum
 * summed over all reactions. The deposition momentum is calculated by subtracting the average momentum from each transportable particle
 * from the available momentum. The list of transportable particles is specified via the list of particle specified in the *a_settings* argument.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_particles           [in]    The list of particles to be transported.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the deposition momentum.
 *
 * @return                          The requested multi-group deposition momentum as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareSingle::multiGroupDepositionMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    std::map<std::string, Transporting::Particle> const &products( a_particles.particles( ) );
    Vector vector = multiGroupAvailableMomentum( a_smr, a_settings, a_temperatureInfo, a_reactionsToExclude );

    for( std::map<std::string, Transporting::Particle>::const_iterator iter = products.begin( ); iter != products.end( ); ++iter ) {
        vector -= multiGroupAverageMomentum( a_smr, a_settings, a_temperatureInfo, iter->first, a_reactionsToExclude );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, gain for the requested particle and label. This is a cross section weighted gain summed over all reactions.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID           [in]    The particle PoPs' id for the whose gain is to be calculated.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the gain.
 *
 * @return                          The requested multi-group gain as a **GIDI::Vector**.
 ***********************************************************************************************************/

Vector ProtareSingle::multiGroupGain( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector( 0 );
    std::string const projectile_ID = projectile( ).ID( );

    for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
        Reaction const *reaction1 = reactionToMultiGroup( a_settings, i1, a_reactionsToExclude );

        if( reaction1 != nullptr ) vector += reaction1->multiGroupGain( a_smr, a_settings, a_temperatureInfo, a_productID, projectile_ID );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 *
 *
 * @return      A list of label, mu cutoff pairs.
 ***********************************************************************************************************/

stringAndDoublePairs ProtareSingle::muCutoffForCoulombPlusNuclearElastic( ) const {

    stringAndDoublePairs muCutoffs;

    for( std::size_t i1 = 0; i1 < m_styles.size( ); ++i1 ) {
        Styles::Base const *style1 = m_styles.get<Styles::Base>( i1 );

        if( style1->moniker( ) == GIDI_CoulombPlusNuclearElasticMuCutoffStyleChars ) {
            Styles::CoulombPlusNuclearElasticMuCutoff const *style2 = static_cast<Styles::CoulombPlusNuclearElasticMuCutoff const *>( style1 );
            
            stringAndDoublePair labelMu( style2->label( ), style2->muCutoff( ) );

            muCutoffs.push_back( std::move( labelMu ) );
        }
    }

    return( muCutoffs );
}

/* *********************************************************************************************************//**
 * Returns the list of DelayedNeutronProduct instances.
 *
 * @return      a_delayedNeutronProducts        The list of delayed neutrons.
 ***********************************************************************************************************/

DelayedNeutronProducts ProtareSingle::delayedNeutronProducts( ) const {

    DelayedNeutronProducts delayedNeutronProducts1;

    if( hasFission( ) ) {
        for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
            Reaction const *reaction1 = m_reactions.get<Reaction>( i1 );

            if( !reaction1->active( ) ) continue;
            if( reaction1->hasFission( ) ) reaction1->delayedNeutronProducts( delayedNeutronProducts1 );
        }
    }

    return( delayedNeutronProducts1 );
}

/* *********************************************************************************************************//**
 * Calls the **incompleteParticles** method for each active reaction in the *reactions* and *orphanProducts* nodes.
 *
 * @param a_settings                [in]    Specifies the requested label.
 * @param a_incompleteParticles     [out]   The list of particles whose **completeParticle** method returns *false*.
 ***********************************************************************************************************/

void ProtareSingle::incompleteParticles( Transporting::Settings const &a_settings, std::set<std::string> &a_incompleteParticles ) const {

    for( std::size_t i1 = 0; i1 < m_reactions.size( ); ++i1 ) {
        Reaction const *reaction1 = m_reactions.get<Reaction>( i1 );

        if( !reaction1->active( ) ) continue;
        reaction1->incompleteParticles( a_settings, a_incompleteParticles );
    }

    for( std::size_t i1 = 0; i1 < m_orphanProducts.size( ); ++i1 ) {
        Reaction const *reaction1 = m_orphanProducts.get<Reaction>( i1 );

        if( !reaction1->active( ) ) continue;
        reaction1->incompleteParticles( a_settings, a_incompleteParticles );
    }
}

/* *********************************************************************************************************//**
 * Write *this* to a file in GNDS/XML format.
 *
 * @param       a_fileName          [in]        Name of file to save XML lines to.
 ***********************************************************************************************************/

void ProtareSingle::saveAs( std::string const &a_fileName ) const {

    GUPI::WriteInfo writeInfo;

    toXMLList( writeInfo, "" );

    std::ofstream fileio;
    fileio.open( a_fileName.c_str( ) );
    for( std::list<std::string>::iterator iter = writeInfo.m_lines.begin( ); iter != writeInfo.m_lines.end( ); ++iter ) {
        fileio << *iter << std::endl;
    }
    fileio.close( );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void ProtareSingle::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string header = LUPI_XML_verionEncoding;
    std::string attributes;

    a_writeInfo.push_back( header );

    attributes  = a_writeInfo.addAttribute( GIDI_projectileChars, projectile( ).ID( ) );
    attributes += a_writeInfo.addAttribute( GIDI_targetChars, GNDS_target( ).ID( ) );
    attributes += a_writeInfo.addAttribute( GIDI_evaluationChars, evaluation( ) );
    attributes += a_writeInfo.addAttribute( GIDI_formatChars, m_formatVersion.format( ) );
    attributes += a_writeInfo.addAttribute( GIDI_projectileFrameChars, frameToString( projectileFrame( ) ) );

    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );

    m_externalFiles.toXMLList( a_writeInfo, indent2 );
    m_styles.toXMLList( a_writeInfo, indent2 );

    std::vector<std::string> pops_XMLList;
    m_internalPoPs.toXMLList( pops_XMLList, indent2 );
    for( std::vector<std::string>::iterator iter = pops_XMLList.begin( ); iter != pops_XMLList.end( ); ++iter ) a_writeInfo.push_back( *iter );

    m_reactions.toXMLList( a_writeInfo, indent2 );
    m_orphanProducts.toXMLList( a_writeInfo, indent2 );
    m_sums.toXMLList( a_writeInfo, indent2 );
    m_fissionComponents.toXMLList( a_writeInfo, indent2 );

    a_writeInfo.addNodeEnder( moniker( ) );
}

/* *********************************************************************************************************//**
 * This method parses the **targetInfo** node in the **evaluated** style node into the *m_targetInfo* member of *this*.
 *
 * @param a_node                [in]    The protare (i.e., reactionSuite) node to be parsed and used to construct a Protare.
 ***********************************************************************************************************/

void ProtareSingle::parseEvaluatedTargetInfo( HAPI::Node const &a_node ) {

    m_targetInfo.parseEvaluatedTargetInfo( a_node );
}

}
