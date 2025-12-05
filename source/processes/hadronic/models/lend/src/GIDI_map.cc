/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "GIDI.hpp"
#include <HAPI.hpp>

static std::string GIDI_basePath( char const *a_path );
static std::string GIDI_basePath( std::string const a_path );
static std::string GIDI_addPaths( std::string const &a_base, std::string const &a_path );

namespace GIDI {

namespace Map {

/* *********************************************************************************************************//**
 * Returns the type of file for *a_filename*.
 *
 * @param a_filename                [in]    Path of the file.
 *
 * @return                                  A **FileType** instance.
 ***********************************************************************************************************/

static FileType fileType( std::string const &a_path ) {

    if( a_path.compare( a_path.size( ) - 2, 2, "h5" ) == 0 ) return( GIDI::FileType::HDF );
    return( GIDI::FileType::XML );
}

/* *********************************************************************************************************//**
 * User data passed to the Map::directory method. It stores the desired projectile, target, library and evalaute infomation
 * as a list of found matches. An empty string for the projectile's id matches all projectiles. 
 * A empty string for the target's id matches all targets. An empty evaluation string matches all evaluations.
 ***********************************************************************************************************/

class MapWalkDirectoryCallbackData {

    public:
        std::string const &m_projectileID;                      /**< The desired projectile's id. */
        std::string const &m_targetID;                          /**< The desired target's id. */
        std::string const &m_library;                           /**< The desired library name. */
        std::string const &m_evaluation;                        /**< The desired evaluation id. */

        std::vector<ProtareBase const *> m_protareEntries;     /**< list of matched protare entries. */

        /* *********************************************************************************************************//**
         *
         * @param a_projectileID        [in]    The projectile's id to match.
         * @param a_targetID            [in]    The target's id to match.
         * @param a_library             [in]    The library to match.
         * @param a_evaluation          [in]    The evaluation to match.
         ***********************************************************************************************************/

        MapWalkDirectoryCallbackData( std::string const &a_projectileID, std::string const &a_targetID, std::string const &a_library, std::string const &a_evaluation ) :
                m_projectileID( a_projectileID ),
                m_targetID( a_targetID ),
                m_library( a_library ),
                m_evaluation( a_evaluation ) {

        }
};

/* *********************************************************************************************************//**
 *
 * @param a_protareEntry        [in]    The protare entry to compare to the data protare parameters in *a_data*.
 * @param a_library             [in]    The library name of the current map file.
 * @param a_data                [in]    A MapWalkDirectoryCallbackData instance.
 * @param a_level               [in]    Nested level of *this* map file. For internal use.
 *
 * @return                              Always returns true.
 ***********************************************************************************************************/

bool MapWalkDirectoryCallback( ProtareBase const *a_protareEntry, std::string const &a_library, void *a_data, LUPI_maybeUnused int a_level ) {

    MapWalkDirectoryCallbackData *mapWalkDirectoryCallbackData = static_cast<MapWalkDirectoryCallbackData *>( a_data );

    if( ( mapWalkDirectoryCallbackData->m_projectileID == "" ) || 
            PoPI::compareSpecialParticleIDs( mapWalkDirectoryCallbackData->m_projectileID, a_protareEntry->projectileID( ) ) ) {
        if( ( mapWalkDirectoryCallbackData->m_targetID == "" ) || ( mapWalkDirectoryCallbackData->m_targetID == a_protareEntry->targetID( ) ) ) {
            if( ( mapWalkDirectoryCallbackData->m_library == "" ) || ( mapWalkDirectoryCallbackData->m_library == a_library ) ) {
                if( ( mapWalkDirectoryCallbackData->m_evaluation == "" ) || ( mapWalkDirectoryCallbackData->m_evaluation == a_protareEntry->evaluation( ) ) ) {
                    mapWalkDirectoryCallbackData->m_protareEntries.push_back( a_protareEntry );
                }
            }
        }
    }
    return( true );
}

/*! \class BaseEntry
 * This is the virtual base class inherited by all map entry classes.
 */

/* *********************************************************************************************************//**
 *
 * @param a_node        [in]    The **HAPI::Node** to be parsed.
 * @param a_basePath    [in]    A path prepended to this entry's path.
 * @param a_parent      [in]    Pointer to the *Map* containing *this*.
 ***********************************************************************************************************/

BaseEntry::BaseEntry( HAPI::Node const &a_node, std::string const &a_basePath, Map const *a_parent ) :
        GUPI::Ancestry( a_node.name( ) ),
        m_name( a_node.name( ) ),
        m_parent( a_parent ),
        m_path( a_node.attribute_as_string( GIDI_pathChars ) ),
        m_cumulativePath( GIDI_addPaths( a_basePath, m_path ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

BaseEntry::~BaseEntry( ) {

}

/* *********************************************************************************************************//**
 * Returns either the *entered*, *cumulative* or *realPath* path for the entry.
 *
 * @param a_form            [in]    The type of path to return.
 * @return                          The requested path.
 ***********************************************************************************************************/

std::string BaseEntry::path( PathForm a_form ) const {

    if( a_form == PathForm::entered ) return( m_path );
    if( a_form == PathForm::cumulative ) return( m_cumulativePath );
    return( LUPI::FileInfo::realPath( m_cumulativePath ) );
}

/* *********************************************************************************************************//**
 * Fills *a_libraries* with the name of all the libraries *this* is contained in. The first library in the list is the 
 * library *this* is defined in and the last is the starting library.
 *
 * @param a_libraries           [out]   The instances that is filled with the library names.
 ***********************************************************************************************************/

void BaseEntry::libraries( std::vector<std::string> &a_libraries ) const {

    parent( )->libraries( a_libraries );
}

/* *********************************************************************************************************//**
 *
 * @param a_node        [in]    The **HAPI::Node** to be parsed to contruct a Import instance.
 * @param a_pops        [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_basePath    [in]    A path prepended to this entry's path.
 * @param a_parent      [in]    Pointer to the *Map* containing *this*.
 ***********************************************************************************************************/

Import::Import( HAPI::Node const &a_node, PoPI::Database const &a_pops, std::string const &a_basePath, Map const *a_parent ) :
        BaseEntry( a_node, a_basePath, a_parent ),
        m_map( nullptr ) {

    m_map = new Map( path( BaseEntry::PathForm::cumulative ), a_pops, a_parent );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Import::~Import( ) {

    delete m_map;
}

/* *********************************************************************************************************//**
 * Returns the Protare entry to the first protare to match *a_projectileID*, *a_targetID*, *a_library* and *a_evaluation*. If
 * *a_evaluation* is an empty string, only *a_projectileID* and *a_targetID* are matched.
 *
 * @param a_projectileID        [in]    The projectile's id to match.
 * @param a_targetID            [in]    The target's id to match.
 * @param a_library             [in]    The library to match.
 * @param a_evaluation          [in]    The evaluation to match.
 *
 * @return                              The const pointer **ProtareBase** for the matched protare.
 ***********************************************************************************************************/

ProtareBase const *Import::findProtareEntry( std::string const &a_projectileID, std::string const &a_targetID,
                std::string const &a_library, std::string const &a_evaluation ) const {

    return( m_map->findProtareEntry( a_projectileID, a_targetID, a_library, a_evaluation ) );
}

/* *********************************************************************************************************//**
 * Returns the list of all Protare entries that match *a_projectileID*, *a_targetID*, *a_library* and *a_evaluation*.
 * The arguments *a_projectileID*, *a_targetID*, *a_library* and *a_evaluation* can be an C++ regex string. An empty
 * string for any of the arguments will match all.
 *
 * @param a_protareEntries      [out]   The list of **ProtareBase** found.
 * @param a_projectileID        [in]    The projectile's id to match.
 * @param a_targetID            [in]    The target's id to match.
 * @param a_library             [in]    The library to match.
 * @param a_evaluation          [in]    The evaluation to match.
 ***********************************************************************************************************/

void Import::findProtareEntries( std::vector<ProtareBase const *> &a_protareEntries, std::regex const &a_projectileID, 
                std::regex const &a_targetID, std::regex const &a_library, std::regex const &a_evaluation ) const {

    return( m_map->findProtareEntries( a_protareEntries, a_projectileID, a_targetID, a_library, a_evaluation ) );
}

/* *********************************************************************************************************//**
 * Returns the path to the first protare to match *a_projectileID*, *a_targetID*, *a_library* and *a_evaluation*. If
 * *a_evaluation* is an empty string, only *a_projectileID* and *a_targetID* are matched.
 *
 * @param a_projectileID        [in]    The projectile's id to match.
 * @param a_targetID            [in]    The target's id to match.
 * @param a_library             [in]    The library to match.
 * @param a_evaluation          [in]    The evaluation to match.
 * @param a_form                [in]    Determines the form of the path returned.
 *
 * @return                              The path to the matched protare.
 ***********************************************************************************************************/

std::string Import::protareFilename( std::string const &a_projectileID, std::string const &a_targetID, std::string const &a_library,
                std::string const &a_evaluation, PathForm a_form ) const {

    return( m_map->protareFilename( a_projectileID, a_targetID, a_library, a_evaluation, a_form ) );
}

/* *********************************************************************************************************//**
 * Returns a list of all evaluations with a match to *a_projectileID* and *a_targetID*.
 *
 * @param a_projectileID        [in]    The projectile's id to match.
 * @param a_targetID            [in]    The target's id to match.
 * @return                              List of evaluations.
 ***********************************************************************************************************/

std::vector<std::string> Import::availableEvaluations( std::string const &a_projectileID, std::string const &a_targetID ) const {

    return( m_map->availableEvaluations( a_projectileID, a_targetID ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Import::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string attributes;

    attributes  = a_writeInfo.addAttribute( GIDI_pathChars, path( ) );
    a_writeInfo.addNodeStarterEnder( a_indent, moniker( ), attributes );
}

/* *********************************************************************************************************//**
 * @param a_node        [in]    The **HAPI::Node** to be parsed.
 * @param a_basePath    [in]    A path prepended to this entry's path.
 * @param a_parent      [in]    Pointer to the *Map* containing *this*.
 ***********************************************************************************************************/

ProtareBase::ProtareBase( HAPI::Node const &a_node, std::string const &a_basePath, Map const *const a_parent ) :
        BaseEntry( a_node, a_basePath, a_parent ),
        m_projectileID( a_node.attribute_as_string( GIDI_projectileChars ) ),
        m_targetID( a_node.attribute_as_string( GIDI_targetChars ) ),
        m_evaluation( a_node.attribute_as_string( GIDI_evaluationChars ) ),
        m_interaction( a_node.attribute_as_string( GIDI_interactionChars ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

ProtareBase::~ProtareBase( ) {

}

/* *********************************************************************************************************//**
 * Returns the library *this* is contained in.
 ***********************************************************************************************************/

std::string const &ProtareBase::library( ) const {

    return( parent( )->library( ) );
}

/* *********************************************************************************************************//**
 * Returns the resolved library *this* is contained in.
 ***********************************************************************************************************/

std::string const &ProtareBase::resolvedLibrary( ) const {

    return( parent( )->resolvedLibrary( ) );
}

/* *********************************************************************************************************//**
 * Compares *a_projectileID*, *a_targetID* and *a_evaluation* to *this* data and returns true if they match
 * and false otherwise. If *a_evaluation* is an empty string, only *a_projectileID* and *a_targetID* are compared.
 *
 * @param a_projectileID        [in]    The projectile's id to match.
 * @param a_targetID            [in]    The target's id to match.
 * @param a_library             [in]    The library to match.
 * @param a_evaluation          [in]    The evaluation to match.
 *
 * @return                              The *this* pointer if *this* matches otherwise a **nullptr**.
 ***********************************************************************************************************/

ProtareBase const *ProtareBase::findProtareEntry( std::string const &a_projectileID, std::string const &a_targetID,
                std::string const &a_library, std::string const &a_evaluation ) const {

    if( !PoPI::compareSpecialParticleIDs( a_projectileID, projectileID( ) ) ) return( nullptr );
    if( a_targetID != targetID( ) ) return( nullptr );
    if( ( a_library != "" ) && ( parent( )->library( ) != a_library ) ) return( nullptr );
    if( ( a_evaluation == "" ) || ( a_evaluation == evaluation( ) ) ) return( this );
    return( nullptr );
}

/* *********************************************************************************************************//**
 * Returns the list of all Protare entries that match *a_projectileID*, *a_targetID*, *a_library* and *a_evaluation*. 
 * The arguments *a_projectileID*, *a_targetID*, *a_library* and *a_evaluation* can be an C++ regex string. An empty
 * string for any of the arguments will match all.
 *
 * @param a_protareEntries      [out]   The list of **ProtareBase** found.
 * @param a_projectileID        [in]    The projectile's id to match.
 * @param a_targetID            [in]    The target's id to match.
 * @param a_library             [in]    The library to match.
 * @param a_evaluation          [in]    The evaluation to match.
 ***********************************************************************************************************/

void ProtareBase::findProtareEntries( std::vector<ProtareBase const *> &a_protareEntries, std::regex const &a_projectileID,
                std::regex const &a_targetID, std::regex const &a_library, std::regex const &a_evaluation ) const {

    if( regex_match( m_projectileID, a_projectileID ) ) {
        if( regex_match( m_targetID, a_targetID ) ) {
            if( regex_match( parent( )->library( ), a_library ) ) {
                if( regex_match( m_evaluation, a_evaluation ) ) a_protareEntries.push_back( this );
            }
        }
    }
}

/* *********************************************************************************************************//**
 * Compares *a_projectileID*, *a_targetID* and *a_evaluation* to *this* data and returns true if they match
 * and false otherwise. If *a_evaluation* is an empty string, only *a_projectileID* and *a_targetID* are compared.
 *
 * @param a_projectileID        [in]    The projectile's id to match.
 * @param a_targetID            [in]    The target's id to match.
 * @param a_evaluation          [in]    The evaluation to match.
 * @return                              true if match and false otherwise.
 ***********************************************************************************************************/

bool ProtareBase::isMatch( std::string const &a_projectileID, std::string const &a_targetID, std::string const &a_evaluation ) const {

    if( !PoPI::compareSpecialParticleIDs( a_projectileID, m_projectileID ) ) return( false );
    if( a_targetID != m_targetID ) return( false );
    if( a_evaluation == "" ) return( true );
    return( a_evaluation == m_evaluation );
}

/* *********************************************************************************************************//**
 *
 * @param a_node        [in]    The **HAPI::Node** to be parsed to contruct a Protare entry instance.
 * @param a_pops        [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_basePath    [in]    A path prepended to this entry's path.
 * @param a_parent      [in]    Pointer to the *Map* containing *this*.
 ***********************************************************************************************************/

Protare::Protare( HAPI::Node const &a_node, PoPI::Database const &a_pops, std::string const &a_basePath, Map const *const a_parent ) :
        ProtareBase( a_node, a_basePath, a_parent ),
        m_isPhotoAtomic( false ) {

    if( interaction( ) == "" ) {                // Some old GNDS 1.10 files do not have an "interaction" attribute.
        setInteraction( GIDI_MapInteractionNuclearChars );
        if( PoPI::IDs::photon == projectileID( ) ) {
            try {
                PoPI::Base const &target( a_pops.get<PoPI::Base>( targetID( ) ) );
                m_isPhotoAtomic = target.isChemicalElement( ); }
            catch (...) {                       // Let's ignore, but user must understand that m_interaction and m_isPhotoAtomic may be wrong.
            }
            if( m_isPhotoAtomic ) setInteraction( GIDI_MapInteractionAtomicChars );
        }
    }
    m_isPhotoAtomic = interaction( ) == GIDI_MapInteractionAtomicChars;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Protare::~Protare( ) {

}

/* *********************************************************************************************************//**
 * Returns a GIDI::Protare instance of the protare reference by the *m_path* member.
 *
 * @param a_construction            [in]    Used to pass user options to the constructor.
 * @param a_pops                    [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_particleSubstitution    [in]    Map of particles to substitute with another particles.
 *
 * @return                          Returns the Protare matching the TNSL protare.
 ***********************************************************************************************************/

GIDI::Protare *Protare::protare( Construction::Settings const &a_construction, PoPI::Database const &a_pops, ParticleSubstitution const &a_particleSubstitution ) const {

    return( protareSingle( a_construction, a_pops, a_particleSubstitution ) );
}

/* *********************************************************************************************************//**
 * Returns a GIDI::ProtareSingle instance of the protare reference by the *m_path* member.
 *
 * @param a_construction            [in]    Used to pass user options to the constructor.
 * @param a_pops                    [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_particleSubstitution    [in]    Map of particles to substitute with another particles.
 *
 * @return                          Returns the Protare matching the TNSL protare.
 ***********************************************************************************************************/

GIDI::ProtareSingle *Protare::protareSingle( Construction::Settings const &a_construction, PoPI::Database const &a_pops, ParticleSubstitution const &a_particleSubstitution ) const {

    std::vector<std::string> libraries1;
    libraries( libraries1 );

    return( new ProtareSingle( a_construction, path( ), fileType( path( ) ), a_pops, a_particleSubstitution, libraries1, interaction( ), true ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Protare::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string attributes;

    attributes  = a_writeInfo.addAttribute( GIDI_projectileChars, projectileID( ) );
    attributes += a_writeInfo.addAttribute( GIDI_targetChars, targetID( ) );
    attributes += a_writeInfo.addAttribute( GIDI_evaluationChars, evaluation( ) );
    attributes += a_writeInfo.addAttribute( GIDI_pathChars, path( PathForm::entered ) );
    attributes += a_writeInfo.addAttribute( GIDI_interactionChars, interaction( ) );
    a_writeInfo.addNodeStarterEnder( a_indent, moniker( ), attributes );
}

/* *********************************************************************************************************//**
 *
 * @param a_node        [in]    The **HAPI::Node** to be parsed to contruct a TNSL entry instance.
 * @param a_pops        [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_basePath    [in]    A path prepended to this entry's path.
 * @param a_parent      [in]    Pointer to the *Map* containing *this*.
 ***********************************************************************************************************/

TNSL::TNSL( HAPI::Node const &a_node, LUPI_maybeUnused PoPI::Database const &a_pops, std::string const &a_basePath, Map const *const a_parent ) :
        ProtareBase( a_node, a_basePath, a_parent ),
        m_standardTarget( a_node.attribute_as_string( GIDI_standardTargetChars ) ),
        m_standardEvaluation( a_node.attribute_as_string( GIDI_standardEvaluationChars ) ) {

    setInteraction( "" );

    HAPI::Node const &protare = a_node.child( GIDI_protareChars );              // Format 0.1 support. This format is deprecated.
    if( !protare.empty() ) {
        m_standardTarget = protare.attribute_as_string( GIDI_targetChars );
        m_standardEvaluation = protare.attribute_as_string( GIDI_evaluationChars );
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

TNSL::~TNSL( ) {

}

/* *********************************************************************************************************//**
 * Returns a GIDI::ProtareTNSL instance of the protare reference by the *m_path* member.
 *
 * @param a_construction            [in]    Used to pass user options to the constructor.
 * @param a_pops                    [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_particleSubstitution    [in]    Map of particles to substitute with another particles.
 *
 * @return                          Returns the Protare matching the TNSL protare.
 ***********************************************************************************************************/

GIDI::Protare *TNSL::protare( Construction::Settings const &a_construction, PoPI::Database const &a_pops, LUPI_maybeUnused ParticleSubstitution const &a_particleSubstitution ) const {

    Map const *map = parent( );

    ParticleSubstitution particleSubstitution;
    std::vector<std::string> libraries1;
    libraries( libraries1 );

    ProtareSingle *protare1 = new ProtareSingle( a_construction, path( ), fileType( path( ) ), a_pops, particleSubstitution, libraries1, GIDI_MapInteractionTNSLChars, false );

    while( true ) {
        Map const *parent = map->parent( );

        if( parent == nullptr ) break;
        map = parent;
    }
    ProtareSingle *protare2 = static_cast<ProtareSingle *>( map->protare( a_construction, a_pops, PoPI::IDs::neutron, m_standardTarget, "", m_standardEvaluation ) );
    if( protare2 == nullptr ) protare2 = static_cast<ProtareSingle *>( map->protare( a_construction, a_pops, PoPI::IDs::neutron, m_standardTarget ) );

    ProtareTNSL *protareTNSL = new ProtareTNSL( a_construction, protare2, protare1 );

    return( protareTNSL );
}

/* *********************************************************************************************************//**
 * Returns a GIDI::ProtareSingle instance of the protare reference by the *m_path* member. Note, this is different from the *protare* method
 * which returns a **GIDI::ProtareTNSL** instance.
 *
 * @param a_construction            [in]    Used to pass user options to the constructor.
 * @param a_pops                    [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_particleSubstitution    [in]    Map of particles to substitute with another particles.
 *
 * @return                          Returns the Protare matching the TNSL protare.
 ***********************************************************************************************************/

GIDI::ProtareSingle *TNSL::protareSingle( Construction::Settings const &a_construction, PoPI::Database const &a_pops, ParticleSubstitution const &a_particleSubstitution ) const {

    std::vector<std::string> libraries1;
    libraries( libraries1 );

    return( new ProtareSingle( a_construction, path( ), fileType( path( ) ), a_pops, a_particleSubstitution, libraries1, GIDI_MapInteractionTNSLChars, false ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void TNSL::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string attributes;

    attributes  = a_writeInfo.addAttribute( GIDI_projectileChars, projectileID( ) );
    attributes += a_writeInfo.addAttribute( GIDI_targetChars, targetID( ) );
    attributes += a_writeInfo.addAttribute( GIDI_evaluationChars, evaluation( ) );
    attributes += a_writeInfo.addAttribute( GIDI_pathChars, path( PathForm::entered ) );
    attributes += a_writeInfo.addAttribute( GIDI_standardTargetChars, standardTarget( ) );
    attributes += a_writeInfo.addAttribute( GIDI_standardEvaluationChars, standardEvaluation( ) );
    a_writeInfo.addNodeStarterEnder( a_indent, moniker( ), attributes );
}
/* *********************************************************************************************************//**
 *
 * @param a_fileName    [in]    The path to the map file to parse to construct a Map instance.
 * @param a_pops        [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_parent      [in]    Pointer to the *Map* containing *this*.
 ***********************************************************************************************************/

Map::Map( std::string const &a_fileName, PoPI::Database const &a_pops, Map const *a_parent ) :
        GUPI::Ancestry( GIDI_mapChars ) {

    initialize( a_fileName, a_pops, a_parent );
}

/* *********************************************************************************************************//**
 *
 * @param a_node        [in]    HAPI::Node corresponding to the map node
 * @param a_fileName    [in]    std::string, the name of the file containing this map.
 * @param a_pops        [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_parent      [in]    Pointer to the *Map* containing *this*.
 ***********************************************************************************************************/

Map::Map( HAPI::Node const &a_node, std::string const &a_fileName, PoPI::Database const &a_pops, Map const *a_parent ) :
        GUPI::Ancestry( a_node.name( ) ) {

    initialize( a_node, a_fileName, a_pops, a_parent );
}

/* *********************************************************************************************************//**
 * This method is called by the fileName constructors, opens the document and calls the other initialize method
 *
 * @param a_fileName    [in]    The path to the map file to parse to construct a Map instance.
 * @param a_pops        [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_parent      [in]    Pointer to the *Map* containing *this*.
 ***********************************************************************************************************/

void Map::initialize( std::string const &a_fileName, PoPI::Database const &a_pops, Map const *a_parent ) {

    HAPI::File *doc = new HAPI::PugiXMLFile( a_fileName.c_str( ), "Map::initialize" );

    HAPI::Node map = doc->first_child( );

    if( strcmp( map.name( ).c_str( ), GIDI_mapChars ) != 0 ) throw Exception( "Invalid map file " + a_fileName );

    initialize( map, a_fileName, a_pops, a_parent );
    delete doc;
}

/* *********************************************************************************************************//**
 * This method is called either by the constructor or by the other initialize method. Does most of the work of parsing
 *
 * @param a_node        [in]    HAPI::Node corresponding to the map node.
 * @param a_fileName    [in]    std::string, the name of the file containing this map
 * @param a_pops        [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_parent      [in]    Pointer to the *Map* containing *this*.
 ***********************************************************************************************************/

void Map::initialize( HAPI::Node const &a_node, std::string const &a_fileName, PoPI::Database const &a_pops, Map const *a_parent ) {

    m_parent = a_parent;
    m_fileName = a_fileName;
    m_realFileName = LUPI::FileInfo::realPath( a_fileName );
    m_projectilesLoaded = false;

    std::string basePath = GIDI_basePath( m_realFileName );

    std::string format = a_node.attribute_as_string( GIDI_formatChars );
    if( ( format == GNDS_formatVersion_2_0_LLNL_4Chars ) || ( format == GIDI_mapFormatVersion_0_2Chars ) ) format = GNDS_formatVersion_2_0Chars;
    if( format != GNDS_formatVersion_2_0Chars ) {
        if( format != GIDI_mapFormatVersion_0_1Chars ) throw Exception( "Unsupported map format" );
    }

    m_library = a_node.attribute_as_string( GIDI_libraryChars );
    for( HAPI::Node child = a_node.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
        if( strcmp( child.name( ).c_str( ), GIDI_importChars ) == 0 ) {
            m_entries.push_back( new Import( child, a_pops, basePath, this ) ); }
        else if( strcmp( child.name( ).c_str( ), GIDI_protareChars ) == 0 ) {
            m_entries.push_back( new Protare( child, a_pops, basePath, this ) ); }
        else if( strcmp( child.name( ).c_str( ), GIDI_TNSLChars ) == 0 ) {
            m_entries.push_back( new TNSL( child, a_pops, basePath, this ) ); }
        else {
            throw Exception( std::string( "Invalid entry '" ) + child.name( ) + std::string( "' in map file " ) + a_fileName );
        }
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Map::~Map( ) {

    for( std::vector<BaseEntry *>::const_iterator iter = m_entries.begin( ); iter < m_entries.end( ); ++iter ) delete *iter;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

std::string const &Map::resolvedLibrary( ) const {

    if( ( m_library == "" ) && ( m_parent != nullptr ) ) return( m_parent->resolvedLibrary( ) );
    return( m_library );
}

/* *********************************************************************************************************//**
 * Fills *a_libraries* with the name of all the libraries *this* is contained in. The first library in the list is the 
 * library *this* is defined in and the last is the starting library.
 *
 * @param a_libraries           [out]   The instances that is filled with the library names.
 ***********************************************************************************************************/

void Map::libraries( std::vector<std::string> &a_libraries ) const {

    a_libraries.push_back( m_library );
    if( m_parent != nullptr ) m_parent->libraries( a_libraries );
}

/* *********************************************************************************************************//**
 * Returns the Protare to the first protare to match *a_projectileID*, *a_targetID* and *a_evaluation*. If
 * *a_evaluation* is an empty string, only *a_projectileID* and *a_targetID* are matched.
 *
 * @param a_projectileID        [in]    The projectile's id to match.
 * @param a_targetID            [in]    The target's id to match.
 * @param a_library             [in]    The library to match.
 * @param a_evaluation          [in]    The evaluation to match.
 *
 * @return                              The const pointer **ProtareBase** for the matched protare.
 ***********************************************************************************************************/

ProtareBase const *Map::findProtareEntry( std::string const &a_projectileID, std::string const &a_targetID, std::string const &a_library,
                std::string const &a_evaluation ) const {

    ProtareBase const *protareEntry = nullptr;

    for( std::vector<BaseEntry *>::const_iterator iter = m_entries.begin( ); iter != m_entries.end( ); ++iter ) {
        protareEntry = (*iter)->findProtareEntry( a_projectileID, a_targetID, a_library, a_evaluation );
        if( protareEntry != nullptr ) break;
    }
    return( protareEntry );
}

/* *********************************************************************************************************//**
 * Returns the list of all Protare entries that match *a_projectileID*, *a_targetID*, *a_library* and *a_evaluation*.
 * The arguments *a_projectileID*, *a_targetID*, *a_library* and *a_evaluation* can be an C++ regex string. An empty
 * string for any of the arguments will match all.
 * 
 * @param a_protareEntries      [out]   The list of **ProtareBase** found.
 * @param a_projectileID        [in]    The projectile's id to match.
 * @param a_targetID            [in]    The target's id to match.
 * @param a_library             [in]    The library to match.
 * @param a_evaluation          [in]    The evaluation to match.
 ***********************************************************************************************************/

void Map::findProtareEntries( std::vector<ProtareBase const *> &a_protareEntries, std::regex const &a_projectileID,
                std::regex const &a_targetID, std::regex const &a_library, std::regex const &a_evaluation ) const {

    for( std::vector<BaseEntry *>::const_iterator iter = m_entries.begin( ); iter != m_entries.end( ); ++iter ) {
        (*iter)->findProtareEntries( a_protareEntries, a_projectileID, a_targetID, a_library, a_evaluation );
    }
}

/* *********************************************************************************************************//**
 * Returns the path to the first protare to match *a_projectileID*, *a_targetID* and *a_evaluation*. If
 * *a_evaluation* is an empty string, only *a_projectileID* and *a_targetID* are matched.
 *
 * @param a_projectileID        [in]    The projectile's id to match.
 * @param a_targetID            [in]    The target's id to match.
 * @param a_library             [in]    The library to match.
 * @param a_evaluation          [in]    The evaluation to match.
 * @param a_form                [in]    Determines the form of the path returned.
 * @return                              The path to the matched protare.
 ***********************************************************************************************************/

std::string Map::protareFilename( std::string const &a_projectileID, std::string const &a_targetID, std::string const &a_library,
                std::string const &a_evaluation, BaseEntry::PathForm a_form ) const {

    ProtareBase const *protareEntry = findProtareEntry( a_projectileID, a_targetID, a_library, a_evaluation );

    if( protareEntry != nullptr ) return( protareEntry->path( a_form ) );
    return( GIDI_emptyFileNameChars );
}

/* *********************************************************************************************************//**
 * Returns *true* if *a_targetID* is a thermal neutron scattering law (TNSL) target id contained this map, including recursion, and *false* otherwise.
 *
 * @param a_targetID            [in]    The target's id.
 *
 * @return                              *true* if target *a_targetID* is present and is a TNSL target, and *false* otherwise.
 ***********************************************************************************************************/

bool Map::isTNSL_target( std::string const &a_targetID ) const {

    ProtareBase const *protareEntry = findProtareEntry( PoPI::IDs::neutron, a_targetID );

    if( protareEntry == nullptr ) return( false );
    if( protareEntry->entryType( ) == EntryType::TNSL ) return( true );
    return( false );
}

/* *********************************************************************************************************//**
 * Returns a list of all evaluations with a match to *a_projectileID* and *a_targetID*.
 *
 * @param a_projectileID        [in]    The projectile's id to match.
 * @param a_targetID            [in]    The target's id to match.
 * @return                              List of evaluations.
 ***********************************************************************************************************/

std::vector<std::string> Map::availableEvaluations( std::string const &a_projectileID, std::string const &a_targetID ) const {

    std::vector<std::string> list;

    for( std::vector<BaseEntry *>::const_iterator iter1 = m_entries.begin( ); iter1 != m_entries.end( ); ++iter1 ) {
        if( (*iter1)->name( ) == GIDI_importChars ) {
            Import *_mapEntry = static_cast<Import *> (*iter1);

            std::vector<std::string> sub_list = _mapEntry->availableEvaluations( a_projectileID, a_targetID );
            for( std::vector<std::string>::const_iterator iter2 = sub_list.begin( ); iter2 != sub_list.end( ); ++iter2 )
                list.push_back( *iter2 ); }
        else {
            ProtareBase *protareEntry = static_cast<ProtareBase *> (*iter1);

            if( protareEntry->isMatch( a_projectileID, a_targetID ) ) list.push_back( protareEntry->evaluation( ) );
        }
    }
    return( list );
}

/* *********************************************************************************************************//**
 * If a protare matching *a_projectileID*, *a_targetID* and *a_evaluation* is found, the Protare constructor is called with
 * its fileName.
 *
 * @param a_construction                [in]    Pass to the Protare constructor.
 * @param a_pops                        [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_projectileID                [in]    The projectile's id to match.
 * @param a_targetID                    [in]    The target's id to match.
 * @param a_library                     [in]    The library to match.
 * @param a_evaluation                  [in]    The evaluation to match.
 * @param a_targetRequiredInGlobalPoPs  [in]    If *true*, the target is required to be in **a_pops**.
 * @param a_ignorePoPs                  [in]    If *true*, no particle is required to be in **a_pops**.
 ***********************************************************************************************************/

GIDI::Protare *Map::protare( Construction::Settings const &a_construction, PoPI::Database const &a_pops, std::string const &a_projectileID, 
                std::string const &a_targetID, std::string const &a_library, std::string const &a_evaluation, LUPI_maybeUnused bool a_targetRequiredInGlobalPoPs, LUPI_maybeUnused bool a_ignorePoPs ) const {

    std::string targetID( a_targetID );
    std::string atomicTargetID;

    ParticleSubstitution particleSubstitution;
    GIDI::Protare *nuclear = nullptr, *atomic = nullptr, *protare;

    if( a_projectileID != PoPI::IDs::neutron ) {                // Check if targetID is for TNSL target. If so, and projectile is not a neutron, need to use standardTarget name.
        ProtareBase const *protareEntry = findProtareEntry( PoPI::IDs::neutron, targetID );

        if( protareEntry != nullptr ) {
            if( protareEntry->entryType( ) == EntryType::TNSL ) targetID = static_cast<TNSL const *>( protareEntry )->standardTarget( );
        }
    }

    if( a_projectileID == PoPI::IDs::photon ) {
        PoPI::ParseIdInfo parseIdInfo( targetID );

        if( a_construction.photoMode( ) != Construction::PhotoMode::nuclearOnly ) {
            atomicTargetID = targetID;                                          // Kludge for 99120 and similar targets.

            if( parseIdInfo.isNuclear( ) ) atomicTargetID = parseIdInfo.symbol( );
            ProtareBase const *protareEntry = findProtareEntry( a_projectileID, atomicTargetID, a_library, a_evaluation );
            if( protareEntry != nullptr ) {
                particleSubstitution.insert( { atomicTargetID, ParticleInfo( targetID, a_pops, a_pops, true ) } );
                atomic = protareEntry->protare( a_construction, a_pops, particleSubstitution );
                particleSubstitution.clear( );
            }
        }
        if( ( a_construction.photoMode( ) != Construction::PhotoMode::atomicOnly ) && ( targetID != atomicTargetID ) ) {
            if( parseIdInfo.isSupported( ) ) {                // Kludge to ignore 99120 and similar targers.
                ProtareBase const *protareEntry = findProtareEntry( a_projectileID, targetID, a_library, a_evaluation );
                if( protareEntry != nullptr ) nuclear = protareEntry->protare( a_construction, a_pops, particleSubstitution );
            }
        } }
    else {
        ProtareBase const *protareEntry = findProtareEntry( a_projectileID, targetID, a_library, a_evaluation );
        if( protareEntry != nullptr ) nuclear = protareEntry->protare( a_construction, a_pops, particleSubstitution );
    }

    if( nuclear == nullptr ) {
        protare = atomic; }
    else if( atomic == nullptr ) {
        protare = nuclear; }
    else {
        ProtareComposite *protareComposite = new ProtareComposite( a_construction );

        protareComposite->setProjectile( nuclear->projectile( ) );
        protareComposite->setTarget( nuclear->target( ) );
        protareComposite->append( nuclear );
        protareComposite->append( atomic );
        protare = protareComposite;
    }

    return( protare );
}

/* *********************************************************************************************************//**
 * Returns a list of all Protare entry's matching the input data.
 *
 * @param a_projectileID        [in]    The projectile's id to match.
 * @param a_targetID            [in]    The target's id to match.
 * @param a_library             [in]    The library to match.
 * @param a_evaluation          [in]    The evaluation to match.
 * @return                              List of all Protare entry's matching input parameters.
 ***********************************************************************************************************/

std::vector<ProtareBase const *> Map::directory( std::string const &a_projectileID, std::string const &a_targetID, std::string const &a_library, 
                std::string const &a_evaluation ) const {

    MapWalkDirectoryCallbackData mapWalkDirectoryCallbackData( a_projectileID, a_targetID, a_library, a_evaluation );

    walk( MapWalkDirectoryCallback, &mapWalkDirectoryCallbackData, 0 );
    return( mapWalkDirectoryCallbackData.m_protareEntries );
}

/* *********************************************************************************************************//**
 * A method to walk a map file. For each Protare entry found, the **a_mapWalkCallBack** function is called with
 * a pointer to the Protare entry and **a_userData** has its arguments.
 *
 * @param           a_mapWalkCallBack       [in]    The callback function.
 * @param           a_userData              [in]    Pointer to user data.
 * @param           a_level                 [in]    Nested level of *this* map file. For internal use.
 *
 * @return          true if no issue is found and false if an issue is found.
 ***********************************************************************************************************/

bool Map::walk( MapWalkCallBack a_mapWalkCallBack, void *a_userData, int a_level ) const {

    for( std::size_t i1 = 0; i1 < size( ); ++i1 ) {
        BaseEntry const *entry = (*this)[i1];

        std::string path = entry->path( BaseEntry::PathForm::cumulative );

        if( entry->name( ) == GIDI_importChars ) {
            Import const *mapEntry = static_cast<Import const *>( entry );
            if( !mapEntry->map( )->walk( a_mapWalkCallBack, a_userData, a_level + 1 ) ) return( true ); }
        else if( ( entry->name( ) == GIDI_protareChars ) || ( entry->name( ) == GIDI_TNSLChars ) ) {
            if( !a_mapWalkCallBack( static_cast<ProtareBase const *>( entry ), m_library, a_userData, a_level ) ) return( true ); }
        else {
            std::cerr << "    ERROR: unknown map entry name: " << entry->name( ) << std::endl;
        }
    }

    return( true );
}

/* *********************************************************************************************************//**
 * Write *this* to a file in GNDS/XML format.
 *
 * @param       a_fileName          [in]        Name of file to save XML lines to.
 ***********************************************************************************************************/

void Map::saveAs( std::string const &a_fileName ) const {

    GUPI::WriteInfo writeInfo;

    toXMLList( writeInfo, "" );

    std::ofstream fileio;
    fileio.open( a_fileName.c_str( ) );
    for( std::list<std::string>::iterator iter = writeInfo.m_lines.begin( ); iter != writeInfo.m_lines.end( ); ++iter ) fileio << *iter << std::endl;
    fileio.close( );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Map::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string header = LUPI_XML_verionEncoding;
    std::string attributes;

    a_writeInfo.push_back( header );

    attributes  = a_writeInfo.addAttribute( GIDI_libraryChars, m_library );
    attributes += a_writeInfo.addAttribute( GIDI_formatChars, GNDS_formatVersion_2_0Chars );

    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );
    for( auto iter = m_entries.begin( ); iter != m_entries.end( ); ++iter ) (*iter)->toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/* *********************************************************************************************************//**
 * Returns the name of the RIS file for *this* map file. All RIS files must have a standard name which is the name of
 * the map file with its extension replaced with ".ris".
 *
 * @return                                      The **std::string** representing the RIS file.
 ***********************************************************************************************************/

std::string Map::RIS_fileName( ) {

    std::size_t found = m_fileName.rfind( '.' );
    std::string RIS_fileName( m_fileName.substr( 0, found ) );

    return( RIS_fileName + ".ris" );
}

/* *********************************************************************************************************//**
 * Returns **true** if the RIS file for *this* map file exists and **false** otherwise.
 *
 * @return                                      A boolean indicating if the RIS file exists or not.
 ***********************************************************************************************************/

bool Map::RIS_fileExist( ) {

    return( LUPI::FileInfo::exists( RIS_fileName( ) ) );
}

/* *********************************************************************************************************//**
 * Load the data from the RIS file if it exists. If it does not exists, no error is reported and the returned
 * **RISI::Projectiles** will be empty.
 *
 * @param       a_energyUnit        [in/out]    The unit desired for threshold energies.
 *
 * @return                                      A const reference to the **RISI::Projectiles** data.
 ***********************************************************************************************************/

RISI::Projectiles const &Map::RIS_load( std::string const &a_energyUnit ) {

    if( !m_projectilesLoaded ) {
        if( RIS_fileExist( ) ) GIDI::RISI::readRIS( RIS_fileName( ), a_energyUnit, m_projectiles );
    }
    m_projectilesLoaded = true;

    return( m_projectiles );
}

/* *********************************************************************************************************//**
 * If *this* has target *a_target* for projectile *a_projectile*, then *a_target* is returned. If it does not have target *a_target*, 
 * for projectile *a_projectile*, then a reasonable substitute is returned if one can be found in *this*. If no reasonable substitute 
 * is found, an empty string is returned.
 * 
 *
 * @param a_pops                [in]    A PoPI::Database instance used to get information about *a_target*..
 * @param a_projectile          [in]    The PoPs id of the projectile.
 * @param a_target              [in]    The PoPs id of the target.
 *
 * @return                              The PoPs id of the replacement target.
 ***********************************************************************************************************/

std::string Map::replacementTarget( PoPI::Database const &a_pops, std::string const &a_projectile, std::string const &a_target ) {

    if( findProtareEntry( a_projectile, a_target ) != nullptr ) return( a_target );

    if( a_pops.exists( a_target ) ) {
        std::string particleID = a_pops.final( a_target );
        if( a_pops.isParticle( particleID ) ) {
            PoPI::Particle const &particle = a_pops.particle( particleID );
            if( particle.hasNucleus( ) ) {
                PoPI::Nuclide const *nuclide = static_cast<PoPI::Nuclide const *>( &particle );
                if( particle.isNucleus( ) ) {
                    PoPI::Nucleus const *nucleus = static_cast<PoPI::Nucleus const *>( &particle );
                    nuclide = nucleus->nuclide( );
                }
                PoPI::Isotope const *isotope = nuclide->isotope( );
                std::string targetRegexString( isotope->chemicalElement( )->symbol( ) + "[0-9]+" );

                std::vector<ProtareBase const *> protareEntries;
                findProtareEntries( protareEntries, std::regex( a_projectile ), std::regex( targetRegexString ) );
                if( protareEntries.size( ) > 0 ) {
                    int offset = 0;
                    std::map<int, std::string> choices;
                    for( auto entryIter = protareEntries.begin( ); entryIter != protareEntries.end( ); ++entryIter ) {
                        PoPI::Nuclide const &nuclide2 = a_pops.get<PoPI::Nuclide const>( (*entryIter)->targetID( ) );
                        int diffA = nuclide->A( ) - nuclide2.A( );
                        offset += diffA;
                        choices[diffA] = nuclide2.ID( );
                    }
                    offset = offset < 0 ? -1 : 1;
                    int diff = 2 * offset;
                    int step = 2;
                    for( int doTwo = 0; doTwo < 2; ++doTwo ) {
                        for( std::size_t index = 0; index < choices.size( ); ++index ) {
                            auto iter = choices.find( diff );
                            if( iter != choices.end( ) ) return( (*iter).second );

                            diff *= -1;
                            iter = choices.find( diff );
                            if( iter != choices.end( ) ) return( (*iter).second );

                            diff = -diff + step * offset;
                        }
                        step = 1;
                    }
                }
            }
        }
    }

    return( "" );
}

}       // End of namespace Map.

}       // End of namespace GIDI.

/* *********************************************************************************************************//**
 * Splits the path at the last path separator (e.g., the '/' charactor on Unix systems) and returns the first (i.e.,
 * directory) part. Returns "." is no '/' is present.
 *
 * @param a_path            The path whose directory is to be returned.
 *
 * @return                  The directory of file **a_path**
 ***********************************************************************************************************/

static std::string GIDI_basePath( char const *a_path ) {

    char *p1, realPath[LUPI_PATH_MAX+1];

    strcpy( realPath, a_path );
    if( ( p1 = strrchr( realPath, '/' ) ) != nullptr ) {
        *p1 = 0; }
    else {
        strcpy( realPath, "." );
    }
    std::string basePath( realPath );
    return( basePath );
}

/* *********************************************************************************************************//**
 * Calls GIDI_basePath( char const *a_path ).
 *
 * @param a_path
 * @return
 ***********************************************************************************************************/

static std::string GIDI_basePath( std::string const a_path ) {

    return( GIDI_basePath( a_path.c_str( ) ) );
}

/* *********************************************************************************************************//**
 * If **a_path** is not an absolute path, prepends **a_path** to it.
 *
 * @param a_base            [in]    Base path to prepend to **a_path**.
 * @param a_path            [in]    Path
 * @return                          Prepend path.
 ***********************************************************************************************************/

static std::string GIDI_addPaths( std::string const &a_base, std::string const &a_path ) {

    std::string path( a_path );

    if( ( a_base.size( ) > 0 ) && ( path[0] != GIDI_FILE_SEPARATOR[0] ) ) path = a_base + GIDI_FILE_SEPARATOR + path;
    return( path );
}
