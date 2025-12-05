/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdexcept>

#include "PoPI.hpp"

#define PoPI_gaugeBosonsChars "gaugeBosons"
#define PoPI_leptonsChars "leptons"
#define PoPI_baryonsChars "baryons"
#define PoPI_unorthodoxesChars "unorthodoxes"

#define MsgSize (8 * 1024)
#ifdef _WIN32
#define __func__ __FUNCTION__
#endif

namespace PoPI {

static void parseAliases( HAPI::Node const &a_node, Database *a_DB );

/*! \class Database
 * The main class for storing **PoPs** data.
 */

/* *********************************************************************************************************//**
 * Database constructor for an initial empty **PoPs** database.
 ***********************************************************************************************************/

Database::Database( ) : 
        m_gaugeBosons( PoPI_gaugeBosonsChars ),
        m_leptons( PoPI_leptonsChars ),
        m_baryons( PoPI_baryonsChars ),
        m_chemicalElements( PoPI_chemicalElementsChars ),
        m_unorthodoxes( PoPI_unorthodoxesChars ) {

}

/* *********************************************************************************************************//**
 * Database constructor for a **PoPs** database with data read from the file *a_fileName*.
 *
 * @param a_fileName                    [in]    The **PoPs** file to read in.
 ***********************************************************************************************************/

Database::Database( std::string const &a_fileName ) : 
        m_gaugeBosons( PoPI_gaugeBosonsChars ),
        m_leptons( PoPI_leptonsChars ),
        m_baryons( PoPI_baryonsChars ),
        m_chemicalElements( PoPI_chemicalElementsChars ),
        m_unorthodoxes( PoPI_unorthodoxesChars ) {

    addFile( a_fileName, false );
}

/* *********************************************************************************************************//**
 * Database constructor for a **PoPs** database with data read from a **HAPI::Node** instance. This method is mainly
 * for internal use.
 *
 * @param a_database                    [in]    A **HAPI::Node** instance containing the data to parse.
 ***********************************************************************************************************/

Database::Database( HAPI::Node const &a_database ) :
        m_gaugeBosons( PoPI_gaugeBosonsChars ),
        m_leptons( PoPI_leptonsChars ),
        m_baryons( PoPI_baryonsChars ),
        m_chemicalElements( PoPI_chemicalElementsChars ),
        m_unorthodoxes( PoPI_unorthodoxesChars ) {

    addDatabase( a_database, false );
}

/* *********************************************************************************************************//**
 * Adds the contents of the file *a_fileName* to *this*.
 *
 * @param a_fileName                    [in]    The **PoPs** file to get data from.
 * @param a_warnIfDuplicate             [in]    This argument is currently not used.
 ***********************************************************************************************************/

void Database::addFile( std::string const &a_fileName, bool a_warnIfDuplicate ) {

    addFile( a_fileName.c_str( ), a_warnIfDuplicate );
}

/* *********************************************************************************************************//**
 * Adds the contents of the file *a_fileName* to *this*.
 *
 * @param a_fileName                    [in]    The **PoPs** file to get data from.
 * @param a_warnIfDuplicate             [in]    This argument is currently not used.
 ***********************************************************************************************************/

void Database::addFile( char const *a_fileName, bool a_warnIfDuplicate ) {

    HAPI::File *doc = new HAPI::PugiXMLFile( a_fileName, "Database::addFile" );
    HAPI::Node database = doc->first_child( );
    addDatabase( database, a_warnIfDuplicate );
    delete doc;
}

/* *********************************************************************************************************//**
 * Adds the contents of the *a_string* to *this*. *a_string* must be an XML string
 * starting with an **PoPs** XML node (i.e., element).
 *
 * @param a_string                      [in]    A *std::string* instance of **PoPs** data in an XML format.
 * @param a_warnIfDuplicate             [in]    This argument is currently not used.
 ***********************************************************************************************************/

void Database::addDatabase( std::string const &a_string, bool a_warnIfDuplicate ) {

    // a_string must contain a complete & well-formed XML document
    pugi::xml_document doc;

    pugi::xml_parse_result result = doc.load_string( a_string.c_str( ) );
    if( result.status != pugi::status_ok ) {
        char Msg[MsgSize+1];

        snprintf( Msg, MsgSize, "ERROR: in file '%s' in method '%s': %s.", __FILE__, __func__, result.description( ) );
        throw Exception( Msg );
    }

    HAPI::PugiXMLNode *database_internal = new HAPI::PugiXMLNode(doc.first_child( ));
    HAPI::Node database(database_internal);
    addDatabase( database, a_warnIfDuplicate );
}

/* *********************************************************************************************************//**
 * Adds the contents of *a_database* to *this*. The top node of *a_database* must be a valid **PoPs** node.
 *
 * @param a_database                    [in]    The **HAPI::Node** node to be added to *this*.
 * @param a_warnIfDuplicate             [in]    This argument is currently not used.
 ***********************************************************************************************************/

void Database::addDatabase( HAPI::Node const &a_database, LUPI_maybeUnused bool a_warnIfDuplicate ) {

    if( a_database.name( ) != PoPI_PoPsChars ) throw Exception( "Node '" + a_database.name( ) + "' is not a 'PoPs' node." );

    LUPI::FormatVersion formatVersion( a_database.attribute( PoPI_formatChars ).value( ) );
    if( !supportedFormat( formatVersion ) ) throw Exception( "Invalid format '" + formatVersion.format( ) + " in file " + a_database.name( ) + "." );
    if( m_formatVersion.format( ) == "" ) m_formatVersion = formatVersion;

    if( m_name == "" ) m_name = a_database.attribute( PoPI_nameChars ).value( );
    if( m_version == "" ) m_version = a_database.attribute( PoPI_versionChars ).value( );

    for( HAPI::Node child = a_database.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
        std::string s_name( child.name( ) );

        if(      s_name == PoPI_gaugeBosonsChars ) {
            m_gaugeBosons.appendFromParentNode( child, this, this ); }
        else if( s_name == PoPI_leptonsChars ) {
            m_leptons.appendFromParentNode( child, this, this ); }
        else if( s_name == PoPI_baryonsChars ) {
            m_baryons.appendFromParentNode( child, this, this ); }
        else if( s_name == PoPI_chemicalElementsChars ) {
            m_chemicalElements.appendFromParentNode( child, this, this ); }
        else if( s_name == PoPI_unorthodoxesChars ) {
            m_unorthodoxes.appendFromParentNode( child, this, this ); }
        else if( s_name == PoPI_aliasesChars ) {
            parseAliases( child, this ); }
        else {
        }
    }

    std::vector<Alias *> unresolvedAliases2;
    for( std::vector<Alias *>::iterator iter = m_unresolvedAliases.begin( ); iter != m_unresolvedAliases.end( ); ++iter ) {
        auto pidIter = m_idsMap.find( (*iter)->pid( ) );            // Locate pid.

        if( pidIter == m_idsMap.end( ) ) {
            unresolvedAliases2.push_back( *iter ); }
        else {
            (*iter)->setPidIndex( pidIter->second );
        }
    }
    m_unresolvedAliases.clear( );
    for( auto iter = unresolvedAliases2.begin( ); iter != unresolvedAliases2.end( ); ++iter ) {
        m_unresolvedAliases.push_back( *iter );
    }
}

/* *********************************************************************************************************//**
 * For internal use only. This method parses a **PoPs** *aliases* node.
 *
 * @param a_node                        [in]    The **HAPI::Node** node to be parsed.
 * @param a_DB                          [in]    The **PoPI::Database** to add the alias data to.
 ***********************************************************************************************************/

static void parseAliases( HAPI::Node const &a_node, Database *a_DB ) {

    for( HAPI::Node child = a_node.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
        std::string name = child.name( );
        Alias *alias = nullptr;

        if( name == PoPI_aliasChars ) {
            alias = new Alias( child, a_DB ); }
        else if( name == PoPI_metaStableChars ) {
            alias = new MetaStable( child, a_DB ); }
        else if( name == PoPI_particleChars ) {         // Needed for GNDS 1.10.
            alias = new Alias( child, a_DB ); }
        else {
            throw Exception( "Node '" + name + "' not supported as a child node of PoPs/aliases." );
        }
        a_DB->addAlias( alias );
    }
}

/* *********************************************************************************************************//**
 * Destructor for a **PoPI::Database** instance.
 ***********************************************************************************************************/

Database::~Database( ) {

    for( std::vector<Alias *>::iterator iter = m_aliases.begin( ); iter != m_aliases.end( ); ++iter ) delete *iter;
}

/* *********************************************************************************************************//**
 * This method returns the list of ids for the aliases with unresolved pids.
 *
 * @return                                      A std::vector<std::string> of the ids of the aliases with unresolved pids.
 ***********************************************************************************************************/

std::vector<std::string> Database::unresolvedAliasIds( ) const {

    std::vector<std::string> ids;

    for( std::vector<Alias *>::const_iterator iter = m_unresolvedAliases.begin( ); iter != m_unresolvedAliases.end( ); ++iter ) {
        ids.push_back( (*iter)->ID( ) );
    }

    return( ids );
}

/* *********************************************************************************************************//**
 * Internally, **PoPI::Database** stores a unique integer (called an index) for each particle in *this*. This method returns the
 * the index for the specified particle.
 *
 * @param a_id                          [in]    The **PoPs** id for the specified particle.
 *
 * @return                                      The internal index for the specified particle.
 ***********************************************************************************************************/

std::size_t Database::operator[]( std::string const &a_id ) const {

    auto iter = m_idsMap.find( a_id );
    if( iter == m_idsMap.end( ) ) {
        std::string errorMessage( "particle '" + a_id + "' not in database -3." );
        throw Exception( errorMessage );
    }

    return( iter->second );
}

/* *********************************************************************************************************//**
 * Returns **true** if the specified index is value and **false** otherwise. This is, if a particle exists within
 * *this* with index *a_index*.
 *
 * @param a_index                       [in]    A particle index to test.
 *
 * @return                                      **true** is the specified index is valid and **false** otherwise.
 ***********************************************************************************************************/

bool Database::exists( std::size_t a_index ) const {

    if( a_index >= m_list.size( ) ) return( false );
    return( true );
}

/* *********************************************************************************************************//**
 * Returns **true** if the specified id exists within *this* and **false** otherwise.
 *
 * @param a_id                          [in]    A particle id to test.
 *
 * @return                                      **true** is the specified id exists in *this* and **false** otherwise.
 ***********************************************************************************************************/

bool Database::exists( std::string const &a_id ) const {

    auto iter = m_idsMap.find( a_id );
    return( iter != m_idsMap.end( ) );
}

/* *********************************************************************************************************//**
 * Returns **true** if the specified intid exists within *this* and **false** otherwise.
 *
 * @param a_intid                       [in]    A particle's intidd to test.
 *
 * @return                                      **true** is the specified intid exists in *this* and **false** otherwise.
 ***********************************************************************************************************/

bool Database::existsIntid( int a_intid ) const {

    return( m_intidsMap.find( a_intid ) != m_intidsMap.end( ) );
}

/* *********************************************************************************************************//**
 * Returns a std::vector of std::string's of all aliases in *this* that resolve to *a_id*.
 *
 * @param a_id                          [in]    A particle's id whose.
 *
 * @return                                      Vector of alias ids.
 ***********************************************************************************************************/

std::vector<std::string> Database::aliasReferences( std::string const &a_id ) {

    std::vector<std::string> ids;

    for( auto aliasIter = m_aliases.begin( ); aliasIter != m_aliases.end( ); ++aliasIter ) {
        if( final( (*aliasIter)->pid( ) ) == a_id ) ids.push_back( (*aliasIter)->ID( ) );
    }

    return( ids );
}

/* *********************************************************************************************************//**
 * This method resolves aliases to return an actual particle specified by *a_id*. That is, if *a_id* is an alias,
 * then its referenced particle is returned. However, if *a_returnAtMetaStableAlias* is **true** and a meta-stable
 * is found while resolving *a_id*, then the meta-stable id will be returned.
 *
 * @param a_id                          [in]    A particle's id whose resolved particle id is requested.
 * @param a_returnAtMetaStableAlias     [in]    If **true**, the resolving will stop if a meta-stable is found.
 *
 * @return                                      The revolved id for *a_id*.
 ***********************************************************************************************************/

std::string Database::final( std::string const &a_id, bool a_returnAtMetaStableAlias ) const {

    std::size_t index( final( (*this)[a_id], a_returnAtMetaStableAlias ) );

    return( m_list[index]->ID( ) );
}

/* *********************************************************************************************************//**
 * This method resolves aliases to return an actual particle specified by *a_index*. That is, if *a_index* is an alias,
 * then its referenced particle is returned. However, if *a_returnAtMetaStableAlias* is **true** and a meta-stable
 * is found while resolving *a_index*, the then meta-stable index will be returned.
 *
 * @param a_index                       [in]    A particle's index whose resolved particle index is requested.
 * @param a_returnAtMetaStableAlias     [in]    If **true**, the resolving will stop if a meta-stable is found.
 *
 * @return                                      The revolved index for *a_index*.
 ***********************************************************************************************************/

std::size_t Database::final( std::size_t a_index, bool a_returnAtMetaStableAlias ) const {

    while( isAlias( a_index ) ) {
        if( a_returnAtMetaStableAlias && isMetaStableAlias( a_index ) ) break;
        a_index = ((Alias *) m_list[a_index])->pidIndex( );
    }
    return( a_index );
}

/* *********************************************************************************************************//**
 * This method returns the chemical element symbol for *a_id* if it is a PoPs id for a chemicalElement, isotope, nuclide,
 * or nucleus object. Otherwise, it returns an empty string. The PoPs id *a_id* must be in *this*.
 *
 * @param a_id                          [in]    A particle's id whose chemical element symbol is requested.
 *
 * @return                                      The string for the chemical element symbol.
 ***********************************************************************************************************/

std::string Database::chemicalElementSymbol( std::string const &a_id ) const {

    std::string symbol1;
    Base const *base = nullptr;

    auto iter = m_idsMap.find( a_id );
    if( iter != m_idsMap.end( ) ) {
        std::string finalId = final( a_id );
        iter = m_idsMap.find( finalId );
        base = m_list[iter->second]; }
    else {
        auto iter2 = m_symbolMap.find( a_id );
        if( iter2 != m_symbolMap.end( ) ) base = m_symbolList[iter2->second];
    }

    if( base != nullptr ) {
        if( base->isNucleus( ) ) base = static_cast<Nucleus const *>( base )->nuclide( );
        if( base->isNuclide( ) ) base = static_cast<Nuclide const *>( base )->isotope( );
        if( base->isIsotope( ) ) base = static_cast<Isotope const *>( base )->chemicalElement( );
        if( base->isChemicalElement( ) ) symbol1 = base->ID( );
    }

    return( symbol1 );
}

/* *********************************************************************************************************//**
 * This method returns the isotope symbol for *a_id* if it is a PoPs id for an isotope, nuclide,
 * or nucleus object. Otherwise, it returns an empty string. The PoPs id *a_id* must be in *this*.
 *
 * @param a_id                          [in]    A particle's id whose isotope symbol is requested.
 *
 * @return                                      The string for the isopte symbol.
 ***********************************************************************************************************/

std::string Database::isotopeSymbol( std::string const &a_id ) const {
    
    std::string symbol1;
    Base const *base = nullptr;

    auto iter = m_idsMap.find( a_id );
    if( iter != m_idsMap.end( ) ) {
        std::string finalId = final( a_id );
        iter = m_idsMap.find( finalId );
        base = m_list[iter->second]; }
    else {
        auto iter2 = m_symbolMap.find( a_id );
        if( iter2 != m_symbolMap.end( ) ) base = m_symbolList[iter2->second];
    }
    
    if( base != nullptr ) {
        if( base->isNucleus( ) ) base = static_cast<Nucleus const *>( base )->nuclide( );
        if( base->isNuclide( ) ) base = static_cast<Nuclide const *>( base )->isotope( );
        if( base->isIsotope( ) ) symbol1 = base->ID( );
    }

    return( symbol1 );
}

/* *********************************************************************************************************//**
 * Returns the intid for particle *a_id* or -1 if *a_id* is not in *this*.
 *
 * @param a_id                          [in]    A particle's id whose intid is requested.
 *
 * @return                                      The intid for *a_id* or -1 if *a_id* not in *this*.
 ***********************************************************************************************************/

int Database::intid( std::string const &a_id ) const {

    int intid2 = -1;

    if( exists( a_id ) ) {
        Base const &base = get<Base const>( a_id );
        intid2 =  base.intid( );
    }

    return( intid2 );
}

/* *********************************************************************************************************//**
 * Returns the intid for particle with index *a_index* or -1 if *a_index* is not in *this*.
 *
 * @param a_index                       [in]    A particle's index whose index is requested.
 *
 * @return                                      The intid for *a_index* or -1 if *a_index* not in *this*.
 ***********************************************************************************************************/

int Database::intid( std::size_t a_index ) const {

    int intid2 = -1;

    if( exists( a_index ) ) {
        Base const &base = get<Base const>( a_index );

        if( base.isParticle( ) ) {
            IDBase const &idBase = static_cast<IDBase const &>( base );
            intid2 =  idBase.intid( );
        }
    }

    return( intid2 );
}

/* *********************************************************************************************************//**
 * Returns the index for particle with intid *a_intid* or -1 if *a_intid* is not in *this*.
 *
 * @param a_index                       [in]    A particle's index whose index is requested.
 *
 * @return                                      The intid for *a_index* or -1 if *a_index* not in *this*.
 ***********************************************************************************************************/

std::size_t Database::indexFromIntid( int a_intid ) const {

    auto iter = m_intidsMap.find( a_intid );
    if( iter == m_intidsMap.end( ) ) {
        throw Exception( "Intid " + LUPI::Misc::argumentsToString( "%d", a_intid ) + " not in pops" );
    }

    return( iter->second );
}

/* *********************************************************************************************************//**
 * This method adds a **PoPI::Base** instance to *this* and returns the unique index for it.
 *
 * @param a_item                        [in]    The **PoPI::Base** instance to add to *this*.
 *
 * @return                                      The index for the added **PoPI::Base** instance.
 ***********************************************************************************************************/

std::size_t Database::add( Base *a_item ) {

    std::size_t index = m_list.size( );

    m_idsMap[a_item->ID( )] = index;
    m_list.push_back( a_item );
    a_item->setIndex( index );

    if( a_item->intid( ) > 0 ) m_intidsMap[a_item->intid( )] = index;

    if( a_item->isAlias( ) ) m_unresolvedAliases.push_back( (Alias *) a_item );
    return( index );
}

/* *********************************************************************************************************//**
 * This method adds a **PoPI::SymbolBase** instance to *this* and returns the unique index for it.
 *
 * @param a_item                        [in]    The **PoPI::SymbolBase** instance to add to *this*.
 *
 * @return                                      The index for the added **PoPI::SymbolBase** instance.
 ***********************************************************************************************************/
    

std::size_t Database::addSymbol( SymbolBase *a_item ) {

    if( a_item->Class( ) == Particle_class::chemicalElement ) return( this->add( a_item ) );

    std::size_t index = m_symbolList.size( );

    m_symbolMap[a_item->symbol( )] = index;
    m_symbolList.push_back( a_item );
    a_item->setIndex( index );

    return( index );
}

/* *********************************************************************************************************//**
 * This method calculates nuclide gamma branching infomation and adds it to *a_nuclideGammaBranchStateInfos*.
 *
 * @param a_nuclideGammaBranchStateInfos [in]   The **NuclideGammaBranchStateInfos** instance to added nuclide gamma branching infomation to.
 * @param a_pops2                               A second PoPs used for storing GRIN added particles.
 * @param a_extraGammaBranchStates              Any additional nuclide needed by GRIN. Currently, one the capture residual.
 ***********************************************************************************************************/

void Database::calculateNuclideGammaBranchStateInfos( NuclideGammaBranchStateInfos &a_nuclideGammaBranchStateInfos, Database const *a_pops2,
                std::vector<std::string> &a_extraGammaBranchStates ) const {


    calculateNuclideGammaBranchStateInfos2( a_nuclideGammaBranchStateInfos );
    if( a_pops2 != nullptr ) {
        a_pops2->calculateNuclideGammaBranchStateInfos2( a_nuclideGammaBranchStateInfos );
        for( auto iter = a_extraGammaBranchStates.begin( ); iter != a_extraGammaBranchStates.end( ); ++iter ) {
            PoPI::Nuclide const &nuclide = a_pops2->get<PoPI::Nuclide>( *iter );
            nuclide.calculateNuclideGammaBranchStateInfos( *a_pops2, a_nuclideGammaBranchStateInfos, true );
        }
    }

    std::vector<NuclideGammaBranchStateInfo *> &nuclideGammaBranchStateInfos = a_nuclideGammaBranchStateInfos.nuclideGammaBranchStateInfos( );
    for( std::size_t i1 = 0; i1 < nuclideGammaBranchStateInfos.size( ); ++i1 ) {
        NuclideGammaBranchStateInfo *nuclideGammaBranchStateInfo = nuclideGammaBranchStateInfos[i1];

        nuclideGammaBranchStateInfo->calculateDerivedData( a_nuclideGammaBranchStateInfos );
    }
}

/* *********************************************************************************************************//**
 * This method calculates nuclide gamma branching infomation and adds it to *a_nuclideGammaBranchStateInfos*.
 *
 * @param a_nuclideGammaBranchStateInfos [in]    The **NuclideGammaBranchStateInfos** instance to added nuclide gamma branching infomation to.
 ***********************************************************************************************************/

void Database::calculateNuclideGammaBranchStateInfos2( NuclideGammaBranchStateInfos &a_nuclideGammaBranchStateInfos ) const {

    for( std::size_t i1 = 0; i1 <  m_chemicalElements.size( ); ++i1 ) {
        ChemicalElement const &chemicalElement = m_chemicalElements[i1];

        chemicalElement.calculateNuclideGammaBranchStateInfos( *this, a_nuclideGammaBranchStateInfos );
    }
}

/* *********************************************************************************************************//**
 * This method returns the mass of the particle/alias with id *a_id*.
 * Currently, *a_unit* is ignored and the mass is returned in unit of amu.
 *
 * @param a_id                          [in]    The PoPs id of the particle/alias.file to write *this* to.
 * @param a_unit                        [in]    The unit of the returned mass.
 ***********************************************************************************************************/

double Database::massValue( std::string const &a_id, std::string const &a_unit ) const {

    Particle const &particle2 = particle( final( a_id ) );

    return( particle2.massValue( a_unit ) );
}

/* *********************************************************************************************************//**
 * Writes an **XML** version of *this* to the file *a_fileName*.
 *
 * @param a_fileName                    [in]    The file to write *this* to.
 ***********************************************************************************************************/

void Database::saveAs( std::string const &a_fileName ) const {

    std::string indent1( "" );
    std::vector<std::string> XMLList;

    XMLList.push_back( "<?xml version=\"1.0\"?>" );
    toXMLList( XMLList, indent1 );

    std::ofstream fileio;
    fileio.open( a_fileName.c_str( ) );
    for( std::vector<std::string>::iterator iter = XMLList.begin( ); iter != XMLList.end( ); ++iter ) {
        fileio << *iter << std::endl;
    }
    fileio.close( );
}

/* *********************************************************************************************************//**
 * Adds the contents of *this* to *a_XMLList* where each item in *a_XMLList* is one line (without linefeeds) to output as an XML representation of *this*.
 *
 * @param a_XMLList                     [in]    The list to add an XML output representation of *this* to.
 * @param a_indent1                     [in]    The amount of indentation to added to each line added to *a_XMLList*.
 ***********************************************************************************************************/

void Database::toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const {

    std::string indent2 = a_indent1 + "  ";
    std::string indent3 = indent2 + "  ";

    std::string header1 = a_indent1 + "<PoPs name=\"" + m_name + "\" version=\"" + m_version + "\" format=\"" + m_formatVersion.format( ) + "\">";
    a_XMLList.push_back( std::move( header1 ) );

    if( m_aliases.size( ) > 0 ) {
        std::string header2 = indent2 + "<" + PoPI_aliasesChars + ">";
        a_XMLList.push_back( std::move( header2 ) );
        for( std::vector<Alias *>::const_iterator iter = m_aliases.begin( ); iter != m_aliases.end( ); ++iter )
            (*iter)->toXMLList( a_XMLList, indent3 );
        appendXMLEnd( a_XMLList, PoPI_aliasesChars );
    }
    m_gaugeBosons.toXMLList( a_XMLList, indent2 );
    m_leptons.toXMLList( a_XMLList, indent2 );
    m_baryons.toXMLList( a_XMLList, indent2 );
    m_unorthodoxes.toXMLList( a_XMLList, indent2 );
    m_chemicalElements.toXMLList( a_XMLList, indent2 );

    appendXMLEnd( a_XMLList, PoPI_PoPsChars );
}

/* *********************************************************************************************************//**
 * Prints a brief outline of the contents of *this*.
 *
 * @param a_printIndices                [in]    If **true**, each particles index is also printed.
 ***********************************************************************************************************/

void Database::print( bool a_printIndices ) {

    for( auto iter = m_idsMap.begin( ); iter != m_idsMap.end( ); ++iter ) {
        std::string label( iter->first );
        std::size_t index = iter->second;
        Base *item = m_list[index];
        std::string is_alias( "" );
        std::string mass( "" );

        if( item->isAlias( ) ) {
            is_alias = " is an alias (final is label = '";
            std::size_t finalIndex = final( index );
            IDBase const &myfinal = get<IDBase>( finalIndex );
            is_alias += std::string( myfinal.ID( ) );
            is_alias += std::string( "')" ); }
        else if( item->isParticle( ) ) {
            Particle *particle = (Particle *) item;

            try {
                double dmass = particle->massValue( "amu" );
                mass = LUPI::Misc::argumentsToString( "  mass = %e amu", dmass ); }
            catch (...) {
                mass = " particle has no mass data.";
            }
        }

        std::cout << iter->first << " (" << item->ID( ) << ") --> ";
        if( a_printIndices ) std::cout << index << " (" << item->index( ) << ")";
        std::cout << is_alias << mass << std::endl;
    }
}

}
