/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <ctype.h>

#include <RISI.hpp>

namespace GIDI {

namespace RISI {

static void readRIS2( std::string const &a_basePath, std::string const &a_fileName, Projectiles &a_projectiles, std::string const &a_energyUnit );

/*! \class Reaction
 * Class to store a reaction for a reaction information summary **RIS**.
 */

/* *********************************************************************************************************//**
 * @param a_effectiveThreshold      [in]    The effective threshold for the reaction.
 * @param a_products                [in]    The list of final products for the reaction.
 * @param a_multiplicities          [in]    The multiplicities for each product in *a_products*.
 * @param a_intermediates           [in]    The list of intermediates products for the reaction.
 * @param a_process                 [in]    The process for the reaction.
 ***********************************************************************************************************/

Reaction::Reaction( double a_effectiveThreshold, std::vector<std::string> const &a_products, std::vector<int> const &a_multiplicities,
                std::vector<std::string> const &a_intermediates, std::string const &a_process, std::string const &reactionLabel,
                std::string const &convarianceFlag ) :
        m_effectiveThreshold( a_effectiveThreshold ),
        m_products( a_products ),
        m_multiplicities( a_multiplicities ),
        m_intermediates( a_intermediates ),
        m_process( a_process ),
        m_reactionLabel( reactionLabel ),
        m_convarianceFlag( convarianceFlag ) {

}

/* *********************************************************************************************************//**
 * 
 *
 * @param  a_projectile         [in]    The **Projectile** instance for the requested projectile.
 * @param  a_level              [in]    The current recursive level.
 * @param  a_maxLevel           [in]    The maximum recursive level requested by the user.
 * @param  a_energyMax          [in]    Only reactions with effective thresholds less than this value are processed.
 * @param  a_products           [in]    The list to add additional products to.
 ***********************************************************************************************************/

void Reaction::products( double a_energyMax, std::set<std::string> &a_products ) const {

    if( m_effectiveThreshold >= a_energyMax ) return;

    for( auto productIter = m_products.begin( ); productIter != m_products.end( ); ++productIter ) a_products.insert( *productIter );
}

/*! \class Protare
 * Class to store a protare for a reaction information summary **RIS**.
 */

/* *********************************************************************************************************//**
 * @param a_projectile              [in]    The PoPs id for the projectile.
 * @param a_target                  [in]    The PoPs id for the target.
 * @param a_evaluation              [in]    The evaluation string for the protare.
 * @param a_protareEnergyUnit       [in]    The unit of energy in the **RIS** file.
 * @param a_requestedEnergyUnit     [in]    The unit of energy specified by the user.
 ***********************************************************************************************************/

Protare::Protare( std::string const &a_projectile, std::string const &a_target, std::string const &a_evaluation,
                std::string const &a_protareEnergyUnit, std::string const &a_requestedEnergyUnit ) :
        m_addMode( 0 ),
        m_projectile( a_projectile ),
        m_target( a_target ),
        m_evaluation( a_evaluation ),
        m_energyConversionFactor( 1.0 ) {

    if( a_protareEnergyUnit != a_requestedEnergyUnit ) {
        if( a_protareEnergyUnit == "eV" ) {
            if( a_requestedEnergyUnit != "MeV" ) throw "RISI::Protare: supported a_requestedEnergyUnit '" + a_requestedEnergyUnit + "'.";
            m_energyConversionFactor = 1e-6; }
        else if( a_protareEnergyUnit == "MeV" ) {
            if( a_requestedEnergyUnit != "eV" ) throw "RISI::Protare: supported a_requestedEnergyUnit '" + a_requestedEnergyUnit + "'.";
            m_energyConversionFactor = 1e6; }
        else {
            throw "RISI::Protare: supported a_protareEnergyUnit '" + a_requestedEnergyUnit + "'.";
        }
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Protare::~Protare( ) {

    for( auto reactionIter = m_reactions.begin( ); reactionIter != m_reactions.end( ); ++reactionIter ) delete *reactionIter;

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

void Protare::Oops( LUPI_maybeUnused std::vector<std::string> const &a_elements ) {

    throw "No mode has been set for adding to the Protare: " + m_projectile + " + " + m_target + ".";
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

void Protare::addAlias( std::vector<std::string> const &a_elements ) {

    m_aliases[a_elements[1]] = a_elements[0];
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

void Protare::addReaction( std::vector<std::string> const &a_elements ) {

    std::vector<std::string> productsString = LUPI::Misc::splitString( a_elements[0], '+', true );
    double effectiveThreshold = m_energyConversionFactor * std::stod( a_elements[1] );
    std::vector<std::string> intermediates = LUPI::Misc::splitString( a_elements[2], ':', true );

    std::string reactionLabel;
    std::string covarianceFlag;

    if( a_elements.size( ) > 5 ) {
        reactionLabel = a_elements[4];
        covarianceFlag = a_elements[5];
    }

    std::vector<std::string> products;
    std::vector<int> multiplicities;

    for( auto iter = productsString.begin( ); iter != productsString.end( ); ++iter ) {
        char *begin = const_cast<char *>( (*iter).c_str( ) ), *end = begin;
        long multiplicity = 1;

        if( isdigit( begin[0] ) ) {
            multiplicity = strtol( begin, &end, 10 );
        }
        products.push_back( end );
        multiplicities.push_back( static_cast<int>( multiplicity ) );
    }

    m_reactions.push_back( new Reaction( effectiveThreshold, products, multiplicities, intermediates, a_elements[3], reactionLabel, covarianceFlag ) );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

void Protare::add( std::vector<std::string> const &a_elements ) {

    if( m_addMode == 0 ) {
        Oops( a_elements ); }
    else if( m_addMode == 1 ) {
        addAlias( a_elements ); }
    else if( m_addMode == 2 ) {
        addReaction( a_elements );
    }
}

/* *********************************************************************************************************//**
 *
 * @param  a_projectile         [in]    The **Projectile** instance for the requested projectile.
 * @param  a_level              [in]    The current recursive level.
 * @param  a_maxLevel           [in]    The maximum recursive level requested by the user.
 * @param  a_energyMax          [in]    Only reactions with effective thresholds less than this value are processed.
 * @param  a_products           [in]    The list to add additional products to.
 ***********************************************************************************************************/

void Protare::products( Projectile const *a_projectile, int a_level, int a_maxLevel, double a_energyMax, std::map<std::string, int> &a_products ) const {

    std::set<std::string> productSet;
    for( auto reactionIter = m_reactions.begin( ); reactionIter != m_reactions.end( ); ++reactionIter )
        (*reactionIter)->products( a_energyMax, productSet );

    for( auto productIter = productSet.begin( ); productIter != productSet.end( ); ++productIter ) {
        a_projectile->products( *productIter, a_level, a_maxLevel, a_energyMax, a_products );
    }
}

/*! \class Target
 * Stores a list of **Protare** instances for a specified target.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Target::~Target( ) {

    for( auto iter = m_protares.begin( ); iter != m_protares.end( ); ++iter ) delete *iter;

}

/* *********************************************************************************************************//**
 * Adds *a_protare* to *this* list of **Protare** instances.
 *
 * @param  a_protare            [in]    The **Protare** instance to add to *this*.
 ***********************************************************************************************************/

void Target::add( Protare *a_protare ) {

    m_protares.push_back( a_protare );
}

/* *********************************************************************************************************//**
 *
 * @param  a_projectile         [in]    The **Projectile** instance for the requested projectile.
 * @param  a_level              [in]    The current recursive level.
 * @param  a_maxLevel           [in]    The maximum recursive level requested by the user.
 * @param  a_energyMax          [in]    Only reactions with effective thresholds less than this value are processed.
 * @param  a_products           [in]    The list to add additional products to.
 ***********************************************************************************************************/

void Target::products( Projectile const *a_projectile, int a_level, int a_maxLevel, double a_energyMax, std::map<std::string, int> &a_products ) const {

    m_protares[0]->products( a_projectile, a_level, a_maxLevel, a_energyMax, a_products );
}

/* *********************************************************************************************************//**
 * Calls the *print* method on each **Projectile** in *this*.
 *
 * @param a_indent      [in]    The **Protare** instance to add to *this*.
 ***********************************************************************************************************/

void Target::print( std::string const &a_indent ) const {

    std::cout << a_indent + m_id << ":";
    for( auto iter = m_protares.begin( ); iter != m_protares.end( ); ++iter ) std::cout << " " << (*iter)->evaluation( );
    std::cout << std::endl;
}

/*! \class Projectile
 * Stores a list of projectiles and their associated **Target** instance.
 */


/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Projectile::~Projectile( ) {

    for( auto iter = m_targets.begin( ); iter != m_targets.end( ); ++iter ) delete (*iter).second;

}

/* *********************************************************************************************************//**
 * Adds *a_protare* to the associated target of *this*.
 *
 * @param  a_protare            [in]    The **Protare** instance to add to *this*.
 ***********************************************************************************************************/

void Projectile::add( Protare *a_protare ) {

    std::string const &target = a_protare->target( );

    auto iter = m_targets.find( target );
    if( iter == m_targets.end( ) ) {
         m_targets[target] = new Target( target );
        iter = m_targets.find( target );
    }

    (*iter).second->add( a_protare );
}

/* *********************************************************************************************************//**
 *
 * @param  a_level              [in]    The current recursive level.
 * @param  a_maxLevel           [in]    The maximum recursive level requested by the user.
 * @param  a_energyMax          [in]    Only reactions with effective thresholds less than this value are processed.
 * @param  a_products           [in]    The list to add additional products to.
 ***********************************************************************************************************/

void Projectile::products( std::string const &a_target, int a_level, int a_maxLevel, double a_energyMax, std::map<std::string, int> &a_products ) const {

    auto protductIter = a_products.find( a_target );
    if( protductIter != a_products.end( ) ) {
        if( a_level > (*protductIter).second ) return; }
    else {
        a_products[a_target] = a_level;             // Adds a_target to a_products.
    }

    if( a_level >= a_maxLevel ) return;

    auto targetIter = m_targets.find( a_target );
    if( targetIter != m_targets.end( ) ) (*targetIter).second->products( this, a_level + 1, a_maxLevel, a_energyMax, a_products );
}

/* *********************************************************************************************************//**
 * Prints *this* id and then calls the *print* method on each **Target** in *this*.
 *
 * @param a_indent      [in]    The **Protare** instance to add to *this*.
 ***********************************************************************************************************/

void Projectile::print( std::string const &a_indent ) const {

    std::cout << a_indent + m_id << std::endl;
    for( auto iter = m_targets.begin( ); iter != m_targets.end( ); ++iter ) (*iter).second->print( a_indent + "  " );
}

/*! \class Projectiles
 * Stores a list of projectiles.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Projectiles::~Projectiles( ) {

    clear( );
}

/* *********************************************************************************************************//**
 * Adds *a_protare* to the associated projectile of *this*.
 *
 * @param  a_protare            [in]    The **Protare** instance to add to *this*.
 ***********************************************************************************************************/

void Projectiles::add( Protare *a_protare ) {

    std::string const &projectile = a_protare->projectile( );

    auto iter = m_projectiles.find( projectile );
    if( iter == m_projectiles.end( ) ) {
         m_projectiles[projectile] = new Projectile( projectile );
        iter = m_projectiles.find( projectile );
    }

    (*iter).second->add( a_protare );
}

/* *********************************************************************************************************//**
 * Clears the contents of the *m_projectiles* member.
 ***********************************************************************************************************/

void Projectiles::clear( ) {

    for( auto iter = m_projectiles.begin( ); iter != m_projectiles.end( ); ++iter ) delete (*iter).second;
    m_projectiles.clear( );
}

/* *********************************************************************************************************//**
 *
 * @param  a_level              [in]    The current recursive level.
 * @param  a_maxLevel           [in]    The maximum recursive level requested by the user.
 * @param  a_energyMax          [in]    Only reactions with effective thresholds less than this value are processed.
 * @param  a_products           [in]    The list to add additional products to.
 ***********************************************************************************************************/

std::vector<std::string> Projectiles::products( std::string const &a_projectile, std::vector<std::string> const &a_seedTargets, int a_maxLevel, 
                    double a_energyMax ) const {

    std::map<std::string, int> productMap;

    auto projectile = m_projectiles.find( a_projectile );
    if( projectile != m_projectiles.end( ) ) {
        for( auto targetIter = a_seedTargets.begin( ); targetIter != a_seedTargets.end( ); ++targetIter )
            (*projectile).second->products( (*targetIter), 0, a_maxLevel, a_energyMax, productMap );
    }

    std::vector<std::string> productList;
    for( auto productIter = productMap.begin( ); productIter != productMap.end( ); ++productIter ) productList.push_back( (*productIter).first );

    return( productList );
}

/* *********************************************************************************************************//**
 * Calls the *print* method on each **Projectile** in *this*.
 *
 * @param a_indent      [in]    The **Protare** instance to add to *this*.
 ***********************************************************************************************************/

void Projectiles::print( std::string const &a_indent ) const {

    for( auto iter = m_projectiles.begin( ); iter != m_projectiles.end( ); ++iter ) (*iter).second->print( a_indent );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

void readRIS( std::string const &a_fileName, std::string const &a_energyUnit, Projectiles &a_projectiles ) {

    a_projectiles.clear( );
    readRIS2( ".", a_fileName, a_projectiles, a_energyUnit );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

static void readRIS2( std::string const &a_basePath, std::string const &a_fileName, Projectiles &a_projectiles, std::string const &a_energyUnit ) {

    std::string errorString;

    std::string fileName( a_fileName );
    if( fileName[0] != '/' ) {                // Only works on Unix like systems.
        fileName = a_basePath + "/" + fileName;
    }

    std::ifstream inputFile;
    inputFile.open( fileName );
    if( !inputFile.good( ) ) throw LUPI::Exception( "Opening RIS file " + a_fileName + " failed." );

    std::string line;
    try {
        if( getline( inputFile, line ) ) {
            std::vector<std::string> elements = LUPI::Misc::splitString( line, ':', true );
            if( elements.size( ) != 2 ) throw LUPI::Exception( "Invalid header line in RIS file '" + fileName + "'" );

            if( elements[0] != "#ris" ) throw LUPI::Exception( "Invalid header tag in RIS file '" + fileName + "'" );

            if( elements[1] != "1.0" ) throw LUPI::Exception( "Invalid header version in RIS file '" + fileName + "'" );

            Protare *protare = nullptr;
            while( getline( inputFile, line ) ) {
                elements = LUPI::Misc::splitString( line, ':', true );
                if( elements.size( ) == 0 ) continue;
                std::string command = elements[0];

                if( command == "#import" ) {
                    readRIS2( LUPI::FileInfo::_dirname( fileName ), elements[1], a_projectiles, a_energyUnit );
                    protare = nullptr; }
                else if( command == "#protare" ) {
                    protare = new Protare(elements[1], elements[2], elements[3], elements[4], a_energyUnit);
                    a_projectiles.add( protare ); }
                else if( command == "#aliases" ) {
                    protare->setAddingAliases( ); }
                else if( command == "#reactions" ) {
                    protare->setAddingReactions( ); }
                else {
                    if( protare == nullptr ) throw LUPI::Exception( "Data without protare defined in RIS file '" + fileName + "'." );
                    protare->add( elements );
                }
            }
        } }
    catch (...) {
        inputFile.close( );
        throw;
    }
}

}           // End of namespace RISI.

}           // End of namespace GIDI.
