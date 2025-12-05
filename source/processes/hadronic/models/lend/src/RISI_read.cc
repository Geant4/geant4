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
#include <iomanip>
#include <fstream>
#include <ctype.h>
#include <algorithm>

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
 * Returns true if the reaction process contains 'fission' (case-insensitive).
 *
 * @return True if this is a fission reaction, false otherwise.
 ***********************************************************************************************************/

bool Reaction::isFission( ) const {

    std::string processLower = m_process;
    std::transform(processLower.begin(), processLower.end(), processLower.begin(), ::tolower);
    return processLower.find("fission") != std::string::npos;
}

/* *********************************************************************************************************//**
 * Returns the multiplicity for the specified product ID.
 *
 * @param a_productId           [in]    The product ID to look up.
 * @return The multiplicity of the specified product, -1 if multiplicity is energy-dependent, or 0 if not found.
 ***********************************************************************************************************/

int Reaction::multiplicity( std::string const &a_productId ) const {

    for( std::size_t index = 0; index < m_products.size( ); ++index ) {
        if( m_products[index] == a_productId ) {
            if( m_multiplicities[index] == 0 ) {
                return -1;
            }
            return m_multiplicities[index];
        }
    }
    return 0;
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

/* *********************************************************************************************************//**
 * This method attempts to print *this* as it appears in a file.
 ***********************************************************************************************************/

void Reaction::printAsRIS_file( int a_labelWidth ) const {

    std::string sep = "";
    std::string productList;
    for( std::size_t index = 0; index < m_products.size( ); ++index ) {
        std::string product = m_products[index];
        if( m_multiplicities[index] != 1 ) product = LUPI::Misc::argumentsToString( "%d", m_multiplicities[index] ) + product;
        productList += sep + product;
        sep = " + ";
    }

    sep = "";
    std::string intermediates;
    for( auto iterIntermediate = m_intermediates.begin( ); iterIntermediate != m_intermediates.end( ); ++iterIntermediate ) {
        intermediates += sep + (*iterIntermediate);
        sep = " ";
    }

    std::string effectiveThresholdStr = LUPI::Misc::doubleToShortestString( m_effectiveThreshold, 15, -3 );
    if( effectiveThresholdStr.find( "e" ) == std::string::npos ) {
        if( effectiveThresholdStr.find( "." ) == std::string::npos ) {
            effectiveThresholdStr += ".0";
        } }
    else {
        std::size_t eMinus = effectiveThresholdStr.find( "e-" );
        if( ( eMinus != std::string::npos ) && ( effectiveThresholdStr.size( ) == eMinus + 3 ) ) {
            effectiveThresholdStr = effectiveThresholdStr.replace( eMinus, 2, "e-0" );
        }
    }
    
    std::cout << "    " << std::left << std::setw( 31 ) << productList
            << " : " << std::left << std::setw( 12 ) << effectiveThresholdStr
            << " : " << std::left << std::setw( 10 ) << intermediates 
            << " : " << std::left << std::setw( 9 ) << m_process 
            << " : " << std::left << std::setw( a_labelWidth ) << m_reactionLabel 
            << " : " << m_convarianceFlag << std::endl;
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
        m_energyUnit( a_protareEnergyUnit ),
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

    std::pair<std::string, std::string> keyName = { a_elements[0], a_elements[1] };
    m_aliases.push_back( keyName );
}

/* *********************************************************************************************************//**
 * Returns true if any reaction in this protare is a fission reaction.
 *
 * @return True if fission reactions are present, false otherwise.
 ***********************************************************************************************************/

bool Protare::fissionPresent( ) const {

    for( auto reactionIter = m_reactions.begin( ); reactionIter != m_reactions.end( ); ++reactionIter ) {
        if( (*reactionIter)->isFission( ) ) return true;
    }
    return false;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

void Protare::addReaction( std::vector<std::string> const &a_elements ) {

    std::vector<std::string> productsString = LUPI::Misc::splitString( a_elements[0], '+', true );
    double effectiveThreshold = m_energyConversionFactor * std::stod( a_elements[1] );
    std::vector<std::string> intermediates = LUPI::Misc::splitString( a_elements[2], ' ', true );

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

/* *********************************************************************************************************//**
 * This method attempts to print *this* as it appears in a file.
 ***********************************************************************************************************/

void Protare::printAsRIS_file( ) const {

    std::size_t labelWidth = 0;
    for( auto iter = m_reactions.begin( ); iter != m_reactions.end( ); ++iter ) {
        if( labelWidth < (*iter)->m_reactionLabel.size( ) ) labelWidth = (*iter)->m_reactionLabel.size( );
    }

    std::cout << "#protare : " << m_projectile << " : " << m_target << " : " << m_evaluation << " : " << m_energyUnit << std::endl;
    if( m_aliases.size( ) > 0 ) {
        std::cout << "#aliases : " << m_aliases.size( ) << std::endl;
        for( auto iter = m_aliases.begin( ); iter != m_aliases.end( ); ++iter ) {
            std::cout << "    " << iter->first << " : " << iter->second << std::endl;
        }
    }
    std::cout << "#reactions : " << m_reactions.size( ) << std::endl;
    for( auto iter = m_reactions.begin( ); iter != m_reactions.end( ); ++iter ) (*iter)->printAsRIS_file( (int) labelWidth );
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
 * Returns true if any protare in this target has fission reactions.
 *
 * @return True if fission reactions are present, false otherwise.
 ***********************************************************************************************************/

bool Target::fissionPresent( ) const {

    for( auto protareIter = m_protares.begin( ); protareIter != m_protares.end( ); ++protareIter ) {
        if( (*protareIter)->fissionPresent( ) ) return true;
    }
    return false;
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

/* *********************************************************************************************************//**
 * Calls the *printAsRIS_file* method on each **Protare** in *this*. This method attempts to print *this* as it appears in a file.
 ***********************************************************************************************************/

void Target::printAsRIS_file( ) const {

    for( auto iter = m_protares.begin( ); iter != m_protares.end( ); ++iter ) (*iter)->printAsRIS_file( );
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
 * Returns true if any target in this projectile has fission reactions.
 *
 * @return True if fission reactions are present, false otherwise.
 ***********************************************************************************************************/

bool Projectile::fissionPresent( std::vector<std::string> targetIds ) const {

    for( auto& targetId : targetIds ) {
        auto iter = m_targets.find( targetId );
        if ( iter == m_targets.end( ) ) {
            throw LUPI::Exception( "Target " + targetId + " missing from .ris file for projectile " + m_id );
        }
        if( (*iter).second->fissionPresent( ) ) return true;
    }
    return false;
}

/* *********************************************************************************************************//**
 * Returns a pointer to the target with the specified name, or nullptr if not found.
 *
 * @param a_targetName     [in]    The name of the target to find.
 * @return A pointer to the Target object, or nullptr if not found.
 ***********************************************************************************************************/

Target const *Projectile::target( std::string const &a_targetName ) const {

    auto targetIter = m_targets.find( a_targetName );
    if( targetIter != m_targets.end( ) ) {
        return (*targetIter).second;
    }
    return nullptr;
}

/* *********************************************************************************************************//**
 * Returns a vector containing the names of all targets available for this projectile.
 *
 * @return A vector of target names.
 ***********************************************************************************************************/

std::vector<std::string> Projectile::targetIds( ) const {

    std::vector<std::string> targetList;
    
    for( auto targetIter = m_targets.begin( ); targetIter != m_targets.end( ); ++targetIter ) {
        targetList.push_back( (*targetIter).first );
    }
    
    return( targetList );
}

/* *********************************************************************************************************//**
 * Populate std::map with product id: max multiplicity for that product that can be created from the given target.
 *
 * @param  a_target             [in]    Target particle id.
 * @param  a_level              [in]    The current recursive level.
 * @param  a_maxLevel           [in]    The maximum recursive level requested by the user.
 * @param  a_energyMax          [in]    Only reactions with effective thresholds less than this value are processed.
 * @param  a_products           [in]    The map to be populated with (product: max multiplicity) pairs.
 ***********************************************************************************************************/

void Projectile::products( std::string const &a_target, int a_level, int a_maxLevel, double a_energyMax, std::map<std::string, int> &a_products ) const {

    auto productIter = a_products.find( a_target );
    if( productIter != a_products.end( ) ) {
        if( a_level < (*productIter).second ) {
            // found a way to make the product in fewer reaction steps
            a_products[a_target] = a_level;
        } else {
            return;
        }
    } else {
        a_products[a_target] = a_level;             // Adds a_target to a_products.
    }

    if( a_level >= a_maxLevel ) return;

    auto targetIter = m_targets.find( a_target );
    if( targetIter != m_targets.end( ) ) (*targetIter).second->products( this, a_level + 1, a_maxLevel, a_energyMax, a_products );
}

/* *********************************************************************************************************//**
 * Filter an initial list of products, returning only those that are available as targets for this projectile
 *
 * @param  a_products           [in]    The list of product IDs to filter.
 * @return A vector of filtered product IDs.
 ***********************************************************************************************************/

std::vector<std::string> Projectile::filterProducts( std::vector<std::string> const &a_productIds ) const {

    std::vector<std::string> filteredProducts;

    for (const auto &productId : a_productIds) {
        if (m_targets.find(productId) != m_targets.end()) {
            filteredProducts.push_back(productId);
        }
    }

    return filteredProducts;
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

/* *********************************************************************************************************//**
 * Calls the *printAsRIS_file* method on each **Target** in *this*. This method attempts to print *this* as it appears in a file.
 *
 ***********************************************************************************************************/

void Projectile::printAsRIS_file( ) const {

    for( auto iter = m_targets.begin( ); iter != m_targets.end( ); ++iter ) (*iter).second->printAsRIS_file( );
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
 * Returns a vector containing the particle ids of all projectiles in this collection.
 *
 * @return A vector of projectile ids.
 ***********************************************************************************************************/

std::vector<std::string> Projectiles::projectileIds( ) const {

    std::vector<std::string> projectileList;
    
    for( auto projectileIter = m_projectiles.begin( ); projectileIter != m_projectiles.end( ); ++projectileIter ) {
        projectileList.push_back( (*projectileIter).first );
    }
    
    return( projectileList );
}

/* *********************************************************************************************************//**
 * Returns a pointer to the projectile with the specified name, or nullptr if not found.
 *
 * @param a_projectile  [in]    The name of the projectile to find.
 * @return A pointer to the Projectile object, or nullptr if not found.
 ***********************************************************************************************************/

Projectile const *Projectiles::projectile( std::string const &a_projectile ) const {

    auto projectileIter = m_projectiles.find( a_projectile );
    if( projectileIter != m_projectiles.end( ) ) {
        return( (*projectileIter).second );
    }
    
    return( nullptr );
}

/* *********************************************************************************************************//**
 * Return a list of products that can be created for the given projectile and initial 'seed' targets.
 *
 * @param  a_projectile         [in]    Projectile particle id (e.g. 'n' or 'photon')
 * @param  a_seedTargets        [in]    List of initial targets, used to determine what products can be created.
 * @param  a_maxLevel           [in]    The maximum recursive level requested by the user.
 * @param  a_energyMax          [in]    Only reactions with effective thresholds less than this value are processed.
 * @param  a_onlyIncludeTargets [in]    Only return products that correspond to protares.
 ***********************************************************************************************************/

std::vector<std::string> Projectiles::products( std::string const &a_projectile, std::vector<std::string> const &a_seedTargets, int a_maxLevel, 
                    double a_energyMax, bool a_onlyIncludeTargets ) const {

    std::map<std::string, int> productMap;

    auto projectile = m_projectiles.find( a_projectile );
    if( projectile != m_projectiles.end( ) ) {
        for( auto targetIter = a_seedTargets.begin( ); targetIter != a_seedTargets.end( ); ++targetIter )
            (*projectile).second->products( (*targetIter), 0, a_maxLevel, a_energyMax, productMap );
    }

    std::vector<std::string> productList;
    for( auto productIter = productMap.begin( ); productIter != productMap.end( ); ++productIter ) productList.push_back( (*productIter).first );

    if( a_onlyIncludeTargets ) {
        productList = (*projectile).second->filterProducts( productList );
    }

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
 * Calls the *printAsRIS_file* method on each **Projectile** in *this*. This method attempts to print *this* as it appears in a file.
 *
 ***********************************************************************************************************/

void Projectiles::printAsRIS_file( ) const {

    std::cout << "#ris : 1.0" << std::endl;
    for( auto iter = m_projectiles.begin( ); iter != m_projectiles.end( ); ++iter ) (*iter).second->printAsRIS_file( );
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
    fileName = LUPI::FileInfo::realPath( fileName );

    std::ifstream inputFile;
    inputFile.open( fileName );
    if( !inputFile.good( ) ) throw LUPI::Exception( "Opening RIS file " + a_fileName + " failed." );

    std::string line;
    try {
        if( getline( inputFile, line ) ) {
            std::vector<std::string> elements = LUPI::Misc::splitString( line, ':', true );
            if( elements.size( ) != 2 ) throw LUPI::Exception( "Invalid header line in RIS file '" + fileName + "'" );

            if( elements[0] != "#ris" ) throw LUPI::Exception( "Invalid header tag in RIS file '" + fileName + "'" );

            std::string formatMajor = elements[1].substr(0, elements[1].find('.'));
            if( formatMajor != "1" ) throw LUPI::Exception( "Invalid header version in RIS file '" + fileName + "'" );

            std::string sep = " : ";
            if( line.find( sep ) == std::string::npos ) sep = ": ";

            Protare *protare = nullptr;
            while( getline( inputFile, line ) ) {
                elements = LUPI::Misc::splitString( line, sep, true );
                if( elements.size( ) == 0 ) continue;
                std::string command = elements[0];

                if( command == "#import" ) {
                    readRIS2( LUPI::FileInfo::_dirname( fileName ), elements[1], a_projectiles, a_energyUnit );
                    protare = nullptr; }
                else if( command == "#protare" ) {
                    protare = new Protare( elements[1], elements[2], elements[3], elements[4], a_energyUnit );
                    a_projectiles.add( protare ); }
                else if( command == "#aliases" ) {
                    if( protare == nullptr ) throw LUPI::Exception( "Aliases without protare defined in RIS file '" + fileName + "'." );
                    protare->setAddingAliases( ); }
                else if( command == "#reactions" ) {
                    if( protare == nullptr ) throw LUPI::Exception( "Reactions without protare defined in RIS file '" + fileName + "'." );
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
