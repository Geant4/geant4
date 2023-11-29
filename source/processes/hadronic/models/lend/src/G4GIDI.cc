//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include <iostream>

#include "G4GIDI.hh"

using namespace std;
using namespace GIDI;

/*
***************************************************************
*/
G4GIDI::G4GIDI( G4int ip, string &dataDirectory ) {

    init( ip );
    addDataDirectory( dataDirectory );
}
/*
***************************************************************
*/
G4GIDI::G4GIDI( G4int ip, list<string> &dataDirectoryList ) {

    init( ip );
    for( auto iter = dataDirectoryList.begin( ); iter != dataDirectoryList.end( ); ++iter )
      addDataDirectory( *iter );
}
/*
***************************************************************
*/
G4GIDI::~G4GIDI( void ) {

    G4GIDI_target *target;
    auto iter = dataDirectories.cbegin();

    while( targets.size( ) > 0 ) {
        target = targets.back( );
        targets.pop_back( );
        delete target;
    } // Loop checking, 11.06.2015, T. Koi

    while( iter != dataDirectories.cend() ) {
        delete *iter;
        dataDirectories.pop_front( );
    }// Loop checking, 11.06.2015, T. Koi
}
/*
***************************************************************
*/
G4int G4GIDI::init( G4int ip ) {

    projectileID = ip;
    if( ip == 0 ) {
        projectile = string( "g" ); }
    else if( ip == 1 ) {
        projectile = string( "n" ); }
    else if( ip == 2 ) {
        projectile = string( "p" ); }
    else if( ip == 3 ) {
        projectile = string( "d" ); }
    else if( ip == 4 ) {
        projectile = string( "t" ); }
    else if( ip == 5 ) {
        projectile = string( "h" ); }
    else if( ip == 6 ) {
        projectile = string( "a" ); }
    else {
        printf( "Invalid projectile ID = %d\n", ip );
        throw 1;
    }
    return( 0 );
}
/*
***************************************************************
*/
G4int G4GIDI::numberOfDataDirectories( void ) {

    return (G4int)dataDirectories.size( );
}
/*
***************************************************************
*/
G4int G4GIDI::addDataDirectory( string &dataDirectory ) {

    for( auto iter = dataDirectories.cbegin( ); iter != dataDirectories.cend( ); ++iter ) {
        if( (*iter)->path( ) == dataDirectory ) return( 0 );
    }

    G4GIDI_map *map = new G4GIDI_map( dataDirectory );
    dataDirectories.push_back( map );

    return( 0 );
}
/*
***************************************************************
*/
G4int G4GIDI::removeDataDirectory( string &dataDirectory ) {

    for( auto iter = dataDirectories.cbegin( ); iter != dataDirectories.cend( ); ++iter ) {
        if( dataDirectory == (*iter)->path( ) ) {
            
        }
    }
    return( 0 );
}
/*
***************************************************************
*/
string G4GIDI::getDataDirectoryAtIndex( G4int index ) {

    unsigned i = (unsigned) index;

    if( index >= 0 ) {
        if( i >= dataDirectories.size( ) ) return( "" );
        for( auto iter = dataDirectories.cbegin( ); iter != dataDirectories.cend( ); ++iter, --index )
          if( index == 0 ) return( (*iter)->fileName( ) );
    }

    return( "" );
}
/*
***************************************************************
*/
vector<string> *G4GIDI::getDataDirectories( void ) {

    std::size_t i = 0;
    vector<string> *v = new vector<string>( numberOfDataDirectories( ) );

    for( auto iter = dataDirectories.cbegin( ); iter != dataDirectories.cend( ); ++iter, ++i )
      (*v)[i] = string( (*iter)->fileName( ) );
    return( v );
}
/*
***************************************************************
*/
G4bool G4GIDI::isThisDataAvailable( string &lib_name, G4int iZ, G4int iA, G4int iM ) {

    G4bool b;
    char *targetName = G4GIDI_Misc_Z_A_m_ToName( iZ, iA, iM );

    if( targetName == nullptr ) return( false );
    string targetSymbol( targetName );
    b = isThisDataAvailable( lib_name, targetSymbol );
    smr_freeMemory( (void **) &targetName );
    return( b );
}
/*
***************************************************************
*/
G4bool G4GIDI::isThisDataAvailable( string &lib_name, string &targetName ) {

    char *path = dataFilename( lib_name, targetName );

    if( path != nullptr ) {
        smr_freeMemory( (void **) &path );
        return( true );
    }
    return( false );
}
/*
***************************************************************
*/
char *G4GIDI::dataFilename( string &lib_name, G4int iZ, G4int iA, G4int iM ) {

    char *targetName = G4GIDI_Misc_Z_A_m_ToName( iZ, iA, iM ), *fileName;

    if( targetName == nullptr ) return( nullptr );
    string targetSymbol( targetName );
    fileName = dataFilename( lib_name, targetSymbol );
    smr_freeMemory( (void **) &targetName );
    return( fileName );
}
/*
***************************************************************
*/
char *G4GIDI::dataFilename( string &lib_name, string &targetSymbol ) {

   char *path;

   for( auto iter = dataDirectories.cbegin( ); iter != dataDirectories.cend( ); ++iter )
      if( ( path = MCGIDI_map_findTarget( nullptr, (*iter)->map, lib_name.c_str(), projectile.c_str( ), targetSymbol.c_str( ) ) ) != nullptr )
         return( path );

   return( nullptr );
}
/*
***************************************************************
*/
vector<string> *G4GIDI::getNamesOfAvailableLibraries( G4int iZ, G4int iA, G4int iM ) {

    char *targetName = G4GIDI_Misc_Z_A_m_ToName( iZ, iA, iM );
    vector<string> *listOfLibraries;

    if( targetName == nullptr ) return( new vector<string>( ) );
    string targetSymbol( targetName );
    listOfLibraries = getNamesOfAvailableLibraries( targetSymbol );
    smr_freeMemory( (void **) &targetName );
    return( listOfLibraries );
}
/*
***************************************************************
*/
vector<string> *G4GIDI::getNamesOfAvailableLibraries( string &targetName ) {

    vector<string> *listOfLibraries = new vector<string>( );

    MCGIDI_map *map;
    MCGIDI_mapEntry *entry;

    for( auto iter = dataDirectories.cbegin( ); iter != dataDirectories.cend( ); ++iter ) {
        map = MCGIDI_map_findAllOfTarget( &((*iter)->smr), (*iter)->map, projectile.c_str( ), targetName.c_str( ) );
        for( entry = MCGIDI_map_getFirstEntry( map ); entry != nullptr; entry = MCGIDI_map_getNextEntry( entry ) ) {
            listOfLibraries->push_back( entry->evaluation );
        }
        MCGIDI_map_free( nullptr, map );
    }
    return( listOfLibraries );
}
/*
***************************************************************
*/
vector<string> *G4GIDI::getNamesOfAvailableTargets( void ) {

    vector<string> *listOfTargets;

    listOfTargets = new vector<string>( );
    if( listOfTargets == nullptr ) return( nullptr );
    for( auto iter_map = dataDirectories.cbegin( ); iter_map != dataDirectories.cend( ); ++iter_map ) {
        if( MCGIDI_map_walkTree( nullptr, (*iter_map)->map, getNamesOfAvailableTargets_walker, (void *) listOfTargets ) != 0 ) {
            delete listOfTargets;
            return( nullptr );
        }
    }
    return( listOfTargets );
}
/*
***************************************************************
*/
G4GIDI_target *G4GIDI::readTarget( string &lib_name, G4int iZ, G4int iA, G4int iM, G4bool bind ) {

    char *targetName = G4GIDI_Misc_Z_A_m_ToName( iZ, iA, iM );
    G4GIDI_target *target;

    if( targetName == nullptr ) return( nullptr );
    string targetSymbol( targetName );
    target = readTarget( lib_name, targetSymbol, bind );
    smr_freeMemory( (void **) &targetName );
    return( target );
}
/*
***************************************************************
*/
G4GIDI_target *G4GIDI::readTarget( string &lib_name, string &targetName, G4bool bind ) {

    for( auto iter_targets = targets.cbegin( ); iter_targets != targets.cend( ); ++iter_targets ) {
        if( (*iter_targets)->name == targetName ) return( nullptr );
    }
    char *path = dataFilename( lib_name, targetName );
    if( path == nullptr ) return( nullptr );

    G4GIDI_target *target = new G4GIDI_target( path );
    if( bind ) targets.push_back( target );
    smr_freeMemory( (void **) &path );
    return( target );
}
/*
***************************************************************
*/
G4GIDI_target *G4GIDI::getAlreadyReadTarget( G4int iZ, G4int iA, G4int iM ) {

    char *targetName = G4GIDI_Misc_Z_A_m_ToName( iZ, iA, iM );
    G4GIDI_target *target;

    if( targetName == nullptr ) return( nullptr );
    string targetSymbol( targetName );
    target = getAlreadyReadTarget( targetSymbol );
    smr_freeMemory( (void **) &targetName );
    return( target );
}
/*
***************************************************************
*/
G4GIDI_target *G4GIDI::getAlreadyReadTarget( string &targetSymbol ) {

    for( auto iter_targets = targets.cbegin( ); iter_targets != targets.cend( ); ++iter_targets ) {
        if( ( (*iter_targets)->name == targetSymbol ) ) return( *iter_targets );
    }
    return( nullptr );
}
/*
***************************************************************
*/
G4int G4GIDI::freeTarget( G4GIDI_target *target ) {

    for( auto iter_targets = targets.cbegin( ); iter_targets != targets.cend( ); ++iter_targets ) {
        if( *iter_targets == target ) {
            targets.erase( iter_targets );
            delete target;
            return( 0 );
        }
    }
    return( 1 );
}
/*
***************************************************************
*/
G4int G4GIDI::freeTarget( G4int iZ, G4int iA, G4int iM ) {

    G4int status;
    char *targetName = G4GIDI_Misc_Z_A_m_ToName( iZ, iA, iM );

    if( targetName == nullptr ) return( 1 );
    string targetSymbol( targetName );
    status = freeTarget( targetSymbol );
    smr_freeMemory( (void **) &targetName );
    return( status );
}
/*
***************************************************************
*/
G4int G4GIDI::freeTarget( string &targetSymbol ) {

    for( auto iter_targets = targets.cbegin( ); iter_targets != targets.cend( ); ++iter_targets ) {
        if( (*iter_targets)->name == targetSymbol ) return( freeTarget( *iter_targets ) );
    }
    return( 1 );
}
/*
***************************************************************
*/
vector<string> *G4GIDI::getListOfReadTargetsNames( void ) {

    vector<string> *listOfTargets;

    listOfTargets = new vector<string>( );
    if( listOfTargets == nullptr ) return( nullptr );
    for( auto iter_targets = targets.cbegin( ); iter_targets != targets.cend( ); ++iter_targets ) {
        listOfTargets->push_back( *(*iter_targets)->getName( ) );
    }
    return( listOfTargets );
}
