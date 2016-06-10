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
/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <iostream>

#include "G4GIDI.hh"

using namespace std;
using namespace GIDI;

/*
***************************************************************
*/
G4GIDI::G4GIDI( int ip, string &dataDirectory ) {

    init( ip );
    addDataDirectory( dataDirectory );
}
/*
***************************************************************
*/
G4GIDI::G4GIDI( int ip, list<string> &dataDirectoryList ) {

    list<string>::iterator iter;

    init( ip );
    for( iter = dataDirectoryList.begin( ); iter != dataDirectoryList.end( ); ++iter ) addDataDirectory( *iter );
}
/*
***************************************************************
*/
G4GIDI::~G4GIDI( void ) {

    G4GIDI_target *target;
    list<G4GIDI_map *>::iterator iter;

    while( targets.size( ) > 0 ) {
        target = targets.back( );
        targets.pop_back( );
        delete target;
    } // Loop checking, 11.06.2015, T. Koi

    while( ( iter = dataDirectories.begin( ) ) != dataDirectories.end( ) ) {
        delete *iter;
        dataDirectories.pop_front( );
    }// Loop checking, 11.06.2015, T. Koi
}
/*
***************************************************************
*/
int G4GIDI::init( int ip ) {

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
int G4GIDI::numberOfDataDirectories( void ) {

    return( dataDirectories.size( ) );
}
/*
***************************************************************
*/
int G4GIDI::addDataDirectory( string &dataDirectory ) {

    list<G4GIDI_map *>::iterator iter;

    for( iter = dataDirectories.begin( ); iter != dataDirectories.end( ); iter++ ) {
        if( (*iter)->path( ) == dataDirectory ) return( 0 );
    }

    G4GIDI_map *map = new G4GIDI_map( dataDirectory );
    dataDirectories.push_back( map );

    return( 0 );
}
/*
***************************************************************
*/
int G4GIDI::removeDataDirectory( string &dataDirectory ) {

    list<G4GIDI_map *>::iterator iter;

    for( iter = dataDirectories.begin( ); iter != dataDirectories.end( ); iter++ ) {
        if( dataDirectory == (*iter)->path( ) ) {
            
        }
    }
    return( 0 );
}
/*
***************************************************************
*/
string G4GIDI::getDataDirectoryAtIndex( int index ) {

    list<G4GIDI_map *>::iterator iter;
    unsigned i = (unsigned) index;

    if( index >= 0 ) {
        if( i >= dataDirectories.size( ) ) return( "" );
        for( iter = dataDirectories.begin( ); iter != dataDirectories.end( ); iter++, index-- ) {
            if( index == 0 ) return( (*iter)->fileName( ) );
        }
    }

    return( "" );
}
/*
***************************************************************
*/
vector<string> *G4GIDI::getDataDirectories( void ) {

    int i = 0;
    list<G4GIDI_map *>::iterator iter;
    vector<string> *v = new vector<string>( numberOfDataDirectories( ) );

    for( iter = dataDirectories.begin( ); iter != dataDirectories.end( ); iter++, i++ ) (*v)[i] = string( (*iter)->fileName( ) );
    return( v );
}
/*
***************************************************************
*/
bool G4GIDI::isThisDataAvailable( string &lib_name, int iZ, int iA, int iM ) {

    bool b;
    char *targetName = G4GIDI_Misc_Z_A_m_ToName( iZ, iA, iM );

    if( targetName == NULL ) return( false );
    string targetSymbol( targetName );
    b = isThisDataAvailable( lib_name, targetSymbol );
    smr_freeMemory( (void **) &targetName );
    return( b );
}
/*
***************************************************************
*/
bool G4GIDI::isThisDataAvailable( string &lib_name, string &targetName ) {

    char *path = dataFilename( lib_name, targetName );

    if( path != NULL ) {
        smr_freeMemory( (void **) &path );
        return( true );
    }
    return( false );
}
/*
***************************************************************
*/
char *G4GIDI::dataFilename( string &lib_name, int iZ, int iA, int iM ) {

    char *targetName = G4GIDI_Misc_Z_A_m_ToName( iZ, iA, iM ), *fileName;

    if( targetName == NULL ) return( NULL );
    string targetSymbol( targetName );
    fileName = dataFilename( lib_name, targetSymbol );
    smr_freeMemory( (void **) &targetName );
    return( fileName );
}
/*
***************************************************************
*/
char *G4GIDI::dataFilename( string &lib_name, string &targetSymbol ) {

/*
    char *path;
    list<G4GIDI_map *>::iterator iter;

    for( iter = dataDirectories.begin( ); iter != dataDirectories.end( ); iter++ ) {
        if( ( path = MCGIDI_map_findTarget( &((*iter)->smr), (*iter)->map, lib_name.c_str( ), projectile.c_str( ), targetSymbol.c_str( ) ) ) != NULL ) {
            return( path );
        }
    }
    return( NULL );
*/
//150121
//
   char *path;
   list<G4GIDI_map *>::iterator iter;

   for( iter = dataDirectories.begin( ); iter != dataDirectories.end( ); iter++ ) {
      if( ( path = MCGIDI_map_findTarget( NULL, (*iter)->map, lib_name.c_str(), projectile.c_str( ), targetSymbol.c_str( ) ) ) != NULL ) {
         return( path );
      }
   }
   return( NULL );
}
/*
***************************************************************
*/
vector<string> *G4GIDI::getNamesOfAvailableLibraries( int iZ, int iA, int iM ) {

    char *targetName = G4GIDI_Misc_Z_A_m_ToName( iZ, iA, iM );
    vector<string> *listOfLibraries;

    if( targetName == NULL ) return( new vector<string>( ) );
    string targetSymbol( targetName );
    listOfLibraries = getNamesOfAvailableLibraries( targetSymbol );
    smr_freeMemory( (void **) &targetName );
    return( listOfLibraries );
}
/*
***************************************************************
*/
vector<string> *G4GIDI::getNamesOfAvailableLibraries( string &targetName ) {

    list<G4GIDI_map *>::iterator iter;
    vector<string> *listOfLibraries = new vector<string>( );

    MCGIDI_map *map;
    MCGIDI_mapEntry *entry;

    for( iter = dataDirectories.begin( ); iter != dataDirectories.end( ); iter++ ) {
        map = MCGIDI_map_findAllOfTarget( &((*iter)->smr), (*iter)->map, projectile.c_str( ), targetName.c_str( ) );
        for( entry = MCGIDI_map_getFirstEntry( map ); entry != NULL; entry = MCGIDI_map_getNextEntry( entry ) ) {
            listOfLibraries->push_back( entry->evaluation );
        }
        MCGIDI_map_free( NULL, map );
    }
    return( listOfLibraries );
}
/*
***************************************************************
*/
vector<string> *G4GIDI::getNamesOfAvailableTargets( void ) {

    vector<string> *listOfTargets;
    list<G4GIDI_map *>::iterator iter_map;

    listOfTargets = new vector<string>( );
    if( listOfTargets == NULL ) return( NULL );
    for( iter_map = dataDirectories.begin( ); iter_map != dataDirectories.end( ); iter_map++ ) {
        if( MCGIDI_map_walkTree( NULL, (*iter_map)->map, getNamesOfAvailableTargets_walker, (void *) listOfTargets ) != 0 ) {
            delete listOfTargets;
            return( NULL );
        }
    }
    return( listOfTargets );
}
/*
***************************************************************
*/
G4GIDI_target *G4GIDI::readTarget( string &lib_name, int iZ, int iA, int iM, bool bind ) {

    char *targetName = G4GIDI_Misc_Z_A_m_ToName( iZ, iA, iM );
    G4GIDI_target *target;

    if( targetName == NULL ) return( NULL );
    string targetSymbol( targetName );
    target = readTarget( lib_name, targetSymbol, bind );
    smr_freeMemory( (void **) &targetName );
    return( target );
}
/*
***************************************************************
*/
G4GIDI_target *G4GIDI::readTarget( string &lib_name, string &targetName, bool bind ) {

    vector<G4GIDI_target *>::iterator iter_targets;

    for( iter_targets = targets.begin( ); iter_targets != targets.end( ); iter_targets++ ) {
        if( (*iter_targets)->name == targetName ) return( NULL );
    }
    char *path = dataFilename( lib_name, targetName );
    if( path == NULL ) return( NULL );

    G4GIDI_target *target = new G4GIDI_target( path );
    if( bind ) targets.push_back( target );
    smr_freeMemory( (void **) &path );
    return( target );
}
/*
***************************************************************
*/
G4GIDI_target *G4GIDI::getAlreadyReadTarget( int iZ, int iA, int iM ) {

    char *targetName = G4GIDI_Misc_Z_A_m_ToName( iZ, iA, iM );
    G4GIDI_target *target;

    if( targetName == NULL ) return( NULL );
    string targetSymbol( targetName );
    target = getAlreadyReadTarget( targetSymbol );
    smr_freeMemory( (void **) &targetName );
    return( target );
}
/*
***************************************************************
*/
G4GIDI_target *G4GIDI::getAlreadyReadTarget( string &targetSymbol ) {

    vector<G4GIDI_target *>::iterator iter_targets;

    for( iter_targets = targets.begin( ); iter_targets != targets.end( ); iter_targets++ ) {
        if( ( (*iter_targets)->name == targetSymbol ) ) return( *iter_targets );
    }
    return( NULL );
}
/*
***************************************************************
*/
int G4GIDI::freeTarget( G4GIDI_target *target ) {

    vector<G4GIDI_target *>::iterator iter_targets;

    for( iter_targets = targets.begin( ); iter_targets != targets.end( ); iter_targets++ ) {
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
int G4GIDI::freeTarget( int iZ, int iA, int iM ) {

    int status;
    char *targetName = G4GIDI_Misc_Z_A_m_ToName( iZ, iA, iM );

    if( targetName == NULL ) return( 1 );
    string targetSymbol( targetName );
    status = freeTarget( targetSymbol );
    smr_freeMemory( (void **) &targetName );
    return( status );
}
/*
***************************************************************
*/
int G4GIDI::freeTarget( string &targetSymbol ) {

    vector<G4GIDI_target *>::iterator iter_targets;

    for( iter_targets = targets.begin( ); iter_targets != targets.end( ); iter_targets++ ) {
        if( (*iter_targets)->name == targetSymbol ) return( freeTarget( *iter_targets ) );
    }
    return( 1 );
}
/*
***************************************************************
*/
vector<string> *G4GIDI::getListOfReadTargetsNames( void ) {

    vector<G4GIDI_target *>::iterator iter_targets;
    vector<string> *listOfTargets;

    listOfTargets = new vector<string>( );
    if( listOfTargets == NULL ) return( NULL );
    for( iter_targets = targets.begin( ); iter_targets != targets.end( ); iter_targets++ ) {
        listOfTargets->push_back( *(*iter_targets)->getName( ) );
    }
    return( listOfTargets );
}
