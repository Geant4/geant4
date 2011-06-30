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
# Copyright (c) 2010, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory 
# Written by Bret R. Beck, beck6@llnl.gov. 
# CODE-461393
# All rights reserved. 
#  
# This file is part of GIDI. For details, see nuclear.llnl.gov. 
# Please also read the "Additional BSD Notice" at nuclear.llnl.gov. 
# 
# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met: 
#
#      1) Redistributions of source code must retain the above copyright notice, 
#         this list of conditions and the disclaimer below.
#      2) Redistributions in binary form must reproduce the above copyright notice, 
#         this list of conditions and the disclaimer (as noted below) in the 
#          documentation and/or other materials provided with the distribution.
#      3) Neither the name of the LLNS/LLNL nor the names of its contributors may be 
#         used to endorse or promote products derived from this software without 
#         specific prior written permission. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT 
# SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS 
# OR SERVICES;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
# AND ON  ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
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
    }

    while( ( iter = dataDirectories.begin( ) ) != dataDirectories.end( ) ) {
        delete *iter;
        dataDirectories.pop_front( );
    }
}
/*
***************************************************************
*/
int G4GIDI::init( int ip ) {

    projectileID = ip;
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
//string G4GIDI::getDataDirectoryAtIndex( int index ) {
string G4GIDI::getDataDirectoryAtIndex( int ) {

#if 0
    list<G4GIDI_map *>::iterator iter;
    unsigned i = (unsigned) index;

    if( i < 0 ) return( "" );
    if( i >= dataDirectories.size( ) ) return( "" );
    for( iter = dataDirectories.begin( ); iter != dataDirectories.end( ); iter++, i-- ) {
        if( i == 0 ) return( (*iter)->fileName( ) );
    }
#endif
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
    xData_free( NULL, targetName );
    return( b );
}
/*
***************************************************************
*/
bool G4GIDI::isThisDataAvailable( string &lib_name, string &targetName ) {

    char *path = dataFilename( lib_name, targetName );

    if( path != NULL ) {
        xData_free( NULL, path );
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
    xData_free( NULL, targetName );
    return( fileName );
}
/*
***************************************************************
*/
char *G4GIDI::dataFilename( string &lib_name, string &targetSymbol ) {

    //char *path, *projectile = "n_1";
    char *path, *projectile = (char*)"n_1";
    list<G4GIDI_map *>::iterator iter;

    for( iter = dataDirectories.begin( ); iter != dataDirectories.end( ); iter++ ) {
        if( ( path = tpia_map_findTarget( &((*iter)->smr), (*iter)->map, lib_name.c_str( ), projectile, targetSymbol.c_str( ) ) ) != NULL ) {
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
    xData_free( NULL, targetName );
    return( listOfLibraries );
}
/*
***************************************************************
*/
vector<string> *G4GIDI::getNamesOfAvailableLibraries( string &targetName ) {

    //char *projectile = "n_1";
    char *projectile = (char*)"n_1";
    list<G4GIDI_map *>::iterator iter;
    vector<string> *listOfLibraries = new vector<string>( );

    tpia_map *map;
    tpia_mapEntry *entry;
    for( iter = dataDirectories.begin( ); iter != dataDirectories.end( ); iter++ ) {
        map = tpia_map_findAllOfTarget( &((*iter)->smr), (*iter)->map, projectile, targetName.c_str( ) );
        for( entry = tpia_map_getFirstEntry( map ); entry != NULL; entry = tpia_map_getNextEntry( entry ) ) {
            listOfLibraries->push_back( entry->evaluation );
        }
        tpia_map_free( NULL, map );
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
        if( tpia_map_walkTree( NULL, (*iter_map)->map, getNamesOfAvailableTargets_walker, (void *) listOfTargets ) != 0 ) {
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
    xData_free( NULL, targetName );
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
    xData_free( NULL, path );
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
    xData_free( NULL, targetName );
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
    xData_free( NULL, targetName );
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
