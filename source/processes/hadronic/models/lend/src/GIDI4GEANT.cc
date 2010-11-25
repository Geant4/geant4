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

#include "GIDI4GEANT.hpp"

/*
***************************************************************
*/
GIDI4GEANT::GIDI4GEANT( int ip, string &dataDirectory ) {

    init( ip );
    addDataDirectory( dataDirectory );
}
/*
***************************************************************
*/
GIDI4GEANT::GIDI4GEANT( int ip, list<string> &dataDirectoryList ) {

    list<string>::iterator iter;

    init( ip );
    for( iter = dataDirectoryList.begin( ); iter != dataDirectoryList.end( ); ++iter ) addDataDirectory( *iter );
}
/*
***************************************************************
*/
GIDI4GEANT::~GIDI4GEANT( void ) {

    GIDI4GEANT_target *target;
    list<GIDI4GEANT_map *>::iterator iter;

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
int GIDI4GEANT::init( int ip ) {

    projectileID = ip;
    return( 0 );
}
/*
***************************************************************
*/
int GIDI4GEANT::numberOfDataDirectories( void ) {

    return( dataDirectories.size( ) );
}
/*
***************************************************************
*/
int GIDI4GEANT::addDataDirectory( string &dataDirectory ) {

    list<GIDI4GEANT_map *>::iterator iter;

    for( iter = dataDirectories.begin( ); iter != dataDirectories.end( ); iter++ ) {
        if( (*iter)->path( ) == dataDirectory ) return( 0 );
    }

    GIDI4GEANT_map *map = new GIDI4GEANT_map( dataDirectory );
    dataDirectories.push_back( map );

    return( 0 );
}
/*
***************************************************************
*/
int GIDI4GEANT::removeDataDirectory( string &dataDirectory ) {

    list<GIDI4GEANT_map *>::iterator iter;

    for( iter = dataDirectories.begin( ); iter != dataDirectories.end( ); iter++ ) {
        if( dataDirectory == (*iter)->path( ) ) {
            
        }
    }
    return( 0 );
}
/*
***************************************************************
*/
//string GIDI4GEANT::getDataDirectoryAtIndex( int index ) {
string GIDI4GEANT::getDataDirectoryAtIndex( int ) {

#if 0
    list<GIDI4GEANT_map *>::iterator iter;
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
vector<string> *GIDI4GEANT::getDataDirectories( void ) {

    int i = 0;
    list<GIDI4GEANT_map *>::iterator iter;
    vector<string> *v = new vector<string>( numberOfDataDirectories( ) );

    for( iter = dataDirectories.begin( ); iter != dataDirectories.end( ); iter++, i++ ) (*v)[i] = string( (*iter)->fileName( ) );
    return( v );
}
/*
***************************************************************
*/
bool GIDI4GEANT::isThisDataAvailable( string &lib_name, int iZ, int iA, int iM ) {

    bool b;
    char *targetName = GIDI4GEANT_Misc_Z_A_m_ToName( iZ, iA, iM );

    if( targetName == NULL ) return( false );
    string targetSymbol( targetName );
    b = isThisDataAvailable( lib_name, targetSymbol );
    xData_free( NULL, targetName );
    return( b );
}
/*
***************************************************************
*/
bool GIDI4GEANT::isThisDataAvailable( string &lib_name, string &targetName ) {

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
char *GIDI4GEANT::dataFilename( string &lib_name, int iZ, int iA, int iM ) {

    char *targetName = GIDI4GEANT_Misc_Z_A_m_ToName( iZ, iA, iM ), *fileName;

    if( targetName == NULL ) return( NULL );
    string targetSymbol( targetName );
    fileName = dataFilename( lib_name, targetSymbol );
    xData_free( NULL, targetName );
    return( fileName );
}
/*
***************************************************************
*/
char *GIDI4GEANT::dataFilename( string &lib_name, string &targetSymbol ) {

    //char *path, *projectile = "n_1";
    char *path, *projectile = (char*)"n_1";
    list<GIDI4GEANT_map *>::iterator iter;

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
vector<string> *GIDI4GEANT::getNamesOfAvailableLibraries( int iZ, int iA, int iM ) {

    char *targetName = GIDI4GEANT_Misc_Z_A_m_ToName( iZ, iA, iM );
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
vector<string> *GIDI4GEANT::getNamesOfAvailableLibraries( string &targetName ) {

    //char *projectile = "n_1";
    char *projectile = (char*)"n_1";
    list<GIDI4GEANT_map *>::iterator iter;
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
vector<string> *GIDI4GEANT::getNamesOfAvailableTargets( void ) {

    vector<string> *listOfTargets;
    list<GIDI4GEANT_map *>::iterator iter_map;

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
GIDI4GEANT_target *GIDI4GEANT::readTarget( string &lib_name, int iZ, int iA, int iM, bool bind ) {

    char *targetName = GIDI4GEANT_Misc_Z_A_m_ToName( iZ, iA, iM );
    GIDI4GEANT_target *target;

    if( targetName == NULL ) return( NULL );
    string targetSymbol( targetName );
    target = readTarget( lib_name, targetSymbol, bind );
    xData_free( NULL, targetName );
    return( target );
}
/*
***************************************************************
*/
GIDI4GEANT_target *GIDI4GEANT::readTarget( string &lib_name, string &targetName, bool bind ) {

    vector<GIDI4GEANT_target *>::iterator iter_targets;

    for( iter_targets = targets.begin( ); iter_targets != targets.end( ); iter_targets++ ) {
        if( (*iter_targets)->name == targetName ) return( NULL );
    }
    char *path = dataFilename( lib_name, targetName );
    if( path == NULL ) return( NULL );

    GIDI4GEANT_target *target = new GIDI4GEANT_target( path );
    if( bind ) targets.push_back( target );
    xData_free( NULL, path );
    return( target );
}
/*
***************************************************************
*/
GIDI4GEANT_target *GIDI4GEANT::getAlreadyReadTarget( int iZ, int iA, int iM ) {

    char *targetName = GIDI4GEANT_Misc_Z_A_m_ToName( iZ, iA, iM );
    GIDI4GEANT_target *target;

    if( targetName == NULL ) return( NULL );
    string targetSymbol( targetName );
    target = getAlreadyReadTarget( targetSymbol );
    xData_free( NULL, targetName );
    return( target );
}
/*
***************************************************************
*/
GIDI4GEANT_target *GIDI4GEANT::getAlreadyReadTarget( string &targetSymbol ) {

    vector<GIDI4GEANT_target *>::iterator iter_targets;

    for( iter_targets = targets.begin( ); iter_targets != targets.end( ); iter_targets++ ) {
        if( ( (*iter_targets)->name == targetSymbol ) ) return( *iter_targets );
    }
    return( NULL );
}
/*
***************************************************************
*/
int GIDI4GEANT::freeTarget( GIDI4GEANT_target *target ) {

    vector<GIDI4GEANT_target *>::iterator iter_targets;

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
int GIDI4GEANT::freeTarget( int iZ, int iA, int iM ) {

    int status;
    char *targetName = GIDI4GEANT_Misc_Z_A_m_ToName( iZ, iA, iM );

    if( targetName == NULL ) return( 1 );
    string targetSymbol( targetName );
    status = freeTarget( targetSymbol );
    xData_free( NULL, targetName );
    return( status );
}
/*
***************************************************************
*/
int GIDI4GEANT::freeTarget( string &targetSymbol ) {

    vector<GIDI4GEANT_target *>::iterator iter_targets;

    for( iter_targets = targets.begin( ); iter_targets != targets.end( ); iter_targets++ ) {
        if( (*iter_targets)->name == targetSymbol ) return( freeTarget( *iter_targets ) );
    }
    return( 1 );
}
/*
***************************************************************
*/
vector<string> *GIDI4GEANT::getListOfReadTargetsNames( void ) {

    vector<GIDI4GEANT_target *>::iterator iter_targets;
    vector<string> *listOfTargets;

    listOfTargets = new vector<string>( );
    if( listOfTargets == NULL ) return( NULL );
    for( iter_targets = targets.begin( ); iter_targets != targets.end( ); iter_targets++ ) {
        listOfTargets->push_back( *(*iter_targets)->getName( ) );
    }
    return( listOfTargets );
}
