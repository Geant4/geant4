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
#include "G4GIDI_map.hh"

using namespace GIDI;

/*
***************************************************************
*/
G4GIDI_map::G4GIDI_map( std::string &dataDirectory ) {

    smr_initialize( &smr, smr_status_Ok, 0 );
    map = MCGIDI_map_readFile( &smr, NULL, dataDirectory.c_str( ) );
    if( !smr_isOk( &smr ) ) {
        smr_print( &smr, 1 );
        throw 1;
    }
}
/*
***************************************************************
*/
G4GIDI_map::~G4GIDI_map( void ) {

    if( map != NULL ) MCGIDI_map_free( NULL, map );
    smr_release( &smr );
}
/*
***************************************************************
*/
std::string G4GIDI_map::fileName( void ) {

    return( map->mapFileName );
}
/*
***************************************************************
*/
std::string G4GIDI_map::path( void ) {

    return( map->path );
}
