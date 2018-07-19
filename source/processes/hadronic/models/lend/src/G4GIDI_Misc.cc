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
#include <string.h>
#include <iostream>
#include <string>
#include <vector>
#include <statusMessageReporting.h>
#include <MCGIDI.h>
#include <MCGIDI_misc.h>
#include "G4GIDI_Misc.hh"
using namespace std;
using namespace GIDI;

/*
***************************************************************
*/
char *G4GIDI_Misc_Z_A_m_ToName( int iZ, int iA, int im ) {

    const char *Z = MCGIDI_misc_ZToSymbol( iZ );
    char S[128], mS[32], *name;

    if( Z == NULL ) return( NULL );
    if( iA == 0 ) {
        if( im != 0 ) return( NULL );
        sprintf( S, "%s_natural", Z ); }
    else {
        sprintf( S, "%s%d", Z, iA );
        if( im != 0 ) {
            //sprintf( mS, "_m%d", im );
            //TK 170509
            //Fix inconsistency of name of excited isomer  between data and code
            sprintf( mS, "m%d", im );
            strcat( S, mS );
        }
    }
    name = (char *) smr_malloc2( NULL, strlen( S ) + 1, 0, "name" );
    if( name != NULL ) strcpy( name, S );
    return( name );
}
/*
***************************************************************
*/
char *G4GIDI_Misc_channelCompound( char *particle1, char *particle2 ) {

    int Z1, A1, m1, Z2, A2, m2, level1, level2;

    if( MCGIDI_miscNameToZAm( NULL, particle1, &Z1, &A1, &m1, &level1 ) ) return( NULL );
    if( MCGIDI_miscNameToZAm( NULL, particle2, &Z2, &A2, &m2, &level2 ) ) return( NULL );
    if( A1 == 0 ) A2 = 0;
    if( A2 == 0 ) A1 = 0;
    return( G4GIDI_Misc_Z_A_m_ToName( Z1 + Z2, A1 + A2, 0 ) );
}
#if 0
/*
***************************************************************
*/
int G4GIDI_Misc_channelProductsCompare( tpia_channel *channel, int nProducts, char **productNames ) {

    int i;
    tpia_product *product;

    if( channel->decayChannel.numberOfProducts != nProducts ) return( 0 );
    for( product = tpia_channel_getFirstProduct( channel ), i = 0; product != NULL; product = tpia_decayChannel_getNextProduct( product ), i++ ) {
        if( strcmp( product->productID->name, productNames[i] ) ) return( 0 );
    }
    return( 1 );
}
#endif
/*
***************************************************************
*/
int getNamesOfAvailableTargets_walker( MCGIDI_mapEntry *entry, int /*level*/, void *userData ) {

    vector<string> *listOfTargets = (vector<string> *) userData;
    vector<string>::iterator iter;

    if( entry->type != MCGIDI_mapEntry_type_target ) return( 0 );
    for( iter = listOfTargets->begin( ); iter != listOfTargets->end( ); iter++ ) {
        if( entry->targetName == iter->c_str( ) ) return( 0 );
    }
    listOfTargets->push_back( entry->targetName );
    return( 0 );
}
