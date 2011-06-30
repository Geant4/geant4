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
#include <string>
#include <vector>
#include <xData.h>
#include <tpia_target.h>
#include <tpia_misc.h>
#include <string.h>
#include "G4GIDI_Misc.hh"

using namespace std;
using namespace GIDI;

/*
***************************************************************
*/
char *G4GIDI_Misc_Z_A_m_ToName( int iZ, int iA, int im ) {

    const char *Z = tpia_misc_ZToSymbol( iZ );
    char S[128], mS[32], *name;

    if( Z == NULL ) return( NULL );
    if( iA == 0 ) {
        if( im != 0 ) return( NULL );
        sprintf( S, "%s_natural", Z ); }
    else {
        sprintf( S, "%s_%d", Z, iA );
        if( im != 0 ) {
            sprintf( mS, "_m%d", im );
            strcat( S, mS );
        }
    }
    name = (char *) xData_malloc2( NULL, strlen( S ) + 1, 0, "name" );
    if( name != NULL ) strcpy( name, S );
    return( name );
}
/*
***************************************************************
*/
char *G4GIDI_Misc_channelCompound( char *particle1, char *particle2 ) {

    int Z1, A1, m1, Z2, A2, m2;

    if( tpia_miscNameToZAm( NULL, particle1, &Z1, &A1, &m1 ) ) return( NULL );
    if( tpia_miscNameToZAm( NULL, particle2, &Z2, &A2, &m2 ) ) return( NULL );
    if( A1 == 0 ) A2 = 0;
    if( A2 == 0 ) A1 = 0;
    return( G4GIDI_Misc_Z_A_m_ToName( Z1 + Z2, A1 + A2, 0 ) );
}
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
/*
***************************************************************
*/
//int getNamesOfAvailableTargets_walker( tpia_mapEntry *entry, int level, void *userData ) {
int getNamesOfAvailableTargets_walker( tpia_mapEntry *entry, int , void *userData ) {

    vector<string> *listOfTargets = (vector<string> *) userData;
    vector<string>::iterator iter;

    if( entry->type != tpia_mapEntry_type_target ) return( 0 );
    for( iter = listOfTargets->begin( ); iter != listOfTargets->end( ); iter++ ) {
        if( entry->targetName == iter->c_str( ) ) return( 0 );
    }
    listOfTargets->push_back( entry->targetName );
    return( 0 );
}
