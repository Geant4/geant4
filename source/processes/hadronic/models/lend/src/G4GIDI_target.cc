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
#include <stdlib.h>

#include "G4GIDI_target.hh"
#include "G4GIDI_mass.hh"
#include "G4GIDI_Misc.hh"

using namespace std;
using namespace GIDI;

/*
***************************************************************
*/
G4GIDI_target::G4GIDI_target( const char *fileName ) {

    init( fileName );
}
/*
***************************************************************
*/
G4GIDI_target::G4GIDI_target( string &fileName ) {

    init( fileName.c_str( ) );
}
/*
***************************************************************
*/
void G4GIDI_target::init( const char *fileName ) {

    int i, j, n, *p, *q;
    tpia_channel *channel;

    smr_initialize( &smr );
    sourceFilename = fileName;
    target = tpia_target_createRead( &smr, fileName );
    if( !smr_isOk( &smr ) ) {
        smr_print( &smr, stderr, 1 );
        throw 1;
    }
    name = target->targetID->name;
    mass = G4GIDI_targetMass( target->targetID->name );
    equalProbableBinSampleMethod = "constant";
    elasticIndices = NULL;
    nElasticIndices = nCaptureIndices = nFissionIndices = nOthersIndices = 0;

    if( ( n = tpia_target_numberOfChannels( &smr, target ) ) > 0 ) {
        p = elasticIndices = (int *) xData_malloc2( NULL, n * sizeof( double ), 1, "elasticIndices" );
        for( i = 0; i < n; i++ ) {      /* Find elastic channel(s). */
            channel = tpia_target_heated_getChannelAtIndex( target->baseHeatedTarget, i );
            if( channel->ENDL_C == 10 ) {
                *(p++) = i;
                nElasticIndices++;
            }
        }
        captureIndices = p;
        for( i = 0; i < n; i++ ) {      /* Find capture channel(s). */
            channel = tpia_target_heated_getChannelAtIndex( target->baseHeatedTarget, i );
            if( channel->ENDL_C == 46 ) {
                *(p++) = i;
                nCaptureIndices++;
            }
        }

        fissionIndices = p;
        for( i = 0; i < n; i++ ) {      /* Find fission channel(s). */
            channel = tpia_target_heated_getChannelAtIndex( target->baseHeatedTarget, i );
            if( channel->fission != NULL ) {
                *(p++) = i;
                nFissionIndices++;
            }
        }
        othersIndices = p;
        for( i = 0; i < n; i++ ) {      /* Find other channel(s). */
            for( j = 0, q = elasticIndices; j < nElasticIndices; j++, q++ ) if( *q == i ) break;
            if( j < nElasticIndices ) continue;
            for( j = 0, q = captureIndices; j < nCaptureIndices; j++, q++ ) if( *q == i ) break;
            if( j < nCaptureIndices ) continue;
            for( j = 0, q = fissionIndices; j < nFissionIndices; j++, q++ ) if( *q == i ) break;
            if( j < nFissionIndices ) continue;
            *p = i;
            p++;
            nOthersIndices++;
        }
    }
}
/*
***************************************************************
*/
G4GIDI_target::~G4GIDI_target( ) {

    tpia_target_free( &smr, target );
    xData_free( &smr, elasticIndices );
    smr_release( &smr );
}
/*
***************************************************************
*/
string *G4GIDI_target::getName( void ) { return( &name ); }
/*
***************************************************************
*/
string *G4GIDI_target::getFilename( void ) { return( &sourceFilename ); }
/*
***************************************************************
*/
int G4GIDI_target::getZ( void ) {
   
    return( target->targetID->Z );
}
/*
***************************************************************
*/
int G4GIDI_target::getA( void ) {
   
    return( target->targetID->A );
}
/*
***************************************************************
*/
int G4GIDI_target::getM( void ) {
   
    return( target->targetID->m );
}
/*
***************************************************************
*/
double G4GIDI_target::getMass( void ) {

    return( mass );
}
/*
***************************************************************
*/
int G4GIDI_target::getTemperatures( double *temperatures ) {

    return( tpia_target_getTemperatures( &smr, target, temperatures ) );
}
/*
***************************************************************
*/
int G4GIDI_target::readTemperature( int index ) {

    return( tpia_target_readHeatedTarget( &smr, target, index, 0 ) );
}
/*
***************************************************************
*/
string G4GIDI_target::getEqualProbableBinSampleMethod( void ) {

    return( equalProbableBinSampleMethod );
}
/*
***************************************************************
*/
int G4GIDI_target::setEqualProbableBinSampleMethod( string method ) {

    if( method == "constant" ) {
        equalProbableBinSampleMethod = "constant"; }
    if( method == "linear" ) {
        equalProbableBinSampleMethod = "linear"; }
    else {
        return( 1 );
    }
    return( 0 );
}
/*
***************************************************************
*/
int G4GIDI_target::getNumberOfChannels( void ) {

    return( tpia_target_numberOfChannels( &smr, target ) );
}
/*
***************************************************************
*/
int G4GIDI_target::getNumberOfProductionChannels( void ) {

    return( tpia_target_numberOfProductionChannels( &smr, target ) );
}
/*
***************************************************************
*/
vector<channelID> *G4GIDI_target::getChannelIDs( void ) { 

    return( getChannelIDs2( target->baseHeatedTarget->channels, tpia_target_numberOfChannels( &smr, target ) ) );
}
/*
***************************************************************
*/
vector<channelID> *G4GIDI_target::getProductionChannelIDs( void ) {

    return( getChannelIDs2( target->baseHeatedTarget->productionChannels, tpia_target_numberOfProductionChannels( &smr, target ) ) );
}
/*
***************************************************************
*/
vector<channelID> *G4GIDI_target::getChannelIDs2( tpia_channel **channels, int n ) {

    int i = 0;
    vector<channelID> *listOfChannels;

    listOfChannels = new vector<channelID>( n );
    for( i = 0; i < n; i++ ) (*listOfChannels)[i].ID = channels[i]->outputChannel;
    return( listOfChannels );
}
/*
***************************************************************
*/
vector<double> *G4GIDI_target::getEnergyGridAtTIndex( int index ) {

    xData_Int i, n;
    double *dEnergyGrid;
    vector<double> *energyGrid;
    vector<double>::iterator iter;

    n = tpia_target_getEnergyGridAtTIndex( &smr, target, index, &dEnergyGrid );
    if( n < 0 ) return( NULL );
    energyGrid = new vector<double>( n );
    for( i = 0, iter = energyGrid->begin( ); i < n; i++, iter++ ) *iter = dEnergyGrid[i];
    return( energyGrid );
}
/*
***************************************************************
*/
double G4GIDI_target::getTotalCrossSectionAtE( double e_in, double temperature ) {

    return( tpia_target_getTotalCrossSectionAtTAndE( NULL, target, temperature, -1, e_in, tpia_crossSectionType_pointwise ) );
}
/*
***************************************************************
*/
double G4GIDI_target::getElasticCrossSectionAtE( double e_in, double temperature ) {

    return( sumChannelCrossSectionAtE( nElasticIndices, elasticIndices, e_in, temperature ) );
}
/*
***************************************************************
*/
double G4GIDI_target::getCaptureCrossSectionAtE( double e_in, double temperature ) {

    return( sumChannelCrossSectionAtE( nCaptureIndices, captureIndices, e_in, temperature ) );
}
/*
***************************************************************
*/
double G4GIDI_target::getFissionCrossSectionAtE( double e_in, double temperature ) {

    return( sumChannelCrossSectionAtE( nFissionIndices, fissionIndices, e_in, temperature ) );
}
/*
***************************************************************
*/
double G4GIDI_target::getOthersCrossSectionAtE( double e_in, double temperature ) {

    return( sumChannelCrossSectionAtE( nOthersIndices, othersIndices, e_in, temperature ) );
}
/*
***************************************************************
*/
double G4GIDI_target::sumChannelCrossSectionAtE( int nIndices, int *indices, double e_in, double temperature ) {

    int i;
    double xsec = 0.;

    for( i = 0; i < nIndices; i++ ) 
        xsec += tpia_target_getIndexChannelCrossSectionAtE( &smr, target, indices[i], temperature, -1, e_in, tpia_crossSectionType_pointwise );
    return( xsec );
}
/*
***************************************************************
*/
int G4GIDI_target::sampleChannelCrossSectionAtE( int nIndices, int *indices, double e_in, double temperature,
        double (*rng)( void * ), void *rngState ) {

    int i;
    double xsec = 0., rxsec = sumChannelCrossSectionAtE( nIndices, indices, e_in, temperature ) * tpia_misc_drng( rng, rngState );

    for( i = 0; i < nIndices - 1; i++ ) {
        xsec += tpia_target_getIndexChannelCrossSectionAtE( &smr, target, indices[i], temperature, -1, e_in, tpia_crossSectionType_pointwise );
        if( xsec >= rxsec ) break;
    }
    return( indices[i] );
}
/*
***************************************************************
*/
//double G4GIDI_target::getElasticFinalState( double e_in, double temperature, double (*rng)( void * ), void *rngState ) {
double G4GIDI_target::getElasticFinalState( double e_in, double , double (*rng)( void * ), void *rngState ) {

    tpia_decaySamplingInfo decaySamplingInfo;
    tpia_channel *channel = tpia_target_heated_getChannelAtIndex_smr( &smr, target->baseHeatedTarget, elasticIndices[0] );
    tpia_product *product;

    decaySamplingInfo.e_in = e_in;
    decaySamplingInfo.isVelocity = 0;
    tpia_frame_setColumn( &smr, &(decaySamplingInfo.frame), 0, tpia_referenceFrame_lab );
    decaySamplingInfo.samplingMethods = &(target->samplingMethods);
    decaySamplingInfo.rng = rng;
    decaySamplingInfo.rngState = rngState;
    product = tpia_decayChannel_getFirstProduct( &(channel->decayChannel) );
    tpia_angular_SampleMu( &smr, &(product->angular), &decaySamplingInfo );

    return( decaySamplingInfo.mu );
}
/*
***************************************************************
*/
vector<G4GIDI_Product> *G4GIDI_target::getCaptureFinalState( double e_in, double temperature, double (*rng)( void * ), void *rngState ) {

    return( getFinalState( nCaptureIndices, captureIndices, e_in, temperature, rng, rngState ) );
}
/*
***************************************************************
*/
vector<G4GIDI_Product> *G4GIDI_target::getFissionFinalState( double e_in, double temperature, double (*rng)( void * ), void *rngState ) {

    return( getFinalState( nFissionIndices, fissionIndices, e_in, temperature, rng, rngState ) );
}
/*
***************************************************************
*/
vector<G4GIDI_Product> *G4GIDI_target::getOthersFinalState( double e_in, double temperature, double (*rng)( void * ), void *rngState ) {

    return( getFinalState( nOthersIndices, othersIndices, e_in, temperature, rng, rngState ) );
}
/*
***************************************************************
*/
vector<G4GIDI_Product> *G4GIDI_target::getFinalState( int nIndices, int *indices, double e_in, double temperature, double (*rng)( void * ), void *rngState ) {

#define nProductsMax 50
    int index = 0, i, n;
    vector<G4GIDI_Product> *products = NULL;
    tpia_decaySamplingInfo decaySamplingInfo;
    tpia_productOutgoingData productDatas[nProductsMax], *productData;

    decaySamplingInfo.e_in = e_in;
    decaySamplingInfo.samplingMethods = &(target->samplingMethods);
    decaySamplingInfo.isVelocity = 0;
    tpia_frame_setColumn( &smr, &(decaySamplingInfo.frame), 0, tpia_referenceFrame_lab );
    decaySamplingInfo.rng = rng;
    decaySamplingInfo.rngState = rngState;

    if( nIndices == 0 ) {
        return( NULL ); }
    else {
        if( nIndices == 1 ) {
            index = indices[0]; }
        else {
            index = sampleChannelCrossSectionAtE( nIndices, indices, e_in, temperature, rng, rngState );
        }
    }
    n = tpia_target_heated_sampleIndexChannelProductsAtE( &smr, target->baseHeatedTarget, index, &decaySamplingInfo, nProductsMax, productDatas );
    if( n > 0 ) {
        if( ( products = new vector<G4GIDI_Product>( n ) ) != NULL ) {
            for( i = 0; i < n; i++ ) {
                productData = &(productDatas[i]);
                (*products)[i].A = productData->productID->A;
                (*products)[i].Z = productData->productID->Z;
                (*products)[i].m = productData->productID->m;
                (*products)[i].kineticEnergy = productData->kineticEnergy;
                (*products)[i].px = productData->px_vx;
                (*products)[i].py = productData->py_vy;
                (*products)[i].pz = productData->pz_vz;
            }
        }
    }

    return( products );
#undef nProductsMax
}
