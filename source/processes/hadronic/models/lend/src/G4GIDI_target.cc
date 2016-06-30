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
#include <stdlib.h>
#include <cmath>

#include "G4GIDI_target.hh"
#include "G4GIDI_mass.hh"
#include "G4GIDI_Misc.hh"

using namespace std;
using namespace GIDI;

/*
***************************************************************
*/
G4GIDI_target::G4GIDI_target( char const *fileName ) {

    init( fileName );
}
/*
***************************************************************
*/
G4GIDI_target::G4GIDI_target( string const &fileName ) {

    init( fileName.c_str( ) );
}
/*
***************************************************************
*/
void G4GIDI_target::init( char const *fileName ) {

    int i, j, n, *p, *q, ir;
    MCGIDI_reaction *reaction;

    smr_initialize( &smr, smr_status_Ok, 1 );
    sourceFilename = fileName;
    target = MCGIDI_target_newRead( &smr, fileName );
    if( !smr_isOk( &smr ) ) {
        smr_print( &smr, 1 );
        throw 1;
    }
    projectilesPOPID = target->projectilePOP->globalPoPsIndex;
    name = target->targetPOP->name;
    mass = G4GIDI_targetMass( target->targetPOP->name );
    equalProbableBinSampleMethod = "constant";
    elasticIndices = NULL;
    nElasticIndices = nCaptureIndices = nFissionIndices = nOthersIndices = 0;

    if( ( n = MCGIDI_target_numberOfReactions( &smr, target ) ) > 0 ) {
        if( ( p = elasticIndices = (int *) smr_malloc2( &smr, n * sizeof( double ), 1, "elasticIndices" ) ) == NULL ) {
            smr_print( &smr, 1 );
            throw 1;
        }
        for( i = 0; i < n; i++ ) {      /* Find elastic channel(s). */
            reaction = MCGIDI_target_heated_getReactionAtIndex( target->baseHeatedTarget, i );
            if( MCGIDI_reaction_getENDF_MTNumber( reaction ) == 2 ) {
                *(p++) = i;
                nElasticIndices++;
            }
        }
        captureIndices = p;
        for( i = 0; i < n; i++ ) {      /* Find capture channel(s). */
            reaction = MCGIDI_target_heated_getReactionAtIndex( target->baseHeatedTarget, i );
            if( MCGIDI_reaction_getENDF_MTNumber( reaction ) == 102 ) {
                *(p++) = i;
                nCaptureIndices++;
            }
        }

        fissionIndices = p;
        for( i = 0; i < n; i++ ) {      /* Find fission channel(s). */
            reaction = MCGIDI_target_heated_getReactionAtIndex( target->baseHeatedTarget, i );
            ir = MCGIDI_reaction_getENDF_MTNumber( reaction );
            if( ( ir != 18 ) && ( ir != 19 ) && ( ir != 20 ) && ( ir != 21 ) && ( ir != 38 ) ) continue;
            *(p++) = i;
            nFissionIndices++;
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
#if 0
printf( "elastic %d: ", nElasticIndices );
for( i = 0; i < nElasticIndices; i++ ) printf( " %d", elasticIndices[i] );
printf( "\ncapture %d: ", nCaptureIndices );
for( i = 0; i < nCaptureIndices; i++ ) printf( " %d", captureIndices[i] );
printf( "\nfission %d: ", nFissionIndices );
for( i = 0; i < nFissionIndices; i++ ) printf( " %d", fissionIndices[i] );
printf( "\nothers %d: ", nOthersIndices );
for( i = 0; i < nOthersIndices; i++ ) printf( " %d", othersIndices[i] );
printf( "\n" );
#endif
    }
}
/*
***************************************************************
*/
G4GIDI_target::~G4GIDI_target( ) {

    MCGIDI_target_free( &smr, target );
    smr_freeMemory( (void **) &elasticIndices );
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
   
    return( target->targetPOP->Z );
}
/*
***************************************************************
*/
int G4GIDI_target::getA( void ) {
   
    return( target->targetPOP->A );
}
/*
***************************************************************
*/
int G4GIDI_target::getM( void ) {
   
    return( target->targetPOP->m );
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

    return( MCGIDI_target_getTemperatures( &smr, target, temperatures ) );
}
/*
***************************************************************
*/
int G4GIDI_target::readTemperature( int index ) {

    return( MCGIDI_target_readHeatedTarget( &smr, target, index ) );
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

    return( MCGIDI_target_numberOfReactions( &smr, target ) );
}
/*
***************************************************************
*/
int G4GIDI_target::getNumberOfProductionChannels( void ) {

    return( MCGIDI_target_numberOfProductionReactions( &smr, target ) );
}
/*
***************************************************************
*/
channelID G4GIDI_target::getChannelsID( int channelIndex ) { 

    MCGIDI_reaction *reaction;

    if( ( reaction = MCGIDI_target_heated_getReactionAtIndex_smr( &smr, target->baseHeatedTarget, channelIndex ) ) == NULL ) {
        smr_print( &smr, 1 );
        throw 1;
    }
    return( string( reaction->outputChannelStr ) );    /* Only works because channelID is defined to be string. */
}
/*
***************************************************************
*/
vector<channelID> *G4GIDI_target::getChannelIDs( void ) { 

    int i, n = MCGIDI_target_numberOfReactions( &smr, target );
    MCGIDI_reaction *reaction;
    vector<channelID> *listOfChannels;

    listOfChannels = new vector<channelID>( n );
    for( i = 0; i < n; i++ ) {
        reaction = MCGIDI_target_heated_getReactionAtIndex( target->baseHeatedTarget, i );
        (*listOfChannels)[i] = reaction->outputChannelStr;
    }
    return( listOfChannels );
}
/*
***************************************************************
*/
vector<channelID> *G4GIDI_target::getProductionChannelIDs( void ) {

    return( NULL );
}
/*
***************************************************************
*/
double G4GIDI_target::getTotalCrossSectionAtE( double e_in, double temperature ) {

    MCGIDI_quantitiesLookupModes mode( projectilesPOPID );

    mode.setProjectileEnergy( e_in );
    mode.setCrossSectionMode( MCGIDI_quantityLookupMode_pointwise );
    mode.setTemperature( temperature );

    return( MCGIDI_target_getTotalCrossSectionAtTAndE( NULL, target, mode, true ) );
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
    MCGIDI_quantitiesLookupModes mode( projectilesPOPID );

    mode.setProjectileEnergy( e_in );
    mode.setCrossSectionMode( MCGIDI_quantityLookupMode_pointwise );
    mode.setTemperature( temperature );

    for( i = 0; i < nIndices; i++ ) 
        xsec += MCGIDI_target_getIndexReactionCrossSectionAtE( &smr, target, indices[i], mode, true );
    return( xsec );
}
/*
***************************************************************
*/
int G4GIDI_target::sampleChannelCrossSectionAtE( int nIndices, int *indices, double e_in, double temperature,
        double (*rng)( void * ), void *rngState ) {

    int i;
    double xsec = 0., rxsec = sumChannelCrossSectionAtE( nIndices, indices, e_in, temperature ) * rng( rngState );
    MCGIDI_quantitiesLookupModes mode( projectilesPOPID );

    mode.setProjectileEnergy( e_in );
    mode.setCrossSectionMode( MCGIDI_quantityLookupMode_pointwise );
    mode.setTemperature( temperature );

    for( i = 0; i < nIndices - 1; i++ ) {
        xsec += MCGIDI_target_getIndexReactionCrossSectionAtE( &smr, target, indices[i], mode, true );
        if( xsec >= rxsec ) break;
    }
    return( indices[i] );
}
/*
***************************************************************
*/
double G4GIDI_target::getElasticFinalState( double e_in, double temperature, double (*rng)( void * ), void *rngState ) {

    MCGIDI_decaySamplingInfo decaySamplingInfo;
    MCGIDI_reaction *reaction = MCGIDI_target_heated_getReactionAtIndex_smr( &smr, target->baseHeatedTarget, elasticIndices[0] );
    MCGIDI_product *product;
    MCGIDI_quantitiesLookupModes mode( projectilesPOPID );

    if( ( product = MCGIDI_outputChannel_getProductAtIndex( &smr, &(reaction->outputChannel), 0 ) ) == NULL ) {
        smr_print( &smr, 1 );
        throw 1;
    }

    mode.setProjectileEnergy( e_in );
    mode.setCrossSectionMode( MCGIDI_quantityLookupMode_pointwise );
    mode.setTemperature( temperature );

    decaySamplingInfo.isVelocity = 0;
    decaySamplingInfo.rng = rng;
    decaySamplingInfo.rngState = rngState;
    if( MCGIDI_product_sampleMu( &smr, product, mode, &decaySamplingInfo ) ) {
        smr_print( &smr, 1 );
        throw 1;
    }

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
vector<G4GIDI_Product> *G4GIDI_target::getFinalState( int nIndices, int *indices, double e_in, double temperature, 
    double (*rng)( void * ), void *rngState ) {

    int index = 0, i, n;
    vector<G4GIDI_Product> *products = NULL;
    MCGIDI_decaySamplingInfo decaySamplingInfo;
    MCGIDI_sampledProductsDatas sampledProductsDatas;
    MCGIDI_sampledProductsData *productData;
    MCGIDI_quantitiesLookupModes mode( projectilesPOPID );

    decaySamplingInfo.isVelocity = 0;
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

    MCGIDI_sampledProducts_initialize( &smr, &sampledProductsDatas, 1000 );
    if( !smr_isOk( &smr ) ) {
        smr_print( &smr, 1 );
        throw 1;
    }

    mode.setProjectileEnergy( e_in );
    mode.setCrossSectionMode( MCGIDI_quantityLookupMode_pointwise );
    mode.setTemperature( temperature );

    n = MCGIDI_target_heated_sampleIndexReactionProductsAtE( &smr, target->baseHeatedTarget, index, mode,
            &decaySamplingInfo, &sampledProductsDatas );
    if( !smr_isOk( &smr ) ) {
        smr_print( &smr, 1 );
        throw 1;
    }
    if( n > 0 ) {
        if( ( products = new vector<G4GIDI_Product>( n ) ) != NULL ) {
            for( i = 0; i < n; i++ ) {
                productData = &(sampledProductsDatas.products[i]);
                (*products)[i].A = productData->pop->A;
                (*products)[i].Z = productData->pop->Z;
                (*products)[i].m = productData->pop->m;
                (*products)[i].kineticEnergy = productData->kineticEnergy;
                (*products)[i].px = productData->px_vx;
                (*products)[i].py = productData->py_vy;
                (*products)[i].pz = productData->pz_vz;
	        (*products)[i].birthTimeSec = productData->birthTimeSec;
            }
        }
    }
    MCGIDI_sampledProducts_release( &smr, &sampledProductsDatas );

    return( products );
}
/*
***************************************************************
*/
double G4GIDI_target::getReactionsThreshold( int index ) {

    return( MCGIDI_target_heated_getReactionsThreshold( &smr, target->baseHeatedTarget, index ) );
}
/*
***************************************************************
*/
double G4GIDI_target::getReactionsDomain( int index, double *EMin, double *EMax ) {

    return( MCGIDI_target_heated_getReactionsDomain( &smr, target->baseHeatedTarget, index, EMin, EMax ) );
}
