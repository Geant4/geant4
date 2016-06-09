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
#include <string.h>
#include <tpia_target.h>
#include <tpia_misc.h>

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

static void _tpia_target_heated_getEnergyGridAndAllocateTotalCrossSections( statusMessageReporting *smr, tpia_target_heated *target, xData_element *element );
static int _tpia_target_heated_releaseAndReturnOne( statusMessageReporting *smr, xData_document *doc, tpia_target_heated *target );
static int _tpia_target_heated_checkElememtsForAccess( statusMessageReporting *smr, xData_document *doc );
static int _tpia_target_heated_checkElememtsForAccess2( statusMessageReporting *smr, xData_element *element );
/*
************************************************************
*/
tpia_target_heated *tpia_target_heated_create( statusMessageReporting *smr ) {

    tpia_target_heated *target;

    //if( ( target = xData_malloc2( smr, sizeof( tpia_target_heated ), 0, "target" ) ) == NULL ) return( NULL );
    if( ( target = (tpia_target_heated*) xData_malloc2( smr, sizeof( tpia_target_heated ), 0, "target" ) ) == NULL ) return( NULL );
    //if( tpia_target_heated_initialize( smr, target ) ) target = xData_free( smr, target );
    if( tpia_target_heated_initialize( smr, target ) ) target = (tpia_target_heated*) xData_free( smr, target );
    return( target );
}
/*
************************************************************
*/
//int tpia_target_heated_initialize( statusMessageReporting *smr, tpia_target_heated *target ) {
int tpia_target_heated_initialize( statusMessageReporting *, tpia_target_heated *target ) {

    memset( target, 0, sizeof( tpia_target_heated ) );
    return( 0 );
}
/*
************************************************************
*/
tpia_target_heated *tpia_target_heated_createRead( statusMessageReporting *smr, const char *fileName, int checkElememtsForAccess ) {

    tpia_target_heated *target;

    if( ( target = tpia_target_heated_create( smr ) ) == NULL ) return( NULL );
    //if( tpia_target_heated_read( smr, target, fileName, checkElememtsForAccess ) != 0 ) target = xData_free( smr, target );
    if( tpia_target_heated_read( smr, target, fileName, checkElememtsForAccess ) != 0 ) target = (tpia_target_heated*) xData_free( smr, target );
    return( target );
}
/*
************************************************************
*/
tpia_target_heated *tpia_target_heated_free( statusMessageReporting *smr, tpia_target_heated *target ) {

    tpia_target_heated_release( smr, target );
    xData_free( smr, target );
    return( NULL );
}
/*
************************************************************
*/
int tpia_target_heated_release( statusMessageReporting *smr, tpia_target_heated *target ) {

    int i;

    //target->path = xData_free( smr, target->path );
    target->path = (char*) xData_free( smr, target->path );
    //target->absPath = xData_free( smr, target->absPath );
    target->absPath = (char*) xData_free( smr, target->absPath );
    target->energyGridLength = 0;
    //target->energyGrid = xData_free( smr, target->energyGrid );
    target->energyGrid = (double*) xData_free( smr, target->energyGrid );
    //target->kerma = xData_free( smr, target->kerma );
    target->kerma = (double*) xData_free( smr, target->kerma );
    //target->totalCrossSectionPointwise.data = xData_free( smr, target->totalCrossSectionPointwise.data );
    target->totalCrossSectionPointwise.data = (double*) xData_free( smr, target->totalCrossSectionPointwise.data );
    //target->totalCrossSectionGrouped.data = xData_free( smr, target->totalCrossSectionGrouped.data );
    target->totalCrossSectionGrouped.data = (double*) xData_free( smr, target->totalCrossSectionGrouped.data );
    xData_releaseAttributionList( smr, &(target->attributes) );
    for( i = 0; i < target->nChannels; i++ ) tpia_channel_free( smr, target->channels[i] );
    target->nChannels = 0;
    //target->channels = xData_free( smr, target->channels );
    target->channels = (tpia_channel**) xData_free( smr, target->channels );
    for( i = 0; i < target->nProductionChannels; i++ ) tpia_channel_free( smr, target->productionChannels[i] );
    target->nProductionChannels = 0;
    //target->productionChannels = xData_free( smr, target->productionChannels );
    target->productionChannels = (tpia_channel**) xData_free( smr, target->productionChannels );
    return( 0 );
}
/*
************************************************************
*/
int tpia_target_heated_read( statusMessageReporting *smr, tpia_target_heated *target, const char *fileName, int checkElememtsForAccess ) {
/*
*   If a target has already been read into this target, user must have called tpia_target_heated_release before calling this routine.
*   Otherwise, there will be memory leaks.
*/
    xData_document *doc = NULL;
    xData_element *element, *channelElement, *channels;
    int nChannels;
    tpia_channel *channel;
    char *name;
    xData_Int i, j;

    tpia_target_heated_initialize( smr, target );
    if( ( target->path = xDataMisc_allocateCopyString2( smr, fileName, "path" ) ) == NULL ) return( 1 );
    if( ( target->absPath = xDataMisc_getAbsPath( smr, fileName ) ) == NULL ) return( _tpia_target_heated_releaseAndReturnOne( smr, doc, target ) );
    if( ( doc = xData_parseReadFile( smr, fileName, NULL, NULL ) ) == NULL ) return( _tpia_target_heated_releaseAndReturnOne( smr, doc, target ) );
    element = xData_getDocumentsElement( doc );
    xData_addToAccessed( smr, element, 1 );
    if( xData_convertAttributeTo_xData_Int( smr, element, "nGroups", &i ) != 0 ) return( _tpia_target_heated_releaseAndReturnOne( smr, doc, target ) );
    target->nGroups = (int) i;
    if( strcmp( element->name, "xTargetHeated" ) != 0 ) {
        tpia_misc_setMessageError_Element( smr, NULL, element, __FILE__, __LINE__, 1, "input file's top element must be xTargetHeated and not %s", 
            element->name ); }
    else {
        xData_copyAttributionList( smr, &(target->attributes), &(element->attributes) );
        if( smr_isOk( smr ) ) target->contents = xData_getAttributesValue( &(target->attributes), "contents" );
        if( ( name = tpia_misc_pointerToAttributeIfAllOk3( smr, target->absPath, 1, &(target->attributes), "projectile" ) ) != NULL )
            target->projectileID = tpia_particle_getInternalID( smr, name );
        if( ( name = tpia_misc_pointerToAttributeIfAllOk3( smr, target->absPath, 1, &(target->attributes), "target" ) ) != NULL )
            target->targetID = tpia_particle_getInternalID( smr, name );
        if( smr_isOk( smr ) ) _tpia_target_heated_getEnergyGridAndAllocateTotalCrossSections( smr, target, element );
        if( smr_isOk( smr ) ) {            /* Get channels. */
            //if( ( channels = xData_getOneElementByTagName( smr, element, "channels", 1 ) ) == NULL ) 
            if( ( channels = xData_getOneElementByTagName( smr, element, (char*)"channels", 1 ) ) == NULL ) 
                return( _tpia_target_heated_releaseAndReturnOne( smr, doc, target ) );
            xData_addToAccessed( smr, channels, 1 );
            if( ( nChannels = xData_numberOfElementsByTagName( smr, channels, "channel" ) ) > 0 ) {
                //if( ( target->channels = xData_malloc2( smr, nChannels * sizeof( tpia_channel * ), 1, "channels" ) ) == NULL ) 
                if( ( target->channels = (tpia_channel**) xData_malloc2( smr, nChannels * sizeof( tpia_channel * ), 1, "channels" ) ) == NULL ) 
                    return( _tpia_target_heated_releaseAndReturnOne( smr, doc, target ) );
                for( channelElement = xData_getFirstElement( channels ); channelElement != NULL; channelElement = xData_getNextElement( channelElement ) ) {
                    if( !strcmp( channelElement->name, "channel" ) ) {
                        if( ( channel = tpia_channel_createGetFromElement( smr, target, channelElement, 1 ) ) == NULL ) break;
                        target->channels[target->nChannels] = channel;
                        target->nChannels++;
                        for( i = channel->crossSectionPointwise.start, j = 0; i < channel->crossSectionPointwise.end; i++, j++ ) 
                            target->totalCrossSectionPointwise.data[i] += channel->crossSectionPointwise.data[j];
                        for( i = channel->crossSectionGrouped.start, j = 0; i < channel->crossSectionGrouped.end; i++, j++ ) 
                            target->totalCrossSectionGrouped.data[i] += channel->crossSectionGrouped.data[j];
                    }
                }
            }
        }
        if( smr_isOk( smr ) ) {            /* Get production channels. */
            //if( ( channels = xData_getOneElementByTagName( smr, element, "productionChannels", 0 ) ) == NULL ) {
            if( ( channels = xData_getOneElementByTagName( smr, element, (char*) "productionChannels", 0 ) ) == NULL ) {
                if( !smr_isOk( smr ) ) return( _tpia_target_heated_releaseAndReturnOne( smr, doc, target ) ); }
            else {
                xData_addToAccessed( smr, channels, 1 );
                if( ( nChannels = xData_numberOfElementsByTagName( smr, channels, "channel" ) ) > 0 ) {
                    //if( ( target->productionChannels = xData_malloc2( smr, nChannels * sizeof( tpia_channel * ), 1, "channels" ) ) != NULL ) {
                    if( ( target->productionChannels = (tpia_channel**) xData_malloc2( smr, nChannels * sizeof( tpia_channel * ), 1, "channels" ) ) != NULL ) {
                        for( channelElement = xData_getFirstElement(channels); channelElement != NULL; channelElement = xData_getNextElement(channelElement) ) {
                            if( !strcmp( channelElement->name, "channel" ) ) {
                                channel = tpia_channel_createGetFromElement( smr, target, channelElement, 1 );
                                if( channel == NULL ) break;
                                target->productionChannels[target->nProductionChannels] = channel;
                                target->nProductionChannels++;
                            }
                        }
                    }
                }
            }
        }
    }
    if( smr_isOk( smr ) && checkElememtsForAccess ) _tpia_target_heated_checkElememtsForAccess( smr, doc );
    xData_parseFree( smr, doc );
    if( !smr_isOk( smr ) ) tpia_target_heated_release( smr, target );
    return( !smr_isOk( smr ) );
}
/*
************************************************************
*/
static void _tpia_target_heated_getEnergyGridAndAllocateTotalCrossSections( statusMessageReporting *smr, tpia_target_heated *target, xData_element *element ) {

    xData_Int i, energyGridIndex, energyGridStart, energyGridEnd, energyGridLength;
    xData_element *energyGrid, *kerma;

    if( !smr_isOk( smr ) ) return;
    //if( ( energyGrid = xData_getOneElementByTagName( smr, element, "energyGrid", 1 ) ) == NULL ) return;
    if( ( energyGrid = xData_getOneElementByTagName( smr, element, (char*) "energyGrid", 1 ) ) == NULL ) return;
    xData_addToAccessed( smr, energyGrid, 1 );
    if( ( energyGrid = xData_getElements_xDataElement( smr, energyGrid ) ) == NULL ) return;
    xData_addToAccessed( smr, energyGrid, 1 );
    xData_getCommonData( smr, energyGrid, &energyGridIndex, &energyGridStart, &energyGridEnd, &energyGridLength );
    if( ( target->energyGrid = xData_1d_x_allocateCopyData( smr, energyGrid ) ) == NULL ) return;
    target->energyGridLength = energyGridLength;
    target->totalCrossSectionPointwise.start = 0;
    target->totalCrossSectionPointwise.end = energyGridLength;
    target->totalCrossSectionPointwise.length = energyGridLength;
    //if( ( target->totalCrossSectionPointwise.data = xData_malloc2( smr, energyGridLength * sizeof( double ), 0, "totalCrossSectionPointwise" ) ) == NULL ) 
    if( ( target->totalCrossSectionPointwise.data = (double*) xData_malloc2( smr, energyGridLength * sizeof( double ), 0, "totalCrossSectionPointwise" ) ) == NULL ) 
        return;
    for( i = 0; i < energyGridLength; i++ ) target->totalCrossSectionPointwise.data[i] = 0.;
    target->totalCrossSectionGrouped.start = 0;
    target->totalCrossSectionGrouped.end = energyGridLength;
    target->totalCrossSectionGrouped.length = energyGridLength;
    //if( ( target->totalCrossSectionGrouped.data = xData_malloc2( smr, target->nGroups * sizeof( double ), 0, "totalCrossSectionGrouped" ) ) == NULL ) return;
    if( ( target->totalCrossSectionGrouped.data = (double*) xData_malloc2( smr, target->nGroups * sizeof( double ), 0, "totalCrossSectionGrouped" ) ) == NULL ) return;
    for( i = 0; i < target->nGroups; i++ ) target->totalCrossSectionGrouped.data[i] = 0.;
    //if( ( kerma = xData_getOneElementByTagName( smr, element, "kerma", 1 ) ) != NULL ) { 
    if( ( kerma = xData_getOneElementByTagName( smr, element, (char*) "kerma", 1 ) ) != NULL ) { 
        xData_addToAccessed( smr, kerma, 1 );
        if( ( kerma = xData_getElements_xDataElement( smr, kerma ) ) == NULL ) return;
        xData_addToAccessed( smr, kerma, 1 );
        if( ( target->kerma = xData_1d_x_allocateCopyData( smr, kerma ) ) == NULL ) return;
    }
}
/*
************************************************************
*/
//int tpia_target_heated_numberOfChannels( statusMessageReporting *smr, tpia_target_heated *target ) {
int tpia_target_heated_numberOfChannels( statusMessageReporting *, tpia_target_heated *target ) {

    return( target->nChannels );
}
/*
************************************************************
*/
//int tpia_target_heated_numberOfProductionChannels( statusMessageReporting *smr, tpia_target_heated *target ) {
int tpia_target_heated_numberOfProductionChannels( statusMessageReporting *, tpia_target_heated *target ) {

    return( target->nProductionChannels );
}
/*
************************************************************
*/
tpia_channel *tpia_target_heated_getChannelAtIndex( tpia_target_heated *target, int index ) {

    tpia_channel *channel = NULL;

    if( ( index >= 0 ) && ( index < target->nChannels ) ) channel = target->channels[index];
    return( channel );
}
/*
************************************************************
*/
tpia_channel *tpia_target_heated_getChannelAtIndex_smr( statusMessageReporting *smr, tpia_target_heated *target, int index ) {

    tpia_channel *channel = tpia_target_heated_getChannelAtIndex( target, index );

    if( channel == NULL ) {
        smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "bad channel index = %d for %s + %s", index, 
            target->projectileID->name, target->targetID->name ); 
    }
    return( channel );
}
/*
************************************************************
*/
tpia_channel *tpia_target_heated_getProductionChannelAtIndex( tpia_target_heated *target, int index ) {

    tpia_channel *channel = NULL;

    if( ( index >= 0 ) && ( index < target->nProductionChannels ) ) channel = target->productionChannels[index];
    return( channel );
}
/*
************************************************************
*/
//xData_Int tpia_target_heated_getEnergyGrid( statusMessageReporting *smr, tpia_target_heated *target, double **energyGrid ) {
xData_Int tpia_target_heated_getEnergyGrid( statusMessageReporting *, tpia_target_heated *target, double **energyGrid ) {

    if( energyGrid != NULL ) *energyGrid = target->energyGrid;
    return( target->energyGridLength );
}
/*
************************************************************
*/
xData_Int tpia_target_heated_getEIndex( tpia_target_heated *target, double e_in ) {

    return( tpia_misc_binarySearch( target->energyGridLength, target->energyGrid, e_in ) );
}
/*
************************************************************
*/
//double tpia_target_heated_getTotalCrossSectionAtE( statusMessageReporting *smr, tpia_target_heated *target, xData_Int iEg, double e_in, 
double tpia_target_heated_getTotalCrossSectionAtE( statusMessageReporting *smr, tpia_target_heated *target, xData_Int , double e_in, 
        int crossSectionType ) {

    double xsec = 0.;

    if( crossSectionType == tpia_crossSectionType_grouped ) {
        xsec = 0; }
    else if( crossSectionType == tpia_crossSectionType_pointwise ) {
        xsec = tpia_misc_getPointwiseCrossSectionAtE( smr, &(target->totalCrossSectionPointwise), target->energyGrid, 
            tpia_target_heated_getEIndex( target, e_in ), e_in );
    }
    return( xsec );
}
/*
************************************************************
*/
double tpia_target_heated_getIndexChannelCrossSectionAtE( statusMessageReporting *smr, tpia_target_heated *target, int index, xData_Int iEg, double e_in, 
        int crossSectionType ) {

    double xsec = 0.;
    tpia_channel *channel = tpia_target_heated_getChannelAtIndex_smr( smr, target, index );

    if( channel != NULL ) xsec = tpia_channel_getCrossSectionAtE( smr, channel, iEg, e_in, crossSectionType );
    return( xsec );
}
/*
************************************************************
*/
int tpia_target_heated_sampleIndexChannelProductsAtE( statusMessageReporting *smr, tpia_target_heated *target, int index, 
        tpia_decaySamplingInfo *decaySamplingInfo, int nProductData, tpia_productOutgoingData *productDatas ) {

    tpia_channel *channel = tpia_target_heated_getChannelAtIndex_smr( smr, target, index );

    if( channel == NULL ) return( -1 );
    return( tpia_decayChannel_sampleProductsAtE( smr, &(channel->decayChannel), decaySamplingInfo, nProductData, productDatas ) );
}
/*
************************************************************
*/
static int _tpia_target_heated_releaseAndReturnOne( statusMessageReporting *smr, xData_document *doc, tpia_target_heated *target ) {

    tpia_target_heated_release( smr, target );
    if( doc != NULL ) xData_parseFree( smr, doc );
    return( 1 );
}
/*
************************************************************
*/
static int _tpia_target_heated_checkElememtsForAccess( statusMessageReporting *smr, xData_document *doc ) {

    xData_element *element = xData_getDocumentsElement( doc );

    _tpia_target_heated_checkElememtsForAccess2( smr, element );
    return( 0 );
}
/*
************************************************************
*/
static int _tpia_target_heated_checkElememtsForAccess2( statusMessageReporting *smr, xData_element *element ) {

    xData_element *child;

    if( xData_getAccessed( smr, element ) != 1 ) printf( "%3d %s\n", xData_getAccessed( smr, element ), element->fullName );
    for( child = xData_getFirstElement( element ); child != NULL; child = xData_getNextElement( child ) ) 
        _tpia_target_heated_checkElememtsForAccess2( smr, child );
    return( 0 );
}

#if defined __cplusplus
}
#endif
