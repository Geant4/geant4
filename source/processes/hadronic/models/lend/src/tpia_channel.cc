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

static int _tpia_channel_getProductData( statusMessageReporting *smr, xData_element *channelElement, tpia_channel *channel );
/*
************************************************************
*/
tpia_channel *tpia_channel_create( statusMessageReporting *smr ) {

    tpia_channel *channel;

    //if( ( channel = xData_malloc2( smr, sizeof( tpia_channel ), 0, "channel" ) ) == NULL ) return( NULL );
    if( ( channel = (tpia_channel*) xData_malloc2( smr, sizeof( tpia_channel ), 0, "channel" ) ) == NULL ) return( NULL );
    if( tpia_channel_initialize( smr, channel ) ) channel = tpia_channel_free( smr, channel );
    return( channel );
}
/*
************************************************************
*/
//int tpia_channel_initialize( statusMessageReporting *smr, tpia_channel *channel ) {
int tpia_channel_initialize( statusMessageReporting *, tpia_channel *channel ) {

    memset( channel, 0, sizeof( tpia_channel ) );
    return( 0 );
}
/*
************************************************************
*/
tpia_channel *tpia_channel_createGetFromElement( statusMessageReporting *smr, tpia_target_heated *target, xData_element *channelElement,
    int pointwiseRequired ) {

    tpia_channel *channel;

    if( ( channel = tpia_channel_create( smr ) ) == NULL ) return( NULL );
    if( tpia_channel_getFromElement( smr, target, channelElement, channel, pointwiseRequired ) != 0 ) channel = tpia_channel_free( smr, channel );
    return( channel );
}
/*
************************************************************
*/
tpia_channel *tpia_channel_free( statusMessageReporting *smr, tpia_channel *channel ) {

    tpia_channel_release( smr, channel );
    xData_free( smr, channel );
    return( NULL );
}
/*
************************************************************
*/
int tpia_channel_release( statusMessageReporting *smr, tpia_channel *channel ) {

    tpia_product *product, *nextProduct;

    xData_releaseAttributionList( smr, &(channel->attributes) );
    //channel->crossSectionPointwise.data = xData_free( smr, channel->crossSectionPointwise.data );
    channel->crossSectionPointwise.data = (double*) xData_free( smr, channel->crossSectionPointwise.data );
    //channel->crossSectionGrouped.data = xData_free( smr, channel->crossSectionGrouped.data );
    channel->crossSectionGrouped.data = (double*) xData_free( smr, channel->crossSectionGrouped.data );
    //channel->availableEnergyGrouped.data = xData_free( smr, channel->availableEnergyGrouped.data );
    channel->availableEnergyGrouped.data = (double*) xData_free( smr, channel->availableEnergyGrouped.data );
    for( product = channel->decayChannel.products; product != NULL; product = nextProduct ) {
        nextProduct = product->next;
        tpia_product_free( smr, product );
    }
    channel->decayChannel.numberOfProducts = 0;
    channel->decayChannel.products = NULL;
    return( 0 );
}
/*
************************************************************
*/
int tpia_channel_getFromElement( statusMessageReporting *smr, tpia_target_heated *target, xData_element *channelElement, tpia_channel *channel, 
    int pointwiseRequired ) {

    xData_Int ll;
    char *p;
    xData_element *element, *pElement, *gElement, *eElement;

    xData_addToAccessed( smr, channelElement, 1 );
    channel->target = target;
    xData_copyAttributionList( smr, &(channel->attributes), &(channelElement->attributes) );
    channel->outputChannel = tpia_misc_pointerToAttributeIfAllOk2(smr, channelElement, 1, &(channel->attributes), "outputChannel" );
    channel->genre = tpia_misc_pointerToAttributeIfAllOk2( smr, channelElement, 1, &(channel->attributes), "genre" );
    channel->QString = tpia_misc_pointerToAttributeIfAllOk2( smr, channelElement, 1, &(channel->attributes), "Q" );
    channel->fission = tpia_misc_pointerToAttributeIfAllOk2( smr, channelElement, 0, &(channel->attributes), "fission" );
    if( smr_isOk( smr ) ) {
        ll = 0;
        if( xData_convertAttributeTo_xData_Int( smr, channelElement, "ENDL_C", &ll ) >= 0 ) channel->ENDL_C = (int) ll;
    }
    if( smr_isOk( smr ) ) {
        ll = 0;
        if( xData_convertAttributeTo_xData_Int( smr, channelElement, "ENDF_MT2", &ll ) >= 0 ) channel->ENDF_MT = (int) ll;
    }
    if( smr_isOk( smr ) ) {
        channel->QIsFloat = 1;
        channel->Q = strtod( channel->QString, &p );      /* Q string may be something like "notApplicable". */
        if( *p != 0 ) {                                     /* In that case set QIsFloat to false. */
            channel->QIsFloat = 0;
            channel->Q = 0.;
        }
        //if( ( element = xData_getOneElementByTagName( smr, channelElement, "crossSection", 1 ) ) != NULL ) {
        if( ( element = xData_getOneElementByTagName( smr, channelElement, (char*) "crossSection", 1 ) ) != NULL ) {
            if( ( tpia_frame_setFromElement( smr, element, 2, &channel->crossSectionFrame ) ) == 0 ) {
                xData_addToAccessed( smr, element, 1 );
                //if( ( pElement = xData_getOneElementByTagName( smr, element, "indexed", 1 ) ) != NULL ) {
                if( ( pElement = xData_getOneElementByTagName( smr, element, (char*) "indexed", 1 ) ) != NULL ) {
                    channel->crossSectionPointwise.data = tpia_misc_get2dxindex_y_data( smr, pElement,
                        &(channel->crossSectionPointwise.start), &(channel->crossSectionPointwise.end), target->energyGrid );
                }
                if( ( gElement = xData_getOneElementByTagName( smr, element, (char*) "grouped", 1 ) ) != NULL ) {
                    tpia_misc_get2d_xShared_yHistogram_data_Grouped( smr, gElement, &(channel->crossSectionGrouped) );
                }
                if( ( channel->crossSectionGrouped.data != NULL ) && ( ( channel->crossSectionPointwise.data != NULL ) || !pointwiseRequired ) ) {
                    if( target->contents == NULL ) {                /* Only supported "crossSection" currently. */
                        if( !tpia_channel_isProduction( smr, channel ) ) {
                            //if( ( eElement = xData_getOneElementByTagName( smr, channelElement, "availableEnergy", 1 ) ) != NULL ) {
                            if( ( eElement = xData_getOneElementByTagName( smr, channelElement, (char*) "availableEnergy", 1 ) ) != NULL ) {
                                xData_addToAccessed( smr, eElement, 1 );
                                //if( ( gElement = xData_getOneElementByTagName( smr, eElement, "grouped", 1 ) ) != NULL ) {
                                if( ( gElement = xData_getOneElementByTagName( smr, eElement, (char*) "grouped", 1 ) ) != NULL ) {
                                    tpia_misc_get2d_xShared_yHistogram_data_Grouped( smr, gElement, &(channel->availableEnergyGrouped) );
                                }
                            }
                        }
                    }
                    if( smr_isOk( smr ) ) _tpia_channel_getProductData( smr, channelElement, channel );
                }
            }
        }
    }
    return( !smr_isOk( smr ) );
}
/*
************************************************************
*/
static int _tpia_channel_getProductData( statusMessageReporting *smr, xData_element *channelElement, tpia_channel *channel ) {

    return( tpia_product_getDecayChannelFromElement( smr, channelElement, channel, NULL, &(channel->decayChannel.products) ) );
}
/*
************************************************************
*/
tpia_product *tpia_channel_getFirstProduct( tpia_channel *channel ) {

    return( tpia_decayChannel_getFirstProduct( &(channel->decayChannel) ) );
}
/*
************************************************************
*/
//tpia_product *tpia_channel_getProductByIndex( statusMessageReporting *smr, tpia_channel *channel, int index ) {
tpia_product *tpia_channel_getProductByIndex( statusMessageReporting *, tpia_channel *channel, int index ) {

    int i = 0;
    tpia_product *p;

    if( index < 0 ) return( NULL );
    for( p = tpia_channel_getFirstProduct( channel ); ( p != NULL ) && ( i < index ); p = tpia_decayChannel_getNextProduct( p ), i++ ) ;
    return( p );
}
/*
************************************************************
*/
//int tpia_channel_numberOfProducts( statusMessageReporting *smr, tpia_channel *channel ) {
int tpia_channel_numberOfProducts( statusMessageReporting *, tpia_channel *channel ) {

    return( channel->decayChannel.numberOfProducts );
}
/*
************************************************************
*/
//int tpia_channel_isProduction( statusMessageReporting *smr, tpia_channel *channel ) {
int tpia_channel_isProduction( statusMessageReporting *, tpia_channel *channel ) {

    return( strcmp( channel->genre, "production" ) == 0 );
}
/*
************************************************************
*/
//double tpia_channel_getCrossSectionAtE( statusMessageReporting *smr, tpia_channel *channel, xData_Int iEg, double e_in,
double tpia_channel_getCrossSectionAtE( statusMessageReporting *smr, tpia_channel *channel, xData_Int , double e_in,
        int crossSectionType ) {

    double xsec = 0.;

    if( crossSectionType == tpia_crossSectionType_grouped ) {
        xsec = 0; }
    else if( crossSectionType == tpia_crossSectionType_pointwise ) {
        xsec = tpia_misc_getPointwiseCrossSectionAtE( smr, &(channel->crossSectionPointwise), channel->target->energyGrid,
            tpia_target_heated_getEIndex( channel->target, e_in ), e_in );
    }
    return( xsec );
}

}
