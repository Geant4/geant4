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
#include <ctype.h>

#include <gString.h>
#include <tpia_target.h>
#include <tpia_misc.h>

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

//static const int tpia_b_unknown = 0,
static const int /*tpia_b_unknown = 0,*/
                 tpia_b_twoBody_angular = tpia_m_angular,
                 tpia_b_twoBody_formFactor = 0,                                       /* ??? */
                 tpia_b_NBody_Legendre = tpia_m_Legendre,
                 tpia_b_NBody_angular_energy = tpia_m_angular | tpia_m_angular_energy,
                 tpia_b_NBody_uncorrelate_Legendre = tpia_m_angular | tpia_m_Legendre,
                 tpia_b_NBody_pairProduction = 0;                                     /* ??? */

const char *tpia_productGenre_unknown = "unknown",
    *tpia_productGenre_twoBody_angular = "twoBody_angular",
    *tpia_productGenre_twoBody_formFactor = "twoBody_formFactor",
    *tpia_productGenre_NBody_Legendre = "NBody_Legendre",
    *tpia_productGenre_NBody_angular_energy = "NBody_angular_energy",
    *tpia_productGenre_NBody_uncorrelate_Legendre = "NBody_uncorrelate_Legendre",
    *tpia_productGenre_NBody_pairProduction = "NBody_pairProduction";

static int _tpia_product_getProductOutgoingData( statusMessageReporting *smr, xData_element *productElement, tpia_product *product );
static int _tpia_product_checkRequiredData( statusMessageReporting *smr, int allowMany, int m, xData_element *productElement, tpia_product *product, char *str );
static int _tpia_product_getDepositionEnergy( statusMessageReporting *smr, xData_element *depositionEnergy, tpia_product *product );
static int _tpia_product_getMultiplicityFromElement( statusMessageReporting *smr, xData_element *data, tpia_product *product );
/*
************************************************************
*/
tpia_product *tpia_product_create( statusMessageReporting *smr ) {

    tpia_product *product;

    //if( ( product = xData_malloc2( smr, sizeof( tpia_product ), 0, "product" ) ) == NULL ) return( NULL );
    if( ( product = (tpia_product*) xData_malloc2( smr, sizeof( tpia_product ), 0, "product" ) ) == NULL ) return( NULL );
    if( tpia_product_initialize( smr, product ) ) product = tpia_product_free( smr, product );
    return( product );
}
/*
************************************************************
*/
int tpia_product_initialize( statusMessageReporting *smr, tpia_product *product ) {

    memset( product, 0, sizeof( tpia_product ) );
    if( tpia_angular_initialize( smr, &(product->angular) ) ) return( 1 );
    if( tpia_Legendre_initialize( smr, &(product->Legendre) ) ) return( 1 );
    return( 0 );
}
/*
************************************************************
*/
tpia_product *tpia_product_createGetFromElement( statusMessageReporting *smr, tpia_channel *channel, tpia_product *parentProduct, 
    xData_element *productElement ) {

    tpia_product *product;

    if( ( product = tpia_product_create( smr ) ) == NULL ) return( NULL );
    if( tpia_product_getFromElement( smr, channel, parentProduct, productElement, product ) != 0 ) product = tpia_product_free( smr, product );
    return( product );
}
/*
************************************************************
*/
tpia_product *tpia_product_free( statusMessageReporting *smr, tpia_product *product ) {

    tpia_product_release( smr, product );
    xData_free( smr, product );
    return( NULL );
}
/*
************************************************************
*/
int tpia_product_release( statusMessageReporting *smr, tpia_product *product) {

    tpia_multiplicity *multiplicity, *multiplicity_next;
    tpia_product *decayProduct, *nextProduct;

    xData_releaseAttributionList( smr, &(product->attributes) );
    //product->depositionEnergyGrouped.data = xData_free( smr, product->depositionEnergyGrouped.data );
    product->depositionEnergyGrouped.data = (double*) xData_free( smr, product->depositionEnergyGrouped.data );

    if( product->multiplicityVsEnergy != NULL ) tpia_multiplicity_free( smr, product->multiplicityVsEnergy );
    for( multiplicity = product->delayedNeutronMultiplicityVsEnergy; multiplicity != NULL; multiplicity = multiplicity_next ) {
        multiplicity_next = multiplicity->next;
        tpia_multiplicity_free( smr, multiplicity );
    }
    tpia_angular_release( smr, &(product->angular) );
    tpia_Legendre_release( smr, &(product->Legendre ) );
    tpia_angularEnergy_release( smr, &(product->angularEnergy) );
    for( decayProduct = product->decayChannel.products; decayProduct != NULL; decayProduct = nextProduct ) {
        nextProduct = decayProduct->next;
        tpia_product_free( smr, decayProduct );
    }
    product->decayChannel.numberOfProducts = 0;
    product->decayChannel.products = NULL;
    return( 0 );
}
/*
************************************************************
*/
int tpia_product_getFromElement( statusMessageReporting *smr, tpia_channel *channel, tpia_product *parentProduct, xData_element *productElement, 
    tpia_product *product ) {

    char const *productGenre;
    char *name, *multiplicity, *e;

    xData_addToAccessed( smr, productElement, 1 );
    product->channel = channel;
    product->parentProduct = parentProduct;
    if( xData_copyAttributionList( smr, &(product->attributes), &(productElement->attributes) ) != 0 ) return( 0 );
    name = tpia_misc_pointerToAttributeIfAllOk2( smr, productElement, 1, &(product->attributes), "particle" );
    if( name != NULL ) {
        product->productID = tpia_particle_getInternalID( smr, name );
        multiplicity = tpia_misc_pointerToAttributeIfAllOk2( smr, productElement, 1, &(product->attributes), "multiplicity" );
        if( multiplicity != NULL ) {
            if( strcmp( multiplicity, "energyDependent" ) && strcmp( multiplicity, "partialProduction" ) ) {    /* Must be an integer. */
                product->multiplicity = strtol( multiplicity, &e, 10 );
                while( isspace( *e ) ) e++;
                if( *e != 0 ) tpia_misc_setMessageError_Element( smr, NULL, productElement, __FILE__, __LINE__, 1, "bad multiplicity = %s", multiplicity );
            }
        }
    }
    if( ( productGenre = tpia_misc_pointerToAttributeIfAllOk2( smr, productElement, 1, &(product->attributes), "genre" ) ) != NULL ) {
        if( strcmp( productGenre, tpia_productGenre_unknown ) == 0 ) {
            product->b_dataRequired = 0;
            product->genre = tpia_productGenre_unknown; }
        else if( strcmp( productGenre, tpia_productGenre_twoBody_angular ) == 0 ) {
            product->b_dataRequired = tpia_b_twoBody_angular;
            product->genre = tpia_productGenre_twoBody_angular; }
        else if( strcmp( productGenre, tpia_productGenre_twoBody_formFactor ) == 0 ) {
            product->b_dataRequired = tpia_b_twoBody_formFactor;
            product->genre = tpia_productGenre_twoBody_formFactor; }
        else if( strcmp( productGenre, tpia_productGenre_NBody_Legendre ) == 0 ) {
            product->b_dataRequired = tpia_b_NBody_Legendre;
            product->genre = tpia_productGenre_NBody_Legendre; }
        else if( strcmp( productGenre, tpia_productGenre_NBody_angular_energy ) == 0 ) {
            product->b_dataRequired = tpia_b_NBody_angular_energy;
            product->genre = tpia_productGenre_NBody_angular_energy; }
        else if( strcmp( productGenre, tpia_productGenre_NBody_uncorrelate_Legendre ) == 0 ) {
            product->b_dataRequired = tpia_b_NBody_uncorrelate_Legendre;
            product->genre = tpia_productGenre_NBody_uncorrelate_Legendre; }
        else if( strcmp( productGenre, tpia_productGenre_NBody_pairProduction ) == 0 ) {
            product->b_dataRequired = tpia_b_NBody_pairProduction;
            product->genre = tpia_productGenre_NBody_pairProduction; }
        else {
            tpia_misc_setMessageError_Element( smr, NULL, productElement, __FILE__, __LINE__, 1, "unsupported product genre = %s", productGenre );
        }
        if( smr_isOk( smr ) ) _tpia_product_getProductOutgoingData( smr, productElement, product );
    }
    return( !smr_isOk( smr ) );
}
/*
************************************************************
*/
static int _tpia_product_getProductOutgoingData( statusMessageReporting *smr, xData_element *productElement, tpia_product *product ) {

    xData_element *data;
    int allowMany = 0;

    for( data = xData_getFirstElement( productElement ); data != NULL; data = xData_getNextElement( data ) ) {
        if( strcmp( data->name, "depositionEnergy" ) == 0 ) {
            //if( _tpia_product_checkRequiredData( smr, allowMany, tpia_m_depositionEnergy, productElement, product, "deposition energy" ) ) return( 1 );
            if( _tpia_product_checkRequiredData( smr, allowMany, tpia_m_depositionEnergy, productElement, product, (char*)"deposition energy" ) ) return( 1 );
            if( _tpia_product_getDepositionEnergy( smr, data, product ) != 0 ) return( 1 ); }
        else if( strcmp( data->name, "multiplicity" ) == 0 ) {
            allowMany = ( product->channel->fission != NULL ) && ( strcmp( product->productID->name, "n_1" ) == 0 );
            //if( _tpia_product_checkRequiredData( smr, allowMany, tpia_m_multiplicity, productElement, product, "multiplicity" ) ) return( 1 );
            if( _tpia_product_checkRequiredData( smr, allowMany, tpia_m_multiplicity, productElement, product, (char*) "multiplicity" ) ) return( 1 );
            if( _tpia_product_getMultiplicityFromElement( smr, data, product ) != 0 ) return( 1 ); }
        else if( strcmp( data->name, "angular" ) == 0 ) {
            //if( _tpia_product_checkRequiredData( smr, allowMany, tpia_m_angular, productElement, product, "angular" ) ) return( 1 );
            if( _tpia_product_checkRequiredData( smr, allowMany, tpia_m_angular, productElement, product, (char*) "angular" ) ) return( 1 );
            if( tpia_angular_getFromElement( smr, data, &(product->angular) ) != 0 ) return( 1 ); }
        else if( strcmp( data->name, "Legendre" ) == 0 ) {
            //if( _tpia_product_checkRequiredData( smr, allowMany, tpia_m_Legendre, productElement, product, "Legendre" ) ) return( 1 );
            if( _tpia_product_checkRequiredData( smr, allowMany, tpia_m_Legendre, productElement, product, (char*) "Legendre" ) ) return( 1 );
            if( tpia_Legendre_getFromElement( smr, data, &(product->Legendre) ) != 0 ) return( 1 ); }
        else if( strcmp( data->name, "angularEnergy" ) == 0 ) {
            if( _tpia_product_checkRequiredData( smr, allowMany, tpia_m_angular_energy, productElement, product, (char*) "angularEnergy" ) ) return( 1 );
            if( tpia_angularEnergy_getFromElement( smr, data, &(product->angularEnergy) ) != 0 ) return( 1 ); }
        else if( strcmp( data->name, "decayChannel" ) == 0 ) {
            xData_addToAccessed( smr, data, 1 );
            if( tpia_product_getDecayChannelFromElement( smr, data, product->channel, product, &(product->decayChannel.products) ) ) return( 1 ); }
        else {
            printf( "   %s\n", data->name );
        }
    }
    if( ( product->b_dataPresent >> tpia_m_commonShift ) != ( product->b_dataRequired >> tpia_m_commonShift ) ) {
        gString gStr;
        int missing = ~product->b_dataPresent & product->b_dataRequired;
        char const *str = "";
        if( gString_initialize( NULL, &gStr, 100, 100 ) == 0 ) {
            if( missing & tpia_m_angular ) gString_addTo( NULL, &gStr, "angular " );
            if( missing & tpia_m_formFactor ) gString_addTo( NULL, &gStr, "formFactor " );
            if( missing & tpia_m_Legendre ) gString_addTo( NULL, &gStr, "Legendre " );
            if( missing & tpia_m_angular_energy ) gString_addTo( NULL, &gStr, "angular_energy " );
            str = gString_string( NULL, &gStr );
        }
        tpia_misc_setMessageError_Element( smr, NULL, productElement, __FILE__, __LINE__, 1, "missing data %s for product %s", str, 
            product->productID->name );
        gString_release( NULL, &gStr );
        return( 1 );
    }
    return( 0 );
}
/*
************************************************************
*/
static int _tpia_product_checkRequiredData(statusMessageReporting *smr, int allowMany, int m, xData_element *productElement, tpia_product *product, char *str) {

    if( !allowMany && ( product->b_dataPresent & m ) ) {
        tpia_misc_setMessageError_Element( smr, NULL, productElement, __FILE__, __LINE__, 1, "multiple %s", str );
        return( 1 );
    }
    if( ( m & ( tpia_m_depositionEnergy | tpia_m_multiplicity | tpia_m_decayChannel ) ) == 0 ) {
        if( ( product->b_dataRequired & m ) == 0 ) {
            tpia_misc_setMessageError_Element( smr, NULL, productElement, __FILE__, __LINE__, 1, "extra product data %s", str );
            return( 1 );
        }
    }
    product->b_dataPresent += m;
    return( 0 );
}
/*
************************************************************
*/
int tpia_product_getDecayChannelFromElement( statusMessageReporting *smr, xData_element *parentElement, tpia_channel *channel, tpia_product *parentProduct,
    tpia_product **priorProductNext ) {

    xData_elementList *list;
    tpia_product *product;
    int i, status = 0;

    list = xData_getElementsByTagName( smr, parentElement, "product" );
    for( i = 0; i < list->n; i++ ) {
        if( ( product = tpia_product_createGetFromElement( smr, channel, parentProduct, list->items[i].element ) ) == NULL ) {
            status = 1;
            break;
        }
        if( parentProduct == NULL ) {
            channel->decayChannel.m1_fullMass_MeV = channel->target->projectileID->fullMass_MeV;
            channel->decayChannel.m2_fullMass_MeV = channel->target->targetID->fullMass_MeV;
            channel->decayChannel.numberOfProducts++; }
        else {
            channel->decayChannel.m1_fullMass_MeV = parentProduct->productID->fullMass_MeV;
            channel->decayChannel.m2_fullMass_MeV = 0.;
            parentProduct->decayChannel.numberOfProducts++;
        }
        *priorProductNext = product;
        priorProductNext = &(product->next);
    }
    xData_freeElementList( smr, list );
    return( status );
}
/*
************************************************************
*/
static int _tpia_product_getDepositionEnergy( statusMessageReporting *smr, xData_element *depositionEnergy, tpia_product *product ) {

    xData_element *data;

    xData_addToAccessed( smr, depositionEnergy, 1 );
    for( data = xData_getFirstElement( depositionEnergy ); data != NULL; data = xData_getNextElement( data ) ) {
        if( strcmp( data->name, "grouped" ) == 0 ) {
            if( tpia_misc_get2d_xShared_yHistogram_data_Grouped( smr, data, &(product->depositionEnergyGrouped) ) ) return( 1 ); }
        else {
            tpia_misc_setMessageError_Element( smr, NULL, depositionEnergy, __FILE__, __LINE__, 1, "unsupported deposition energy type = %s", data->name );
            return( 1 );
        }
    }
    return( 0 );
}
/*
************************************************************
*/
static int _tpia_product_getMultiplicityFromElement( statusMessageReporting *smr, xData_element *data, tpia_product *product ) {

    tpia_multiplicity *multiplicity, *prior, *current;
    const char *timeScale;
    int isDelayedNeutrons;
    double dTimeScale;

    if( tpia_multiplicity_getTimeScaleFromElement( smr, data, &timeScale, &isDelayedNeutrons, &dTimeScale ) ) return( 1 );
    if( ( isDelayedNeutrons == 0 ) && ( product->multiplicityVsEnergy != NULL ) ) {
        tpia_misc_setMessageError_Element( smr, NULL, data, __FILE__, __LINE__, 1, "extra product multiplicity data" );
        return( 1 );
    }
    if( ( multiplicity = tpia_multiplicity_createGetFromElement( smr, data, product->channel->target->nGroups ) ) == NULL ) return( 1 );
    if( isDelayedNeutrons == 0 ) {
        product->multiplicityVsEnergy = multiplicity; }
    else {
        if( product->delayedNeutronMultiplicityVsEnergy == NULL ) {
            product->delayedNeutronMultiplicityVsEnergy = multiplicity; }
        else {
            if( product->delayedNeutronMultiplicityVsEnergy->timeScale > multiplicity->timeScale ) {
                multiplicity->next = product->delayedNeutronMultiplicityVsEnergy;
                product->delayedNeutronMultiplicityVsEnergy = multiplicity; }
            else {
                for( current = product->delayedNeutronMultiplicityVsEnergy->next, prior = product->delayedNeutronMultiplicityVsEnergy; current != NULL; 
                    current = current->next ) {
                    if( current->timeScale > multiplicity->timeScale ) {
                        multiplicity->next = current;
                        prior->next = multiplicity;
                        break;
                    }
                    prior = current;
                }
                if( current == NULL ) prior->next = multiplicity;
            }
        }
    }
    return( 0 );
}
/*
************************************************************
*/
//long tpia_product_dataRequired( statusMessageReporting *smr, tpia_product *product ) {
long tpia_product_dataRequired( statusMessageReporting *, tpia_product *product ) {

    return( product->b_dataRequired );
}
/*
************************************************************
*/
tpia_product *tpia_product_getFirstProduct( tpia_product *product ) {

    return( tpia_decayChannel_getFirstProduct( &(product->decayChannel) ) );
}
/*
************************************************************
*/
//tpia_product *tpia_product_getProductByIndex( statusMessageReporting *smr, tpia_product *product, int index ) {
tpia_product *tpia_product_getProductByIndex( statusMessageReporting *, tpia_product *product, int index ) {

    int i = 0;
    tpia_product *p;

    if( index < 0 ) return( NULL );
    for( p = tpia_product_getFirstProduct( product ); ( p != NULL ) && ( i < index ); p = tpia_decayChannel_getNextProduct( p ), i++ ) ;
    return( p );
}
/*
************************************************************
*/
//int tpia_product_doesDecay( statusMessageReporting *smr, tpia_product *product ) {
int tpia_product_doesDecay( statusMessageReporting *, tpia_product *product ) {

    return( product->decayChannel.products != NULL );
}
/*
************************************************************
*/
//int tpia_product_numberOfProducts( statusMessageReporting *smr, tpia_product *product ) {
int tpia_product_numberOfProducts( statusMessageReporting *, tpia_product *product ) {

    return( product->decayChannel.numberOfProducts );
}
/*
************************************************************
*/
//int tpia_product_isDataPresent( statusMessageReporting *smr, tpia_product *product, int b_data ) {
int tpia_product_isDataPresent( statusMessageReporting *, tpia_product *product, int b_data ) {

    return( product->b_dataPresent && b_data );
}
/*
************************************************************
*/
//int tpia_product_sampleMultiplicity( statusMessageReporting *smr, tpia_product *product, double e_in, double r ) {
int tpia_product_sampleMultiplicity( statusMessageReporting *, tpia_product *product, double e_in, double r ) {

    int i, multiplicity;
    tpia_multiplicity *multiplicityVsEnergy = product->multiplicityVsEnergy;
    double *p = multiplicityVsEnergy->pointwise, dMult;

    if( e_in <= p[0] ) {
        dMult = p[1]; }
    else if( e_in >= p[2 * ( multiplicityVsEnergy->numberOfPointwise - 1 )] ) {
        dMult = p[2 * multiplicityVsEnergy->numberOfPointwise - 1]; }
    else {
        for( i = 0; i < multiplicityVsEnergy->numberOfPointwise - 1; i++, p += 2 ) if( e_in < p[2] ) break;
        dMult = ( e_in - p[0] ) / ( p[2] - p[0] );
        dMult = dMult * p[3] + ( 1. - dMult ) * p[1];
    }
    multiplicity = (int) dMult;
    if( r < ( dMult - multiplicity ) ) multiplicity++;

    return( multiplicity );
}

#if defined __cplusplus
}
#endif
