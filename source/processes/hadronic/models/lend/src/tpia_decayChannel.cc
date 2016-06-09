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
#ifdef WIN32
   #define _USE_MATH_DEFINES
#endif
#include <cmath>
#include "tpia_target.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

/*
************************************************************
*/
tpia_product *tpia_decayChannel_getFirstProduct( tpia_decayChannel *decayChannel ) {

    return( decayChannel->products );
}
/*
************************************************************
*/
tpia_product *tpia_decayChannel_getNextProduct( tpia_product *product ) {

    return( product->next );
}
/*
************************************************************
*/
int tpia_decayChannel_sampleProductsAtE( statusMessageReporting *smr, tpia_decayChannel *decayChannel, tpia_decaySamplingInfo *decaySamplingInfo,
        int nProductData, tpia_productOutgoingData *productDatas ) {

    int i, n = 0, multiplicity, secondTwoBody = 0, labFrame = tpia_referenceFrame_lab;
    tpia_product *product, *nextProduct;
    double phi, p;

    if( nProductData < decayChannel->numberOfProducts ) {
        smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "nProductData = %d < decayChannel->numberOfProducts = %d", nProductData,
            decayChannel->numberOfProducts );
        return( -1 );
    }
    for( i = 0, product = tpia_decayChannel_getFirstProduct( decayChannel ); product != NULL; i++, product = tpia_decayChannel_getNextProduct( product ) ) {
        nextProduct = tpia_decayChannel_getNextProduct( product );
        if( !secondTwoBody ) {
            if( ( multiplicity = product->multiplicity ) == 0 ) multiplicity = tpia_product_sampleMultiplicity( smr, product, decaySamplingInfo->e_in,
                tpia_misc_drng( decaySamplingInfo->rng, decaySamplingInfo->rngState ) );
            while( multiplicity > 0 ) {
                if( n >= nProductData ) break;          /* This needs work. */
                multiplicity--;
                decaySamplingInfo->genre = product->genre;
                decaySamplingInfo->productID = product->productID;
                decaySamplingInfo->mu = 0;
                decaySamplingInfo->Ep = 0;
                productDatas[n].genre = product->genre;
                productDatas[n].isVelocity = decaySamplingInfo->isVelocity;
                tpia_frame_setColumns( smr, &(productDatas[n].frame), 1, &labFrame );
                productDatas[n].productID = product->productID;
                productDatas[n].decayChannel = &(product->decayChannel);
                if( strcmp( product->genre, "twoBody_angular" ) == 0 ) {
                    secondTwoBody = 1;
                    productDatas[n+1].productID = nextProduct->productID;
                    productDatas[n].genre = product->genre;
                    tpia_angular_SampleMu( smr, &(product->angular), decaySamplingInfo );   /* Need to test for success. */
                    if( smr_isOk( smr ) ) {
                        phi = 2. * M_PI * tpia_misc_drng( decaySamplingInfo->rng, decaySamplingInfo->rngState );
                        productDatas[n].isVelocity = decaySamplingInfo->isVelocity;
                        productDatas[n].frame = decaySamplingInfo->frame;
                        tpia_kinetics_2BodyReaction( smr, decayChannel, decaySamplingInfo->e_in, decaySamplingInfo->mu, phi, &(productDatas[n]) );
                    } }
                else if( strcmp( product->genre, "NBody_Legendre" ) == 0 ) {
                    tpia_Legendre_SampleEp( smr, &(product->Legendre), 1, decaySamplingInfo ); }
                else if( strcmp( product->genre, "NBody_angular_energy" ) == 0 ) {
                    tpia_angular_SampleMu( smr, &(product->angular), decaySamplingInfo );   /* Need to test for success. */
                    tpia_angularEnergy_SampleEp( smr, &(product->angularEnergy), decaySamplingInfo ); }
                else if( strcmp( product->genre, "NBody_uncorrelate_Legendre" ) == 0 ) {
                    tpia_angular_SampleMu( smr, &(product->angular), decaySamplingInfo );   /* Need to test for success. */
                    tpia_Legendre_SampleEp( smr, &(product->Legendre), 0, decaySamplingInfo ); }
                else if( strcmp( product->genre, "unknown" ) == 0 ) {
                    }
                else {
                    printf( "Unknown spectral data form product name = %s, genre = %s\n", product->productID->name, product->genre );
                }
                if( !smr_isOk( smr ) ) return( -1 );
                if( secondTwoBody ) {
                    n++;
                    productDatas[n].productID = nextProduct->productID;
                    productDatas[n].genre = nextProduct->genre; }
                else {
                    productDatas[n].kineticEnergy = decaySamplingInfo->Ep;
                    p = std::sqrt( decaySamplingInfo->Ep * ( decaySamplingInfo->Ep + 2. * product->productID->fullMass_MeV ) );
                    productDatas[n].pz_vz = p * decaySamplingInfo->mu;
                    p = std::sqrt( 1. - decaySamplingInfo->mu * decaySamplingInfo->mu ) * p;
                    phi = 2. * M_PI * tpia_misc_drng( decaySamplingInfo->rng, decaySamplingInfo->rngState );
                    productDatas[n].px_vx = p * std::sin( phi );
                    productDatas[n].py_vy = p * std::cos( phi );
                }
                n++;
            }
        }
    }
    return( n );
}

#if defined __cplusplus
}
#endif
