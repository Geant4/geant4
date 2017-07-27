/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <string.h>
#include <cmath>

#include "MCGIDI.h"
#include "MCGIDI_misc.h"

#if defined __cplusplus
#include "G4Log.hh"
namespace GIDI {
using namespace GIDI;
#endif

/*
************************************************************
*/
MCGIDI_outputChannel *MCGIDI_outputChannel_new( statusMessageReporting *smr ) {

    MCGIDI_outputChannel *outputChannel;

    if( ( outputChannel = (MCGIDI_outputChannel *) smr_malloc2( smr, sizeof( MCGIDI_outputChannel ), 0, "outputChannel" ) ) == NULL ) return( NULL );
    if( MCGIDI_outputChannel_initialize( smr, outputChannel ) ) outputChannel = MCGIDI_outputChannel_free( smr, outputChannel );
    return( outputChannel );
}
/*
************************************************************
*/
int MCGIDI_outputChannel_initialize( statusMessageReporting * /*smr*/, MCGIDI_outputChannel *outputChannel ) {

    memset( outputChannel, 0, sizeof( MCGIDI_outputChannel ) );
    return( 0 );
}
/*
************************************************************
*/
MCGIDI_outputChannel *MCGIDI_outputChannel_free( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel ) {

    MCGIDI_outputChannel_release( smr, outputChannel );
    smr_freeMemory( (void **) &outputChannel );
    return( NULL );
}
/*
************************************************************
*/
int MCGIDI_outputChannel_release( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel ) {

    int i;

    for( i = 0; i < outputChannel->numberOfProducts; i++ ) MCGIDI_product_release( smr, &(outputChannel->products[i]) );
    smr_freeMemory( (void **) &(outputChannel->products) );
    MCGIDI_outputChannel_initialize( smr, outputChannel );

    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_outputChannel_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_POPs *pops, MCGIDI_outputChannel *outputChannel,
        MCGIDI_reaction *reaction, MCGIDI_product *parent ) {

    int n, delayedNeutronIndex = 0;
    char const *genre, *Q;
    xDataTOM_element *child;

    MCGIDI_outputChannel_initialize( smr, outputChannel );

    outputChannel->reaction = reaction;
    outputChannel->parent = parent;
    if( ( genre = xDataTOM_getAttributesValueInElement( element, "genre" ) ) == NULL ) goto err;
    if( ( parent != NULL ) && ( strcmp( genre, "NBody" ) ) ) {
        smr_setReportError2( smr, smr_unknownID, 1, "decay channel's genre can only be 'uncorreclated' (a.k.a. 'NBody') and not '%s'", genre );
        goto err;
    }
    if( strcmp( genre, "twoBody" ) == 0 ) {
        outputChannel->genre = MCGIDI_channelGenre_twoBody_e; }
    else if( strcmp( genre, "NBody" ) == 0 ) {
        outputChannel->genre = MCGIDI_channelGenre_uncorrelated_e; }
    else if( strcmp( genre, "sumOfRemainingOutputChannels" ) == 0 ) {
        outputChannel->genre = MCGIDI_channelGenre_sumOfRemaining_e; }
    else {
        smr_setReportError2( smr, smr_unknownID, 1, "unsupported genre = '%s'", genre );
        goto err;
    }
    if( ( Q = xDataTOM_getAttributesValueInElement( element, "Q" ) ) == NULL ) goto err;
    outputChannel->QIsFloat = !MCGIDI_misc_PQUStringToDoubleInUnitOf( smr, Q, "MeV", &(outputChannel->Q) );

    if( ( n = xDataTOM_numberOfElementsByName( smr, element, "product" ) ) == 0 ) {
        smr_setReportError2p( smr, smr_unknownID, 1, "outputChannel does not have any products" );
        goto err;
    }
    if( ( outputChannel->products = (MCGIDI_product *) smr_malloc2( smr, n * sizeof( MCGIDI_product ), 0, "outputChannel->products" ) ) == NULL ) goto err;

    for( child = xDataTOME_getFirstElement( element ); child != NULL; child = xDataTOME_getNextElement( child ) ) {
        if( strcmp( child->name, "product" ) == 0 ) {
            if( MCGIDI_product_parseFromTOM( smr, child, outputChannel, pops, &(outputChannel->products[outputChannel->numberOfProducts]),
                &delayedNeutronIndex ) ) goto err;
            outputChannel->numberOfProducts++; }
        else if( strcmp( child->name, "fissionEnergyReleased" ) == 0 ) {     /* ????????? Need to support. */
            continue; }
        else {
            printf( "outputChannel child not currently supported = %s\n", child->name );
        }
    }
    if( outputChannel->genre == MCGIDI_channelGenre_twoBody_e ) {
        double projectileMass_MeV, targetMass_MeV, productMass_MeV, residualMass_MeV;

        projectileMass_MeV = MCGIDI_reaction_getProjectileMass_MeV( smr, reaction );
        targetMass_MeV = MCGIDI_reaction_getTargetMass_MeV( smr, reaction );
        productMass_MeV = MCGIDI_product_getMass_MeV( smr, &(outputChannel->products[0]) );
        residualMass_MeV = MCGIDI_product_getMass_MeV( smr, &(outputChannel->products[1]) );
        MCGIDI_product_setTwoBodyMasses( smr, &(outputChannel->products[0]), projectileMass_MeV, targetMass_MeV, productMass_MeV, residualMass_MeV );
    }

    return( 0 );

err:
    MCGIDI_outputChannel_release( smr, outputChannel );
    return( 1 );
}
/*
************************************************************
*/
int MCGIDI_outputChannel_numberOfProducts( MCGIDI_outputChannel *outputChannel ) {

    return( outputChannel->numberOfProducts );
}
/*
************************************************************
*/
MCGIDI_product *MCGIDI_outputChannel_getProductAtIndex( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel, int i ) {

    if( ( i < 0 ) || ( i >= outputChannel->numberOfProducts ) ) {
        smr_setReportError2( smr, smr_unknownID, 1, "bad product index = %d: outputChannel as only %d products", i, outputChannel->numberOfProducts );
        return( NULL );
    }
    return( &(outputChannel->products[i]) );
}
/*
************************************************************
*/
int MCGIDI_outputChannel_getDomain( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel, double *EMin, double *EMax ) {

    if( outputChannel->reaction != NULL ) return( MCGIDI_reaction_getDomain( smr, outputChannel->reaction, EMin, EMax ) );
    return( MCGIDI_product_getDomain( smr, outputChannel->parent, EMin, EMax ) );
}
/*
************************************************************
*/
MCGIDI_target_heated *MCGIDI_outputChannel_getTargetHeated( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel ) {

    if( outputChannel->reaction != NULL ) return( MCGIDI_reaction_getTargetHeated( smr, outputChannel->reaction ) );
    return( MCGIDI_product_getTargetHeated( smr, outputChannel->parent ) );
}
/*
************************************************************
*/
double MCGIDI_outputChannel_getProjectileMass_MeV( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel ) {

    if( outputChannel->reaction != NULL ) return( MCGIDI_reaction_getProjectileMass_MeV( smr, outputChannel->reaction ) );
    return( MCGIDI_product_getProjectileMass_MeV( smr, outputChannel->parent ) );
}
/*
************************************************************
*/
double MCGIDI_outputChannel_getTargetMass_MeV( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel ) {

    if( outputChannel->reaction != NULL ) return( MCGIDI_reaction_getTargetMass_MeV( smr, outputChannel->reaction ) );
    return( MCGIDI_product_getTargetMass_MeV( smr, outputChannel->parent ) );
}
/*
************************************************************
*/
double MCGIDI_outputChannel_getQ_MeV( statusMessageReporting * /*smr*/, MCGIDI_outputChannel *outputChannel, double /*e_in*/ ) {

    return( outputChannel->Q );
}
/*
************************************************************
*/
double MCGIDI_outputChannel_getFinalQ( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel, double e_in ) {

    int iProduct;
    double Q = outputChannel->Q;
    MCGIDI_product *product;

    for( iProduct = 0; iProduct < outputChannel->numberOfProducts; iProduct++ ) {
        product = &(outputChannel->products[iProduct]);
        if( product->decayChannel.genre != MCGIDI_channelGenre_undefined_e ) Q += MCGIDI_outputChannel_getFinalQ( smr, &(product->decayChannel), e_in );
        if( !smr_isOk( smr ) ) break;
    }
    return( Q );
}
/*
************************************************************
*/
int MCGIDI_outputChannel_sampleProductsAtE( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel, MCGIDI_quantitiesLookupModes &modes,
        MCGIDI_decaySamplingInfo *decaySamplingInfo, MCGIDI_sampledProductsDatas *productDatas, double *masses_ ) {

    int i1, multiplicity, secondTwoBody = 0, isDecayChannel = ( outputChannel->reaction == NULL );
    double e_in = modes.getProjectileEnergy( );
    MCGIDI_product *product;
    double phi, p, masses[3];
    MCGIDI_distribution *distribution;
    MCGIDI_sampledProductsData productData[2];

    if( isDecayChannel ) {
        masses[0] = masses_[0];              /* More work may be needed here. */
        masses[1] = masses_[1]; }
    else {
        masses[0] = MCGIDI_reaction_getProjectileMass_MeV( smr, outputChannel->reaction );
        masses[1] = MCGIDI_reaction_getTargetMass_MeV( smr, outputChannel->reaction );
    }

    for( i1 = 0; i1 < outputChannel->numberOfProducts; i1++ ) {
        product = &(outputChannel->products[i1]);
        if( product->decayChannel.genre != MCGIDI_channelGenre_undefined_e ) {
            if( MCGIDI_outputChannel_sampleProductsAtE( smr, &(product->decayChannel), modes, decaySamplingInfo, productDatas, masses ) < 0 ) return( -1 ); }
        else {
            distribution = &(product->distribution);
            if( distribution->type == MCGIDI_distributionType_none_e ) continue;
            if( !secondTwoBody ) {
                if( ( multiplicity = product->multiplicity ) == 0 ) multiplicity = MCGIDI_product_sampleMultiplicity( smr, product, e_in,
                    decaySamplingInfo->rng( decaySamplingInfo->rngState ) );
                while( multiplicity > 0 ) { 

                    multiplicity--;
                    decaySamplingInfo->pop = product->pop;
                    decaySamplingInfo->mu = 0;
                    decaySamplingInfo->Ep = 0;
                    productData[0].isVelocity = decaySamplingInfo->isVelocity;
                    productData[0].pop = product->pop;
                    productData[0].delayedNeutronIndex = product->delayedNeutronIndex;
                    productData[0].delayedNeutronRate = product->delayedNeutronRate;
                    productData[0].birthTimeSec = 0;
                    if( product->delayedNeutronRate > 0 ) {
                        productData[0].birthTimeSec = -G4Log( decaySamplingInfo->rng( decaySamplingInfo->rngState ) ) / product->delayedNeutronRate;
                    }

                    switch( outputChannel->genre ) {
                    case MCGIDI_channelGenre_twoBody_e :
                        secondTwoBody = 1;
                        MCGIDI_angular_sampleMu( smr, distribution->angular, modes, decaySamplingInfo );
                        if( smr_isOk( smr ) ) {
                            phi = 2. * M_PI * decaySamplingInfo->rng( decaySamplingInfo->rngState );
                            MCGIDI_kinetics_2BodyReaction( smr, distribution->angular, e_in, decaySamplingInfo->mu, phi, productData );
                            if( !smr_isOk( smr ) ) return( -1 );
                            productData[1].pop = product[1].pop;
                            productData[1].delayedNeutronIndex = product[1].delayedNeutronIndex;
                            productData[1].delayedNeutronRate = product->delayedNeutronRate;
                            productData[1].birthTimeSec = 0;
                            MCGIDI_sampledProducts_addProduct( smr, productDatas, productData );
                            if( !smr_isOk( smr ) ) return( -1 );
                            MCGIDI_sampledProducts_addProduct( smr, productDatas, &(productData[1]) );
                            if( !smr_isOk( smr ) ) return( -1 );
                        }
                        break;
                    case MCGIDI_channelGenre_uncorrelated_e :
                    case MCGIDI_channelGenre_sumOfRemaining_e :
                        masses[2] = MCGIDI_product_getMass_MeV( smr, product );
                        switch( distribution->type ) {
                        case MCGIDI_distributionType_uncorrelated_e :
                            MCGIDI_uncorrelated_sampleDistribution( smr, distribution, modes, decaySamplingInfo );
                            break;
                        case MCGIDI_distributionType_energyAngular_e :
                            MCGIDI_energyAngular_sampleDistribution( smr, distribution, modes, decaySamplingInfo );
                            break;
                        case MCGIDI_distributionType_KalbachMann_e :
                            MCGIDI_KalbachMann_sampleEp( smr, distribution->KalbachMann, modes, decaySamplingInfo );
                            break;
                        case MCGIDI_distributionType_angularEnergy_e :
                            MCGIDI_angularEnergy_sampleDistribution( smr, distribution->angularEnergy, modes, decaySamplingInfo );
                            break;
                        default :
                            printf( "Unknown spectral data form product name = %s, channel genre = %d\n", product->pop->name, outputChannel->genre );
                            break;
                        }
                        break;
                    case MCGIDI_channelGenre_undefined_e :
                        printf( "Channel is undefined\n" );
			break;
                    case MCGIDI_channelGenre_twoBodyDecay_e :
                        printf( "Channel is twoBodyDecay\n" );
			break;
                    case MCGIDI_channelGenre_uncorrelatedDecay_e :
                        printf( "Channel is uncorrelatedDecay\n" );
			break;
                    default :
                        printf( "Unsupported channel genre = %d\n", outputChannel->genre );
			break;
                    }
                    if( !smr_isOk( smr ) ) return( -1 );
                    if( !secondTwoBody ) {
                        if( decaySamplingInfo->frame == xDataTOM_frame_centerOfMass ) {
                            if( MCGIDI_kinetics_COM2Lab( smr, modes, decaySamplingInfo, masses ) != 0 ) return( -1 );
                        }
                        productData[0].kineticEnergy = decaySamplingInfo->Ep;
                        p = std::sqrt( decaySamplingInfo->Ep * ( decaySamplingInfo->Ep + 2. * product->pop->mass_MeV ) );
                        if( productData[0].isVelocity ) p *= MCGIDI_speedOfLight_cm_sec / std::sqrt( p * p + product->pop->mass_MeV * product->pop->mass_MeV );
                        productData[0].pz_vz = p * decaySamplingInfo->mu;
                        p = std::sqrt( 1. - decaySamplingInfo->mu * decaySamplingInfo->mu ) * p;
                        phi = 2. * M_PI * decaySamplingInfo->rng( decaySamplingInfo->rngState );
                        productData[0].px_vx = p * std::sin( phi );
                        productData[0].py_vy = p * std::cos( phi );
                        MCGIDI_sampledProducts_addProduct( smr, productDatas, productData );
                        if( !smr_isOk( smr ) ) return( -1 );
                    }
                } // Loop checking, 11.06.2015, T. Koi
            }
        }
    }
    return( productDatas->numberOfProducts );
}

#if defined __cplusplus
}
#endif

