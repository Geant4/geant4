/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <string.h>
#include <cmath>

#include "MCGIDI.h"
#include "MCGIDI_misc.h"
#include "MCGIDI_fromTOM.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

typedef struct polynomialCallbackArgs_s {
    int length;
    double energyFactor;
    double *coefficients;
} polynomialCallbackArgs;

static int MCGIDI_product_parsePiecewiseMultiplicity( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_product *product );
static ptwXYPoints *MCGIDI_product_parsePolynomialMultiplicity( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_product *product );
static int MCGIDI_product_parseWeightedReferenceMultiplicityFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_product *product, 
    ptwXYPoints **multiplicityVsEnergy, ptwXYPoints **norms );
static double MCGIDI_product_evaluatePolynomial( double x, polynomialCallbackArgs *args );
/*
************************************************************
*/
MCGIDI_product *MCGIDI_product_new( statusMessageReporting *smr ) {

    MCGIDI_product *product;

    if( ( product = (MCGIDI_product *) smr_malloc2( smr, sizeof( MCGIDI_product ), 0, "product" ) ) == NULL ) return( NULL );
    if( MCGIDI_product_initialize( smr, product ) ) product = MCGIDI_product_free( smr, product );
    return( product );
}
/*
************************************************************
*/
int MCGIDI_product_initialize( statusMessageReporting * /*smr*/, MCGIDI_product *product ) {

    memset( product, 0, sizeof( MCGIDI_product ) );
    product->delayedNeutronIndex = -1;
    return( 0 );
}
/*
************************************************************
*/
MCGIDI_product *MCGIDI_product_free( statusMessageReporting *smr, MCGIDI_product *product ) {

    MCGIDI_product_release( smr, product );
    smr_freeMemory( (void **) &product );
    return( NULL );
}
/*
************************************************************
*/
int MCGIDI_product_release( statusMessageReporting *smr, MCGIDI_product *product ) {

    int i;

    if( product->label != NULL ) smr_freeMemory( (void **) &(product->label) );

    if( product->multiplicityVsEnergy != NULL ) ptwXY_free( product->multiplicityVsEnergy );
    if( product->piecewiseMultiplicities != NULL ) {
        for( i = 0; i < product->numberOfPiecewiseMultiplicities; i++ ) ptwXY_free( product->piecewiseMultiplicities[i] );
        smr_freeMemory( (void **) &(product->piecewiseMultiplicities) );
    }
    if( product->norms != NULL ) ptwXY_free( product->norms );

    MCGIDI_distribution_release( smr, &(product->distribution) );
    MCGIDI_outputChannel_release( smr, &(product->decayChannel) );

    MCGIDI_product_initialize( smr, product );
    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_product_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_outputChannel *outputChannel,
        MCGIDI_POPs *pops, MCGIDI_product *product, int *delayedNeutronIndex ) {

    char const *name, *label, *delayedNeutron, *multiplicityStr, *multiplicityUnits[2] = { "MeV", "" };
    xDataTOM_element *multiplicity, *multiplicityTOM, *decayChannelElement;
    nfu_status status;
    ptwXYPoints *multiplicityVsEnergy = NULL, *norms1 = NULL, *norms2 = NULL;

    MCGIDI_product_initialize( smr, product );

    product->outputChannel = outputChannel;
    if( ( name = xDataTOM_getAttributesValueInElement( element, "name" ) ) == NULL ) goto err;
    if( ( product->pop = MCGIDI_POPs_findParticle( pops, name ) ) == NULL ) {
        smr_setReportError2( smr, smr_unknownID, 1, "product '%s' not found in pops", name );
        goto err;
    }
    if( ( label = xDataTOM_getAttributesValueInElement( element, "label" ) ) != NULL ) {
        if( ( product->label = smr_allocateCopyString2( smr, label, "product->label" ) ) == NULL ) goto err;
    }

    if( ( delayedNeutron = xDataTOM_getAttributesValueInElement( element, "emissionMode" ) ) != NULL ) {
        if( strcmp( delayedNeutron, "delayed" ) == 0 ) {
            if( ( delayedNeutron = xDataTOM_getAttributesValueInElement( element, "decayRate" ) ) == NULL ) {
                goto err;
            }
            if( MCGIDI_misc_PQUStringToDoubleInUnitOf( smr, delayedNeutron, "1/s", &(product->delayedNeutronRate) ) != 0 ) goto err;
            product->delayedNeutronIndex = *delayedNeutronIndex;
            (*delayedNeutronIndex)++;
        }
    }

    if( ( multiplicityStr = xDataTOM_getAttributesValueInElement( element, "multiplicity" ) ) == NULL ) goto err;
    if( xDataTOME_convertAttributeToInteger( NULL, element, "multiplicity", &(product->multiplicity) ) ) {
        if( strcmp( multiplicityStr, "energyDependent" ) ) {
            smr_setReportError2( smr, smr_unknownID, 1, "invalid multiplicity '%s' for product '%s'", multiplicityStr, name );
            goto err;
        }
        if( ( multiplicity = xDataTOME_getOneElementByName( smr, element, "multiplicity", 1 ) ) == NULL ) goto err;
        if( ( multiplicityTOM = xDataTOME_getOneElementByName( NULL, multiplicity, "weightedReference", 0 ) ) != NULL ) {
            if( MCGIDI_product_parseWeightedReferenceMultiplicityFromTOM( smr, multiplicityTOM, product, &multiplicityVsEnergy, &norms1 ) ) goto err; }
        else if( ( multiplicityTOM = xDataTOME_getOneElementByName( NULL, multiplicity, "piecewise", 0 ) ) != NULL ) {
            if( MCGIDI_product_parsePiecewiseMultiplicity( smr, multiplicityTOM, product ) ) goto err; }
        else if( ( multiplicityTOM = xDataTOME_getOneElementByName( NULL, multiplicity, "polynomial", 0 ) ) != NULL ) {
            if( ( multiplicityVsEnergy = MCGIDI_product_parsePolynomialMultiplicity( smr, multiplicityTOM, product ) ) == NULL ) goto err; }
        else {
/* ??????? Need to check interpolation. */
            if( ( multiplicityTOM = xDataTOME_getOneElementByName( smr, multiplicity, "pointwise", 1 ) ) == NULL ) goto err;
            if( ( multiplicityVsEnergy = MCGIDI_misc_dataFromElement2ptwXYPointsInUnitsOf( smr, multiplicityTOM, multiplicityUnits ) ) == NULL ) goto err;
        }
    }

    if( strcmp( product->pop->name, "gamma" ) == 0 ) {
        if( ( norms2 = ptwXY_new( ptwXY_interpolationLinLin, NULL, 2., 1e-3, 200, 10, &status, 0 ) ) == NULL ) {
            smr_setReportError2( smr, smr_unknownID, 1, "ptwXY_new err = %d: %s\n", status, nfu_statusMessage( status ) );
            goto err;
        }
    }
    if( MCGIDI_distribution_parseFromTOM( smr, element, product, pops, norms2 ) ) goto err;
    if( norms2 != NULL ) {
        if( ptwXY_length( norms2 ) < 2 ) {
            norms2 = ptwXY_free( norms2 ); }
        else {
            if( ptwXY_simpleCoalescePoints( norms2 ) != nfu_Okay ) goto err;
            if( ( ptwXY_getYMin( norms2 ) > 0.99 ) && ( ptwXY_getYMax( norms2 ) < 1.01 ) ) norms2 = ptwXY_free( norms2 );
        }
    }
    if( ( norms1 != NULL ) && ( norms2 != NULL ) ) {
        smr_setReportError2p( smr, smr_unknownID, 1, "norm1 and norm2 are both not NULL" );
        goto err;
    }

    product->multiplicityVsEnergy = multiplicityVsEnergy;
    product->norms = norms1;
    if( norms2 != NULL ) product->norms = norms2;

    if( ( decayChannelElement = xDataTOME_getOneElementByName( NULL, element, "decayChannel", 0 ) ) != NULL ) {
        if( MCGIDI_outputChannel_parseFromTOM( smr, decayChannelElement, pops, &(product->decayChannel), NULL, product ) ) goto err;
    }

    return( 0 );

err:
    if( multiplicityVsEnergy != NULL ) ptwXY_free( multiplicityVsEnergy );
    if( norms1 != NULL ) ptwXY_free( norms1 );
    if( norms2 != NULL ) ptwXY_free( norms2 );
    MCGIDI_product_release( smr, product );
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_product_parsePiecewiseMultiplicity( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_product *product ) {

    int i;
    xDataTOM_XYs *XYs;
    xDataTOM_regionsXYs *regionsXYs = (xDataTOM_regionsXYs *) element->xDataInfo.data;
    ptwXYPoints *multiplicityVsEnergy;
    char const *multiplicityUnits[2] = { "MeV", "" };

    if( ( product->piecewiseMultiplicities = (ptwXYPoints **) smr_malloc2( smr, regionsXYs->length * sizeof( ptwXYPoints * ), 1, "piecewiseMultiplicities" ) ) == NULL ) return( 1 );

    for( i = 0; i < regionsXYs->length; i++ ) {
/* ??????? Need to check interpolation. */
        XYs = &(regionsXYs->XYs[i]);
        if( ( multiplicityVsEnergy = MCGIDI_misc_dataFromXYs2ptwXYPointsInUnitsOf( smr, XYs, ptwXY_interpolationLinLin, multiplicityUnits ) 
            ) == NULL ) return( 1 );
        product->piecewiseMultiplicities[i] = multiplicityVsEnergy;
        product->numberOfPiecewiseMultiplicities++;
    }

    return( 0 );
}
/*
************************************************************
*/
static ptwXYPoints *MCGIDI_product_parsePolynomialMultiplicity( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_product *product ) {

    int length;
    double *coefficients;
    char const *energyUnit;
    ptwXYPoints *ptwXY = NULL;
    nfu_status status;
    double EMin, EMax;
    polynomialCallbackArgs args;

    if( MCGIDI_product_getDomain( smr, product, &EMin, &EMax ) ) goto err;

    length = xDataTOM_polynomial_getDataFromXDataInfo( (xDataTOM_xDataInfo *) &(element->xDataInfo), &coefficients );
    if( ( ptwXY = ptwXY_new( ptwXY_interpolationLinLin, NULL, 2., 1e-3, length, 10, &status, 0 ) ) == NULL ) {
        smr_setReportError2( smr, smr_unknownID, 1, "ptwXY_new err = %d: %s\n", status, nfu_statusMessage( status ) );
        goto err;
    }

    if( ( energyUnit = xDataTOM_axes_getUnit( smr, &(element->xDataInfo.axes), 0 ) ) == NULL ) goto err;
    args.energyFactor = MCGIDI_misc_getUnitConversionFactor( smr, energyUnit, "MeV" );
    if( !smr_isOk( smr ) ) goto err;

    args.length = length;
    args.coefficients = coefficients;
    ptwXY_setValueAtX( ptwXY, EMin, MCGIDI_product_evaluatePolynomial( EMin, &args ) );
    ptwXY_setValueAtX( ptwXY, EMax, MCGIDI_product_evaluatePolynomial( EMax, &args ) );
    if( length > 2 ) {          /* ?????????????? This needs work. */
        int i, n = 4 * length;
        double E = EMin, dE = ( EMax - EMin ) / n;

        for( i = 1; i < n; i++ ) {
            E += dE;
            ptwXY_setValueAtX( ptwXY, E, MCGIDI_product_evaluatePolynomial( E, &args ) );
        }
    }
    return( ptwXY );

err:
    if( ptwXY != NULL ) ptwXY_free( ptwXY );
    return( NULL );
}
/*
************************************************************
*/
static double MCGIDI_product_evaluatePolynomial( double x, polynomialCallbackArgs *args ) {

    int i;
    double value = 0.;

    x /= args->energyFactor;
    for( i = args->length; i > 0; i-- ) value = value * x + args->coefficients[i-1];

    return( value );
}
/*
************************************************************
*/
static int MCGIDI_product_parseWeightedReferenceMultiplicityFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_product * /*product*/, 
        ptwXYPoints **multiplicityVsEnergy, ptwXYPoints **norms ) {

    char const *link, *energyInMWUnits[2] = { "MeV", "" };
    xDataTOM_element *reference, *productTOM, *multiplicity, *weights, *pointwise;

    if( ( reference = xDataTOME_getOneElementByName( smr, element, "reference", 1 ) ) == NULL ) goto err;
    if( ( link = xDataTOM_getAttributesValueInElement( reference, "xlink:href" ) ) == NULL ) goto err;
    if( ( productTOM = xDataTOM_getLinksElement( smr, reference, link ) ) == NULL ) goto err;
    if( ( multiplicity = xDataTOME_getOneElementByName( smr, productTOM, "multiplicity", 1 ) ) == NULL ) goto err;
                        /* Currently, only pointwise supported. */
    if( ( pointwise = xDataTOME_getOneElementByName( smr, multiplicity, "pointwise", 1 ) ) == NULL ) goto err;
    if( ( *multiplicityVsEnergy = MCGIDI_misc_dataFromElement2ptwXYPointsInUnitsOf( smr, pointwise, energyInMWUnits ) ) == NULL ) goto err;

    if( ( weights = xDataTOME_getOneElementByName( smr, element, "weights", 1 ) ) == NULL ) goto err;
    if( ( pointwise = xDataTOME_getOneElementByName( smr, weights, "pointwise", 1 ) ) == NULL ) goto err;
    if( ( *norms = MCGIDI_misc_dataFromElement2ptwXYPointsInUnitsOf( smr, pointwise, energyInMWUnits ) ) == NULL ) goto err;

    return( 0 );

err:
    if( *multiplicityVsEnergy != NULL ) *multiplicityVsEnergy = ptwXY_free( *multiplicityVsEnergy );
    if( *norms != NULL ) *norms = ptwXY_free( *norms );
    return( 1 );
}
/*
************************************************************
*/
int MCGIDI_product_getDomain( statusMessageReporting *smr, MCGIDI_product *product, double *EMin, double *EMax ) {

    return( MCGIDI_outputChannel_getDomain( smr, product->outputChannel, EMin, EMax ) );
}
/*
************************************************************
*/
int MCGIDI_product_setTwoBodyMasses( statusMessageReporting *smr, MCGIDI_product *product, double projectileMass_MeV, double targetMass_MeV,
    double productMass_MeV, double residualMass_MeV ) {

    return( MCGIDI_angular_setTwoBodyMasses( smr, product->distribution.angular, projectileMass_MeV, targetMass_MeV, productMass_MeV, residualMass_MeV ) );
}
/*
************************************************************
*/
double MCGIDI_product_getMass_MeV( statusMessageReporting * /*smr*/, MCGIDI_product *product ) {

    return( MCGIDI_POP_getMass_MeV( product->pop ) );
}
/*
************************************************************
*/
MCGIDI_target_heated *MCGIDI_product_getTargetHeated( statusMessageReporting *smr, MCGIDI_product *product ) {

    return( MCGIDI_outputChannel_getTargetHeated( smr, product->outputChannel ) );
}
/*
************************************************************
*/
double MCGIDI_product_getProjectileMass_MeV( statusMessageReporting *smr, MCGIDI_product *product ) {

    return( MCGIDI_outputChannel_getProjectileMass_MeV( smr, product->outputChannel ) );
}
/*
************************************************************
*/
double MCGIDI_product_getTargetMass_MeV( statusMessageReporting *smr, MCGIDI_product *product ) {

    return( MCGIDI_outputChannel_getTargetMass_MeV( smr, product->outputChannel ) );
}
/*
************************************************************
*/
int MCGIDI_product_sampleMultiplicity( statusMessageReporting * /*smr*/, MCGIDI_product *product, double e_in, double r ) {

    int i, multiplicity;
    double y, norm = 1.0;
    ptwXYPoints *ptwXY = product->multiplicityVsEnergy;

    if( product->piecewiseMultiplicities != NULL ) {
        for( i = 0; i < product->numberOfPiecewiseMultiplicities - 1; i++ ) {
            if( e_in < ptwXY_getXMax( product->piecewiseMultiplicities[i] ) ) break;
        }
        ptwXY = product->piecewiseMultiplicities[i];
    }
    y = MCGIDI_sampling_ptwXY_getValueAtX( ptwXY, e_in );
    if( product->norms != NULL ) norm = MCGIDI_sampling_ptwXY_getValueAtX( product->norms, e_in );
    y *= norm;
    multiplicity = (int) y;
    if( r < ( y - multiplicity ) ) multiplicity++;

    return( multiplicity );
}
/*
************************************************************
*/
int MCGIDI_product_sampleMu( statusMessageReporting *smr, MCGIDI_product *product, MCGIDI_quantitiesLookupModes &modes, 
        MCGIDI_decaySamplingInfo *decaySamplingInfo ) {

    if( product->distribution.type != MCGIDI_distributionType_angular_e ) {
        smr_setReportError2( smr, smr_unknownID, 1, "product distribution is not angular: type = %d", product->distribution.type );
        return( 1 );
    }
    return( MCGIDI_angular_sampleMu( smr, product->distribution.angular, modes, decaySamplingInfo ) );
}


/*
************************************************************
*/
int MCGIDI_sampledProducts_initialize( statusMessageReporting *smr, MCGIDI_sampledProductsDatas *sampledProductsDatas, int incrementSize ) {

    if( incrementSize < 10 ) incrementSize = 10;
    sampledProductsDatas->numberOfProducts = 0;
    sampledProductsDatas->numberAllocated = 0;
    sampledProductsDatas->incrementSize = incrementSize;
    sampledProductsDatas->products = NULL;
    return( MCGIDI_sampledProducts_remalloc( smr, sampledProductsDatas ) );
}
/*
************************************************************
*/
int MCGIDI_sampledProducts_release( statusMessageReporting * /*smr*/, MCGIDI_sampledProductsDatas *sampledProductsDatas ) {

    smr_freeMemory( (void **) &(sampledProductsDatas->products) );
    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_sampledProducts_remalloc( statusMessageReporting *smr, MCGIDI_sampledProductsDatas *sampledProductsDatas ) {

    int size = sampledProductsDatas->numberAllocated + sampledProductsDatas->incrementSize;

    if( ( sampledProductsDatas->products = (MCGIDI_sampledProductsData *) smr_realloc2( smr, sampledProductsDatas->products, 
        size * sizeof( MCGIDI_sampledProductsData ), "products" ) ) != NULL ) {
        sampledProductsDatas->numberAllocated = size;
        return( 0 );
    }
    sampledProductsDatas->numberOfProducts = 0;
    sampledProductsDatas->numberAllocated = 0;
    return( 1 );
}
/*
************************************************************
*/
int MCGIDI_sampledProducts_addProduct( statusMessageReporting *smr, MCGIDI_sampledProductsDatas *sampledProductsDatas, MCGIDI_sampledProductsData *sampledProductsData ) {

    if( sampledProductsDatas->numberOfProducts == sampledProductsDatas->numberAllocated ) {
        if( ( MCGIDI_sampledProducts_remalloc( smr, sampledProductsDatas ) ) != 0 ) return( 1 );
    }
    sampledProductsDatas->products[sampledProductsDatas->numberOfProducts] = *sampledProductsData;
    sampledProductsDatas->numberOfProducts++;
    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_sampledProducts_number( MCGIDI_sampledProductsDatas *sampledProductsDatas ) {

    return( sampledProductsDatas->numberOfProducts );
}
/*
************************************************************
*/
MCGIDI_sampledProductsData *MCGIDI_sampledProducts_getProductAtIndex( MCGIDI_sampledProductsDatas *sampledProductsDatas, int index ) {

    if( index < 0 ) return( NULL );
    if( index >= sampledProductsDatas->numberOfProducts ) return( NULL );
    return( &(sampledProductsDatas->products[index]) );
}

#if defined __cplusplus
}
#endif

