/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <string.h>
#include <cmath>

#include "MCGIDI_fromTOM.h"
#include "MCGIDI_misc.h"
#include "MCGIDI_private.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

/*
************************************************************
*/
MCGIDI_angular *MCGIDI_angular_new( statusMessageReporting *smr ) {

    MCGIDI_angular *angular;

    if( ( angular = (MCGIDI_angular *) smr_malloc2( smr, sizeof( MCGIDI_angular ), 0, "angular" ) ) == NULL ) return( NULL );
    if( MCGIDI_angular_initialize( smr, angular ) ) angular = MCGIDI_angular_free( smr, angular );
    return( angular );
}
/*
************************************************************
*/
int MCGIDI_angular_initialize( statusMessageReporting * /*smr*/, MCGIDI_angular *angular ) {

    memset( angular, 0, sizeof( MCGIDI_angular ) );
    return( 0 );
}
/*
************************************************************
*/
MCGIDI_angular *MCGIDI_angular_free( statusMessageReporting *smr, MCGIDI_angular *angular ) {

    MCGIDI_angular_release( smr, angular );
    smr_freeMemory( (void **) &angular );
    return( NULL );
}
/*
************************************************************
*/
int MCGIDI_angular_release( statusMessageReporting *smr, MCGIDI_angular *angular ) {


    MCGIDI_sampling_pdfsOfXGivenW_release( smr, &(angular->dists) );

    MCGIDI_angular_initialize( smr, angular );
    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_angular_setTwoBodyMasses( statusMessageReporting * /*smr*/, MCGIDI_angular *angular, double projectileMass_MeV, double targetMass_MeV, 
    double productMass_MeV, double residualMass_MeV ) {

    if( angular == NULL ) return( 0 );         /* ???????? This needs work. Happens when first product of a two-body reaction as no distribution. */
    angular->projectileMass_MeV = projectileMass_MeV;
    angular->targetMass_MeV = targetMass_MeV;
    angular->productMass_MeV = productMass_MeV;
    angular->residualMass_MeV = residualMass_MeV;
    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_angular_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_distribution *distribution, ptwXYPoints *norms ) {

    MCGIDI_angular *angular = NULL;
    xDataTOM_element *angularElement, *linearElement, *frameElement = NULL;
    char const *nativeData;
    ptwXYPoints *pdfXY = NULL;
    ptwXPoints *cdfX = NULL;
    ptwXY_interpolation interpolationXY, interpolationWY;

    if( ( angularElement = xDataTOME_getOneElementByName( smr, element, "angular", 1 ) ) == NULL ) goto err;
    if( ( angular = MCGIDI_angular_new( smr ) ) == NULL ) goto err;

    if( ( nativeData = xDataTOM_getAttributesValueInElement( angularElement, "nativeData" ) ) == NULL ) goto err;
    if( strcmp( nativeData, "isotropic" ) == 0 ) {
        if( ( frameElement = xDataTOME_getOneElementByName( smr, angularElement, "isotropic", 1 ) ) == NULL ) {
            smr_setReportError2( smr, smr_unknownID, 1, "angular type missing for nativeData = '%s'", nativeData );
            goto err;
        }
        angular->type = MCGIDI_angularType_isotropic; }
    else if( strcmp( nativeData, "recoil" ) == 0 ) {            /* BRB. Needs work to get referenced product data?????? */
        angular->type = MCGIDI_angularType_recoil; }
    else {
        int i, j, n;
        double norm, energyFactor;
        nfu_status status;
        xDataTOM_XYs *XYs;
        xDataTOM_W_XYs *W_XYs;
        ptwXYPoint *point;
        MCGIDI_pdfsOfXGivenW *dists = &(angular->dists);
        MCGIDI_pdfOfX *dist;
        char const *energyUnit, *multiplicityProbabilityUnits[2] = { "", "" };

        if( ( linearElement = xDataTOME_getOneElementByName( NULL, angularElement, "linear", 0 ) ) == NULL ) {
            if( ( linearElement = xDataTOME_getOneElementByName( smr, angularElement, "pointwise", 1 ) ) == NULL ) {
                smr_setReportError2( smr, smr_unknownID, 1, "unsupported angular type: nativeData = '%s'", nativeData );
                goto err;
            }
        }
        frameElement = linearElement;

        if( MCGIDI_fromTOM_interpolation( smr, linearElement, 0, &interpolationWY ) ) goto err;
        if( MCGIDI_fromTOM_interpolation( smr, linearElement, 1, &interpolationXY ) ) goto err;
        dists->interpolationWY = interpolationWY;
        dists->interpolationXY = interpolationXY;

        if( ( W_XYs = (xDataTOM_W_XYs *) xDataTOME_getXDataIfID( smr, linearElement, "W_XYs" ) ) == NULL ) goto err;
        if( ( dists->Ws = (double *) smr_malloc2( smr, W_XYs->length * sizeof( double ), 1, "dists->Ws" ) ) == NULL ) goto err;
        if( ( dists->dist = (MCGIDI_pdfOfX *) smr_malloc2( smr, W_XYs->length * sizeof( MCGIDI_pdfOfX ), 0, "dists->dist" ) ) == NULL ) goto err;

        energyUnit = xDataTOM_subAxes_getUnit( smr, &(W_XYs->subAxes), 0 );
        if( !smr_isOk( smr ) ) goto err;
        energyFactor = MCGIDI_misc_getUnitConversionFactor( smr, energyUnit, "MeV" );
        if( !smr_isOk( smr ) ) goto err;

        for( i = 0; i < W_XYs->length; i++ ) {
            XYs = &(W_XYs->XYs[i]);
            dist = &(dists->dist[i]);
            dists->Ws[i] = XYs->value * energyFactor;
            if( ( pdfXY = MCGIDI_misc_dataFromXYs2ptwXYPointsInUnitsOf( smr, XYs, interpolationXY, multiplicityProbabilityUnits ) ) == NULL ) goto err;
            if( ptwXY_simpleCoalescePoints( pdfXY ) != nfu_Okay ) goto err;
            dist->numberOfXs = n = (int) ptwXY_length( pdfXY );

            if( ( dist->Xs = (double *) smr_malloc2( smr, 3 * n * sizeof( double ), 0, "dist->Xs" ) ) == NULL ) goto err;
            dists->numberOfWs++;
            dist->pdf = &(dist->Xs[n]);
            dist->cdf = &(dist->pdf[n]);

            for( j = 0; j < n; j++ ) {
                point = ptwXY_getPointAtIndex_Unsafely( pdfXY, j );
                dist->Xs[j] = point->x;
                dist->pdf[j] = point->y;
            }

            if( ( cdfX = ptwXY_runningIntegral( pdfXY, &status ) ) == NULL ) {
                smr_setReportError2( smr, smr_unknownID, 1, "ptwXY_runningIntegral err = %d: %s\n", status, nfu_statusMessage( status ) );
                goto err;
            }
            norm = ptwX_getPointAtIndex_Unsafely( cdfX, n - 1 );
            if( norms != NULL ) {
                ptwXY_setValueAtX( norms, XYs->value, norm ); }
            else if( std::fabs( 1. - norm ) > 0.99 ) {
                smr_setReportError2( smr, smr_unknownID, 1, "bad norm = %e for angular.linear data", norm );
                goto err;
            }
            for( j = 0; j < n; j++ ) dist->cdf[j] = ptwX_getPointAtIndex_Unsafely( cdfX, j ) / norm;
            for( j = 0; j < n; j++ ) dist->pdf[j] /= norm;
            pdfXY = ptwXY_free( pdfXY );
            cdfX = ptwX_free( cdfX );
        }
        angular->type = MCGIDI_angularType_linear;
    }

    if( frameElement != NULL ) {
        if( ( angular->frame = MCGIDI_misc_getProductFrame( smr, frameElement ) ) == xDataTOM_frame_invalid ) goto err;
    }

    distribution->angular = angular;
    distribution->type = MCGIDI_distributionType_angular_e;

    return( 0 );

err:
    if( pdfXY != NULL ) ptwXY_free( pdfXY );
    if( cdfX != NULL ) cdfX = ptwX_free( cdfX );
    if ( angular != NULL ) MCGIDI_angular_free( smr, angular );
    return( 1 );
}
/*
************************************************************
*/
int MCGIDI_angular_sampleMu( statusMessageReporting *smr, MCGIDI_angular *angular, MCGIDI_quantitiesLookupModes &modes, 
        MCGIDI_decaySamplingInfo *decaySamplingInfo ) {

    double randomMu = decaySamplingInfo->rng( decaySamplingInfo->rngState );
    MCGIDI_pdfsOfXGivenW_sampled sampled;

    switch( angular->type ) {
    case MCGIDI_angularType_isotropic :
        decaySamplingInfo->frame = angular->frame;
        decaySamplingInfo->mu = 1. - 2. * decaySamplingInfo->rng( decaySamplingInfo->rngState );
        break;
    case MCGIDI_angularType_linear :
        decaySamplingInfo->frame = angular->frame;
        sampled.smr = smr;
        sampled.w = modes.getProjectileEnergy( );
        MCGIDI_sampling_sampleX_from_pdfsOfXGivenW( &(angular->dists), &sampled, randomMu );
        decaySamplingInfo->mu = sampled.x;
        break;
    case MCGIDI_angularType_recoil :
    default :
        smr_setReportError2( smr, smr_unknownID, 1, "angular type = %d not supported", angular->type );
    }
    return( !smr_isOk( smr ) );
}

#if defined __cplusplus
}
#endif

