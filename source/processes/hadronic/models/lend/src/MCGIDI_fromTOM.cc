/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
/*
#include <unistd.h>
#include <ctype.h>
*/

#include "MCGIDI_fromTOM.h"
#include "MCGIDI_misc.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

static int MCGIDI_fromTOM_pdfOfXGivenW( statusMessageReporting *smr, ptwXYPoints *pdfXY, MCGIDI_pdfsOfXGivenW *dists, int i, double *norm );
/*
************************************************************
*/
int MCGIDI_fromTOM_pdfsOfXGivenW( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_pdfsOfXGivenW *dists, ptwXYPoints *norms,
    char const *toUnits[3] ) {

    int i;
    double norm, wUnitFactor;
    char const *wFromUnit, *toUnitsXY[2] = { toUnits[1], toUnits[2] };
    xDataTOM_XYs *XYs;
    xDataTOM_W_XYs *W_XYs; 
    ptwXYPoints *pdfXY = NULL;
    ptwXY_interpolation interpolationXY, interpolationWY;

    wFromUnit = xDataTOM_axes_getUnit( smr, &(element->xDataInfo.axes), 0 );
    if( !smr_isOk( smr ) ) goto err;
    wUnitFactor = MCGIDI_misc_getUnitConversionFactor( smr, wFromUnit, toUnits[0] );
    if( !smr_isOk( smr ) ) goto err;

    if( MCGIDI_fromTOM_interpolation( smr, element, 0, &interpolationWY ) ) goto err;
    if( MCGIDI_fromTOM_interpolation( smr, element, 1, &interpolationXY ) ) goto err;
    dists->interpolationWY = interpolationWY;
    dists->interpolationXY = interpolationXY;
    if( norms != NULL ) {
        if( interpolationWY == ptwXY_interpolationOther ) {
            smr_setReportError2p( smr, smr_unknownID, 1, "interpolationWY ptwXY_interpolationOther not supported" );
            goto err;
        }
    }

    W_XYs = (xDataTOM_W_XYs *) xDataTOME_getXDataIfID( smr, element, "W_XYs" );
    if( ( dists->Ws = (double *) smr_malloc2( smr, W_XYs->length * sizeof( double ), 1, "dists->Ws" ) ) == NULL ) goto err;
    if( ( dists->dist = (MCGIDI_pdfOfX *) smr_malloc2( smr, W_XYs->length * sizeof( MCGIDI_pdfOfX ), 0, "dists->dist" ) ) == NULL ) goto err;

    for( i = 0; i < W_XYs->length; i++ ) { 
        XYs = &(W_XYs->XYs[i]);
        dists->Ws[i] = wUnitFactor * XYs->value;
        if( ( pdfXY =  MCGIDI_misc_dataFromXYs2ptwXYPointsInUnitsOf( smr, XYs, interpolationXY, toUnitsXY ) ) == NULL ) goto err;
        if( MCGIDI_fromTOM_pdfOfXGivenW( smr, pdfXY, dists, i, &norm ) ) goto err;
        if( norms != NULL ) {
            ptwXY_setValueAtX( norms, XYs->value, norm ); }
        else if( std::fabs( 1. - norm ) > 0.99 ) {
            smr_setReportError2( smr, smr_unknownID, 1, "bad norm = %e for data", norm );
            goto err;
        }
        ptwXY_free( pdfXY );
        pdfXY = NULL;
    }

    return( 0 );

err:
    if( pdfXY != NULL ) ptwXY_free( pdfXY );
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_fromTOM_pdfOfXGivenW( statusMessageReporting *smr, ptwXYPoints *pdfXY, MCGIDI_pdfsOfXGivenW *dists, int i, double *norm ) {

    if( MCGIDI_fromTOM_pdfOfX( smr, pdfXY, &(dists->dist[i]), norm ) ) return( 1 );
    dists->numberOfWs++;
    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_fromTOM_pdfOfX( statusMessageReporting *smr, ptwXYPoints *pdfXY, MCGIDI_pdfOfX *dist, double *norm ) {

    int j1, n1 = (int) ptwXY_length( pdfXY );
    nfu_status status; 
    ptwXPoints *cdfX = NULL;
    ptwXYPoint *point; 

    dist->numberOfXs = 0;
    dist->Xs = NULL;
    if( ptwXY_simpleCoalescePoints( pdfXY ) != nfu_Okay ) goto err;

    if( ( dist->Xs = (double *) smr_malloc2( smr, 3 * n1 * sizeof( double ), 0, "dist->Xs" ) ) == NULL ) goto err;
    dist->pdf = &(dist->Xs[n1]);
    dist->cdf = &(dist->pdf[n1]);

    for( j1 = 0; j1 < n1; j1++ ) { 
            point = ptwXY_getPointAtIndex_Unsafely( pdfXY, j1 );
            dist->Xs[j1] = point->x;
            dist->pdf[j1] = point->y;
    }

    if( ( cdfX = ptwXY_runningIntegral( pdfXY, &status ) ) == NULL ) {
            smr_setReportError2( smr, smr_unknownID, 1, "ptwXY_runningIntegral err = %d: %s\n", status, nfu_statusMessage( status ) );
            goto err;
    }
    *norm = ptwX_getPointAtIndex_Unsafely( cdfX, n1 - 1 );
    if( *norm == 0. ) {             /* Should only happend for gammas. */
        double inv_norm, sum = 0;

        inv_norm = 1.0 / ( dist->Xs[n1-1] - dist->Xs[0] );
        for( j1 = 0; j1 < n1; ++j1 ) {
            if( j1 > 0 ) sum += dist->Xs[j1] - dist->Xs[j1-1];
            dist->pdf[j1] = 1;
            dist->cdf[j1] = sum * inv_norm;
        }
        dist->cdf[n1-1] = 1.; }
    else {
        for( j1 = 0; j1 < n1; j1++ ) dist->cdf[j1] = ptwX_getPointAtIndex_Unsafely( cdfX, j1 ) / *norm; 
        for( j1 = 0; j1 < n1; j1++ ) dist->pdf[j1] /= *norm;
    }
    ptwX_free( cdfX ); 

    dist->numberOfXs = n1;
    return( 0 );

err:
    if( dist->Xs != NULL ) smr_freeMemory( (void **) &(dist->Xs) );
    if( cdfX != NULL ) ptwX_free( cdfX ); 
    return( 1 );
}
/*
************************************************************
*/
int MCGIDI_fromTOM_interpolation( statusMessageReporting *smr, xDataTOM_element *element, int index, ptwXY_interpolation *interpolation ) {

    enum xDataTOM_interpolationFlag independent, dependent;
    enum xDataTOM_interpolationQualifier qualifier;

    if( xDataTOME_getInterpolation( smr, element, index, &independent, &dependent, &qualifier ) ) return( 1 );

    *interpolation = ptwXY_interpolationOther;

    if( dependent == xDataTOM_interpolationFlag_flat ) {
        *interpolation = ptwXY_interpolationFlat; }
    else if( independent == xDataTOM_interpolationFlag_linear ) {
        if( dependent == xDataTOM_interpolationFlag_linear ) {
            *interpolation = ptwXY_interpolationLinLin; }
        else if( dependent == xDataTOM_interpolationFlag_log ) {
            *interpolation = ptwXY_interpolationLinLog;
        } }
    else if( independent == xDataTOM_interpolationFlag_log ) {
        if( dependent == xDataTOM_interpolationFlag_linear ) {
            *interpolation = ptwXY_interpolationLogLin; }
        else if( dependent == xDataTOM_interpolationFlag_log ) {
            *interpolation = ptwXY_interpolationLogLog;
        }
    }

    return( 0 );
}

#if defined __cplusplus
}
#endif
