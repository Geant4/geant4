/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <stdlib.h>
#include <cmath>

#include "MCGIDI.h"

#if defined __cplusplus
#include "G4Log.hh"
#include "G4Pow.hh"
namespace GIDI {
using namespace GIDI;
#endif

/*
************************************************************
*/
int MCGIDI_sampling_pdfsOfXGivenW_initialize( statusMessageReporting * /*smr*/, MCGIDI_pdfsOfXGivenW *dists ) {

    memset( dists, 0, sizeof( MCGIDI_pdfsOfXGivenW ) );
    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_sampling_pdfsOfXGivenW_release( statusMessageReporting *smr, MCGIDI_pdfsOfXGivenW *dists ) {

    int i;

    for( i = 0; i < dists->numberOfWs; i++ ) MCGIDI_sampling_pdfsOfX_release( smr, &(dists->dist[i]) );
    smr_freeMemory( (void **) &(dists->Ws) );
    smr_freeMemory( (void **) &(dists->dist) );
    MCGIDI_sampling_pdfsOfXGivenW_initialize( smr, dists );
    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_sampling_pdfsOfX_release( statusMessageReporting * /*smr*/, MCGIDI_pdfOfX *dist ) {

    smr_freeMemory( (void **) &(dist->Xs) );
    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_sampling_sampleX_from_pdfsOfXGivenW( MCGIDI_pdfsOfXGivenW *dists, MCGIDI_pdfsOfXGivenW_sampled *sampled, double rngValue ) {

    int iW, iX1;

    sampled->interpolationWY = dists->interpolationWY;
    sampled->interpolationXY = dists->interpolationXY;
    iW = sampled->iW = MCGIDI_misc_binarySearch( dists->numberOfWs, dists->Ws, sampled->w );
    sampled->frac = 1;

    if( iW == -2 ) {            /* w < first value of Ws. */
        return( MCGIDI_sampling_sampleX_from_pdfOfX( dists->dist, sampled, rngValue ) ); }
    else if( iW == -1 ) {       /* w > last value of Ws. */
        return( MCGIDI_sampling_sampleX_from_pdfOfX( &(dists->dist[dists->numberOfWs-1]), sampled, rngValue ) ); }
    else {
        if( MCGIDI_sampling_sampleX_from_pdfOfX( &(dists->dist[iW]), sampled, rngValue ) ) return( 1 );
        if( dists->interpolationWY != ptwXY_interpolationFlat ) {    // ptwXY_interpolationOther was not allowed at startup.
            double xSampled = sampled->x, frac = 1.;

            iX1 = sampled->iX1;
            if( MCGIDI_sampling_sampleX_from_pdfOfX( &(dists->dist[iW+1]), sampled, rngValue ) ) return( 1 );

            if( dists->interpolationWY == ptwXY_interpolationLinLin ) {
                frac = ( dists->Ws[iW+1] - sampled->w ) / ( dists->Ws[iW+1] - dists->Ws[iW] );
                sampled->x = frac * xSampled + ( 1 - frac ) * sampled->x; }
            else if( dists->interpolationWY == ptwXY_interpolationLogLin ) {
                frac = G4Log( dists->Ws[iW+1] / sampled->w ) / G4Log( dists->Ws[iW+1] / dists->Ws[iW] );
                sampled->x = frac * xSampled + ( 1 - frac ) * sampled->x; }
            else if( dists->interpolationWY == ptwXY_interpolationLinLog ) {
                frac = ( dists->Ws[iW+1] - sampled->w ) / ( dists->Ws[iW+1] - dists->Ws[iW] );
                sampled->x = xSampled * G4Pow::GetInstance()->powA( sampled->x / xSampled, frac ); }
            else if( dists->interpolationWY == ptwXY_interpolationLogLog ) {
                frac = G4Log( dists->Ws[iW+1] / sampled->w ) / G4Log( dists->Ws[iW+1] / dists->Ws[iW] );
                sampled->x = xSampled * G4Pow::GetInstance()->powA( sampled->x / xSampled, frac ); }
            else {  // This should never happen.
                smr_setReportError2( sampled->smr, smr_unknownID, 1, "bad interpolation = %d\n", dists->interpolationWY );
                return( 1 );
            }

            sampled->iX2 = sampled->iX1;
            sampled->iX1 = iX1;
            sampled->frac = frac;
        }
    }

    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_sampling_sampleX_from_pdfOfX( MCGIDI_pdfOfX *dist, MCGIDI_pdfsOfXGivenW_sampled *sampled, double rngValue ) {

    int iX;
    double d1, d2, frac;

    iX = sampled->iX1 = MCGIDI_misc_binarySearch( dist->numberOfXs, dist->cdf, rngValue );

    if( iX < 0 ) {                 /* This should never happen. */
        smr_setReportError2( sampled->smr, smr_unknownID, 1, "bad iX = %d\n", iX );
        sampled->x = dist->Xs[0];
        return( 1 );
    }
    if( sampled->interpolationXY == ptwXY_interpolationFlat ) {
        frac = ( dist->cdf[iX+1] - rngValue ) / ( dist->cdf[iX+1] - dist->cdf[iX] );
        sampled->x = frac * dist->Xs[iX] + ( 1 - frac ) * dist->Xs[iX+1]; }
    else {
        double s1 = dist->pdf[iX+1] - dist->pdf[iX];

        if( s1 == 0. ) {
            if( dist->pdf[iX] == 0 ) {
                sampled->x = dist->Xs[iX];
                if( iX == 0 ) sampled->x = dist->Xs[1]; }
            else {
                frac = ( dist->cdf[iX+1] - rngValue ) / ( dist->cdf[iX+1] - dist->cdf[iX] );
                sampled->x = frac * dist->Xs[iX] + ( 1 - frac ) * dist->Xs[iX+1];
            } }
        else {
            s1 = s1 / ( dist->Xs[iX+1] - dist->Xs[iX] );
            d1 = rngValue - dist->cdf[iX];
            d2 = dist->cdf[iX+1] - rngValue;
            if( d2 > d1 ) {     /* Closer to iX. */
                sampled->x = dist->Xs[iX] + ( std::sqrt( dist->pdf[iX] * dist->pdf[iX] + 2. * s1 * d1 ) - dist->pdf[iX] ) / s1; }
            else {              /* Closer to iX + 1. */
                sampled->x = dist->Xs[iX+1] - ( dist->pdf[iX+1] - std::sqrt( dist->pdf[iX+1] * dist->pdf[iX+1] - 2. * s1 * d2 ) ) / s1;
            }
        }
    }

    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_sampling_doubleDistribution( statusMessageReporting *smr, MCGIDI_pdfsOfXGivenW *pdfOfWGivenV, MCGIDI_pdfsOfXGivenW *pdfOfXGivenVAndW,  
        MCGIDI_quantitiesLookupModes &modes, MCGIDI_decaySamplingInfo *decaySamplingInfo ) {

    int iV;
    double e_in = modes.getProjectileEnergy( );
    double randomW = decaySamplingInfo->rng( decaySamplingInfo->rngState ), randomX = decaySamplingInfo->rng( decaySamplingInfo->rngState );
    MCGIDI_pdfsOfXGivenW_sampled sampledX, sampledW;
    ptwXY_interpolation interpolationWY = pdfOfWGivenV->interpolationWY; 

    sampledX.smr = smr;
    sampledW.smr = smr;
    sampledW.interpolationXY = pdfOfWGivenV->interpolationXY;
    iV = MCGIDI_misc_binarySearch( pdfOfWGivenV->numberOfWs, pdfOfWGivenV->Ws, e_in );
    if( iV < 0 ) {
        interpolationWY = ptwXY_interpolationFlat;
        if( iV == -2 ) {
            iV = 0; }
        else {
            iV = pdfOfWGivenV->numberOfWs - 1;
        }
        e_in = pdfOfWGivenV->Ws[iV];
    }
        
    MCGIDI_sampling_sampleX_from_pdfOfX( &(pdfOfWGivenV->dist[iV]), &sampledW, randomW );
    sampledX.w = sampledW.x;
    MCGIDI_sampling_sampleX_from_pdfsOfXGivenW( &(pdfOfXGivenVAndW[iV]), &sampledX, randomX );

    if( interpolationWY != ptwXY_interpolationFlat ) {
        double x = sampledX.x, w = sampledW.x, Vs[3] = { e_in, pdfOfWGivenV->Ws[iV], pdfOfWGivenV->Ws[iV+1] };

        MCGIDI_sampling_sampleX_from_pdfOfX( &(pdfOfWGivenV->dist[iV+1]), &sampledW, randomW );
        sampledX.w = sampledW.x;
        MCGIDI_sampling_sampleX_from_pdfsOfXGivenW( &(pdfOfXGivenVAndW[iV+1]), &sampledX, randomX );

        MCGIDI_sampling_interpolationValues( smr, interpolationWY, Vs, w, sampledW.x, &sampledW.x );
        MCGIDI_sampling_interpolationValues( smr, interpolationWY, Vs, x, sampledX.x, &sampledX.x );
    }
    
    decaySamplingInfo->mu = sampledW.x;
    decaySamplingInfo->Ep = sampledX.x;

    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_sampling_interpolationValues( statusMessageReporting *smr, ptwXY_interpolation interpolation, double *ws, double y1, double y2, double *y ) {

    double frac;

    if( interpolation == ptwXY_interpolationLinLin ) {
        frac = ( ws[2] - ws[0] ) / ( ws[2] - ws[1] );
        *y = frac * y1 + ( 1 - frac ) * y2; }
    else if( interpolation == ptwXY_interpolationLogLin ) {
        frac = G4Log( ws[2] / ws[0] ) / G4Log( ws[2] / ws[1] );
        *y = frac * y1 + ( 1 - frac ) * y2; }
    else if( interpolation == ptwXY_interpolationLinLog ) {
        frac = ( ws[2] - ws[0] ) / ( ws[2] - ws[1] );
        *y = y1 * G4Pow::GetInstance()->powA( y2 / y1, frac ); }
    else if( interpolation == ptwXY_interpolationLogLog ) {
        frac = G4Log( ws[2] / ws[0] ) / G4Log( ws[2] / ws[1] );
        *y = y2 * G4Pow::GetInstance()->powA( y2 / y1, frac ); }
    else {  // This should never happen.
        smr_setReportError2( smr, smr_unknownID, 1, "bad interpolation = %d\n", interpolation );
        return( 1 );
    }
    return( 0 );
}
/*
************************************************************
*/
double MCGIDI_sampling_ptwXY_getValueAtX( ptwXYPoints *ptwXY, double x1 ) {

    double y1;

    if( ptwXY_getValueAtX( ptwXY, x1, &y1 ) == nfu_XOutsideDomain ) {
        if( x1 < ptwXY_getXMin( ptwXY ) ) {
            ptwXY_getValueAtX( ptwXY, ptwXY_getXMin( ptwXY ), &y1 ); }
        else {
            ptwXY_getValueAtX( ptwXY, ptwXY_getXMax( ptwXY ), &y1 );
        }
    }
    return( y1 );
}

#if defined __cplusplus
}
#endif

