/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <string.h>

#include "MCGIDI_fromTOM.h"
#include "MCGIDI.h"
#include "MCGIDI_misc.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

static int MCGIDI_energyAngular_linear_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_distribution *distribution );
/*
************************************************************
*/
MCGIDI_energyAngular *MCGIDI_energyAngular_new( statusMessageReporting *smr ) {

    MCGIDI_energyAngular *energyAngular;

    if( ( energyAngular = (MCGIDI_energyAngular *) smr_malloc2( smr, sizeof( MCGIDI_energyAngular ), 0, "energyAngular" ) ) == NULL ) return( NULL );
    if( MCGIDI_energyAngular_initialize( smr, energyAngular ) ) energyAngular = MCGIDI_energyAngular_free( smr, energyAngular );
    return( energyAngular );
}
/*
************************************************************
*/
int MCGIDI_energyAngular_initialize( statusMessageReporting * /*smr*/, MCGIDI_energyAngular *energyAngular ) {

    memset( energyAngular, 0, sizeof( MCGIDI_energyAngular ) );
    return( 0 );
}
/*
************************************************************
*/
MCGIDI_energyAngular *MCGIDI_energyAngular_free( statusMessageReporting *smr, MCGIDI_energyAngular *energyAngular ) {

    MCGIDI_energyAngular_release( smr, energyAngular );
    smr_freeMemory( (void **) &energyAngular );
    return( NULL );
}
/*
************************************************************
*/
int MCGIDI_energyAngular_release( statusMessageReporting *smr, MCGIDI_energyAngular *energyAngular ) {

    int i;

    for( i = 0; i < energyAngular->pdfOfEpGivenE.numberOfWs; i++ ) {
        MCGIDI_sampling_pdfsOfXGivenW_release( smr, &(energyAngular->pdfOfMuGivenEAndEp[i]) );
    }
    smr_freeMemory( (void **) &(energyAngular->pdfOfMuGivenEAndEp) );
    MCGIDI_sampling_pdfsOfXGivenW_release( smr, &(energyAngular->pdfOfEpGivenE) );
    MCGIDI_energyAngular_initialize( smr, energyAngular );

    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_energyAngular_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_distribution *distribution ) {

    xDataTOM_element *energyAngularElement;
    char const *nativeData;

    if( ( energyAngularElement = xDataTOME_getOneElementByName( smr, element, "energyAngular", 1 ) ) == NULL ) goto err;

    if( ( nativeData = xDataTOM_getAttributesValueInElement( energyAngularElement, "nativeData" ) ) == NULL ) goto err;
    if( strcmp( nativeData, "KalbachMann" ) == 0 ) {
        return( MCGIDI_KalbachMann_parseFromTOM( smr, energyAngularElement, distribution ) ); }
    else if( strcmp( nativeData, "linear" ) == 0 ) {
        return( MCGIDI_energyAngular_linear_parseFromTOM( smr, energyAngularElement, distribution ) ); }
    else {
        smr_setReportError2( smr, smr_unknownID, 1, "energyAngular nativeData = '%s' not supported", nativeData );
        goto err;
    }

    return( 0 );

err:
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_energyAngular_linear_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_distribution *distribution ) {

    int iV, iW;
    double y, norm, energyInFactor, energyOutFactor;
    char const *energyUnit, *multiplicityProbabilityUnits[2] = { "", "1/MeV" };
    xDataTOM_element *linear;
    ptwXY_interpolation interpolationXY, interpolationWY, interpolationVY;
    xDataTOM_XYs *XYs;
    xDataTOM_W_XYs *W_XYs;
    xDataTOM_V_W_XYs *V_W_XYs;
    MCGIDI_pdfsOfXGivenW *pdfOfEpGivenE, *pdfOfMuGivenEAndEp = NULL, *pdfOfMuGivenEAndEp2 = NULL;
    MCGIDI_energyAngular *energyAngular = NULL;
    ptwXYPoints *pdfXY1 = NULL, *pdfXY2 = NULL;
    nfu_status status;

    if( ( linear = xDataTOME_getOneElementByName( smr, element, "linear", 1 ) ) == NULL ) goto err;

    if( MCGIDI_fromTOM_interpolation( smr, linear, 0, &interpolationVY ) ) goto err;
    if( MCGIDI_fromTOM_interpolation( smr, linear, 1, &interpolationWY ) ) goto err;
    if( MCGIDI_fromTOM_interpolation( smr, linear, 2, &interpolationXY ) ) goto err;

    if( ( energyAngular = MCGIDI_energyAngular_new( smr ) ) == NULL ) goto err;
    if( ( energyAngular->frame = MCGIDI_misc_getProductFrame( smr, linear ) ) == xDataTOM_frame_invalid ) goto err;

    pdfOfEpGivenE = &(energyAngular->pdfOfEpGivenE);
    pdfOfEpGivenE->interpolationWY = interpolationVY;
    pdfOfEpGivenE->interpolationXY = interpolationWY;

    if( ( V_W_XYs = (xDataTOM_V_W_XYs *) xDataTOME_getXDataIfID( smr, linear, "V_W_XYs" ) ) == NULL ) goto err;
    if( ( pdfOfEpGivenE->Ws = (double *) smr_malloc2( smr, V_W_XYs->length * sizeof( double ), 1, "pdfOfEpGivenE->Ws" ) ) == NULL ) goto err;
    if( ( pdfOfEpGivenE->dist = (MCGIDI_pdfOfX *) smr_malloc2( smr, V_W_XYs->length * sizeof( MCGIDI_pdfOfX ), 0, "pdfOfEpGivenE->dist" ) ) == NULL ) goto err;
    if( ( pdfOfMuGivenEAndEp = (MCGIDI_pdfsOfXGivenW *) smr_malloc2( smr, V_W_XYs->length * sizeof( MCGIDI_pdfsOfXGivenW ), 1, "pdfOfMuGivenEAndEp" ) ) == NULL ) goto err;

    energyUnit = xDataTOM_subAxes_getUnit( smr, &(V_W_XYs->subAxes), 0 );
    if( !smr_isOk( smr ) ) goto err;
    energyInFactor = MCGIDI_misc_getUnitConversionFactor( smr, energyUnit, "MeV" );
    if( !smr_isOk( smr ) ) goto err;

    energyUnit = xDataTOM_subAxes_getUnit( smr, &(V_W_XYs->subAxes), 1 );
    if( !smr_isOk( smr ) ) goto err;
    energyOutFactor = MCGIDI_misc_getUnitConversionFactor( smr, energyUnit, "MeV" );
    if( !smr_isOk( smr ) ) goto err;

    for( iV = 0; iV < V_W_XYs->length; iV++ ) {
        W_XYs = &(V_W_XYs->W_XYs[iV]);
        pdfOfMuGivenEAndEp2 = &(pdfOfMuGivenEAndEp[iV]);
        pdfOfMuGivenEAndEp2->interpolationWY = interpolationWY;
        pdfOfMuGivenEAndEp2->interpolationXY = interpolationXY;
        if( ( pdfXY2 = ptwXY_new( interpolationWY, NULL, 2., 1e-6, W_XYs->length, 10, &status, 0 ) ) == NULL ) goto errA;
        if( ( pdfOfMuGivenEAndEp2->Ws = (double *) smr_malloc2( smr, W_XYs->length * sizeof( double ), 1, "pdfOfMuGivenEAndEp2->Ws" ) ) == NULL ) goto err;
        if( ( pdfOfMuGivenEAndEp2->dist = (MCGIDI_pdfOfX *) smr_malloc2( smr, W_XYs->length * sizeof( MCGIDI_pdfOfX ), 0, "pdfOfMuGivenEAndEp2->dist" ) ) == NULL ) goto err;
        for( iW = 0; iW < W_XYs->length; iW++ ) {
            XYs = &(W_XYs->XYs[iW]);
            if( ( pdfXY1 =  MCGIDI_misc_dataFromXYs2ptwXYPointsInUnitsOf( smr, XYs, interpolationXY, multiplicityProbabilityUnits ) ) == NULL ) goto err;
            y = ptwXY_integrateDomain( pdfXY1, &status );
            if( ( status = ptwXY_setValueAtX( pdfXY2, energyOutFactor * XYs->value, y ) ) != nfu_Okay ) goto errA;
            if( status != nfu_Okay ) goto errA;

            if( y == 0 ) {
                if( ( status = ptwXY_add_double( pdfXY1, 0.5 ) ) != nfu_Okay ) goto errA;
            }
            pdfOfMuGivenEAndEp2->Ws[iW] = energyOutFactor * XYs->value;
            if( MCGIDI_fromTOM_pdfOfX( smr, pdfXY1, &(pdfOfMuGivenEAndEp2->dist[iW]), &norm ) ) goto err;
            pdfOfMuGivenEAndEp2->numberOfWs++;

            pdfXY1 = ptwXY_free( pdfXY1 );
        }
        pdfOfEpGivenE->Ws[iV] = energyInFactor * W_XYs->value;
        if( MCGIDI_fromTOM_pdfOfX( smr, pdfXY2, &(pdfOfEpGivenE->dist[iV]), &norm ) ) goto err;
        pdfOfEpGivenE->numberOfWs++;

        pdfXY2 = ptwXY_free( pdfXY2 );
    }
    energyAngular->pdfOfMuGivenEAndEp = pdfOfMuGivenEAndEp;
    distribution->energyAngular = energyAngular;
    distribution->type = MCGIDI_distributionType_energyAngular_e;

    return( 0 );

errA:
    smr_setReportError2( smr, smr_unknownID, 1, "ptwXY_integrateDomain err = %d: %s\n", status, nfu_statusMessage( status ) );
err:
    if( pdfXY1 != NULL ) ptwXY_free( pdfXY1 );
    if( pdfXY2 != NULL ) ptwXY_free( pdfXY2 );
    if( energyAngular != NULL ) MCGIDI_energyAngular_free( smr, energyAngular );
/* ????????? Need to free pdfOfMuGivenEAndEp, now may be handled by MCGIDI_energyAngular_free. Need to check. */
    return( 1 );
}
/*
************************************************************
*/
int MCGIDI_energyAngular_sampleDistribution( statusMessageReporting *smr, MCGIDI_distribution *distribution, MCGIDI_quantitiesLookupModes &modes,
        MCGIDI_decaySamplingInfo *decaySamplingInfo ) {

    double Ep;
    MCGIDI_energyAngular *energyAngular = distribution->energyAngular;

    MCGIDI_sampling_doubleDistribution( smr, &(energyAngular->pdfOfEpGivenE), energyAngular->pdfOfMuGivenEAndEp, modes, decaySamplingInfo );
    Ep = decaySamplingInfo->mu;
    decaySamplingInfo->mu = decaySamplingInfo->Ep;
    decaySamplingInfo->Ep = Ep;
    decaySamplingInfo->frame = energyAngular->frame;

    return( 0 );
}

#if defined __cplusplus
}
#endif

