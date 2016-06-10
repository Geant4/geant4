/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <string.h>

#include "MCGIDI.h"
#include "MCGIDI_fromTOM.h"
#include "MCGIDI_misc.h"
#include "MCGIDI_private.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

static int MCGIDI_angularEnergy_parsePointwiseFromTOM( statusMessageReporting *smr, xDataTOM_element *pointwise, MCGIDI_distribution *distribution );
/*
************************************************************
*/
MCGIDI_angularEnergy *MCGIDI_angularEnergy_new( statusMessageReporting *smr ) {

    MCGIDI_angularEnergy *angularEnergy;

    if( ( angularEnergy = (MCGIDI_angularEnergy *) smr_malloc2( smr, sizeof( MCGIDI_angularEnergy ), 0, "angularEnergy" ) ) == NULL ) return( NULL );
    if( MCGIDI_angularEnergy_initialize( smr, angularEnergy ) ) angularEnergy = MCGIDI_angularEnergy_free( smr, angularEnergy );
    return( angularEnergy );
}
/*
************************************************************
*/
int MCGIDI_angularEnergy_initialize( statusMessageReporting * /*smr*/, MCGIDI_angularEnergy *angularEnergy ) {

    memset( angularEnergy, 0, sizeof( MCGIDI_angularEnergy ) );
    return( 0 );
}
/*
************************************************************
*/
MCGIDI_angularEnergy *MCGIDI_angularEnergy_free( statusMessageReporting *smr, MCGIDI_angularEnergy *angularEnergy ) {

    MCGIDI_angularEnergy_release( smr, angularEnergy );
    smr_freeMemory( (void **) &angularEnergy );
    return( NULL );
}
/*
************************************************************
*/
int MCGIDI_angularEnergy_release( statusMessageReporting *smr, MCGIDI_angularEnergy *angularEnergy ) {

    int i;

    for( i = 0; i < angularEnergy->pdfOfMuGivenE.numberOfWs; i++ ) {
        MCGIDI_sampling_pdfsOfXGivenW_release( smr, &(angularEnergy->pdfOfEpGivenEAndMu[i]) );
    }
    smr_freeMemory( (void **) &(angularEnergy->pdfOfEpGivenEAndMu) );
    MCGIDI_sampling_pdfsOfXGivenW_release( smr, &(angularEnergy->pdfOfMuGivenE) );
    MCGIDI_angularEnergy_initialize( smr, angularEnergy );

    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_angularEnergy_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_distribution *distribution ) {

    xDataTOM_element *angularEnergyElement, *pointwise = NULL;
    char const *nativeData;

    if( ( angularEnergyElement = xDataTOME_getOneElementByName( smr, element, "angularEnergy", 1 ) ) == NULL ) goto err;

    if( ( nativeData = xDataTOM_getAttributesValueInElement( angularEnergyElement, "nativeData" ) ) == NULL ) goto err;
    if( strcmp( nativeData, "pointwise" ) == 0 ) {
        if( ( pointwise = xDataTOME_getOneElementByName( smr, angularEnergyElement, "pointwise", 1 ) ) == NULL ) goto err; }
    else if( strcmp( nativeData, "linear" ) == 0 ) {
        if( ( pointwise = xDataTOME_getOneElementByName( smr, angularEnergyElement, "linear", 1 ) ) == NULL ) goto err; }
    else {
        smr_setReportError2( smr, smr_unknownID, 1, "angularEnergy nativeData = '%s' not supported", nativeData );
        goto err;
    }
    if( pointwise != NULL ) return( MCGIDI_angularEnergy_parsePointwiseFromTOM( smr, pointwise, distribution ) );

    return( 0 );

err:
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_angularEnergy_parsePointwiseFromTOM( statusMessageReporting *smr, xDataTOM_element *pointwise, MCGIDI_distribution *distribution ) {

    int iV, iW;
    double y, norm, energyInFactor;
    char const *energyUnit, *energyOutProbabilityUnits[2] = { "MeV", "1/MeV" };
    MCGIDI_angularEnergy *angularEnergy = NULL;
    ptwXY_interpolation interpolationXY, interpolationWY, interpolationVY;
    xDataTOM_XYs *XYs;
    xDataTOM_W_XYs *W_XYs;
    xDataTOM_V_W_XYs *V_W_XYs;
    MCGIDI_pdfsOfXGivenW *pdfOfMuGivenE, *pdfOfEpGivenEAndMu = NULL, *pdfOfEpGivenEAndMu2 = NULL;
    ptwXYPoints *pdfXY1 = NULL, *pdfXY2 = NULL;
    nfu_status status;

    if( MCGIDI_fromTOM_interpolation( smr, pointwise, 0, &interpolationVY ) ) goto err;
    if( MCGIDI_fromTOM_interpolation( smr, pointwise, 1, &interpolationWY ) ) goto err;
    if( MCGIDI_fromTOM_interpolation( smr, pointwise, 2, &interpolationXY ) ) goto err;
    if( ( angularEnergy = MCGIDI_angularEnergy_new( smr ) ) == NULL ) goto err;

    if( ( angularEnergy->frame = MCGIDI_misc_getProductFrame( smr, pointwise ) ) == xDataTOM_frame_invalid ) goto err;

    pdfOfMuGivenE = &(angularEnergy->pdfOfMuGivenE);
    pdfOfMuGivenE->interpolationWY = interpolationVY;
    pdfOfMuGivenE->interpolationXY = interpolationWY;

    if( ( V_W_XYs = (xDataTOM_V_W_XYs *) xDataTOME_getXDataIfID( smr, pointwise, "V_W_XYs" ) ) == NULL ) goto err;
    if( ( pdfOfMuGivenE->Ws = (double *) smr_malloc2( smr, V_W_XYs->length * sizeof( double ), 1, "pdfOfMuGivenE->Ws" ) ) == NULL ) goto err;
    if( ( pdfOfMuGivenE->dist = (MCGIDI_pdfOfX *) smr_malloc2( smr, V_W_XYs->length * sizeof( MCGIDI_pdfOfX ), 0, "pdfOfMuGivenE->dist" ) ) == NULL ) goto err;
    if( ( pdfOfEpGivenEAndMu = (MCGIDI_pdfsOfXGivenW *) smr_malloc2( smr, V_W_XYs->length * sizeof( MCGIDI_pdfsOfXGivenW ), 1, "pdfOfEpGivenEAndMu" ) ) == NULL ) goto err;

    energyUnit = xDataTOM_subAxes_getUnit( smr, &(V_W_XYs->subAxes), 0 );
    if( !smr_isOk( smr ) ) goto err;
    energyInFactor = MCGIDI_misc_getUnitConversionFactor( smr, energyUnit, "MeV" );
    if( !smr_isOk( smr ) ) goto err;

    for( iV = 0; iV < V_W_XYs->length; iV++ ) {
        W_XYs = &(V_W_XYs->W_XYs[iV]);
        pdfOfEpGivenEAndMu2 = &(pdfOfEpGivenEAndMu[iV]);
        pdfOfEpGivenEAndMu2->interpolationWY = interpolationWY;
        pdfOfEpGivenEAndMu2->interpolationXY = interpolationXY;
        if( ( pdfXY2 = ptwXY_new( interpolationWY, NULL, 2., 1e-6, W_XYs->length, 10, &status, 0 ) ) == NULL ) goto errA;
        if( ( pdfOfEpGivenEAndMu2->Ws = (double *) smr_malloc2( smr, W_XYs->length * sizeof( double ), 1, "pdfOfEpGivenEAndMu2->Ws" ) ) == NULL ) goto err;
        if( ( pdfOfEpGivenEAndMu2->dist = (MCGIDI_pdfOfX *) smr_malloc2( smr, W_XYs->length * sizeof( MCGIDI_pdfOfX ), 0, "pdfOfEpGivenEAndMu2->dist" ) ) == NULL ) goto err;
        for( iW = 0; iW < W_XYs->length; iW++ ) {
            XYs = &(W_XYs->XYs[iW]);
            if( ( pdfXY1 =  MCGIDI_misc_dataFromXYs2ptwXYPointsInUnitsOf( smr, XYs, interpolationXY, energyOutProbabilityUnits ) ) == NULL ) goto err;
            y = ptwXY_integrateDomain( pdfXY1, &status );
            if( ( status = ptwXY_setValueAtX( pdfXY2, XYs->value, y ) ) != nfu_Okay ) goto errA;
            if( status != nfu_Okay ) goto errA;

            if( y == 0 ) {
                if( ( status = ptwXY_add_double( pdfXY1, 0.5 ) ) != nfu_Okay ) goto errA;
            }
            pdfOfEpGivenEAndMu2->Ws[iW] = XYs->value;
            if( MCGIDI_fromTOM_pdfOfX( smr, pdfXY1, &(pdfOfEpGivenEAndMu2->dist[iW]), &norm ) ) goto err;
            pdfOfEpGivenEAndMu2->numberOfWs++;

            pdfXY1 = ptwXY_free( pdfXY1 );
        }
        pdfOfMuGivenE->Ws[iV] = energyInFactor * W_XYs->value;
        if( MCGIDI_fromTOM_pdfOfX( smr, pdfXY2, &(pdfOfMuGivenE->dist[iV]), &norm ) ) goto err;
        pdfOfMuGivenE->numberOfWs++;

        pdfXY2 = ptwXY_free( pdfXY2 );
    }

    angularEnergy->pdfOfEpGivenEAndMu = pdfOfEpGivenEAndMu;
    distribution->angularEnergy = angularEnergy;
    distribution->type = MCGIDI_distributionType_angularEnergy_e;

    return( 0 );

errA:
    smr_setReportError2( smr, smr_unknownID, 1, "ptwXY_integrateDomain err = %d: %s\n", status, nfu_statusMessage( status ) );
err:
    if( pdfXY1 != NULL ) ptwXY_free( pdfXY1 );
    if( pdfXY2 != NULL ) ptwXY_free( pdfXY2 );
/* Need to free pdfOfEpGivenEAndMu. */
    if( angularEnergy != NULL ) MCGIDI_angularEnergy_free( smr, angularEnergy );
    return( 1 );
}
/*
************************************************************
*/
int MCGIDI_angularEnergy_sampleDistribution( statusMessageReporting *smr, MCGIDI_angularEnergy *angularEnergy, MCGIDI_quantitiesLookupModes &modes,
        MCGIDI_decaySamplingInfo *decaySamplingInfo ) {

    int status = MCGIDI_sampling_doubleDistribution( smr, &(angularEnergy->pdfOfMuGivenE), angularEnergy->pdfOfEpGivenEAndMu, modes, decaySamplingInfo );

    decaySamplingInfo->frame = angularEnergy->frame;
    return( status );
}

#if defined __cplusplus
}
#endif

