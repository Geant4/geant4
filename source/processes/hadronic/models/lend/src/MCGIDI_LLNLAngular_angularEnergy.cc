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

static int MCGIDI_LLNL_angularEnergy_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_distribution *distribution );
static int MCGIDI_LLNL_angularEnergy_parsePointwiseFromTOM( statusMessageReporting *smr, xDataTOM_element *pointwise, MCGIDI_distribution *distribution );
/*
************************************************************
*/
int MCGIDI_LLNLAngular_angularEnergy_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_distribution *distribution ) {

    xDataTOM_element *angularEnergyElement;

    if( ( angularEnergyElement = xDataTOME_getOneElementByName( smr, element, "LLNLAngular_angularEnergy", 1 ) ) == NULL ) return( 1 );

    if( MCGIDI_angular_parseFromTOM( smr, angularEnergyElement, distribution, NULL ) ) goto err;
    if( MCGIDI_LLNL_angularEnergy_parseFromTOM( smr, angularEnergyElement, distribution ) ) goto err;

    return( 0 );

err:
    if( distribution->angular ) distribution->angular = MCGIDI_angular_free( smr, distribution->angular );
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_LLNL_angularEnergy_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_distribution *distribution ) {

    xDataTOM_element *angularEnergyElement, *pointwise = NULL;
    char const *nativeData;

    if( ( angularEnergyElement = xDataTOME_getOneElementByName( smr, element, "LLNLAngularEnergy", 1 ) ) == NULL ) goto err;

    if( ( nativeData = xDataTOM_getAttributesValueInElement( angularEnergyElement, "nativeData" ) ) == NULL ) goto err;
    if( strcmp( nativeData, "pointwise" ) == 0 ) {
        if( ( pointwise = xDataTOME_getOneElementByName( smr, angularEnergyElement, "pointwise", 1 ) ) == NULL ) goto err; }
    else if( strcmp( nativeData, "linear" ) == 0 ) {
        if( ( pointwise = xDataTOME_getOneElementByName( smr, angularEnergyElement, "linear", 1 ) ) == NULL ) goto err; }
    else {
        smr_setReportError2( smr, smr_unknownID, 1, "angularEnergy nativeData = '%s' not supported", nativeData );
        goto err;
    }
    if( pointwise != NULL ) return( MCGIDI_LLNL_angularEnergy_parsePointwiseFromTOM( smr, pointwise, distribution ) );

    return( 0 );

err:
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_LLNL_angularEnergy_parsePointwiseFromTOM( statusMessageReporting *smr, xDataTOM_element *pointwise, MCGIDI_distribution *distribution ) {

    int iV = 0, iW;
    double y1, norm/*, energyInFactor*/;
    char const /*energyUnit,*/ *energyOutProbabilityUnits[2] = { "MeV", "1/MeV" };
    MCGIDI_angularEnergy *angularEnergy = NULL;
    ptwXY_interpolation interpolationXY, interpolationWY, interpolationVY;
    xDataTOM_XYs *XYs;
    xDataTOM_W_XYs *W_XYs;
    xDataTOM_V_W_XYs *V_W_XYs;
    MCGIDI_pdfsOfXGivenW *pdfOfEpGivenEAndMu = NULL, *pdfOfEpGivenEAndMu2 = NULL;
    ptwXYPoints *pdfXY1 = NULL;
    nfu_status status;
    enum xDataTOM_frame frame;

    if( ( frame = MCGIDI_misc_getProductFrame( smr, pointwise ) ) == xDataTOM_frame_invalid ) goto err;
    if( MCGIDI_fromTOM_interpolation( smr, pointwise, 0, &interpolationVY ) ) goto err;
    if( MCGIDI_fromTOM_interpolation( smr, pointwise, 1, &interpolationWY ) ) goto err;
    if( MCGIDI_fromTOM_interpolation( smr, pointwise, 2, &interpolationXY ) ) goto err;

    if( ( V_W_XYs = (xDataTOM_V_W_XYs *) xDataTOME_getXDataIfID( smr, pointwise, "V_W_XYs" ) ) == NULL ) goto err;
    /*energyUnit = xDataTOM_subAxes_getUnit( smr, &(V_W_XYs->subAxes), 0 );*/
    if( !smr_isOk( smr ) ) goto err;
    /*energyInFactor = MCGIDI_misc_getUnitConversionFactor( smr, energyUnit, "MeV" );*/
    if( !smr_isOk( smr ) ) goto err;

    if( ( pdfOfEpGivenEAndMu = (MCGIDI_pdfsOfXGivenW *) smr_malloc2( smr, V_W_XYs->length * sizeof( MCGIDI_pdfsOfXGivenW ), 1, "pdfOfEpGivenEAndMu" ) ) == NULL ) goto err;
    for( iV = 0; iV < V_W_XYs->length; iV++ ) {
        W_XYs = &(V_W_XYs->W_XYs[iV]);
        pdfOfEpGivenEAndMu2 = &(pdfOfEpGivenEAndMu[iV]);
        pdfOfEpGivenEAndMu2->Ws = NULL;
        pdfOfEpGivenEAndMu2->dist = NULL;

        pdfOfEpGivenEAndMu2->interpolationWY = interpolationWY;
        pdfOfEpGivenEAndMu2->interpolationXY = interpolationXY;
        if( ( pdfOfEpGivenEAndMu2->Ws = (double *) smr_malloc2( smr, W_XYs->length * sizeof( double ), 1, "pdfOfEpGivenEAndMu2->Ws" ) ) == NULL ) goto err;
        if( ( pdfOfEpGivenEAndMu2->dist = (MCGIDI_pdfOfX *) smr_malloc2( smr, W_XYs->length * sizeof( MCGIDI_pdfOfX ), 0, "pdfOfEpGivenEAndMu2->dist" ) ) == NULL ) goto err;

        for( iW = 0; iW < W_XYs->length; iW++ ) {
            XYs = &(W_XYs->XYs[iW]);
            if( ( pdfXY1 =  MCGIDI_misc_dataFromXYs2ptwXYPointsInUnitsOf( smr, XYs, interpolationXY, energyOutProbabilityUnits ) ) == NULL ) goto err;
            y1 = ptwXY_integrateDomain( pdfXY1, &status );
            if( status != nfu_Okay ) goto errA;

            if( y1 == 0 ) {
                if( ( status = ptwXY_add_double( pdfXY1, 0.5 ) ) != nfu_Okay ) goto errA;
            }
            pdfOfEpGivenEAndMu2->Ws[iW] = XYs->value;
            if( MCGIDI_fromTOM_pdfOfX( smr, pdfXY1, &(pdfOfEpGivenEAndMu2->dist[iW]), &norm ) ) goto err;
            pdfOfEpGivenEAndMu2->numberOfWs++;

            pdfXY1 = ptwXY_free( pdfXY1 );
        }
        pdfOfEpGivenEAndMu2 = NULL;
    }

    if( ( angularEnergy = MCGIDI_angularEnergy_new( smr ) ) == NULL ) goto err;
    angularEnergy->frame = frame;

    angularEnergy->pdfOfMuGivenE.numberOfWs = distribution->angular->dists.numberOfWs;
    angularEnergy->pdfOfMuGivenE.interpolationWY = distribution->angular->dists.interpolationWY;
    angularEnergy->pdfOfMuGivenE.interpolationXY = distribution->angular->dists.interpolationXY;
    angularEnergy->pdfOfMuGivenE.Ws = distribution->angular->dists.Ws;
    angularEnergy->pdfOfMuGivenE.dist = distribution->angular->dists.dist;
    smr_freeMemory( (void **) &(distribution->angular) );
    distribution->angular = NULL;

    angularEnergy->pdfOfEpGivenEAndMu = pdfOfEpGivenEAndMu;
    distribution->angularEnergy = angularEnergy;
    distribution->type = MCGIDI_distributionType_angularEnergy_e;

    return( 0 );

errA:
    smr_setReportError2( smr, smr_unknownID, 1, "ptwXY_integrateDomain err = %d: %s\n", status, nfu_statusMessage( status ) );
err:
    if( pdfXY1 != NULL ) ptwXY_free( pdfXY1 );
    if( pdfOfEpGivenEAndMu2 != NULL ) MCGIDI_sampling_pdfsOfXGivenW_release( smr, pdfOfEpGivenEAndMu2 );
    if( pdfOfEpGivenEAndMu != NULL ) {
        for( ; iV > 0; iV-- ) MCGIDI_sampling_pdfsOfXGivenW_release( smr, &(pdfOfEpGivenEAndMu[iV]) );
        smr_freeMemory( (void **) &pdfOfEpGivenEAndMu );
    }
    return( 1 );
}

#if defined __cplusplus
}
#endif
