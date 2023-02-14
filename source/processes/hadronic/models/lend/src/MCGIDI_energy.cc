/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <string.h>
#include <cmath>

#ifdef WIN32
#define M_PI 3.141592653589793238463
#endif

#include "MCGIDI_fromTOM.h"
#include "MCGIDI_misc.h"
#include "MCGIDI_private.h"
#include <nf_specialFunctions.h>

#if defined __cplusplus
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"
namespace GIDI {
using namespace GIDI;
#endif

static int MCGIDI_energy_parseWeightFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_energyWeightedFunctional *weightedFunctional );
static int MCGIDI_energy_parseWeightedFunctionalsFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_energy *energy );
static int MCGIDI_energy_parseGeneralEvaporationFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_energy *energy );
static int MCGIDI_energy_parseEvaporationFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_energy *energy );
static int MCGIDI_energy_parseWattFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_energy *energy );
static int MCGIDI_energy_parseSimpleMaxwellianFissionFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_energy *energy );
static int MCGIDI_energy_parseMadlandNixFromTOM( statusMessageReporting *smr, xDataTOM_element *functional, MCGIDI_energy *energy );
static nfu_status MCGIDI_energy_parseMadlandNixFromTOM_callback( double x, double *y, void *argList );
static double MCGIDI_energy_parseMadlandNixFromTOM_callback_g( double Ep, double EFL, double T_M, nfu_status *status );
static int MCGIDI_energy_parseNBodyPhaseSpaceFromTOM( statusMessageReporting *smr, xDataTOM_element *functional, MCGIDI_energy *energy,
        MCGIDI_distribution *distribution );

static int MCGIDI_energy_sampleSimpleMaxwellianFission( statusMessageReporting *smr, double e_in_U_theta, MCGIDI_decaySamplingInfo *decaySamplingInfo );
static int MCGIDI_energy_sampleEvaporation( statusMessageReporting *smr, double e_in_U_theta, MCGIDI_decaySamplingInfo *decaySamplingInfo );
static int MCGIDI_energy_sampleWatt( statusMessageReporting *smr, double e_in_U, double Watt_a, double Watt_b, MCGIDI_decaySamplingInfo *decaySamplingInfo );
static int MCGIDI_energy_sampleWeightedFunctional( statusMessageReporting *smr, MCGIDI_energy *energy, 
    MCGIDI_quantitiesLookupModes &modes, MCGIDI_decaySamplingInfo *decaySamplingInfo );
static nfu_status MCGIDI_energy_NBodyPhaseSpacePDF_callback( double x, double *y, void *argList );
/*
************************************************************
*/
MCGIDI_energy *MCGIDI_energy_new( statusMessageReporting *smr ) {

    MCGIDI_energy *energy;

    if( ( energy = (MCGIDI_energy *) smr_malloc2( smr, sizeof( MCGIDI_energy ), 0, "energy" ) ) == NULL ) return( NULL );
    if( MCGIDI_energy_initialize( smr, energy ) ) energy = MCGIDI_energy_free( smr, energy );
    return( energy );
}
/*
************************************************************
*/
int MCGIDI_energy_initialize( statusMessageReporting * /*smr*/, MCGIDI_energy *energy ) {

    memset( energy, 0, sizeof( MCGIDI_energy ) );
    return( 0 );
}
/*
************************************************************
*/
MCGIDI_energy *MCGIDI_energy_free( statusMessageReporting *smr, MCGIDI_energy *energy ) {

    MCGIDI_energy_release( smr, energy );
    smr_freeMemory( (void **) &energy );
    return( NULL );
}
/*
************************************************************
*/
int MCGIDI_energy_release( statusMessageReporting *smr, MCGIDI_energy *energy ) {

    int i;

    MCGIDI_sampling_pdfsOfXGivenW_release( smr, &(energy->dists) );
    if( energy->theta ) energy->theta = ptwXY_free( energy->theta );
    if( energy->Watt_a ) energy->Watt_a = ptwXY_free( energy->Watt_a );
    if( energy->Watt_b ) energy->Watt_b = ptwXY_free( energy->Watt_b );
    if( ( energy->type == MCGIDI_energyType_generalEvaporation ) || ( energy->type == MCGIDI_energyType_NBodyPhaseSpace ) ) {
        MCGIDI_sampling_pdfsOfX_release( smr, &(energy->g) ); }
    else if( energy->type == MCGIDI_energyType_weightedFunctional ) {
        for( i = 0; i < energy->weightedFunctionals.numberOfWeights; i++ ) {
            ptwXY_free( energy->weightedFunctionals.weightedFunctional[i].weight );
            MCGIDI_energy_free( smr, energy->weightedFunctionals.weightedFunctional[i].energy );
        }
    }

    MCGIDI_energy_initialize( smr, energy );
    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_energy_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_distribution *distribution, ptwXYPoints *norms,
        enum MCGIDI_energyType energyType, double gammaEnergy_MeV ) {

    MCGIDI_energy *energy = NULL;
    xDataTOM_element *energyElement, *linearElement, *functional, *frameElement;
    char const *nativeData;
    double projectileMass_MeV, targetMass_MeV;

    if( ( energy = MCGIDI_energy_new( smr ) ) == NULL ) goto err;

    projectileMass_MeV = MCGIDI_product_getProjectileMass_MeV( smr, distribution->product );
    targetMass_MeV = MCGIDI_product_getTargetMass_MeV( smr, distribution->product );
    energy->e_inCOMFactor = targetMass_MeV / ( projectileMass_MeV + targetMass_MeV );

    if( ( energyType == MCGIDI_energyType_primaryGamma ) || ( energyType == MCGIDI_energyType_discreteGamma ) ) {
        energy->type = energyType;
        energy->gammaEnergy_MeV = gammaEnergy_MeV;
        energy->frame = xDataTOM_frame_lab;            /* BRB. This should not be hardwired?????? Probably needs to be changed in GND also. */
        if( energyType == MCGIDI_energyType_primaryGamma ) energy->primaryGammaMassFactor = energy->e_inCOMFactor; }
    else {
        if( ( energyElement = xDataTOME_getOneElementByName( smr, element, "energy", 1 ) ) == NULL ) goto err;
        if( ( nativeData = xDataTOM_getAttributesValueInElement( energyElement, "nativeData" ) ) == NULL ) goto err;
        if( ( linearElement = xDataTOME_getOneElementByName( NULL, energyElement, "linear", 0 ) ) == NULL ) 
            linearElement = xDataTOME_getOneElementByName( NULL, energyElement, "pointwise", 0 );
        if( linearElement == NULL ) {
            if( ( functional = xDataTOME_getOneElementByName( NULL, energyElement, "generalEvaporation", 0 ) ) != NULL ) {
                if( MCGIDI_energy_parseGeneralEvaporationFromTOM( smr, functional, energy ) ) goto err; }
            else if( ( functional = xDataTOME_getOneElementByName( NULL, energyElement, "simpleMaxwellianFission", 0 ) ) != NULL ) {
                if( MCGIDI_energy_parseSimpleMaxwellianFissionFromTOM( smr, functional, energy ) ) goto err; }
            else if( ( functional = xDataTOME_getOneElementByName( NULL, energyElement, "evaporation", 0 ) ) != NULL ) {
                if( MCGIDI_energy_parseEvaporationFromTOM( smr, functional, energy ) ) goto err; }
            else if( ( functional = xDataTOME_getOneElementByName( NULL, energyElement, "Watt", 0 ) ) != NULL ) {
                if( MCGIDI_energy_parseWattFromTOM( smr, functional, energy ) ) goto err; }
            else if( ( functional = xDataTOME_getOneElementByName( NULL, energyElement, "MadlandNix", 0 ) ) != NULL ) {
                if( MCGIDI_energy_parseMadlandNixFromTOM( smr, functional, energy ) ) goto err; }
            else if( ( functional = xDataTOME_getOneElementByName( NULL, energyElement, "NBodyPhaseSpace", 0 ) ) != NULL ) {
                if( MCGIDI_energy_parseNBodyPhaseSpaceFromTOM( smr, functional, energy, distribution ) ) goto err; }
            else if( ( functional = xDataTOME_getOneElementByName( NULL, energyElement, "weightedFunctionals", 0 ) ) != NULL ) {
                if( MCGIDI_energy_parseWeightedFunctionalsFromTOM( smr, functional, energy ) ) goto err; }
            else {
                smr_setReportError2( smr, smr_unknownID, 1, "unsupported energy type: nativeData = '%s'", nativeData );
                goto err;
            }
            frameElement = functional; }
        else {
            char const *toUnits[3] = { "MeV", "MeV", "1/MeV" };

            frameElement = linearElement;
            if( MCGIDI_fromTOM_pdfsOfXGivenW( smr, linearElement, &(energy->dists), norms, toUnits ) ) goto err;
            energy->type = MCGIDI_energyType_linear;
        }
        if( ( energy->frame = MCGIDI_misc_getProductFrame( smr, frameElement ) ) == xDataTOM_frame_invalid ) goto err;
    }
    distribution->energy = energy;

    return( 0 );

err:
    if( energy != NULL ) MCGIDI_energy_free( smr, energy );
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_energy_parseWeightedFunctionalsFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_energy *energy ) {

    int i;
    xDataTOM_element *child;

    for( i = 0, child = xDataTOME_getFirstElement( element ); child != NULL; i++, child = xDataTOME_getNextElement( child ) ) {
        if( strcmp( child->name, "weighted" ) ) goto err;
        if( MCGIDI_energy_parseWeightFromTOM( smr, child, &(energy->weightedFunctionals.weightedFunctional[i]) ) ) goto err;
        energy->weightedFunctionals.numberOfWeights++;
    }
    energy->type = MCGIDI_energyType_weightedFunctional;
    return( 0 );

err:
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_energy_parseWeightFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_energyWeightedFunctional *weightedFunctional ) {

    xDataTOM_element *child;
    MCGIDI_energy *energy = NULL;
    ptwXYPoints *weight = NULL;
    char const *toUnits[2] = { "MeV", "" };

    if( ( energy = MCGIDI_energy_new( smr ) ) == NULL ) goto err;
    for( child = xDataTOME_getFirstElement( element ); child != NULL; child = xDataTOME_getNextElement( child ) ) {
        if( strcmp( child->name, "weight" ) == 0 ) {
            if( ( weight = MCGIDI_misc_dataFromElement2ptwXYPointsInUnitsOf( smr, child, toUnits ) ) == NULL ) goto err; }
        else if( strcmp( child->name, "evaporation" ) == 0 ) {
            if( MCGIDI_energy_parseEvaporationFromTOM( smr, child, energy ) ) goto err; }
        else {
            smr_setReportError2( smr, smr_unknownID, 1, "unsupported energy type = '%s' in weighted functional", child->name );
            goto err;
        }
    }
    weightedFunctional->weight = weight;
    weightedFunctional->energy = energy;
    return( 0 );

err:
    if( weight != NULL ) ptwXY_free( weight );
    if( energy != NULL ) MCGIDI_energy_free( smr, energy );
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_energy_parseGeneralEvaporationFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_energy *energy ) {

    double norm;
    xDataTOM_element *thetaTOM, *gTOM;
    ptwXYPoints *theta = NULL, *g = NULL;
    char const *toUnits[2] = { "MeV", "MeV" };

    if( ( thetaTOM = xDataTOME_getOneElementByName( smr, element, "theta", 1 ) ) == NULL ) goto err;
    if( ( theta = MCGIDI_misc_dataFromElement2ptwXYPointsInUnitsOf( smr, thetaTOM, toUnits ) ) == NULL ) goto err;

    if( ( gTOM = xDataTOME_getOneElementByName( smr, element, "g", 1 ) ) == NULL ) goto err;
    toUnits[0] = "";
    toUnits[1] = "";
    if( ( g = MCGIDI_misc_dataFromElement2ptwXYPointsInUnitsOf( smr, gTOM, toUnits ) ) == NULL ) goto err;
    if( MCGIDI_fromTOM_pdfOfX( smr, g, &(energy->g), &norm ) ) goto err;
    energy->gInterpolation = ptwXY_getInterpolation( g );
    g = ptwXY_free( g );
    if( std::fabs( 1. - norm ) > 0.001 ) printf( "bad norm = %e\n", norm );

    energy->type = MCGIDI_energyType_generalEvaporation;
    energy->theta = theta;
    return( 0 );

err:
    if( theta != NULL ) ptwXY_free( theta );
    if( g != NULL ) ptwXY_free( g );
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_energy_parseSimpleMaxwellianFissionFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_energy *energy ) {

    char const *U, *toUnits[2] = { "MeV", "MeV" };
    xDataTOM_element *thetaTOM;

    if( ( U = xDataTOM_getAttributesValueInElement( element, "U" ) ) == NULL ) {
        smr_setReportError2( smr, smr_unknownID, 1, "functional form '%s' missing 'U' attribute", element->name );
        goto err;
    }
    if( MCGIDI_misc_PQUStringToDoubleInUnitOf( smr, U, "MeV", &(energy->U) ) != 0 ) goto err;
    if( ( thetaTOM = xDataTOME_getOneElementByName( smr, element, "theta", 1 ) ) == NULL ) goto err;
    if( ( energy->theta = MCGIDI_misc_dataFromElement2ptwXYPointsInUnitsOf( smr, thetaTOM, toUnits ) ) == NULL ) goto err;
    energy->type = MCGIDI_energyType_simpleMaxwellianFission;
    return( 0 );

err:
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_energy_parseEvaporationFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_energy *energy ) {

    char const *U, *toUnits[2] = { "MeV", "MeV" };
    xDataTOM_element *thetaTOM;

    if( ( U = xDataTOM_getAttributesValueInElement( element, "U" ) ) == NULL ) {
        smr_setReportError2( smr, smr_unknownID, 1, "functional form '%s' missing 'U' attribute", element->name );
        goto err;
    }
    if( MCGIDI_misc_PQUStringToDoubleInUnitOf( smr, U, "MeV", &(energy->U) ) != 0 ) goto err;
    if( ( thetaTOM = xDataTOME_getOneElementByName( smr, element, "theta", 1 ) ) == NULL ) goto err;
    if( ( energy->theta = MCGIDI_misc_dataFromElement2ptwXYPointsInUnitsOf( smr, thetaTOM, toUnits ) ) == NULL ) goto err;
    energy->type = MCGIDI_energyType_evaporation;
    return( 0 );

err:
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_energy_parseWattFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_energy *energy ) {

    char const *U, *toUnits[2] = { "MeV", "MeV" };
    xDataTOM_element *aOrBTOM;

    if( ( U = xDataTOM_getAttributesValueInElement( element, "U" ) ) == NULL ) {
        smr_setReportError2( smr, smr_unknownID, 1, "functional form '%s' missing 'U' attribute", element->name );
        goto err;
    }
    if( MCGIDI_misc_PQUStringToDoubleInUnitOf( smr, U, "MeV", &(energy->U) ) != 0 ) goto err;

    if( ( aOrBTOM = xDataTOME_getOneElementByName( smr, element, "a", 1 ) ) == NULL ) goto err;
    if( ( energy->Watt_a = MCGIDI_misc_dataFromElement2ptwXYPointsInUnitsOf( smr, aOrBTOM, toUnits ) ) == NULL ) goto err;

    toUnits[1] = "1/MeV";
    if( ( aOrBTOM = xDataTOME_getOneElementByName( smr, element, "b", 1 ) ) == NULL ) goto err;
    if( ( energy->Watt_b = MCGIDI_misc_dataFromElement2ptwXYPointsInUnitsOf( smr, aOrBTOM, toUnits ) ) == NULL ) goto err;

    energy->type = MCGIDI_energyType_Watt;
    return( 0 );

err:
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_energy_parseMadlandNixFromTOM( statusMessageReporting *smr, xDataTOM_element *functional, MCGIDI_energy *energy ) {

    int iE, length, nXs, i1, n;
    double E=0., T_M=0., EFL=0., EFH=0., argList[3] = { 0., 0., 0. },
           xs[] = { 1e-5, 1e-3, 1e-1, 1e1, 1e3, 1e5, 3e7 }, norm;
    ptwXYPoints *ptwXY_TM = NULL, *pdfXY = NULL;
    ptwXYPoint *point;
    ptwXPoints *cdfX = NULL;
    nfu_status status = nfu_Okay;
    xDataTOM_element *TM_TOM;
    xDataTOM_XYs *XYs;
    MCGIDI_pdfsOfXGivenW *dists = &(energy->dists);
    MCGIDI_pdfOfX *dist;
    char const *EF, *TMUnits[2] = { "MeV", "MeV" };

    nXs = sizeof( xs ) / sizeof( xs[0] );

    if( ( EF = xDataTOM_getAttributesValueInElement( functional, "EFL" ) ) == NULL ) {
        smr_setReportError2( smr, smr_unknownID, 1, "MadlandNix '%s' missing 'EFL' attribute", functional->name );
        goto err;
    }
    if( MCGIDI_misc_PQUStringToDoubleInUnitOf( smr, EF, TMUnits[0], &EFL ) != 0 ) goto err;
    argList[0] = EFL;

    if( ( EF = xDataTOM_getAttributesValueInElement( functional, "EFH" ) ) == NULL ) {
        smr_setReportError2( smr, smr_unknownID, 1, "MadlandNix '%s' missing 'EFH' attribute", functional->name );
        goto err;
    }
    if( MCGIDI_misc_PQUStringToDoubleInUnitOf( smr, EF, TMUnits[0], &EFH ) != 0 ) goto err;
    argList[1] = EFH;

    if( ( TM_TOM = xDataTOME_getOneElementByName( smr, functional, "T_M", 1 ) ) == NULL ) goto err;
    if( ( XYs = (xDataTOM_XYs *) xDataTOME_getXDataIfID( smr, TM_TOM, "XYs" ) ) == NULL ) goto err;
    if( ( ptwXY_TM = MCGIDI_misc_dataFromXYs2ptwXYPointsInUnitsOf( smr, XYs, ptwXY_interpolationLinLin, TMUnits ) ) == NULL ) goto err;

    length = (int) ptwXY_length( ptwXY_TM );
    dists->interpolationWY = ptwXY_interpolationLinLin;
    dists->interpolationXY = ptwXY_interpolationLinLin;     /* Ignoring what the data says as it is probably wrong. */
    if( ( dists->Ws = (double *) smr_malloc2( smr, length * sizeof( double ), 1, "dists->Ws" ) ) == NULL ) goto err;
    if( ( dists->dist = (MCGIDI_pdfOfX *) smr_malloc2( smr, length * sizeof( MCGIDI_pdfOfX ), 0, "dists->dist" ) ) == NULL ) goto err;

    for( iE = 0; iE < length; iE++ ) {
        ptwXY_getXYPairAtIndex( ptwXY_TM, iE, &E, &T_M );
        argList[2] = T_M;
        dist = &(dists->dist[iE]);
        dists->Ws[iE] = E;

        if( ( pdfXY = ptwXY_createFromFunction( nXs, xs, (ptwXY_createFromFunction_callback) MCGIDI_energy_parseMadlandNixFromTOM_callback, 
            (void *) argList, 1e-3, 0, 12, &status ) ) == NULL ) goto err;
        if( ( status = ptwXY_normalize( pdfXY ) ) != nfu_Okay ) {
            smr_setReportError2( smr, smr_unknownID, 1, "ptwXY_normalize err = %d: %s\n", status, nfu_statusMessage( status ) );
            goto err;
        }

        if( ptwXY_simpleCoalescePoints( pdfXY ) != nfu_Okay ) goto err;
        dist->numberOfXs = n = (int) ptwXY_length( pdfXY );

        if( ( dist->Xs = (double *) smr_malloc2( smr, 3 * n * sizeof( double ), 0, "dist->Xs" ) ) == NULL ) goto err;
        dists->numberOfWs++;
        dist->pdf = &(dist->Xs[n]);
        dist->cdf = &(dist->pdf[n]);

        for( i1 = 0; i1 < n; i1++ ) {
            point = ptwXY_getPointAtIndex_Unsafely( pdfXY, i1 );
            dist->Xs[i1] = point->x;
            dist->pdf[i1] = point->y;
        }

        if( ( cdfX = ptwXY_runningIntegral( pdfXY, &status ) ) == NULL ) {
            smr_setReportError2( smr, smr_unknownID, 1, "ptwXY_runningIntegral err = %d: %s\n", status, nfu_statusMessage( status ) );
            goto err;
        }

        norm = ptwX_getPointAtIndex_Unsafely( cdfX, n - 1 );
        for( i1 = 0; i1 < n; i1++ ) dist->cdf[i1] = ptwX_getPointAtIndex_Unsafely( cdfX, i1 ) / norm;
        for( i1 = 0; i1 < n; i1++ ) dist->pdf[i1] /= norm;
        pdfXY = ptwXY_free( pdfXY );
        cdfX = ptwX_free( cdfX );
    }

    energy->type = MCGIDI_energyType_MadlandNix;

    ptwXY_free( ptwXY_TM );
    return( 0 );

err:
    if( ptwXY_TM != NULL ) ptwXY_free( ptwXY_TM );
    if( pdfXY != NULL ) ptwXY_free( pdfXY );
    if( cdfX != NULL ) cdfX = ptwX_free( cdfX );

    return( 1 );
}
/*
************************************************************
*/
static nfu_status MCGIDI_energy_parseMadlandNixFromTOM_callback( double Ep, double *y, void *argList ) {

    double *parameters = (double *) argList, EFL, EFH, T_M;
    nfu_status status = nfu_Okay;

    EFL = parameters[0];
    EFH = parameters[1];
    T_M = parameters[2];
    *y = MCGIDI_energy_parseMadlandNixFromTOM_callback_g( Ep, EFL, T_M, &status );
    if( status == nfu_Okay ) *y += MCGIDI_energy_parseMadlandNixFromTOM_callback_g( Ep, EFH, T_M, &status );
    *y *= 0.5;
    return( status );
}
/*
************************************************************
*/
static double MCGIDI_energy_parseMadlandNixFromTOM_callback_g( double Ep, double E_F, double T_M, nfu_status *status ) {

    double u1, u2, E1, E2 = 0., gamma1 = 0., gamma2 = 0., signG = 1;

    u1 = std::sqrt( Ep ) - std::sqrt( E_F );
    u1 *= u1 / T_M;
    u2 = std::sqrt( Ep ) + std::sqrt( E_F );
    u2 *= u2 / T_M;
    E1 = 0;                      /* u1^3/2 * E1 is zero for u1 = 0. but E1 is infinity, whence, the next test. */
    if( u1 != 0 ) E1 = nf_exponentialIntegral( 1, u1, status );
    if( *status == nfu_Okay ) E2 = nf_exponentialIntegral( 1, u2, status );
    if( *status != nfu_Okay ) return( 0. );
    if( u1 > 2. ) {
        signG = -1;
        gamma1 = nf_incompleteGammaFunctionComplementary( 1.5, u1, status );
        if( *status == nfu_Okay ) gamma2 = nf_incompleteGammaFunctionComplementary( 1.5, u2, status ); }
    else {
        gamma1 = nf_incompleteGammaFunction( 1.5, u1, status );
        if( *status == nfu_Okay ) gamma2 = nf_incompleteGammaFunction( 1.5, u2, status );
    }
    if( *status != nfu_Okay ) return( 0. );
    return( ( u2 * std::sqrt( u2 ) * E2 - u1 * std::sqrt( u1 ) * E1 + signG * ( gamma2 - gamma1 ) ) / ( 3 * std::sqrt( E_F * T_M ) ) );
}
/*
************************************************************
*/
static int MCGIDI_energy_parseNBodyPhaseSpaceFromTOM( statusMessageReporting *smr, xDataTOM_element *functional, MCGIDI_energy *energy,
        MCGIDI_distribution *distribution ) {

    int argList[1];
    double xs[2] = { 0.0, 1.0 }, productMass_MeV, norm;
    ptwXYPoints *pdf = NULL;
    nfu_status status;
    char const *mass;

    if( xDataTOME_convertAttributeToInteger( NULL, functional, "numberOfProducts", &(energy->NBodyPhaseSpace.numberOfProducts) ) != 0 ) goto err;
    if( ( mass = xDataTOM_getAttributesValueInElement( functional, "mass" ) ) == NULL ) {
        smr_setReportError2( smr, smr_unknownID, 1, "functional form '%s' missing 'mass' attribute", functional->name );
        goto err;
    }
    if( MCGIDI_misc_PQUStringToDouble( smr, mass, "amu", MCGIDI_AMU2MeV, &(energy->NBodyPhaseSpace.mass) ) ) goto err;
    argList[0] = energy->NBodyPhaseSpace.numberOfProducts;
    if( ( pdf = ptwXY_createFromFunction( 2, xs, MCGIDI_energy_NBodyPhaseSpacePDF_callback, (void *) argList, 1e-3, 0, 16, &status ) ) == NULL ) {
        smr_setReportError2( smr, smr_unknownID, 1, "creating NBodyPhaseSpace pdf failed with ptwXY_createFromFunction error = %d (%s)", 
            status, nfu_statusMessage( status ) );
        goto err;
    }
    if( MCGIDI_fromTOM_pdfOfX( smr, pdf, &(energy->g), &norm ) ) goto err;
    productMass_MeV = MCGIDI_product_getMass_MeV( smr, distribution->product );
    if( !smr_isOk( smr ) ) goto err;
    energy->NBodyPhaseSpace.massFactor = ( 1. - productMass_MeV / ( MCGIDI_AMU2MeV * energy->NBodyPhaseSpace.mass ) ); /* ??????? Hardwired MCGIDI_AMU2MeV */
    energy->NBodyPhaseSpace.Q_MeV = MCGIDI_outputChannel_getQ_MeV( smr, distribution->product->outputChannel, 0. );
    if( !smr_isOk( smr ) ) goto err;

    ptwXY_free( pdf );
    energy->type = MCGIDI_energyType_NBodyPhaseSpace;

    return( 0 );

err:
    if( pdf != NULL ) ptwXY_free( pdf );
    return( 1 );
}
/*
************************************************************
*/
static nfu_status MCGIDI_energy_NBodyPhaseSpacePDF_callback( double x, double *y, void *argList ) {

    int numberOfProducts = *((int *) argList);
    double e = 0.5 * ( 3 * numberOfProducts - 8 );

    *y = std::sqrt( x ) * G4Pow::GetInstance()->powA( 1.0 - x, e );
    return( nfu_Okay );
}
/*
************************************************************
*/
int MCGIDI_energy_sampleEnergy( statusMessageReporting *smr, MCGIDI_energy *energy, MCGIDI_quantitiesLookupModes &modes,
        MCGIDI_decaySamplingInfo *decaySamplingInfo ) {
/*
*   This function must be called before angular sampling as it sets the frame but does not test it.
*/
    double theta, randomEp, Watt_a, Watt_b, e_in = modes.getProjectileEnergy( );
    MCGIDI_pdfsOfXGivenW_sampled sampled;

    decaySamplingInfo->frame = energy->frame;
    switch( energy->type ) {
    case MCGIDI_energyType_primaryGamma :
        decaySamplingInfo->Ep = energy->gammaEnergy_MeV + e_in * energy->primaryGammaMassFactor;
        break;
    case MCGIDI_energyType_discreteGamma :
        decaySamplingInfo->Ep = energy->gammaEnergy_MeV;
        break;
    case MCGIDI_energyType_linear :
        randomEp = decaySamplingInfo->rng( decaySamplingInfo->rngState );
        sampled.smr = smr;
        sampled.w = e_in;
        MCGIDI_sampling_sampleX_from_pdfsOfXGivenW( &(energy->dists), &sampled, randomEp );
        decaySamplingInfo->Ep = sampled.x;
        break;
    case MCGIDI_energyType_generalEvaporation :
        sampled.interpolationXY = energy->gInterpolation;
        MCGIDI_sampling_sampleX_from_pdfOfX( &(energy->g), &sampled, decaySamplingInfo->rng( decaySamplingInfo->rngState ) );
        theta = MCGIDI_sampling_ptwXY_getValueAtX( energy->theta, e_in );
        decaySamplingInfo->Ep = theta * sampled.x;
        break;
    case MCGIDI_energyType_simpleMaxwellianFission :
        theta = MCGIDI_sampling_ptwXY_getValueAtX( energy->theta, e_in );
        MCGIDI_energy_sampleSimpleMaxwellianFission( smr, ( e_in - energy->U ) / theta, decaySamplingInfo );
        decaySamplingInfo->Ep *= theta;
        break;
    case MCGIDI_energyType_evaporation :
        theta = MCGIDI_sampling_ptwXY_getValueAtX( energy->theta, e_in );
        MCGIDI_energy_sampleEvaporation( smr, ( e_in - energy->U ) / theta, decaySamplingInfo );
        decaySamplingInfo->Ep *= theta;
        break;
    case MCGIDI_energyType_Watt :
        Watt_a = MCGIDI_sampling_ptwXY_getValueAtX( energy->Watt_a, e_in );
        Watt_b = MCGIDI_sampling_ptwXY_getValueAtX( energy->Watt_b, e_in );
        MCGIDI_energy_sampleWatt( smr, e_in - energy->U, Watt_a, Watt_b, decaySamplingInfo );
        break;
    case MCGIDI_energyType_MadlandNix :
        MCGIDI_sampling_sampleX_from_pdfsOfXGivenW( &(energy->dists), &sampled, decaySamplingInfo->rng( decaySamplingInfo->rngState ) );
        decaySamplingInfo->Ep = sampled.x;
        break;
    case MCGIDI_energyType_NBodyPhaseSpace :
        MCGIDI_sampling_sampleX_from_pdfOfX( &(energy->g), &sampled, decaySamplingInfo->rng( decaySamplingInfo->rngState ) );
        decaySamplingInfo->Ep = ( energy->e_inCOMFactor * e_in + energy->NBodyPhaseSpace.Q_MeV ) * energy->NBodyPhaseSpace.massFactor * sampled.x;
        break;
    case MCGIDI_energyType_weightedFunctional :
        MCGIDI_energy_sampleWeightedFunctional( smr, energy, modes, decaySamplingInfo );
        break;
    default :
        smr_setReportError2( smr, smr_unknownID, 1, "energy type = %d not supported", energy->type );
    }

    return( !smr_isOk( smr ) );
}
/*
************************************************************
*/
static int MCGIDI_energy_sampleSimpleMaxwellianFission( statusMessageReporting * /*smr*/, double e_in_U_theta, MCGIDI_decaySamplingInfo *decaySamplingInfo ) {

    int i1;
    double a = e_in_U_theta, b, c, x, norm_a, xMin = 0., xMax = a, sqrt_x, sqrt_pi_2 = std::sqrt( M_PI ) / 2.;

    sqrt_x = std::sqrt( a );
    norm_a = sqrt_pi_2 * erf( sqrt_x ) - sqrt_x * G4Exp( -a );
    b = norm_a * decaySamplingInfo->rng( decaySamplingInfo->rngState );
    for( i1 = 0; i1 < 16; i1++ ) {
        x = 0.5 * ( xMin + xMax );
        sqrt_x = std::sqrt( x );
        c = sqrt_pi_2 * erf( sqrt_x ) - sqrt_x * G4Exp( -x );
        if( b < c ) {
            xMax = x; }
        else {
            xMin = x;
        }
    }
        /* To order e, the correct x is x + e where e = 1 + ( 1 - b * exp( x ) ) / x. */
    decaySamplingInfo->Ep = x;

    return( 0 );
}
/*
************************************************************
*/
static int MCGIDI_energy_sampleEvaporation( statusMessageReporting * /*smr*/, double e_in_U_theta, MCGIDI_decaySamplingInfo *decaySamplingInfo ) {

    int i1;
    double a = e_in_U_theta, b, c, x, norm_a, xMin = 0., xMax = a;

    norm_a = 1 - ( 1 + a ) * G4Exp( -a );
    b = 1. - norm_a * decaySamplingInfo->rng( decaySamplingInfo->rngState );
    for( i1 = 0; i1 < 16; i1++ ) {
        x = 0.5 * ( xMin + xMax );
        c = ( 1 + x ) * G4Exp( -x );
        if( b > c ) {
            xMax = x; }
        else {
            xMin = x;
        }
    }
        /* To order e, the correct x is x + e where e = 1 + ( 1 - b * std::exp( x ) ) / x. */
    decaySamplingInfo->Ep = x;

    return( 0 );
}
/*
************************************************************
*/
static int MCGIDI_energy_sampleWatt( statusMessageReporting * /*smr*/, double e_in_U, double Watt_a, double Watt_b, MCGIDI_decaySamplingInfo *decaySamplingInfo ) {
/*
*   From MCAPM via Sample Watt Spectrum as in TART ( Kalos algorithm ).
*/
    double WattMin = 0., WattMax = e_in_U, x, y, z, energyOut = 0., rand1, rand2;

    x = 1. + ( Watt_b / ( 8. * Watt_a ) );
    y = ( x + std::sqrt( x * x - 1. ) ) / Watt_a;
    z = Watt_a * y - 1.;
   G4int icounter=0;
    G4int icounter_max=1024;
    do
    {
      icounter++;
      if ( icounter > icounter_max ) {
         G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
         break;
      }
        rand1 = -G4Log( decaySamplingInfo->rng( decaySamplingInfo->rngState ) );
        rand2 = -G4Log( decaySamplingInfo->rng( decaySamplingInfo->rngState ) );
        energyOut = y * rand1;
    }
    while( ( ( rand2 - z * ( rand1 + 1. ) ) * ( rand2 - z * ( rand1 + 1. ) ) > Watt_b * y * rand1 ) || ( energyOut < WattMin ) || ( energyOut > WattMax ) ); // Loop checking, 11.06.2015, T. Koi
    decaySamplingInfo->Ep = energyOut;

    return( 0 );
}
/*
************************************************************
*/
static int MCGIDI_energy_sampleWeightedFunctional( statusMessageReporting *smr, MCGIDI_energy *energy, 
        MCGIDI_quantitiesLookupModes &modes, MCGIDI_decaySamplingInfo *decaySamplingInfo ) {
/*
c   This routine assumes that the weights sum to 1.
*/
    int iW;
    double rW = decaySamplingInfo->rng( decaySamplingInfo->rngState ), cumulativeW = 0., weight;
    MCGIDI_energyWeightedFunctional *weightedFunctional = NULL;

    for( iW = 0; iW < energy->weightedFunctionals.numberOfWeights; iW++ ) {
        weightedFunctional = &(energy->weightedFunctionals.weightedFunctional[iW]);
        weight = MCGIDI_sampling_ptwXY_getValueAtX( weightedFunctional->weight, modes.getProjectileEnergy( ) );
        cumulativeW += weight;
        if( cumulativeW >= rW ) break;
    }
    return( MCGIDI_energy_sampleEnergy( smr, weightedFunctional->energy, modes, decaySamplingInfo ) );
}

#if defined __cplusplus
}
#endif

