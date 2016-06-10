/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <string.h>
#include <cmath>

#if defined __cplusplus
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"
namespace GIDI {
using namespace GIDI;
#endif
const double C1 = 0.04, C2 = 1.8e-6/*, C3 = 6.7e-7*/;
/*
const double Et1 = 130., Et3 = 41.;
*/
#if defined __cplusplus
}
#endif

#include "MCGIDI_fromTOM.h"
#include "MCGIDI_misc.h"
#include "MCGIDI_private.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif


static int MCGIDI_KalbachMann_parseFromTOM2( statusMessageReporting *smr, int dataPerEout, int index, xDataTOM_KalbachMannCoefficients *coefficientsXData, 
        double energyInFactor, double energyOutFactor, MCGIDI_KalbachMann *KalbachMann );
static double MCGIDI_KalbachMann_S_a_or_b( double Z_AB, double N_AB, double Z_C, double N_C, double I_ab );
/*
************************************************************
*/
MCGIDI_KalbachMann *MCGIDI_KalbachMann_new( statusMessageReporting *smr, ptwXY_interpolation interpolationWY, 
        ptwXY_interpolation interpolationXY ) {

    MCGIDI_KalbachMann *KalbachMann;

    if( ( KalbachMann = (MCGIDI_KalbachMann *) smr_malloc2( smr, sizeof( MCGIDI_KalbachMann ), 0, "KalbachMann" ) ) == NULL ) return( NULL );
    if( MCGIDI_KalbachMann_initialize( smr, KalbachMann, interpolationWY, interpolationXY ) ) KalbachMann = MCGIDI_KalbachMann_free( smr, KalbachMann );
    return( KalbachMann );
}
/*
************************************************************
*/
int MCGIDI_KalbachMann_initialize( statusMessageReporting * /*smr*/, MCGIDI_KalbachMann *KalbachMann, ptwXY_interpolation interpolationWY, ptwXY_interpolation interpolationXY ) {

    memset( KalbachMann, 0, sizeof( MCGIDI_KalbachMann ) );
    KalbachMann->dists.interpolationWY = interpolationWY;
    KalbachMann->dists.interpolationXY = interpolationXY;
    return( 0 );
}
/*
************************************************************
*/
MCGIDI_KalbachMann *MCGIDI_KalbachMann_free( statusMessageReporting *smr, MCGIDI_KalbachMann *KalbachMann ) {

    MCGIDI_KalbachMann_release( smr, KalbachMann );
    smr_freeMemory( (void **) &KalbachMann );
    return( NULL );
}
/*
************************************************************
*/
int MCGIDI_KalbachMann_release( statusMessageReporting *smr, MCGIDI_KalbachMann *KalbachMann ) {

    int i;
    MCGIDI_pdfsOfXGivenW *dists = &(KalbachMann->dists);

    for( i = 0; i < dists->numberOfWs; i++ ) {
        smr_freeMemory( (void **) &(KalbachMann->ras[i].rs) );
        smr_freeMemory( (void **) &(dists->dist[i].Xs) );
    }
    smr_freeMemory( (void **) &(KalbachMann->ras) );
    smr_freeMemory( (void **) &(dists->Ws) );
    smr_freeMemory( (void **) &(dists->dist) );

    MCGIDI_KalbachMann_initialize( smr, KalbachMann, ptwXY_interpolationLinLin, ptwXY_interpolationLinLin );
    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_KalbachMann_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_distribution *distribution ) {

    MCGIDI_KalbachMann *KalbachMann = NULL;
    xDataTOM_element *KalbachMannElement;
    int index, dataPerEout = 3;
    double energyInFactor, energyOutFactor;
    xDataTOM_xDataInfo *xDataInfo;
    xDataTOM_KalbachMann *KalbachMannXData;
    ptwXY_interpolation interpolationXY, interpolationWY;
    char const *energyFromUnit, *energyToUnit = "MeV";

    MCGIDI_POP *productPOP = distribution->product->pop;
    double productZ = productPOP->Z, productA = productPOP->A, productN = productA - productZ;
    MCGIDI_target_heated *targetHeated = MCGIDI_product_getTargetHeated( smr, distribution->product );
    MCGIDI_POP *projectilePOP = MCGIDI_target_heated_getPOPForProjectile( smr, targetHeated );
    double projectileZ = projectilePOP->Z, projectileA = projectilePOP->A, projectileN = projectileA - projectileZ;
    MCGIDI_POP *targetPOP = MCGIDI_target_heated_getPOPForTarget( smr, targetHeated );
    double targetZ = targetPOP->Z, targetA = targetPOP->A, targetN = targetA - targetZ;
    double Ia = 0., Ib = 0., Ma = -1, mb = -1;

    if( ( targetA == 0 ) && ( targetZ == 6 ) ) {    /* Special case for C_000 evaluation. */
        targetN = 6;
        targetA = 12;
    }
    if( ( KalbachMannElement = xDataTOME_getOneElementByName( smr, element, "KalbachMann", 1 ) ) == NULL ) goto err;

    if( MCGIDI_fromTOM_interpolation( smr, KalbachMannElement, 0, &interpolationWY ) ) goto err;
    if( MCGIDI_fromTOM_interpolation( smr, KalbachMannElement, 1, &interpolationXY ) ) goto err;

    xDataInfo = &(KalbachMannElement->xDataInfo);
    KalbachMannXData = (xDataTOM_KalbachMann *) xDataInfo->data;
    if( KalbachMannXData->type == xDataTOM_KalbachMannType_fra ) dataPerEout = 4;

    energyFromUnit = xDataTOM_axes_getUnit( smr, &(xDataInfo->axes), 0 );
    if( !smr_isOk( smr ) ) goto err;
    energyInFactor = MCGIDI_misc_getUnitConversionFactor( smr, energyFromUnit, energyToUnit );
    if( !smr_isOk( smr ) ) goto err;

    energyFromUnit = xDataTOM_axes_getUnit( smr, &(xDataInfo->axes), 1 );
    if( !smr_isOk( smr ) ) goto err;
    energyOutFactor = MCGIDI_misc_getUnitConversionFactor( smr, energyFromUnit, energyToUnit );
    if( !smr_isOk( smr ) ) goto err;

    if( ( KalbachMann = distribution->KalbachMann = MCGIDI_KalbachMann_new( smr, interpolationWY, interpolationXY ) ) == NULL ) goto err;

/*
    double productMass MCGIDI_product_getMass_MeV( smr, distribution->product ), residualMass;
*/
    KalbachMann->energyToMeVFactor = MCGIDI_misc_getUnitConversionFactor( smr, energyToUnit, "MeV" );
    KalbachMann->massFactor = (double) productZ + productN;     /* This is not correct as masses are needed not Z and N. */
    KalbachMann->massFactor /= projectileN + projectileZ + targetZ + targetN - productZ + productN;
    KalbachMann->massFactor += 1.;

    if( projectileZ == 0 ) {
        if( projectileN == 1 ) Ma = 1; }
    else if( projectileZ == 1 ) {
        if( projectileN == 1 ) {
            Ma = 1; }
        else if( projectileN == 2 ) {
            Ia = 2.22;
            Ma = 1; } }
    else if( projectileZ == 2 ) {
        if( projectileN == 2 ) {
            Ia = 28.3;
            Ma = 0;
        }
    }
        
    if( productZ == 0 ) {
        if( productN == 1 ) mb = 0.5; }
    else if( productZ == 1 ) {
        if( productN == 1 ) {
            mb = 1; }
        else if( productN == 2 ) {
            Ia = 2.22;
            mb = 1; }
        else if( productN == 3 ) {
            Ib = 8.48;
            mb = 1; } }
    else if( productZ == 2 ) {
        if( productN == 1 ) {
            Ib = 7.72;
            mb = 1; }
        else if( productN == 2 ) {
            Ib = 28.3;
            mb = 2;
        }
    }
        
    KalbachMann->Ma = Ma;
    KalbachMann->mb = mb;

    KalbachMann->Sa = MCGIDI_KalbachMann_S_a_or_b( targetZ, targetN, targetZ + projectileZ, targetN + projectileN, Ia );
    KalbachMann->Sb = MCGIDI_KalbachMann_S_a_or_b( projectileZ + targetZ - productZ, projectileN + targetN - productN, 
        targetZ + projectileZ, targetN + projectileN, Ib );

    KalbachMann->dists.numberOfWs = 0;
    if( ( KalbachMann->dists.Ws = (double *) smr_malloc2( smr, KalbachMannXData->numberOfEnergies * sizeof( double ), 0, "KalbachMann->dists->Ws" ) ) == NULL ) goto err;
    if( ( KalbachMann->dists.dist = (MCGIDI_pdfOfX *) smr_malloc2( smr, KalbachMannXData->numberOfEnergies * sizeof( MCGIDI_pdfOfX ), 0, "KalbachMann->dists->dist" ) ) == NULL ) goto err;
    if( ( KalbachMann->ras = (MCGIDI_KalbachMann_ras *) smr_malloc2( smr, KalbachMannXData->numberOfEnergies * sizeof( MCGIDI_KalbachMann_ras ), 0, "KalbachMann->ras" ) ) == NULL ) goto err;

    for( index = 0; index < KalbachMannXData->numberOfEnergies; index++ ) {
        if( MCGIDI_KalbachMann_parseFromTOM2( smr, dataPerEout, index, &(KalbachMannXData->coefficients[index]), 
            energyInFactor, energyOutFactor, KalbachMann ) ) goto err;
    }

    if( ( KalbachMann->frame = MCGIDI_misc_getProductFrame( smr, KalbachMannElement ) ) == xDataTOM_frame_invalid ) goto err;
    distribution->type = MCGIDI_distributionType_KalbachMann_e;

    return( 0 );

err:
    if( KalbachMann != NULL ) MCGIDI_KalbachMann_free( smr, KalbachMann );
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_KalbachMann_parseFromTOM2( statusMessageReporting *smr, int dataPerEout, int index, xDataTOM_KalbachMannCoefficients *coefficientsXData, 
        double energyInFactor, double energyOutFactor, MCGIDI_KalbachMann *KalbachMann ) {

    int i, j, n = coefficientsXData->length / dataPerEout;
    MCGIDI_pdfsOfXGivenW *dists = &(KalbachMann->dists);
    MCGIDI_pdfOfX *dist = &(dists->dist[index]);
    double norm, *p, *rs = NULL, *as_ = NULL, *Xs = NULL, *pdf, *cdf;
    nfu_status status;
    ptwXYPoints *pdfXY = NULL;
    ptwXYPoint *point;
    ptwXPoints *cdfX = NULL;
    char const *ptwFunc = "";

    if( ( Xs = (double *) smr_malloc2( smr, 3 * n * sizeof( double ), 0, "Xs" ) ) == NULL ) goto err;
    pdf = &(Xs[n]);
    cdf = &(pdf[n]);

    if( ( rs = (double *) smr_malloc2( smr, ( dataPerEout - 2 ) * n * sizeof( double ), 0, "rs" ) ) == NULL ) goto err;
    if( dataPerEout == 4 ) as_ = &(rs[n]);

    ptwFunc = "ptwXY_new";
    if( ( pdfXY = ptwXY_new( KalbachMann->dists.interpolationXY, NULL, 2., 1e-3, n, 10, &status, 0 ) ) == NULL ) goto errXY;

    ptwFunc = "ptwXY_setXYPairAtIndex";
    for( i = 0, p = coefficientsXData->coefficients; i < n; i++, p += dataPerEout ) {
        if( ( status = ptwXY_setValueAtX( pdfXY, p[0], p[1] ) ) != nfu_Okay ) goto errXY;
        rs[i] = p[2];
        if( dataPerEout == 4 ) as_[i] = p[3];
    }

    for( j = 0; j < n; j++ ) {
        point = ptwXY_getPointAtIndex_Unsafely( pdfXY, j );
        Xs[j] = energyOutFactor * point->x;
        pdf[j] = point->y / energyOutFactor;
    }

    ptwFunc = "ptwXY_runningIntegral";
    if( ( cdfX = ptwXY_runningIntegral( pdfXY, &status ) ) == NULL ) goto errXY;
    norm = ptwX_getPointAtIndex_Unsafely( cdfX, n - 1 );
    if( std::fabs( 1. - norm ) > 0.99 ) {
        smr_setReportError2( smr, smr_unknownID, 1, "bad norm = %e for angular.linear data", norm );
        goto err;
    }
    for( j = 0; j < n; j++ ) cdf[j] = ptwX_getPointAtIndex_Unsafely( cdfX, j ) / norm;
    for( j = 0; j < n; j++ ) pdf[j] /= norm;

    dists->numberOfWs++;
    dists->Ws[index] = energyInFactor * coefficientsXData->value;
    dist->numberOfXs = n;
    dist->Xs = Xs;
    dist->pdf = pdf;
    dist->cdf = cdf;
    KalbachMann->ras[index].rs = rs;
    KalbachMann->ras[index].as = as_;

    pdfXY = ptwXY_free( pdfXY );
    cdfX = ptwX_free( cdfX );
    return( 0 );
    
errXY:
    smr_setReportError2( smr, smr_unknownID, 1, "%s error = %d: %s\n", ptwFunc, status, nfu_statusMessage( status ) );

err:
    if( Xs != NULL ) smr_freeMemory( (void **) &Xs);
    if( rs != NULL ) smr_freeMemory( (void **) &rs);
    if( pdfXY != NULL ) ptwXY_free( pdfXY );
    if( cdfX != NULL ) cdfX = ptwX_free( cdfX );
    return( 1 );
}
/*
************************************************************
*/
static double MCGIDI_KalbachMann_S_a_or_b( double Z_AB, double N_AB, double Z_C, double N_C, double I_ab ) {

    double A_AB = Z_AB + N_AB, A_C = Z_C + N_C;
    double invA_AB_third = 1.0 / G4Pow::GetInstance()->A13( A_AB ), invA_C_third = 1.0 /  G4Pow::GetInstance()->A13 ( A_C );
    double NZA_AB = ( N_AB - Z_AB ) * ( N_AB - Z_AB ) / A_AB, NZA_C = ( N_C - Z_C ) * ( N_C - Z_C ) / A_C, S;

    S =   15.68 * ( A_C - A_AB )                                             - 28.07 * ( NZA_C - NZA_AB )
        - 18.56 * ( A_C * invA_C_third - A_AB * invA_AB_third )              + 33.22 * ( NZA_C * invA_C_third - NZA_AB * invA_AB_third )
        - 0.717 * ( Z_C * Z_C * invA_C_third - Z_AB * Z_AB * invA_AB_third ) + 1.211 * ( Z_C * Z_C / A_C - Z_AB * Z_AB / A_AB )
        - I_ab;
    return( S );
}
/*
************************************************************
*/
int MCGIDI_KalbachMann_sampleEp( statusMessageReporting *smr, MCGIDI_KalbachMann *KalbachMann, MCGIDI_quantitiesLookupModes &modes, 
        MCGIDI_decaySamplingInfo *decaySamplingInfo ) {

    double Epl, Epu, Ep, r, r2, rl, ru, a, a2, al, au, mu, randomEp = decaySamplingInfo->rng( decaySamplingInfo->rngState );
    MCGIDI_pdfsOfXGivenW_sampled sampled;
    MCGIDI_pdfsOfXGivenW *dists = &(KalbachMann->dists);
    ptwXY_interpolation interpolationWY;

    sampled.smr = smr;
    sampled.w = modes.getProjectileEnergy( );
    MCGIDI_sampling_sampleX_from_pdfsOfXGivenW( dists, &sampled, randomEp );

    interpolationWY = sampled.interpolationWY;
    if( sampled.iW < 0 ) {
        interpolationWY = ptwXY_interpolationFlat;
        if( sampled.iW == -2 ) {        /* ???????????? This should probably report a warning. */
            sampled.iW = 0; }
        else if( sampled.iW == -1 ) {
            sampled.iW = dists->numberOfWs - 1;
        }
    }

    Ep = sampled.x;                                                 /* Sampled Ep. */
    if( sampled.interpolationXY == ptwXY_interpolationFlat ) {      /* Now sample r. */
        r = KalbachMann->ras[sampled.iW].rs[sampled.iX1]; }
    else {
        Epl = dists->dist[sampled.iW].Xs[sampled.iX1];
        Epu = dists->dist[sampled.iW].Xs[sampled.iX1+1];
        rl = KalbachMann->ras[sampled.iW].rs[sampled.iX1];
        ru = KalbachMann->ras[sampled.iW].rs[sampled.iX1+1];
        r = ( ru - rl ) / ( Epu - Epl ) * ( Ep - Epl ) + rl;
    }
    if( interpolationWY == ptwXY_interpolationLinLin ) {
        if( sampled.interpolationXY == ptwXY_interpolationFlat ) {
            r2 = KalbachMann->ras[sampled.iW+1].rs[sampled.iX2]; }
        else {
            Epl = dists->dist[sampled.iW+1].Xs[sampled.iX2];
            Epu = dists->dist[sampled.iW+1].Xs[sampled.iX2+1];
            rl = KalbachMann->ras[sampled.iW+1].rs[sampled.iX2];
            ru = KalbachMann->ras[sampled.iW+1].rs[sampled.iX2+1];
            r2 = ( ru - rl ) / ( Epu - Epl ) * ( Ep - Epl ) + rl;
        }
        r = sampled.frac * r + ( 1. - sampled.frac ) * r2;
    }

    if( KalbachMann->ras[0].as == NULL ) {                          /* Now determine a. */
        double X1, X3_2;
        double eb = KalbachMann->massFactor * KalbachMann->energyToMeVFactor * Ep + KalbachMann->Sb;

        X1 = eb;             /* Not valid for ea > Et1. */
        X3_2 = eb * eb;       /* Not valid for ea > Et3. */
        a = X1 * ( C1 + C2 * X1 * X1 ) + C2 * KalbachMann->Ma * KalbachMann->mb * X3_2 * X3_2; }
    else {
        if( sampled.interpolationXY == ptwXY_interpolationFlat ) {
            a = KalbachMann->ras[sampled.iW].as[sampled.iX1]; }
        else {
            Epl = dists->dist[sampled.iW].Xs[sampled.iX1];
            Epu = dists->dist[sampled.iW].Xs[sampled.iX1+1];
            al = KalbachMann->ras[sampled.iW].as[sampled.iX1];
            au = KalbachMann->ras[sampled.iW].as[sampled.iX1+1];
            a = ( au - al ) / ( Epu - Epl ) * ( Ep - Epl ) + al;
        }
        a2 = 0.;
        if( interpolationWY == ptwXY_interpolationLinLin ) {
            if( sampled.interpolationXY == ptwXY_interpolationFlat ) {
                a2 = KalbachMann->ras[sampled.iW+1].as[sampled.iX2]; }
            else {
                Epl = dists->dist[sampled.iW+1].Xs[sampled.iX2];
                Epu = dists->dist[sampled.iW+1].Xs[sampled.iX2+1];
                al = KalbachMann->ras[sampled.iW+1].as[sampled.iX2];
                au = KalbachMann->ras[sampled.iW+1].as[sampled.iX2+1];
                a2 = ( au - al ) / ( Epu - Epl ) * ( Ep - Epl ) + al;
            }
        }
        a = sampled.frac * a + ( 1. - sampled.frac ) * a2;
    }

        /* In the following: Cosh[ a mu ] + r Sinh[ a mu ] = ( 1 - r ) Cosh[ a mu ] + r ( Cosh[ a mu ] + Sinh[ a mu ] ). */
    if( decaySamplingInfo->rng( decaySamplingInfo->rngState ) >= r ) {   /* Sample the '( 1 - r ) Cosh[ a mu ]' term. */
        double T = ( 2. * decaySamplingInfo->rng( decaySamplingInfo->rngState ) - 1. ) * std::sinh( a );

        mu = G4Log( T + std::sqrt( T * T + 1. ) ) / a; }
    else {                                                              /* Sample the 'r ( Cosh[ a mu ] + Sinh[ a mu ]' term. */
        double rng1 = decaySamplingInfo->rng( decaySamplingInfo->rngState ), exp_a = G4Exp( a );

        mu = G4Log( rng1 * exp_a + ( 1. - rng1 ) / exp_a ) / a;
    }
    if( mu < -1 ) {
        mu = -1;}
    else if( mu >  1 ) {
        mu = 1;
    }

    decaySamplingInfo->frame = KalbachMann->frame;
    decaySamplingInfo->Ep = Ep;
    decaySamplingInfo->mu = mu;
    return( !smr_isOk( smr ) );
}

#if defined __cplusplus
}
#endif
