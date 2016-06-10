/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <float.h>

#include "nf_integration.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

typedef struct nf_GnG_adaptiveQuadrature_info_s {
    nfu_status status;
    nf_Legendre_GaussianQuadrature_callback integrandFunction;
    void *argList;
    nf_GnG_adaptiveQuadrature_callback quadratureFunction;
    double estimate;
    int evaluations, maxDepth, maxDepthReached;
} nf_GnG_adaptiveQuadrature_info;

static double initialPoints[] = { 0.2311, 0.4860, 0.6068, 0.8913, 0.9501 };
static int numberOfInitialPoints = sizeof( initialPoints ) / sizeof( initialPoints[0] );

static double nf_GnG_adaptiveQuadrature2( nf_GnG_adaptiveQuadrature_info *adaptiveQuadrature_info, double currentIntrgral, double x1, double x2, int depth );
/*
============================================================
*/
nfu_status nf_GnG_adaptiveQuadrature( nf_GnG_adaptiveQuadrature_callback quadratureFunction, nf_Legendre_GaussianQuadrature_callback integrandFunction, 
    void *argList, double x1, double x2, int maxDepth, double tolerance, double *integral, long *evaluations ) {
/*
*   See W. Gander and W. Gautschi, "Adaptive quadrature--revisited", BIT 40 (2000), 84-101.
*/
    int i1;
    double estimate = 0., y1, integral_, coarse;
    nfu_status status = nfu_Okay;
    nf_GnG_adaptiveQuadrature_info adaptiveQuadrature_info = { nfu_Okay, integrandFunction, argList, quadratureFunction, 0., 0, maxDepth, 0 };

    *integral = 0.;
    *evaluations = 0;
    if( x1 == x2 ) return( nfu_Okay );

    if( tolerance < 10 * DBL_EPSILON ) tolerance = 10 * DBL_EPSILON;
    if( maxDepth > nf_GnG_adaptiveQuadrature_MaxMaxDepth ) maxDepth = nf_GnG_adaptiveQuadrature_MaxMaxDepth;

    for( i1 = 0; i1 < numberOfInitialPoints; i1++ ) {
        if( ( status = integrandFunction( x1 + ( x2 - x1 ) * initialPoints[i1], &y1, argList ) ) != nfu_Okay ) return( status );
        estimate += y1;
    }
    if( ( status = quadratureFunction( integrandFunction, argList, x1, x2, &integral_ ) ) != nfu_Okay ) return( status );
    estimate = 0.5 * ( estimate * ( x2 - x1 ) / numberOfInitialPoints + integral_ );
    if( estimate == 0. ) estimate = x2 - x1;
    adaptiveQuadrature_info.estimate = tolerance * estimate / DBL_EPSILON;

    if( ( status = quadratureFunction( integrandFunction, argList, x1, x2, &coarse ) ) != nfu_Okay ) return( status );
    integral_ = nf_GnG_adaptiveQuadrature2( &adaptiveQuadrature_info, coarse, x1, x2, 0 );

    for( i1 = 0; i1 < 2; i1++ ) {       /* Estimate may be off by more than a factor of 10. Iterate at most 2 times. */
        if( integral_ == 0. ) break;
        y1 = integral_ / estimate;
        if( ( y1 > 0.1 ) && ( y1 < 10. ) ) break;

        estimate = integral_;
        adaptiveQuadrature_info.estimate = tolerance * integral_ / DBL_EPSILON;
        *evaluations += adaptiveQuadrature_info.evaluations;
        adaptiveQuadrature_info.evaluations = 0;
        integral_ = nf_GnG_adaptiveQuadrature2( &adaptiveQuadrature_info, integral_, x1, x2, 0 );
    }

    *evaluations += adaptiveQuadrature_info.evaluations;
    if( adaptiveQuadrature_info.status == nfu_Okay ) *integral = integral_;
    return( adaptiveQuadrature_info.status );
}
/*
============================================================
*/
static double nf_GnG_adaptiveQuadrature2( nf_GnG_adaptiveQuadrature_info *adaptiveQuadrature_info, double coarse, double x1, double x2, int depth ) {

    double xm, intregral1, intregral2, fine, extrapolate;

    if( adaptiveQuadrature_info->status != nfu_Okay ) return( 0. );
    if( x1 == x2 ) return( 0. );

    adaptiveQuadrature_info->evaluations++;
    depth++;
    if( depth > adaptiveQuadrature_info->maxDepthReached ) adaptiveQuadrature_info->maxDepthReached = depth;

    xm = 0.5 * ( x1 + x2 );
    if( ( adaptiveQuadrature_info->status = adaptiveQuadrature_info->quadratureFunction( adaptiveQuadrature_info->integrandFunction, 
        adaptiveQuadrature_info->argList, x1, xm, &intregral1 ) ) != nfu_Okay ) return( 0. );
    if( ( adaptiveQuadrature_info->status = adaptiveQuadrature_info->quadratureFunction( adaptiveQuadrature_info->integrandFunction, 
        adaptiveQuadrature_info->argList, xm, x2, &intregral2 ) ) != nfu_Okay ) return( 0. );
    fine = intregral1 + intregral2;
    extrapolate = ( 16. * fine - coarse ) / 15.;
    if( extrapolate != 0 ) {
        if( adaptiveQuadrature_info->estimate + ( extrapolate - fine ) == adaptiveQuadrature_info->estimate ) return( fine );
    }
    if( depth > adaptiveQuadrature_info->maxDepth ) return( fine );
    return( nf_GnG_adaptiveQuadrature2( adaptiveQuadrature_info, intregral1, x1, xm, depth ) + 
            nf_GnG_adaptiveQuadrature2( adaptiveQuadrature_info, intregral2, xm, x2, depth ) );
}

#if defined __cplusplus
}
#endif
