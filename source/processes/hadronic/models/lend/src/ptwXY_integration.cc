/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <cmath>
#include <float.h>

#include "ptwXY.h"
#include <nf_Legendre.h>
#include <nf_integration.h>

#if defined __cplusplus
#include "G4Log.hh"
#include "G4Pow.hh"
namespace GIDI {
using namespace GIDI;
#endif

typedef struct ptwXY_integrateWithFunctionInfo_s {
    int degree;
    ptwXY_createFromFunction_callback func;
    void *argList;
    ptwXY_interpolation interpolation;
    double x1, x2, y1, y2;
} ptwXY_integrateWithFunctionInfo;

static nfu_status ptwXY_integrateWithFunction2( nf_Legendre_GaussianQuadrature_callback integrandFunction, void *argList, double x1,
        double x2, double *integral );
static nfu_status ptwXY_integrateWithFunction3( double x, double *y, void *argList );
/*
************************************************************
*/
nfu_status ptwXY_f_integrate( ptwXY_interpolation interpolation, double x1, double y1, double x2, double y2, double *value ) {

    nfu_status status = nfu_Okay;
    double r;

    *value = 0.;
    switch( interpolation ) {
    case ptwXY_interpolationLinLin :                            /* x linear, y linear */
        *value = 0.5 * ( y1 + y2 ) * ( x2 - x1 );
        break;
    case ptwXY_interpolationLinLog :                            /* x linear, y log */
        if( ( y1 <= 0. ) || ( y2 <= 0. ) ) {
            status = nfu_badIntegrationInput; }
        else {
            r = y2 / y1;
            if( std::fabs( r - 1. ) < 1e-4 ) {
                r = r - 1.;
                *value = y1 * ( x2 - x1 ) / ( 1. + r * ( -0.5 + r * ( 1. / 3. + r * ( -0.25 + .2 * r ) ) ) ); }
            else {
                *value = ( y2 - y1 ) * ( x2 - x1 ) / G4Log( r );
            }
        }
        break;
    case ptwXY_interpolationLogLin :                            /* x log, y linear */
        if( ( x1 <= 0. ) || ( x2 <= 0. ) ) {
            status = nfu_badIntegrationInput; }
        else {
            r = x2 / x1;
            if( std::fabs( r - 1. ) < 1e-4 ) {
                r = r - 1.;
                r = r * ( -0.5 + r * ( 1. / 3. + r * ( -0.25 + .2 * r ) ) );
                *value = x1 * ( y2 - y1 ) * r / ( 1. + r ) + y2 * ( x2 - x1 ); }
            else {
                *value = ( y1 - y2 ) * ( x2 - x1 ) / G4Log( r ) + x2 * y2 - x1 * y1;
            }
        }
        break;
    case ptwXY_interpolationLogLog :                            /* x log, y log */
        if( ( x1 <= 0. ) || ( x2 <= 0. ) || ( y1 <= 0. ) || ( y2 <= 0. ) ) {
            status = nfu_badIntegrationInput; }
        else {
            int i, n;
            double a, z, lx, ly, s, f;

            r = y2 / y1;
            if( std::fabs( r - 1. ) < 1e-4 ) {
                ly = ( y2 - y1 ) / y1;
                ly = ly * ( 1. + ly * ( -0.5 + ly * ( 1. / 3. - 0.25 * ly ) ) ); }
            else {
                ly = G4Log( r );
            }
            r = x2 / x1;
            if( std::fabs( r - 1. ) < 1e-4 ) {
                lx = ( x2 - x1 ) / x1;
                lx = lx * ( 1 + lx * ( -0.5 + lx * ( 1. / 3. - 0.25 * lx ) ) ); }
            else {
                lx = G4Log( r );
            }
            a = ly / lx;
            if( std::fabs( r - 1. ) < 1e-3 ) {
                z = ( x2 - x1 ) / x1;
                n = (int) a;
                if( n > 10 ) n = 12;
                if( n < 4 ) n = 6;
                a = a - n + 1;
                f = n + 1.;
                for( i = 0, s = 0.; i < n; i++, a++, f-- ) s = ( 1. + s ) * a * z / f;
                *value = y1 * ( x2 - x1 ) * ( 1. + s ); }
            else {
                *value = y1 * x1 * ( G4Pow::GetInstance()->powA( r, a + 1. ) - 1. ) / ( a + 1. );
            }
        }
        break;
    case ptwXY_interpolationFlat :                            /* x ?, y flat */
        *value = y1 * ( x2 - x1 );
        break;
    case ptwXY_interpolationOther :
        status = nfu_otherInterpolation;
    }
    return( status );
}
/*
************************************************************
*/
double ptwXY_integrate( ptwXYPoints *ptwXY, double xMin, double xMax, nfu_status *status ) {

    int64_t i, n = ptwXY->length;
    double sum = 0., dSum, x, y, x1, x2, y1, y2, _sign = 1.;
    ptwXYPoint *point;

    if( ( *status = ptwXY->status ) != nfu_Okay ) return( 0. );
    *status = nfu_otherInterpolation;
    if( ptwXY->interpolation == ptwXY_interpolationOther ) return( 0. );

    if( xMax < xMin ) {
        x = xMin;
        xMin = xMax;
        xMax = x;
        _sign = -1.;
    }
    if( n < 2 ) return( 0. );

    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY ) ) != nfu_Okay ) return( 0. );
    for( i = 0, point = ptwXY->points; i < n; i++, point++ ) {
        if( point->x >= xMin ) break;
    }
    if( i == n ) return( 0. );
    x2 = point->x;
    y2 = point->y;
    if( i > 0 ) {
        if( x2 > xMin ) {
            x1 = point[-1].x;
            y1 = point[-1].y;
            if( ( *status = ptwXY_interpolatePoint( ptwXY->interpolation, xMin, &y, x1, y1, x2, y2 ) ) != nfu_Okay ) return( 0. );
            if( x2 > xMax ) {
                double yMax;

                if( ( *status = ptwXY_interpolatePoint( ptwXY->interpolation, xMax, &yMax, x1, y1, x2, y2 ) ) != nfu_Okay ) return( 0. );
                if( ( *status = ptwXY_f_integrate( ptwXY->interpolation, xMin, y, xMax, yMax, &sum ) ) != nfu_Okay ) return( 0. );
                return( sum ); }
            else {
                if( ( *status = ptwXY_f_integrate( ptwXY->interpolation, xMin, y, x2, y2, &sum ) ) != nfu_Okay ) return( 0. );
            }
        }
    }
    i++;
    point++;
    for( ; i < n; i++, point++ ) {
        x1 = x2;
        y1 = y2;
        x2 = point->x;
        y2 = point->y;
        if( x2 > xMax ) {
            if( ( *status = ptwXY_interpolatePoint( ptwXY->interpolation, xMax, &y, x1, y1, x2, y2 ) ) != nfu_Okay ) return( 0. );
            if( ( *status = ptwXY_f_integrate( ptwXY->interpolation, x1, y1, xMax, y, &dSum ) ) != nfu_Okay ) return( 0. );
            sum += dSum;
            break;
        }
        if( ( *status = ptwXY_f_integrate( ptwXY->interpolation, x1, y1, x2, y2, &dSum ) ) != nfu_Okay ) return( 0. );
        sum += dSum;
    }

    return( _sign * sum );
}
/*
************************************************************
*/
double ptwXY_integrateDomain( ptwXYPoints *ptwXY, nfu_status *status ) {

    if( ( *status = ptwXY->status ) != nfu_Okay ) return( 0. );
    if( ptwXY->length > 0 ) return( ptwXY_integrate( ptwXY, ptwXY_getXMin( ptwXY ), ptwXY_getXMax( ptwXY ), status ) );
    return( 0. );
}
/*
************************************************************
*/
nfu_status ptwXY_normalize( ptwXYPoints *ptwXY ) {
/*
*   This function assumes ptwXY_integrateDomain checks status and coalesces the points.
*/

    int64_t i;
    nfu_status status; 
    double sum = ptwXY_integrateDomain( ptwXY, &status );

    if( status != nfu_Okay ) return( status );
    if( sum == 0. ) {
        status = nfu_badNorm; }
    else {
        for( i = 0; i < ptwXY->length; i++ ) ptwXY->points[i].y /= sum;
    }
    return( status );
}
/*
************************************************************
*/
double ptwXY_integrateDomainWithWeight_x( ptwXYPoints *ptwXY, nfu_status *status ) {

    if( ( *status = ptwXY->status ) != nfu_Okay ) return( 0. );
    if( ptwXY->length < 2 ) return( 0. );
    return( ptwXY_integrateWithWeight_x( ptwXY, ptwXY_getXMin( ptwXY ), ptwXY_getXMax( ptwXY ), status ) );
}
/*
************************************************************
*/
double ptwXY_integrateWithWeight_x( ptwXYPoints *ptwXY, double xMin, double xMax, nfu_status *status ) {

    int64_t i, n = ptwXY->length;
    double sum = 0., x, y, x1, x2, y1, y2, _sign = 1.;
    ptwXYPoint *point;

    if( ( *status = ptwXY->status ) != nfu_Okay ) return( 0. );
    *status = nfu_unsupportedInterpolation;
    if( ( ptwXY->interpolation != ptwXY_interpolationLinLin ) && 
        ( ptwXY->interpolation != ptwXY_interpolationFlat ) ) return( 0. );

    if( n < 2 ) return( 0. );
    if( xMax < xMin ) {
        x = xMin;
        xMin = xMax;
        xMax = x;
        _sign = -1.;
    }

    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY ) ) != nfu_Okay ) return( 0. );
    for( i = 0, point = ptwXY->points; i < n; ++i, ++point ) {
        if( point->x >= xMin ) break;
    }
    if( i == n ) return( 0. );
    x2 = point->x;
    y2 = point->y;
    if( i > 0 ) {
        if( x2 > xMin ) {
            if( ( *status = ptwXY_interpolatePoint( ptwXY->interpolation, xMin, &y, point[-1].x, point[-1].y, x2, y2 ) ) != nfu_Okay ) return( 0. );
            x2 = xMin;
            y2 = y;
            --i;
            --point;
        }
    }
    ++i;
    ++point;
    for( ; i < n; ++i, ++point ) {
        x1 = x2;
        y1 = y2;
        x2 = point->x;
        y2 = point->y;
        if( x2 > xMax ) {
            if( ( *status = ptwXY_interpolatePoint( ptwXY->interpolation, xMax, &y, x1, y1, x2, y2 ) ) != nfu_Okay ) return( 0. );
            x2 = xMax;
            y2 = y;
        }
        switch( ptwXY->interpolation ) {
        case ptwXY_interpolationFlat :
            sum += ( x2 - x1 ) * y1 * 3 * ( x1 + x2 );
            break;
        case ptwXY_interpolationLinLin :
            sum += ( x2 - x1 ) * ( y1 * ( 2 * x1 + x2 ) + y2 * ( x1 + 2 * x2 ) );
            break;
        default :       /* Only to stop compilers from complaining. */
            break;
        }
        if( x2 == xMax ) break;
    }

    return( _sign * sum / 6 );
}
/*
************************************************************
*/
double ptwXY_integrateDomainWithWeight_sqrt_x( ptwXYPoints *ptwXY, nfu_status *status ) {

    if( ( *status = ptwXY->status ) != nfu_Okay ) return( 0. );
    if( ptwXY->length < 2 ) return( 0. );
    return( ptwXY_integrateWithWeight_sqrt_x( ptwXY, ptwXY_getXMin( ptwXY ), ptwXY_getXMax( ptwXY ), status ) );
}
/*
************************************************************
*/
double ptwXY_integrateWithWeight_sqrt_x( ptwXYPoints *ptwXY, double xMin, double xMax, nfu_status *status ) {

    int64_t i, n = ptwXY->length;
    double sum = 0., x, y, x1, x2, y1, y2, _sign = 1., sqrt_x1, sqrt_x2, inv_apb, c;
    ptwXYPoint *point;

    if( ( *status = ptwXY->status ) != nfu_Okay ) return( 0. );
    *status = nfu_unsupportedInterpolation;
    if( ( ptwXY->interpolation != ptwXY_interpolationLinLin ) &&
        ( ptwXY->interpolation != ptwXY_interpolationFlat ) ) return( 0. );

    if( n < 2 ) return( 0. );
    if( xMax < xMin ) {
        x = xMin;
        xMin = xMax;
        xMax = x;
        _sign = -1.;
    }

    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY ) ) != nfu_Okay ) return( 0. );
    for( i = 0, point = ptwXY->points; i < n; ++i, ++point ) {
        if( point->x >= xMin ) break;
    }
    if( i == n ) return( 0. );
    x2 = point->x;
    y2 = point->y;
    if( i > 0 ) {
        if( x2 > xMin ) {
            if( ( *status = ptwXY_interpolatePoint( ptwXY->interpolation, xMin, &y, point[-1].x, point[-1].y, x2, y2 ) ) != nfu_Okay ) return( 0. );
            x2 = xMin;
            y2 = y;
            --i;
            --point;
        }
    }
    ++i;
    ++point;
    sqrt_x2 = std::sqrt( x2 );
    for( ; i < n; ++i, ++point ) {
        x1 = x2;
        y1 = y2;
        sqrt_x1 = sqrt_x2;
        x2 = point->x;
        y2 = point->y;
        if( x2 > xMax ) {
            if( ( *status = ptwXY_interpolatePoint( ptwXY->interpolation, xMax, &y, x1, y1, x2, y2 ) ) != nfu_Okay ) return( 0. );
            x2 = xMax;
            y2 = y;
        }
        sqrt_x2 = std::sqrt( x2 );
        inv_apb = sqrt_x1 + sqrt_x2;
        c = 2. * ( sqrt_x1 * sqrt_x2 + x1 + x2 );
        switch( ptwXY->interpolation ) {
        case ptwXY_interpolationFlat :
            sum += ( sqrt_x2 - sqrt_x1 ) * y1 * 2.5 * c;
            break;
        case ptwXY_interpolationLinLin :
            sum += ( sqrt_x2 - sqrt_x1 ) * ( y1 * ( c + x1 * ( 1. + sqrt_x2 / inv_apb ) ) + y2 * ( c + x2 * ( 1. + sqrt_x1 / inv_apb ) ) );
            break;
        default :       /* Only to stop compilers from complaining. */
            break;
        }
        if( x2 == xMax ) break;
    }

    return( 2. / 15. * _sign * sum );
}
/*
************************************************************
*/
ptwXPoints *ptwXY_groupOneFunction( ptwXYPoints *ptwXY, ptwXPoints *groupBoundaries, ptwXY_group_normType normType, ptwXPoints *ptwX_norm, nfu_status *status ) {

    int64_t i, igs, ngs;
    double x1, y1, x2, y2, y2p, xg1, xg2, sum;
    ptwXYPoints *f;
    ptwXPoints *groupedData = NULL;

    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY ) ) != nfu_Okay ) return( NULL );
    if( ( *status = groupBoundaries->status ) != nfu_Okay ) return( NULL );
    *status = nfu_otherInterpolation;
    if( ptwXY->interpolation == ptwXY_interpolationOther ) return( NULL );

    ngs = ptwX_length( groupBoundaries ) - 1;
    if( normType == ptwXY_group_normType_norm ) {
        if( ptwX_norm == NULL ) {
            *status = nfu_badNorm;
            return( NULL );
        }
        *status = ptwX_norm->status;
        if( ptwX_norm->status != nfu_Okay ) return( NULL );
        if( ptwX_length( ptwX_norm ) != ngs ) {
            *status = nfu_badNorm;
            return( NULL );
        }
    }

    if( ( f = ptwXY_intersectionWith_ptwX( ptwXY, groupBoundaries, status ) ) == NULL ) return( NULL );
    if( f->length == 0 ) return( ptwX_createLine( ngs, ngs, 0, 0, status ) );

    if( ( groupedData = ptwX_new( ngs, status ) ) == NULL ) goto err;
    xg1 = groupBoundaries->points[0];
    x1 = f->points[0].x;
    y1 = f->points[0].y;
    for( igs = 0, i = 1; igs < ngs; igs++ ) {
        xg2 = groupBoundaries->points[igs+1];
        sum = 0;
        if( xg2 > x1 ) {
            for( ; i < f->length; i++, x1 = x2, y1 = y2 ) {
                x2 = f->points[i].x;
                if( x2 > xg2 ) break;
                y2p = y2 = f->points[i].y;
                if( f->interpolation == ptwXY_interpolationFlat ) y2p = y1;
                sum += ( y1 + y2p ) * ( x2 - x1 );
            }
        }
        if( sum != 0. ) {
            if( normType == ptwXY_group_normType_dx ) {
                sum /= ( xg2 - xg1 ); }
            else if( normType == ptwXY_group_normType_norm ) {
                if( ptwX_norm->points[igs] == 0. ) {
                    *status = nfu_divByZero;
                    goto err;
                }
                sum /= ptwX_norm->points[igs];
            }
        }
        groupedData->points[igs] = 0.5 * sum;
        groupedData->length++;
        xg1 = xg2;
    }

    ptwXY_free( f );
    return( groupedData );

err:
    ptwXY_free( f );
    if( groupedData != NULL ) ptwX_free( groupedData );
    return( NULL );
}
/*
************************************************************
*/
ptwXPoints *ptwXY_groupTwoFunctions( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, ptwXPoints *groupBoundaries, ptwXY_group_normType normType, 
        ptwXPoints *ptwX_norm, nfu_status *status ) {

    int64_t i, igs, ngs;
    double x1, fy1, gy1, x2, fy2, gy2, fy2p, gy2p, xg1, xg2, sum;
    ptwXYPoints *f = NULL, *ff, *g = NULL, *gg = NULL;
    ptwXPoints *groupedData = NULL;

    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY1 ) ) != nfu_Okay ) return( NULL );
    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY2 ) ) != nfu_Okay ) return( NULL );
    if( ( *status = groupBoundaries->status ) != nfu_Okay ) return( NULL );
    *status = nfu_otherInterpolation;
    if( ptwXY1->interpolation == ptwXY_interpolationOther ) return( NULL );
    if( ptwXY2->interpolation == ptwXY_interpolationOther ) return( NULL );

    ngs = ptwX_length( groupBoundaries ) - 1;
    if( normType == ptwXY_group_normType_norm ) {
        if( ptwX_norm == NULL ) {
            *status = nfu_badNorm;
            return( NULL );
        }
        if( ( *status = ptwX_norm->status ) != nfu_Okay ) return( NULL );
        if( ptwX_length( ptwX_norm ) != ngs ) {
            *status = nfu_badNorm;
            return( NULL );
        }
    }

    if( ( ff = ptwXY_intersectionWith_ptwX( ptwXY1, groupBoundaries, status ) ) == NULL ) return( NULL );
    if( ( gg = ptwXY_intersectionWith_ptwX( ptwXY2, groupBoundaries, status ) ) == NULL ) goto err;
    if( ( ff->length == 0 ) || ( gg->length == 0 ) ) {
        ptwXY_free( ff );
        ptwXY_free( gg );
        return( ptwX_createLine( ngs, ngs, 0, 0, status ) );
    }

    if( ( *status = ptwXY_tweakDomainsToMutualify( ff, gg, 4, 0 ) ) != nfu_Okay ) goto err;
    if( ( f = ptwXY_union( ff, gg, status, ptwXY_union_fill ) ) == NULL ) goto err;
    if( ( g = ptwXY_union( gg, f, status, ptwXY_union_fill ) ) == NULL ) goto err;

    if( ( groupedData = ptwX_new( ngs, status ) ) == NULL ) goto err;
    xg1 = groupBoundaries->points[0];
    x1 = f->points[0].x;
    fy1 = f->points[0].y;
    gy1 = g->points[0].y;
    for( igs = 0, i = 1; igs < ngs; igs++ ) {
        xg2 = groupBoundaries->points[igs+1];
        sum = 0;
        if( xg2 > x1 ) {
            for( ; i < f->length; i++, x1 = x2, fy1 = fy2, gy1 = gy2 ) {
                x2 = f->points[i].x;
                if( x2 > xg2 ) break;
                fy2p = fy2 = f->points[i].y;
                if( f->interpolation == ptwXY_interpolationFlat ) fy2p = fy1;
                gy2p = gy2 = g->points[i].y;
                if( g->interpolation == ptwXY_interpolationFlat ) gy2p = gy1;
                sum += ( ( fy1 + fy2p ) * ( gy1 + gy2p ) + fy1 * gy1 + fy2p * gy2p ) * ( x2 - x1 );
            }
        }
        if( sum != 0. ) {
            if( normType == ptwXY_group_normType_dx ) {
                sum /= ( xg2 - xg1 ); }
            else if( normType == ptwXY_group_normType_norm ) {
                if( ptwX_norm->points[igs] == 0. ) {
                    *status = nfu_divByZero;
                    goto err;
                }
                sum /= ptwX_norm->points[igs];
            }
        }
        groupedData->points[igs] = sum / 6.;
        groupedData->length++;
        xg1 = xg2;
    }

    ptwXY_free( f );
    ptwXY_free( g );
    ptwXY_free( ff );
    ptwXY_free( gg );
    return( groupedData );

err:
    ptwXY_free( ff );
    if( gg != NULL ) ptwXY_free( gg );
    // Coverity #63063
    if( f != NULL ) ptwXY_free( f );
    if( g != NULL ) ptwXY_free( g );
    if( groupedData != NULL ) ptwX_free( groupedData );
    return( NULL );
}
/*
************************************************************
*/
ptwXPoints *ptwXY_groupThreeFunctions( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, ptwXYPoints *ptwXY3, ptwXPoints *groupBoundaries, 
        ptwXY_group_normType normType, ptwXPoints *ptwX_norm, nfu_status *status ) {

    int64_t i, igs, ngs;
    double x1, fy1, gy1, hy1, x2, fy2, gy2, hy2, fy2p, gy2p, hy2p, xg1, xg2, sum;
    ptwXYPoints *f = NULL, *ff, *fff = NULL, *g = NULL, *gg = NULL, *h = NULL, *hh = NULL;
    ptwXPoints *groupedData = NULL;

    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY1 ) ) != nfu_Okay ) return( NULL );
    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY2 ) ) != nfu_Okay ) return( NULL );
    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY3 ) ) != nfu_Okay ) return( NULL );
    if( ( *status = groupBoundaries->status ) != nfu_Okay ) return( NULL );
    *status = nfu_otherInterpolation;
    if( ptwXY1->interpolation == ptwXY_interpolationOther ) return( NULL );
    if( ptwXY2->interpolation == ptwXY_interpolationOther ) return( NULL );
    if( ptwXY3->interpolation == ptwXY_interpolationOther ) return( NULL );

    ngs = ptwX_length( groupBoundaries ) - 1;
    if( normType == ptwXY_group_normType_norm ) {
        if( ptwX_norm == NULL ) {
            *status = nfu_badNorm;
            return( NULL );
        }
        if( ( *status = ptwX_norm->status ) != nfu_Okay ) return( NULL );
        if( ptwX_length( ptwX_norm ) != ngs ) {
            *status = nfu_badNorm;
            return( NULL );
        }
    }

    if( ( ff = ptwXY_intersectionWith_ptwX( ptwXY1, groupBoundaries, status ) ) == NULL ) return( NULL );
    if( ( gg = ptwXY_intersectionWith_ptwX( ptwXY2, groupBoundaries, status ) ) == NULL ) goto err;
    if( ( hh = ptwXY_intersectionWith_ptwX( ptwXY3, groupBoundaries, status ) ) == NULL ) goto err;
    if( ( ff->length == 0 ) || ( gg->length == 0 ) || ( hh->length == 0 ) ) return( ptwX_createLine( ngs, ngs, 0, 0, status ) );

    if( ( *status = ptwXY_tweakDomainsToMutualify( ff, gg, 4, 0 ) ) != nfu_Okay ) goto err;
    if( ( *status = ptwXY_tweakDomainsToMutualify( ff, hh, 4, 0 ) ) != nfu_Okay ) goto err;
    if( ( *status = ptwXY_tweakDomainsToMutualify( gg, hh, 4, 0 ) ) != nfu_Okay ) goto err;
    if( ( fff = ptwXY_union(  ff,  gg, status, ptwXY_union_fill ) ) == NULL ) goto err;
    if( (   h = ptwXY_union(  hh, fff, status, ptwXY_union_fill ) ) == NULL ) goto err;
    if( (   f = ptwXY_union( fff,   h, status, ptwXY_union_fill ) ) == NULL ) goto err;
    if( (   g = ptwXY_union(  gg,   h, status, ptwXY_union_fill ) ) == NULL ) goto err;

    if( ( groupedData = ptwX_new( ngs, status ) ) == NULL ) goto err;
    xg1 = groupBoundaries->points[0];
    x1 = f->points[0].x;
    fy1 = f->points[0].y;
    gy1 = g->points[0].y;
    hy1 = h->points[0].y;
    for( igs = 0, i = 1; igs < ngs; igs++ ) {
        xg2 = groupBoundaries->points[igs+1];
        sum = 0;
        if( xg2 > x1 ) {
            for( ; i < f->length; i++, x1 = x2, fy1 = fy2, gy1 = gy2, hy1 = hy2 ) {
                x2 = f->points[i].x;
                if( x2 > xg2 ) break;
                fy2p = fy2 = f->points[i].y;
                if( f->interpolation == ptwXY_interpolationFlat ) fy2p = fy1;
                gy2p = gy2 = g->points[i].y;
                if( g->interpolation == ptwXY_interpolationFlat ) gy2p = gy1;
                hy2p = hy2 = h->points[i].y;
                if( h->interpolation == ptwXY_interpolationFlat ) hy2p = hy1;
                sum += ( ( fy1 + fy2p ) * ( gy1 + gy2p ) * ( hy1 + hy2p ) + 2 * fy1 * gy1 * hy1 + 2 * fy2p * gy2p * hy2p ) * ( x2 - x1 );
            }
        }
        if( sum != 0. ) {
            if( normType == ptwXY_group_normType_dx ) {
                sum /= ( xg2 - xg1 ); }
            else if( normType == ptwXY_group_normType_norm ) {
                if( ptwX_norm->points[igs] == 0. ) {
                    *status = nfu_divByZero;
                    goto err;
                }
                sum /= ptwX_norm->points[igs];
            }
        }
        groupedData->points[igs] = sum / 12.;
        groupedData->length++;
        xg1 = xg2;
    }

    ptwXY_free( f );
    ptwXY_free( g );
    ptwXY_free( h );
    ptwXY_free( ff );
    ptwXY_free( gg );
    ptwXY_free( hh );
    ptwXY_free( fff );
    return( groupedData );

err:
    ptwXY_free( ff );
    if( fff != NULL ) ptwXY_free( fff );
    if( gg != NULL ) ptwXY_free( gg );
    if( hh != NULL ) ptwXY_free( hh );
    if( f != NULL ) ptwXY_free( f );
    if( g != NULL ) ptwXY_free( g );
    if( h != NULL ) ptwXY_free( h );
    if( groupedData != NULL ) ptwX_free( groupedData );
    return( NULL );
}
/*
************************************************************
*/
ptwXPoints *ptwXY_runningIntegral( ptwXYPoints *ptwXY, nfu_status *status ) {

    int i;
    ptwXPoints *runningIntegral = NULL;
    double integral = 0., sum;

    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY ) ) != nfu_Okay ) return( NULL );
    if( ( runningIntegral = ptwX_new( ptwXY->length, status ) ) == NULL ) goto err;

    if( ( *status = ptwX_setPointAtIndex( runningIntegral, 0, 0. ) ) != nfu_Okay ) goto err;
    for( i = 1; i < ptwXY->length; i++ ) {
        if( ( *status = ptwXY_f_integrate( ptwXY->interpolation, ptwXY->points[i-1].x, ptwXY->points[i-1].y, 
            ptwXY->points[i].x, ptwXY->points[i].y, &sum ) ) != nfu_Okay ) goto err;
        integral += sum;
        if( ( *status = ptwX_setPointAtIndex( runningIntegral, i, integral ) ) != nfu_Okay ) goto err;
    }
    return( runningIntegral );

err:
    if( runningIntegral != NULL ) ptwX_free( runningIntegral );
    return( NULL );
}
/*
************************************************************
*/
double ptwXY_integrateWithFunction( ptwXYPoints *ptwXY, ptwXY_createFromFunction_callback func, void *argList,
        double xMin, double xMax, int degree, int recursionLimit, double tolerance, nfu_status *status ) {

    int64_t i1, i2, n1 = ptwXY->length;
    long evaluations;
    double integral = 0., integral_, sign = -1., xa, xb;
    ptwXY_integrateWithFunctionInfo integrateWithFunctionInfo;
    ptwXYPoint *point;

    if( ( *status = ptwXY->status ) != nfu_Okay ) return( 0. );

    if( xMin == xMax ) return( 0. );
    if( n1 < 2 ) return( 0. );

    ptwXY_simpleCoalescePoints( ptwXY );

    if( xMin > xMax ) {
        sign = xMin;
        xMin = xMax;
        xMax = sign;
        sign = -1.;
    }
    if( xMin >= ptwXY->points[n1-1].x ) return( 0. );
    if( xMax <= ptwXY->points[0].x ) return( 0. );

    for( i1 = 0; i1 < ( n1 - 1 ); i1++ ) {
        if( ptwXY->points[i1+1].x > xMin ) break;
    }
    for( i2 = n1 - 1; i2 > i1; i2-- ) {
        if( ptwXY->points[i2-1].x < xMax ) break;
    }
    point = &(ptwXY->points[i1]);

    integrateWithFunctionInfo.degree = degree;
    integrateWithFunctionInfo.func = func;
    integrateWithFunctionInfo.argList = argList;
    integrateWithFunctionInfo.interpolation = ptwXY->interpolation;
    integrateWithFunctionInfo.x2 = point->x;
    integrateWithFunctionInfo.y2 = point->y;

    xa = xMin;
    for( ; i1 < i2; i1++ ) {
        integrateWithFunctionInfo.x1 = integrateWithFunctionInfo.x2;
        integrateWithFunctionInfo.y1 = integrateWithFunctionInfo.y2;
        ++point;
        integrateWithFunctionInfo.x2 = point->x;
        integrateWithFunctionInfo.y2 = point->y;
        xb = point->x;
        if( xb > xMax ) xb = xMax;
        *status = nf_GnG_adaptiveQuadrature( ptwXY_integrateWithFunction2, ptwXY_integrateWithFunction3, &integrateWithFunctionInfo,
            xa, xb, recursionLimit, tolerance, &integral_, &evaluations );
        if( *status != nfu_Okay ) return( 0. );
        integral += integral_;
        xa = xb;
    }

    return( integral );
}
/*
************************************************************
*/
static nfu_status ptwXY_integrateWithFunction2( nf_Legendre_GaussianQuadrature_callback integrandFunction, void *argList, double x1,
        double x2, double *integral ) {

    ptwXY_integrateWithFunctionInfo *integrateWithFunctionInfo = (ptwXY_integrateWithFunctionInfo *) argList;
    nfu_status status;

    status = nf_Legendre_GaussianQuadrature( integrateWithFunctionInfo->degree, x1, x2, integrandFunction, argList, integral );
    return( status );
}
/*
************************************************************
*/
static nfu_status ptwXY_integrateWithFunction3( double x, double *y, void *argList ) {

    double yf;
    ptwXY_integrateWithFunctionInfo *integrateWithFunctionInfo = (ptwXY_integrateWithFunctionInfo *) argList;
    nfu_status status;

    if( ( status = ptwXY_interpolatePoint( integrateWithFunctionInfo->interpolation, x, &yf, 
            integrateWithFunctionInfo->x1, integrateWithFunctionInfo->y1, 
            integrateWithFunctionInfo->x2, integrateWithFunctionInfo->y2 ) ) == nfu_Okay ) {
        status = integrateWithFunctionInfo->func( x, y, integrateWithFunctionInfo->argList );
        *y *= yf;
    }
    return( status );
}

#if defined __cplusplus
}
#endif
