/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdlib.h>
#include <cmath>
#include <float.h>

#include "ptwXY.h"

#if defined __cplusplus
#include <cmath>
#include "G4Exp.hh"
#include "G4Log.hh"
namespace GIDI {
using namespace GIDI;
#endif

static nfu_status ptwXY_createGaussianCenteredSigma1_2( ptwXYPoints *ptwXY, double x1, double y1, double x2, double y2, int addX1Point );
/*
************************************************************
*/
ptwXPoints *ptwXY_getXArray( ptwXYPoints *ptwXY, nfu_status *status ) {

    int64_t i, n;
    ptwXPoints *xArray;

    if( ( *status = ptwXY->status ) != nfu_Okay ) return( NULL );
    n = ptwXY->length;

    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY ) ) != nfu_Okay ) return( NULL );
    if( ( xArray = ptwX_new( n, status ) ) == NULL ) return( NULL );
    for( i = 0; i < n; i++ ) xArray->points[i] = ptwXY->points[i].x;
    xArray->length = n;

    return( xArray );
}
/*
************************************************************
*/
nfu_status ptwXY_dullEdges( ptwXYPoints *ptwXY, double lowerEps, double upperEps, int positiveXOnly ) {

#define minEps 5e-16

    nfu_status status;
    double xm, xp, dx, y, x1, y1, x2, y2, sign;
    ptwXYPoint *p;

/* This routine can only be used for linear interpolation for the y-axes since for log interpolation, y cannot be 0. 
This needs to be fixed and documented. */
    if( ( status = ptwXY->status ) != nfu_Okay ) return( status );
    if( ptwXY->interpolation == ptwXY_interpolationFlat ) return( nfu_invalidInterpolation );
    if( ptwXY->interpolation == ptwXY_interpolationOther ) return( nfu_otherInterpolation );

    if( ptwXY->length < 2 ) return( nfu_Okay );

    if( lowerEps != 0. ) {
        if( std::fabs( lowerEps ) < minEps ) {
            sign = 1;
            if( lowerEps < 0. ) sign = -1;
            lowerEps = sign * minEps;
        }

        p = ptwXY_getPointAtIndex_Unsafely( ptwXY, 0 );
        x1 = p->x;
        y1 = p->y;
        p = ptwXY_getPointAtIndex_Unsafely( ptwXY, 1 );
        x2 = p->x;
        y2 = p->y;

        if( y1 != 0. ) {
            dx = std::fabs( x1 * lowerEps );
            if( x1 == 0 ) dx = std::fabs( lowerEps );
            xm = x1 - dx;
            xp = x1 + dx;
            if( ( xp + dx ) < x2 ) {
                if( ( status = ptwXY_getValueAtX( ptwXY, xp, &y ) ) != nfu_Okay ) return( status );
                if( ( status = ptwXY_setValueAtX( ptwXY, xp,  y ) ) != nfu_Okay ) return( status ); }
            else {
                xp = x2;
                y = y2;
            }
            if( lowerEps > 0 ) {
                if( ( status = ptwXY_setValueAtX( ptwXY, x1, 0. ) ) != nfu_Okay ) return( status ); }
            else {
                if( ( xm < 0. ) && ( x1 >= 0. ) && positiveXOnly ) {
                    if( ( status = ptwXY_setValueAtX( ptwXY, x1, 0. ) ) != nfu_Okay ) return( status ); }
                else {
                    if( ( status = ptwXY_setValueAtX( ptwXY, xm, 0. ) ) != nfu_Okay ) return( status );
                    if( ( status = ptwXY_interpolatePoint( ptwXY->interpolation, x1, &y, xm, 0., xp, y )  ) != nfu_Okay ) return( status );
                    if( ( status = ptwXY_setValueAtX( ptwXY, x1, y ) ) != nfu_Okay ) return( status );
                }
            }
        }
    }

    if( upperEps != 0. ) {
        if( std::fabs( upperEps ) < minEps ) {
            sign = 1;
            if( upperEps < 0. ) sign = -1;
            upperEps = sign * minEps;
        }

        p = ptwXY_getPointAtIndex_Unsafely( ptwXY, ptwXY->length - 2 );
        x1 = p->x;
        y1 = p->y;
        p = ptwXY_getPointAtIndex_Unsafely( ptwXY, ptwXY->length - 1 );
        x2 = p->x;
        y2 = p->y;

        if( y2 != 0. ) {
            dx = std::fabs( x2 * upperEps );
            if( x2 == 0 ) dx = std::fabs( upperEps );
            xm = x2 - dx;
            xp = x2 + dx;
            if( ( xm - dx ) > x1 ) {
                if( ( status = ptwXY_getValueAtX( ptwXY, xm, &y ) ) != nfu_Okay ) return( status );
                if( ( status = ptwXY_setValueAtX( ptwXY, xm,  y ) ) != nfu_Okay ) return( status ); }
            else {
                xm = x1;
                y = y1;
            }
            if( upperEps < 0 ) {
                if( ( status = ptwXY_setValueAtX( ptwXY, x2, 0. ) ) != nfu_Okay ) return( status ); }
            else {
                if( ( status = ptwXY_setValueAtX( ptwXY, xp, 0. ) ) != nfu_Okay ) return( status );
                if( ( status = ptwXY_interpolatePoint( ptwXY->interpolation, x2, &y, xm, y, xp, 0. )  ) != nfu_Okay ) return( status );
                if( ( status = ptwXY_setValueAtX( ptwXY, x2, y ) ) != nfu_Okay ) return( status );
            }
        }
    }

    return( ptwXY->status );

#undef minEps
}
/*
************************************************************
*/
nfu_status ptwXY_mergeClosePoints( ptwXYPoints *ptwXY, double epsilon ) {

    int64_t i, i1, j, k, n = ptwXY->length;
    double x, y;
    ptwXYPoint *p1, *p2;

    if( n < 2 ) return( ptwXY->status );
    if( epsilon < 4 * DBL_EPSILON ) epsilon = 4 * DBL_EPSILON;
    if( ptwXY_simpleCoalescePoints( ptwXY ) != nfu_Okay ) return( ptwXY->status );

    p2 = ptwXY->points;
    x = p2->x;
    for( i1 = 1, p2++; i1 < ( n - 1 ); i1++, p2++ ) {                 /* The first point shall remain the first point and all points close to it are deleted. */
        if( ( p2->x - x ) > 0.5 * epsilon * ( std::fabs( x ) + std::fabs( p2->x ) ) ) break;
    }
    if( i1 != 1 ) {
        for( i = i1, p1 = &(ptwXY->points[1]); i < n; i++, p1++, p2++ ) *p1 = *p2;
        n = ptwXY->length = ptwXY->length - i1 + 1;
    }

    p1 = &(ptwXY->points[n-1]);
    x = p1->x;
    for( i1 = n - 2, p1--; i1 > 0; i1--, p1-- ) {            /* The last point shall remain the last point and all points close to it are deleted. */
        if( x - p1->x > 0.5 * epsilon * ( std::fabs( x ) + std::fabs( p1->x ) ) ) break;
    }
    if( i1 != ( n - 2 ) ) {
        ptwXY->points[i1 + 1] = ptwXY->points[n - 1];
        n = ptwXY->length = i1 + 2;
    }

    for( i = 1; i < n - 1; i++ ) {
        p1 = &(ptwXY->points[i]);
        x = p1->x;
        y = p1->y;
        for( j = i + 1, p2 = &(ptwXY->points[i+1]); j < n - 1; j++, p2++ ) {
            if( ( p2->x - p1->x ) > 0.5 * epsilon * ( std::fabs( p2->x ) + std::fabs( p1->x ) ) ) break;
            x += p2->x;
            y += p2->y;
        }
        if( ( k = ( j - i ) ) > 1 ) {
            p1->x = x / k;
            p1->y = y / k;
            for( p1 = &(ptwXY->points[i+1]); j < n; j++, p1++, p2++ ) *p1 = *p2;
            n -= ( k - 1 );
        }
    }
    ptwXY->length = n;

    return( ptwXY->status );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_intersectionWith_ptwX( ptwXYPoints *ptwXY, ptwXPoints *ptwX, nfu_status *status ) {

    int64_t i, i1, i2, lengthX = ptwX_length( ptwX );
    double x, y, xMin, xMax;
    ptwXYPoints *n = NULL;

    if( ( *status = ptwXY->status ) != nfu_Okay ) return( NULL );
    if( ( *status = ptwX->status ) != nfu_Okay ) return( NULL );
    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY ) ) != nfu_Okay ) goto Err;
    *status = nfu_otherInterpolation;
    if( ptwXY->interpolation == ptwXY_interpolationOther ) return( NULL );
    if( ( n = ptwXY_clone( ptwXY, status ) ) == NULL ) return( NULL );
    if( ptwXY->length == 0 ) return( n );
    xMin = ptwXY->points[0].x;
    xMax = ptwXY->points[ptwXY->length - 1].x;

    if( ( xMin >= ptwX->points[lengthX-1] ) || ( xMax <= ptwX->points[0] ) ) {  /* No overlap. */
        n->length = 0;
        return( n );
    }

    for( i = 0; i < lengthX; i++ ) {        /* Fill in ptwXY at x-points in ptwX. */
        x = ptwX->points[i];
        if( x <= xMin ) continue;
        if( x >= xMax ) break;
        if( ( *status = ptwXY_getValueAtX( ptwXY, x, &y ) ) != nfu_Okay ) goto Err;
        if( ( *status = ptwXY_setValueAtX( n, x, y ) ) != nfu_Okay ) goto Err;
    }
    if( ( *status = ptwXY_simpleCoalescePoints( n ) ) != nfu_Okay ) goto Err;

    i1 = 0;
    i2 = n->length - 1;
    if( lengthX > 0 ) {
        x = ptwX->points[0];
        if( x > n->points[i1].x ) {
            for( ; i1 < n->length; i1++ ) {
                if( n->points[i1].x == x ) break;
            }
        }

        x = ptwX->points[lengthX - 1];
        if( x < n->points[i2].x ) {
            for( ; i2 > i1; i2-- ) {
                if( n->points[i2].x == x ) break;
            }
        }
    }
    i2++;

    if( i1 != 0 ) {
        for( i = i1; i < i2; i++ ) n->points[i - i1] = n->points[i];
    }
    n->length = i2 - i1;

    return( n );

Err:
     ptwXY_free( n );
     return( NULL );
}
/*
************************************************************
*/
nfu_status ptwXY_areDomainsMutual( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 ) {

    nfu_status status;
    int64_t n1 = ptwXY1->length, n2 = ptwXY2->length;
    ptwXYPoint *xy1, *xy2;

    if( ( status = ptwXY1->status ) != nfu_Okay ) return( status );
    if( ( status = ptwXY2->status ) != nfu_Okay ) return( status );
    if( n1 == 0 ) return( nfu_empty );
    if( n2 == 0 ) return( nfu_empty );
    if( n1 < 2 ) { 
        status = nfu_tooFewPoints; }
    else if( n2 < 2 ) { 
        status = nfu_tooFewPoints; }
    else {
        xy1 = ptwXY_getPointAtIndex_Unsafely( ptwXY1, 0 );
        xy2 = ptwXY_getPointAtIndex_Unsafely( ptwXY2, 0 );
        if( xy1->x < xy2->x ) {
            if( xy2->y != 0. ) status = nfu_domainsNotMutual; }
        else if( xy1->x > xy2->x ) {
            if( xy1->y != 0. ) status = nfu_domainsNotMutual;
        }

        if( status == nfu_Okay ) {
            xy1 = ptwXY_getPointAtIndex_Unsafely( ptwXY1, n1 - 1 );
            xy2 = ptwXY_getPointAtIndex_Unsafely( ptwXY2, n2 - 1 );
            if( xy1->x < xy2->x ) {
                if( xy1->y != 0. ) status = nfu_domainsNotMutual; }
            else if( xy1->x > xy2->x ) {
                if( xy2->y != 0. ) status = nfu_domainsNotMutual;
            }
        }
    }
    return( status );
}
/*
************************************************************
*/
nfu_status ptwXY_tweakDomainsToMutualify( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, int epsilonFactor, double epsilon ) {

    nfu_status status;
    int64_t n1 = ptwXY1->length, n2 = ptwXY2->length;
    double sum, diff;
    ptwXYPoint *xy1, *xy2;

    epsilon = std::fabs( epsilon ) + std::fabs( epsilonFactor * DBL_EPSILON );

    if( ( status = ptwXY1->status ) != nfu_Okay ) return( status );
    if( ( status = ptwXY2->status ) != nfu_Okay ) return( status );
    if( n1 == 0 ) return( nfu_empty );
    if( n2 == 0 ) return( nfu_empty );
    if( n1 < 2 ) { 
        status = nfu_tooFewPoints; }
    else if( n2 < 2 ) { 
        status = nfu_tooFewPoints; }
    else {
        xy1 = ptwXY_getPointAtIndex_Unsafely( ptwXY1, 0 );
        xy2 = ptwXY_getPointAtIndex_Unsafely( ptwXY2, 0 );
        if( xy1->x < xy2->x ) {
            if( xy2->y != 0. ) {
                sum = std::fabs( xy1->x ) + std::fabs( xy2->x );
                diff = std::fabs( xy2->x - xy1->x );
                if( diff > epsilon * sum ) {
                    status = nfu_domainsNotMutual; }
                else {
                    xy1->x = xy2->x;
                }
            } }
        else if( xy1->x > xy2->x ) {
            if( xy1->y != 0. ) {
                sum = std::fabs( xy1->x ) + std::fabs( xy2->x );
                diff = std::fabs( xy2->x - xy1->x );
                if( diff > epsilon * sum ) {
                    status = nfu_domainsNotMutual; }
                else {
                    xy2->x = xy1->x;
                }
            }
        }

        if( status == nfu_Okay ) {
            xy1 = ptwXY_getPointAtIndex_Unsafely( ptwXY1, n1 - 1 );
            xy2 = ptwXY_getPointAtIndex_Unsafely( ptwXY2, n2 - 1 );
            if( xy1->x < xy2->x ) {
                if( xy1->y != 0. ) {
                    sum = std::fabs( xy1->x ) + std::fabs( xy2->x );
                    diff = std::fabs( xy2->x - xy1->x );
                    if( diff > epsilon * sum ) {
                        status = nfu_domainsNotMutual; }
                    else {
                        xy2->x = xy1->x;
                    }
                } }
            else if( xy1->x > xy2->x ) {
                if( xy2->y != 0. ) {
                    sum = std::fabs( xy1->x ) + std::fabs( xy2->x );
                    diff = std::fabs( xy2->x - xy1->x );
                    if( diff > epsilon * sum ) {
                        status = nfu_domainsNotMutual; }
                    else {
                        xy1->x = xy2->x;
                    }
                }
            }
        }
    }
    return( status );
}
/*
************************************************************
*/
nfu_status ptwXY_mutualifyDomains( ptwXYPoints *ptwXY1, double lowerEps1, double upperEps1, int positiveXOnly1,
                                          ptwXYPoints *ptwXY2, double lowerEps2, double upperEps2, int positiveXOnly2 ) {

    nfu_status status;
    int64_t n1 = ptwXY1->length, n2 = ptwXY2->length;
    ptwXYPoint *xy1, *xy2;

    switch( status = ptwXY_areDomainsMutual( ptwXY1, ptwXY2 ) ) {
    case nfu_Okay :
    case nfu_empty :
        return( nfu_Okay );
    case nfu_domainsNotMutual :
        break;
    default :
        return( status );
    }

    if( ptwXY1->interpolation == ptwXY_interpolationOther ) return( nfu_otherInterpolation );
    if( ptwXY2->interpolation == ptwXY_interpolationOther ) return( nfu_otherInterpolation );
    if( ptwXY1->interpolation == ptwXY_interpolationFlat ) return( nfu_invalidInterpolation );
    if( ptwXY2->interpolation == ptwXY_interpolationFlat ) return( nfu_invalidInterpolation );

    xy1 = ptwXY_getPointAtIndex_Unsafely( ptwXY1, 0 );
    xy2 = ptwXY_getPointAtIndex_Unsafely( ptwXY2, 0 );
    if( xy1->x < xy2->x ) {
        lowerEps1 = 0.;
        if( xy2->y == 0. ) lowerEps2 = 0.; }
    else if( xy1->x > xy2->x ) {
        lowerEps2 = 0.;
        if( xy1->y == 0. ) lowerEps1 = 0.; }
    else {
        lowerEps1 = lowerEps2 = 0.;
    }

    xy1 = ptwXY_getPointAtIndex_Unsafely( ptwXY1, n1 - 1 );
    xy2 = ptwXY_getPointAtIndex_Unsafely( ptwXY2, n2 - 1 );
    if( xy1->x < xy2->x ) {
        upperEps2 = 0.;
        if( xy1->y == 0. ) upperEps1 = 0.; }
    else if( xy1->x > xy2->x ) {
        upperEps1 = 0.;
        if( xy2->y == 0. ) upperEps2 = 0.; }
    else {
        upperEps1 = upperEps2 = 0.;
    }

    if( ( lowerEps1 != 0. ) || ( upperEps1 != 0. ) ) 
        if( ( status = ptwXY_dullEdges( ptwXY1, lowerEps1, upperEps1, positiveXOnly1 ) ) != nfu_Okay ) return( status );
    if( ( lowerEps2 != 0. ) || ( upperEps2 != 0. ) ) 
        if( ( status = ptwXY_dullEdges( ptwXY2, lowerEps2, upperEps2, positiveXOnly2 ) ) != nfu_Okay ) return( status );

    return( status );
}
/*
************************************************************
*/
nfu_status ptwXY_copyToC_XY( ptwXYPoints *ptwXY, int64_t index1, int64_t index2, int64_t allocatedSize, int64_t *numberOfPoints, double *xys ) {

    int64_t i;
    double *d = xys;
    nfu_status status;
    ptwXYPoint *pointFrom;

    if( ptwXY->status != nfu_Okay ) return( ptwXY->status );
    if( ( status = ptwXY_simpleCoalescePoints( ptwXY ) ) != nfu_Okay ) return( status );

    if( index1 < 0 ) index1 = 0;
    if( index2 > ptwXY->length ) index2 = ptwXY->length;
    if( index2 < index1 ) index2 = index1;
    *numberOfPoints = index2 - index1;
    if( allocatedSize < ( index2 - index1 ) ) return( nfu_insufficientMemory );

    for( i = index1, pointFrom = ptwXY->points; i < index2; i++, pointFrom++ ) {
        *(d++) = pointFrom->x;
        *(d++) = pointFrom->y;
    }

    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwXY_valueTo_ptwXAndY( ptwXYPoints *ptwXY, double **xs, double **ys ) {

    int64_t i1, length = ptwXY_length( ptwXY );
    double *xps, *yps;
    ptwXYPoint *pointFrom;
    nfu_status status;

    if( ptwXY->status != nfu_Okay ) return( ptwXY->status );
    if( ( status = ptwXY_simpleCoalescePoints( ptwXY ) ) != nfu_Okay ) return( status );

    if( ( *xs = (double *) malloc( length * sizeof( double ) ) ) == NULL ) return( nfu_mallocError );
    if( ( *ys = (double *) malloc( length * sizeof( double ) ) ) == NULL ) {
        free( *xs );
        *xs = NULL;
        return( nfu_mallocError );
    }

    for( i1 = 0, pointFrom = ptwXY->points, xps = *xs, yps = *ys; i1 < length; ++i1, ++pointFrom, ++xps, ++yps ) {
        *xps = pointFrom->x;
        *yps = pointFrom->y;
    }

    return( nfu_Okay );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_valueTo_ptwXY( double x1, double x2, double y, nfu_status *status ) {

    ptwXYPoints *n;

    *status = nfu_XNotAscending;
    if( x1 >= x2 ) return( NULL );
    *status = nfu_Okay;
    if( ( n = ptwXY_new( ptwXY_interpolationLinLin, NULL, ptwXY_maxBiSectionMax, ptwXY_minAccuracy, 2, 0, status, 0 ) ) == NULL ) return( NULL );
    ptwXY_setValueAtX( n, x1, y );
    ptwXY_setValueAtX( n, x2, y );
    return( n );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_createGaussianCenteredSigma1( double accuracy, nfu_status *status ) {

    int64_t i, n;
    ptwXYPoint *pm, *pp;
    double x1, y1, x2, y2, accuracy2, yMin = 1e-10;
    ptwXYPoints *gaussian;

    if( accuracy < 1e-5 ) accuracy = 1e-5;
    if( accuracy > 1e-1 ) accuracy = 1e-1;
    if( ( gaussian = ptwXY_new( ptwXY_interpolationLinLin, NULL, 1., accuracy, 200, 100, status, 0 ) ) == NULL ) return( NULL );
    accuracy2 = accuracy = gaussian->accuracy;
    if( accuracy2 > 5e-3 ) accuracy2 = 5e-3;

    x1 = -std::sqrt( -2. * G4Log( yMin ) );
    y1 = yMin;
    x2 = -5.2;
    y2 = G4Exp( -0.5 * x2 * x2 );
    if( ( *status = ptwXY_setValueAtX( gaussian, x1, y1 ) ) != nfu_Okay ) goto Err;
    gaussian->accuracy = 20 * accuracy2;
    if( ( *status = ptwXY_createGaussianCenteredSigma1_2( gaussian, x1, y1, x2, y2, 1 ) ) != nfu_Okay ) goto Err;
    x1 = x2;
    y1 = y2;
    x2 = -4.;
    y2 = G4Exp( -0.5 * x2 * x2 );
    gaussian->accuracy = 5 * accuracy2;
    if( ( *status = ptwXY_createGaussianCenteredSigma1_2( gaussian, x1, y1, x2, y2, 1 ) ) != nfu_Okay ) goto Err;
    x1 = x2;
    y1 = y2;
    x2 = -1;
    y2 = G4Exp( -0.5 * x2 * x2 );
    gaussian->accuracy = accuracy;
    if( ( *status = ptwXY_createGaussianCenteredSigma1_2( gaussian, x1, y1, x2, y2, 1 ) ) != nfu_Okay ) goto Err;
    x1 = x2;
    y1 = y2;
    x2 =  0;
    y2 = G4Exp( -0.5 * x2 * x2 );
    if( ( *status = ptwXY_createGaussianCenteredSigma1_2( gaussian, x1, y1, x2, y2, 1 ) ) != nfu_Okay ) goto Err;

    n = gaussian->length;
    if( ( *status = ptwXY_coalescePoints( gaussian, 2 * n + 1, NULL, 0 ) ) != nfu_Okay ) goto Err;
    if( ( *status = ptwXY_setValueAtX( gaussian, 0., 1. ) ) != nfu_Okay ) goto Err;
    pp = &(gaussian->points[gaussian->length]);
    for( i = 0, pm = pp - 2; i < n; i++, pp++, pm-- ) {
        *pp = *pm;
        pp->x *= -1;
    }
    gaussian->length = 2 * n + 1;

    return( gaussian );

Err:
    ptwXY_free( gaussian );
    return( NULL );
}
/*
************************************************************
*/
static nfu_status ptwXY_createGaussianCenteredSigma1_2( ptwXYPoints *ptwXY, double x1, double y1, double x2, double y2, int addX1Point ) {

    nfu_status status = nfu_Okay;
    int morePoints = 0;
    double x = 0.5 * ( x1 + x2 );
    double y = G4Exp( -x * x / 2 ), yMin = ( y1 * ( x2 - x ) + y2 * ( x - x1 ) ) / ( x2 - x1 );

    if( std::fabs( y - yMin ) > y * ptwXY->accuracy ) morePoints = 1;
    if( morePoints && ( status = ptwXY_createGaussianCenteredSigma1_2( ptwXY, x, y, x2, y2, 0 ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_setValueAtX( ptwXY, x, y ) ) != nfu_Okay ) return( status );
    if( morePoints && ( status = ptwXY_createGaussianCenteredSigma1_2( ptwXY, x1, y1, x, y, 0 ) ) != nfu_Okay ) return( status );
    if( addX1Point ) status = ptwXY_setValueAtX( ptwXY, x1, y1 );
    return( status );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_createGaussian( double accuracy, double xCenter, double sigma, double amplitude, double xMin, double xMax, 
        double /*dullEps*/, nfu_status *status ) {

    int64_t i;
    ptwXYPoints *gaussian, *sliced;
    ptwXYPoint *point;

    if( ( gaussian = ptwXY_createGaussianCenteredSigma1( accuracy, status ) ) == NULL ) return( NULL );
    for( i = 0, point = gaussian->points; i < gaussian->length; i++, point++ ) {
        point->x = point->x * sigma + xCenter;
        point->y *= amplitude;
    }
    if( ( gaussian->points[0].x < xMin ) || ( gaussian->points[gaussian->length - 1].x > xMax ) ) {
        if( ( sliced = ptwXY_xSlice( gaussian, xMin, xMax, 10, 1, status ) ) == NULL ) goto Err;
        ptwXY_free( gaussian );
        gaussian = sliced;
    }

    return( gaussian );

Err:
    ptwXY_free( gaussian );
    return( NULL );
}

#if defined __cplusplus
}
#endif
