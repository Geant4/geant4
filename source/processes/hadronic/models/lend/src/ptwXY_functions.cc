/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <cmath>

#include "ptwXY.h"

#if defined __cplusplus
#include "G4Exp.hh"
#include "G4Pow.hh"
namespace GIDI {
using namespace GIDI;
#endif

static nfu_status ptwXY_pow_callback( ptwXYPoint *point, void *argList );
static nfu_status ptwXY_exp_s( ptwXYPoints *ptwXY, double x1, double y1, double z1, double x2, double y2, double z2, int level );
static nfu_status ptwXY_convolution2( ptwXYPoints *f1, ptwXYPoints *f2, double y, double yMin, double *c );
static nfu_status ptwXY_convolution3( ptwXYPoints *convolute, ptwXYPoints *f1, ptwXYPoints *f2, double y1, double c1, double y2, double c2, double yMin );
/*
************************************************************
*/
nfu_status ptwXY_pow( ptwXYPoints *ptwXY, double v ) { 

    return( ptwXY_applyFunction( ptwXY, ptwXY_pow_callback, (void *) &v, 0 ) );
}
/*
************************************************************
*/
static nfu_status ptwXY_pow_callback( ptwXYPoint *point, void *argList ) {

    nfu_status status = nfu_Okay;
    double *v = (double *) argList;

    point->y = G4Pow::GetInstance()->powA( point->y, *v );
    /* ???? Need to test for valid y-value. */
    return( status );
}
/*
************************************************************
*/
nfu_status ptwXY_exp( ptwXYPoints *ptwXY, double a ) { 

    int64_t i, length;
    nfu_status status;
    double x1, y1, z1, x2, y2, z2;

    length = ptwXY->length;
    if( length < 1 ) return( ptwXY->status );
    if( ptwXY->interpolation == ptwXY_interpolationFlat ) return( nfu_invalidInterpolation );
    if( ptwXY->interpolation == ptwXY_interpolationOther ) return( nfu_otherInterpolation );

    if( ( status = ptwXY_simpleCoalescePoints( ptwXY ) ) != nfu_Okay ) return( status );
    x2 = ptwXY->points[length-1].x;
    y2 = a * ptwXY->points[length-1].y;
    z2 = ptwXY->points[length-1].y = G4Exp( y2 );
    for( i = length - 2; i >= 0; i-- ) {
        x1 = ptwXY->points[i].x;
        y1 = a * ptwXY->points[i].y;
        z1 = ptwXY->points[i].y = G4Exp( y1 );
        if( ( status = ptwXY_exp_s( ptwXY, x1, y1, z1, x2, y2, z2, 0 ) ) != nfu_Okay ) return( status );
        x2 = x1;
        y2 = y1;
    }
    return( status );
}
/*
************************************************************
*/
static nfu_status ptwXY_exp_s( ptwXYPoints *ptwXY, double x1, double y1, double z1, double x2, double y2, double z2, int level ) { 

    nfu_status status;
    double x, y, dx, dy, z, zp, s;

    if( ( x1 == x2 ) || ( y1 == y2 ) ) return( nfu_Okay );
    if( level >= ptwXY->biSectionMax ) return( nfu_Okay );
    level++;
    dx = x2 - x1;
    dy = y2 - y1;
    s = dy / dx;
    x = 1. / s + x2 - z2 * dx / ( z2 - z1 );
    y = ( y1 * ( x2 - x ) + y2 * ( x - x1 ) ) / dx;
    z = z1 * G4Exp( 1 - dy / ( G4Exp( dy ) - 1 ) );
    zp = ( z2 - z1 ) / ( y2 - y1 );

    if( std::fabs( z - zp ) < std::fabs( z * ptwXY->accuracy ) ) return( nfu_Okay );
    if( ( status = ptwXY_setValueAtX( ptwXY, x, z ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_exp_s( ptwXY, x, y, z, x2, y2, z2, level ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_exp_s( ptwXY, x1, y1, z1, x, y, z, level ) ) != nfu_Okay ) return( status );
    return( status );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_convolution( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status, int mode ) {
/*
*   Currently, only supports linear-linear interpolation.
*
*   This function calculates        c(y) = integral dx f1(x) * f2(y-x)
*
*/
    int64_t i1, i2, n1, n2, n;
    ptwXYPoints *f1 = ptwXY1, *f2 = ptwXY2, *convolute;
    double accuracy = ptwXY1->accuracy, yMin, yMax, c, y, dy;

    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY1 ) ) != nfu_Okay ) return( NULL );
    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY2 ) ) != nfu_Okay ) return( NULL );

    *status = nfu_unsupportedInterpolation;
    if( ( ptwXY1->interpolation != ptwXY_interpolationLinLin ) || ( ptwXY2->interpolation != ptwXY_interpolationLinLin ) ) return( NULL );
    *status = nfu_Okay;

    n1 = f1->length;
    n2 = f2->length;

    if( ( n1 == 0 ) || ( n2 == 0 ) ) {
        convolute = ptwXY_new( ptwXY_interpolationLinLin, NULL, 1., accuracy, 0, 0, status, 0 );
        return( convolute );
    }

    if( ( n1 == 1 ) || ( n2 == 1 ) ) {
        *status = nfu_tooFewPoints;
        return( NULL );
    }

    if( accuracy < ptwXY2->accuracy ) accuracy = ptwXY2->accuracy;
    n = n1 * n2;
    if( mode == 0 ) {
        mode = 1;
        if( n > 1000 ) mode = -1;
    }
    if( n > 100000 ) mode = -1;
    if( ( convolute = ptwXY_new( ptwXY_interpolationLinLin, NULL, 1., accuracy, 400, 40, status, 0 ) ) == NULL ) return( NULL );

    yMin = f1->points[0].x + f2->points[0].x;
    yMax = f1->points[n1 - 1].x + f2->points[n2 - 1].x;

    if( ( *status = ptwXY_setValueAtX( convolute, yMin, 0. ) ) != nfu_Okay ) goto Err;

    if( mode < 0 ) {
        dy = ( yMax - yMin ) / 400;
        for( y = yMin + dy; y < yMax; y += dy ) {
            if( ( *status = ptwXY_convolution2( f1, f2, y, yMin, &c ) ) != nfu_Okay ) goto Err;
            if( ( *status = ptwXY_setValueAtX( convolute, y, c ) ) != nfu_Okay ) goto Err;
        } }
    else {
        for( i1 = 0; i1 < n1; i1++ ) {
            for( i2 = 0; i2 < n2; i2++ ) {
                y = yMin + ( f1->points[i1].x - f1->points[0].x ) + ( f2->points[i2].x - f2->points[0].x );
                if( y <= yMin ) continue;
                if( y >= yMax ) continue;
                if( ( *status = ptwXY_convolution2( f1, f2, y, yMin, &c ) ) != nfu_Okay ) goto Err;
                if( ( *status = ptwXY_setValueAtX( convolute, y, c ) ) != nfu_Okay ) goto Err;
            }
        }
    }
    if( ( *status = ptwXY_setValueAtX( convolute, yMax, 0. ) ) != nfu_Okay ) goto Err;
    if( ( *status = ptwXY_simpleCoalescePoints( convolute ) ) != nfu_Okay ) goto Err;
    for( i1 = convolute->length - 1; i1 > 0; i1-- ) {
        if( ( *status = ptwXY_convolution3( convolute, f1, f2, convolute->points[i1 - 1].x, convolute->points[i1 - 1].y, 
            convolute->points[i1].x, convolute->points[i1].y, yMin ) ) != nfu_Okay ) goto Err;
    }

    return( convolute );

Err:
    ptwXY_free( convolute );
    return( NULL );
}
/*
************************************************************
*/
static nfu_status ptwXY_convolution2( ptwXYPoints *f1, ptwXYPoints *f2, double y, double yMin, double *c ) {

    int64_t i1 = 0, i2 = 0, n1 = f1->length, n2 = f2->length, mode;
    double dx1, dx2, x1MinP, x1Min, x2Max;
    double f1x1 = 0,  f1y1 = 0,  f1x2 = 0,  f1y2 = 0,  f2x1 = 0,  f2y1 = 0,  f2x2 = 0,  f2y2 = 0;
    double f1x1p, f1y1p, f1x2p, f1y2p, f2x1p, f2y1p, f2x2p, f2y2p;
    ptwXY_lessEqualGreaterX legx;
    ptwXYOverflowPoint lessThanEqualXPoint, greaterThanXPoint;
    nfu_status status;

    x2Max = f2->points[0].x + ( y - yMin );
    if( x2Max > f2->points[n2 - 1].x ) x2Max = f2->points[n2 - 1].x;
    x1Min = f1->points[0].x;
    x1MinP = y - f2->points[n2 - 1].x;
    if( x1Min < x1MinP ) x1Min = x1MinP;
    *c = 0.;

    switch( legx = ptwXY_getPointsAroundX( f1, x1Min, &lessThanEqualXPoint, &greaterThanXPoint ) ) {
    case ptwXY_lessEqualGreaterX_empty :                                /* These three should not happen. */
    case ptwXY_lessEqualGreaterX_lessThan :
    case ptwXY_lessEqualGreaterX_greater :
        return( nfu_Okay );
    case ptwXY_lessEqualGreaterX_equal :
    case ptwXY_lessEqualGreaterX_between :
        i1 = lessThanEqualXPoint.index;
        f1x1 = f1->points[i1].x;
        f1y1p = f1y1 = f1->points[i1].y;
        i1++;
        if( i1 == n1 ) return( nfu_Okay );
        f1x2 = f1->points[i1].x;
        f1y2 = f1->points[i1].y;
        if( legx == ptwXY_lessEqualGreaterX_between ) {
            if( ( status = ptwXY_interpolatePoint( f1->interpolation, x1Min, &f1y1p, f1x1, f1y1, f1x2, f1y2 ) ) != nfu_Okay ) return( status );
        }
        break;
    }

    switch( legx = ptwXY_getPointsAroundX( f2, x2Max, &lessThanEqualXPoint, &greaterThanXPoint ) ) {
    case ptwXY_lessEqualGreaterX_empty :                                /* These three should not happen. */
    case ptwXY_lessEqualGreaterX_lessThan :
    case ptwXY_lessEqualGreaterX_greater :
        return( nfu_Okay );
    case ptwXY_lessEqualGreaterX_equal :
    case ptwXY_lessEqualGreaterX_between :
        i2 = lessThanEqualXPoint.index;
        if( i2 < f2->length - 1 ) i2++;
        f2x2 = f2->points[i2].x;
        f2y2p = f2y2 = f2->points[i2].y;
        i2--;
        f2x1 = f2->points[i2].x;
        f2y1 = f2->points[i2].y;
        if( legx == ptwXY_lessEqualGreaterX_between ) {
            if( ( status = ptwXY_interpolatePoint( f2->interpolation, x2Max, &f2y2p, f2x1, f2y1, f2x2, f2y2 ) ) != nfu_Okay ) return( status );
        }
        break;
    }

    f1x1p = x1Min;
    f2x2p = x2Max;
    f1y2p = f1y2;
    f2y1p = f2y1;
    while( ( i1 < n1 ) && ( i2 >= 0 ) ) { // Loop checking, 11.06.2015, T. Koi
        dx1 = f1x2  - f1x1p;
        dx2 = f2x2p - f2x1;
        mode = 2;
        if( i1 < n1 ) {
            if( dx1 < dx2 ) mode = 1;
        }
        if( mode == 1 ) {                                               /* Point in f1 is limiting dx step size. */
            f2x1p = f2x2p - dx1;
            if( f2x1p < f2->points[i2].x ) {                            /* Round off issue may cause this. */
                f2x1p = f2x2;
                f2y1p = f2y2; }
            else {
                if( ( status = ptwXY_interpolatePoint( f2->interpolation, f2x1p, &f2y1p, f2x1, f2y1, f2x2, f2y2 ) ) != nfu_Okay ) return( status );
            }
            *c += ( ( f1y1p + f1y2p ) * ( f2y1p + f2y2p ) + f1y1p * f2y2p + f1y2p * f2y1p ) * dx1; /* Note the reversing of f2y1p and f2y2p. */
            i1++;
            if( i1 == n1 ) break;
            f1x1p = f1x1 = f1x2;
            f1y1p = f1y1 = f1y2;
            f1x2 = f1->points[i1].x;
            f1y2p = f1y2 = f1->points[i1].y;
            f2x2p = f2x1p;
            f2y2p = f2y1p;
            f2y1p = f2y1; }
        else {
            f1x2p = f1x1p + dx2;
            if( ( f1x2p > f1->points[i1].x ) || ( dx1 == dx2 ) ) {      /* Round off issue may cause first test to trip. */
                f1x2p = f1x2;
                f1y2p = f1y2; }
            else {
                if( ( status = ptwXY_interpolatePoint( f1->interpolation, f1x2p, &f1y2p, f1x1, f1y1, f1x2, f1y2 ) ) != nfu_Okay ) return( status );
            }
            *c += ( ( f1y1p + f1y2p ) * ( f2y1p + f2y2p ) + f1y1p * f2y2p + f1y2p * f2y1p ) * dx2;  /* Note the reversing of f2y1p and f2y2p. */
            if( i2 == 0 ) break;
            i2--;
            f2x2p = f2x2 = f2x1;
            f2y2p = f2y2 = f2y1;
            f2x1 = f2->points[i2].x;
            f2y1p = f2y1 = f2->points[i2].y;
            f1x1p = f1x2p;
            if( dx1 == dx2 ) {
                f1x1p = f1x1 = f1x2;
                f1y1p = f1y1 = f1y2;
                i1++;
                f1x2 = f1->points[i1].x;
                f1y2p = f1y2 = f1->points[i1].y; }
            else {
                f1y1p = f1y2p;
                f1y2p = f1->points[i1].y;
            }
        }
    }
    *c /= 6.;
    return( nfu_Okay );
}
/*
************************************************************
*/
static nfu_status ptwXY_convolution3( ptwXYPoints *convolute, ptwXYPoints *f1, ptwXYPoints *f2, double y1, double c1, double y2, double c2, double yMin ) {

    nfu_status status;
    double yMid = 0.5 * ( y1 + y2 ), cMid = 0.5 * ( c1 + c2 ), c;

    if( ( y2 - yMid ) <= 1e-5 * ( ptwXY_getXMax( convolute ) - ptwXY_getXMin( convolute ) ) ) return( nfu_Okay );
    if( ( status = ptwXY_convolution2( f1, f2, yMid, yMin, &c ) ) != nfu_Okay ) return( status );
    if( std::fabs( c - cMid ) <= convolute->accuracy * 0.5 * ( std::fabs( c ) + std::fabs( cMid ) ) ) return( nfu_Okay );
    if( ( status = ptwXY_setValueAtX( convolute, yMid, c ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_convolution3( convolute, f1, f2, y1, c1, yMid, c, yMin ) ) != nfu_Okay ) return( status );
    return( ptwXY_convolution3( convolute, f1, f2, yMid, c, y2, c2, yMin ) );
}

#if defined __cplusplus
}
#endif
