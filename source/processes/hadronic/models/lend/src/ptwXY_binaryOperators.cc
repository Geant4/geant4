/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <cmath>
#include <float.h>

#include "ptwXY.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

static double ptwXY_mod2( double v, double m, int pythonMod );
static nfu_status ptwXY_mul2_s_ptwXY( ptwXYPoints *n, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, double x1, double y1, double x2, double y2, int level );
static nfu_status ptwXY_div_s_ptwXY( ptwXYPoints *n, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, double x1, double y1, double x2, double y2, 
    int level, int isNAN1, int isNAN2 );
static ptwXYPoints *ptwXY_div_ptwXY_forFlats( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status, int safeDivide );
static nfu_status ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXYPoints *ptwXY1, double x, double *y );
/*
************************************************************
*/
nfu_status ptwXY_slopeOffset( ptwXYPoints *ptwXY, double slope, double offset ) { 

    int64_t i, nonOverflowLength = ptwXY_getNonOverflowLength( ptwXY );
    ptwXYPoint *p;
    ptwXYOverflowPoint *o, *overflowHeader = &(ptwXY->overflowHeader);

    if( ptwXY->status != nfu_Okay ) return( ptwXY->status );

    for( i = 0, p = ptwXY->points; i < nonOverflowLength; i++, p++ ) p->y = slope * p->y + offset;
    for( o = overflowHeader->next; o != overflowHeader; o = o->next ) o->point.y = slope * o->point.y + offset; 
    return( ptwXY->status );
}   
/*
************************************************************
*/
nfu_status ptwXY_add_double( ptwXYPoints *ptwXY, double value ) { return( ptwXY_slopeOffset( ptwXY, 1., value ) ); }
nfu_status ptwXY_sub_doubleFrom( ptwXYPoints *ptwXY, double value ) { return( ptwXY_slopeOffset( ptwXY,  1., -value ) ); }
nfu_status ptwXY_sub_fromDouble( ptwXYPoints *ptwXY, double value ) { return( ptwXY_slopeOffset( ptwXY, -1.,  value ) ); }
nfu_status ptwXY_mul_double( ptwXYPoints *ptwXY, double value ) { return( ptwXY_slopeOffset( ptwXY, value, 0. ) ); }
nfu_status ptwXY_div_doubleFrom( ptwXYPoints *ptwXY, double value ) { 

    if( value == 0. ) {
        ptwXY->status = nfu_divByZero; }
    else {
        ptwXY_slopeOffset( ptwXY, 1. / value, 0. ); 
    }
    return( ptwXY->status );
}
nfu_status ptwXY_div_fromDouble( ptwXYPoints *ptwXY, double value ) {
/*
*   This does not do any infilling and it should?????????
*/

    int64_t i, nonOverflowLength = ptwXY_getNonOverflowLength( ptwXY );
    ptwXYPoint *p;
    ptwXYOverflowPoint *o, *overflowHeader = &(ptwXY->overflowHeader);

    if( ptwXY->status != nfu_Okay ) return( ptwXY->status );
    if( ptwXY->interpolation == ptwXY_interpolationOther ) return( nfu_otherInterpolation );

    for( i = 0, p = ptwXY->points; i < nonOverflowLength; i++, p++ ) if( p->y == 0. ) ptwXY->status = nfu_divByZero;
    for( o = overflowHeader->next; o != overflowHeader; o = o->next ) if( o->point.y == 0. ) ptwXY->status = nfu_divByZero;
    if( ptwXY->status != nfu_divByZero ) {
        for( i = 0, p = ptwXY->points; i < nonOverflowLength; i++, p++ ) p->y = value / p->y;
        for( o = overflowHeader->next; o != overflowHeader; o = o->next ) o->point.y = value / o->point.y; 
    }
    return( ptwXY->status );
}
/*
************************************************************
*/
nfu_status ptwXY_mod( ptwXYPoints *ptwXY, double m, int pythonMod ) { 

    int64_t i, nonOverflowLength = ptwXY_getNonOverflowLength( ptwXY );
    ptwXYPoint *p;
    ptwXYOverflowPoint *o, *overflowHeader = &(ptwXY->overflowHeader);

    if( ptwXY->status != nfu_Okay ) return( ptwXY->status );
    if( m == 0 ) return( ptwXY->status = nfu_divByZero );

    for( i = 0, p = ptwXY->points; i < nonOverflowLength; i++, p++ ) p->y = ptwXY_mod2( p->y, m, pythonMod );
    for( o = overflowHeader->next; o != overflowHeader; o = o->next ) o->point.y = ptwXY_mod2( o->point.y, m, pythonMod );
    return( ptwXY->status );
}
/*
************************************************************
*/
static double ptwXY_mod2( double v, double m, int pythonMod ) {

    double r = std::fmod( std::fabs( v ), std::fabs( m ) );

    if( pythonMod ) {
        if( ( v * m ) < 0. ) r = std::fabs( m ) - std::fabs( r );
        if( m < 0. ) r *= -1.; }
    else {
        if( v < 0. ) r *= -1.;
    }

    return( r );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_binary_ptwXY( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, double v1, double v2, double v1v2, nfu_status *status ) {

    int64_t i;
    int unionOptions = ptwXY_union_fill | ptwXY_union_mergeClosePoints;
    double y;
    ptwXYPoints *n;
    ptwXYPoint *p;

    *status = nfu_otherInterpolation;
    if( ptwXY1->interpolation == ptwXY_interpolationOther ) return( NULL );
    if( ptwXY2->interpolation == ptwXY_interpolationOther ) return( NULL );
    if( ( *status = ptwXY_areDomainsMutual( ptwXY1, ptwXY2 ) ) != nfu_Okay ) return( NULL );
    if( ( ptwXY1->interpolation == ptwXY_interpolationFlat ) || ( ptwXY2->interpolation == ptwXY_interpolationFlat ) ) {
        *status = nfu_invalidInterpolation;
        if( ( ptwXY1->interpolation != ptwXY2->interpolation ) ) return( NULL );
    }
    if( ( n = ptwXY_union( ptwXY1, ptwXY2, status, unionOptions ) ) != NULL ) {
        for( i = 0, p = n->points; i < n->length; i++, p++ ) {
            if( ( *status = ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXY2, p->x, &y ) ) != nfu_Okay ) goto Err;
            p->y = v1 * p->y + v2 * y + v1v2 * y * p->y;
        }
    }
    return( n );
Err:
    if( n ) ptwXY_free( n );
    return( NULL );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_add_ptwXY( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status ) {

    ptwXYPoints *sum;

    if( ptwXY1->length == 0 ) {
        sum = ptwXY_clone( ptwXY2, status ); }
    else if( ptwXY2->length == 0 ) {
        sum = ptwXY_clone( ptwXY1, status ); }
    else {
        sum = ptwXY_binary_ptwXY( ptwXY1, ptwXY2, 1., 1., 0., status );
    }
    return( sum );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_sub_ptwXY( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status ) {

    ptwXYPoints *diff;

    if( ptwXY1->length == 0 ) {
        diff = ptwXY_clone( ptwXY2, status );
        if( ( *status = ptwXY_neg( diff ) ) != nfu_Okay ) diff = ptwXY_free( diff ); }
    else if( ptwXY2->length == 0 ) {
        diff = ptwXY_clone( ptwXY1, status ); }
    else {
        diff = ptwXY_binary_ptwXY( ptwXY1, ptwXY2, 1., -1., 0., status );
    }
    return( diff );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_mul_ptwXY( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status ) {

    ptwXYPoints *mul;

    if( ptwXY1->length == 0 ) {
        mul = ptwXY_clone( ptwXY1, status ); }
    else if( ptwXY2->length == 0 ) {
        mul = ptwXY_clone( ptwXY2, status ); }
    else {
        mul = ptwXY_binary_ptwXY( ptwXY1, ptwXY2, 0., 0., 1., status );
    }
    return( mul );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_mul2_ptwXY( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status ) {

    int64_t i, length;
    ptwXYPoints *n = NULL;
    int found;
    double x1, y1, x2, y2, u1, u2, v1, v2, xz1 = 0, xz2 = 0, x;

    *status = nfu_otherInterpolation;
    if( ptwXY1->interpolation == ptwXY_interpolationOther ) return( NULL );
    if( ptwXY2->interpolation == ptwXY_interpolationOther ) return( NULL );
    if( ( n = ptwXY_mul_ptwXY( ptwXY1, ptwXY2, status ) ) == NULL ) return( n );
    if( ptwXY1->interpolation == ptwXY_interpolationFlat ) return( n );
    if( ptwXY2->interpolation == ptwXY_interpolationFlat ) return( n );
    length = n->length - 1;
    if( length > 0 ) {
        x2 = n->points[length].x;
        for( i = length - 1; i >= 0; i-- ) {             /* Find and add y zeros not currently in n's. */
            x1 = n->points[i].x;
            if( ( *status = ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXY1, x1, &u1 ) ) != nfu_Okay ) goto Err;
            if( ( *status = ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXY1, x2, &u2 ) ) != nfu_Okay ) goto Err;
            if( ( *status = ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXY2, x1, &v1 ) ) != nfu_Okay ) goto Err;
            if( ( *status = ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXY2, x2, &v2 ) ) != nfu_Okay ) goto Err;
            found = 0;
            if( u1 * u2 < 0 ) {
                xz1 = ( u1 * x2 - u2 * x1 ) / ( u1 - u2 );
                if( ( *status = ptwXY_setValueAtX( n, xz1, 0. ) ) != nfu_Okay ) goto Err;
                found = 1;
            }
            if( v1 * v2 < 0 ) {
                xz2 = ( v1 * x2 - v2 * x1 ) / ( v1 - v2 );
                if( ( *status = ptwXY_setValueAtX( n, xz2, 0. ) ) != nfu_Okay ) goto Err;
                found += 1;
            }
            if( found > 1 ) {
                x = 0.5 * ( xz1 + xz2 );
                if( ( *status = ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXY1, x, &u1 ) ) != nfu_Okay ) goto Err;
                if( ( *status = ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXY2, x, &v1 ) ) != nfu_Okay ) goto Err;
                if( ( *status = ptwXY_setValueAtX( n, x, u1 * v1 ) ) != nfu_Okay ) goto Err;
            }
            x2 = x1;
        }

        if( ( *status = ptwXY_simpleCoalescePoints( n ) ) != nfu_Okay ) goto Err;
        length = n->length;
        x2 = n->points[n->length-1].x;
        y2 = n->points[n->length-1].y;
        for( i = n->length - 2; i >= 0; i-- ) {             /* Make interpolation fit accuracy. Work backwards so new points will not mess up loop. */
            x1 = n->points[i].x;
            y1 = n->points[i].y;
            if( ( *status = ptwXY_mul2_s_ptwXY( n, ptwXY1, ptwXY2, x1, y1, x2, y2, 0 ) ) != nfu_Okay ) goto Err;
            x2 = x1;
            y2 = y1;
        }
        ptwXY_update_biSectionMax( n, (double) length );
    }
    return( n );

Err:
    if( n ) ptwXY_free( n );
    return( NULL );
}
/*
************************************************************
*/
static nfu_status ptwXY_mul2_s_ptwXY( ptwXYPoints *n, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, double x1, double y1, double x2, double y2, int level ) {

    nfu_status status;
    double u1, u2, v1, v2, x, y, yp, dx, a1, a2;

    if( ( x2 - x1 ) < ClosestAllowXFactor * DBL_EPSILON * ( std::fabs( x1 ) + std::fabs( x2 ) ) ) return( nfu_Okay );
    if( level >= n->biSectionMax ) return( nfu_Okay );
    level++;
    if( ( status = ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXY1, x1, &u1 ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXY1, x2, &u2 ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXY2, x1, &v1 ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXY2, x2, &v2 ) ) != nfu_Okay ) return( status );
    if( ( u1 == u2 ) || ( v1 == v2 ) ) return( nfu_Okay );
    a1 = u1 * v1;
    if( y1 == 0 ) a1 = 0.;                                  /* Fix rounding problem. */
    a2 = u2 * v2;
    if( y2 == 0 ) a2 = 0.;                                  /* Fix rounding problem. */
    if( ( a1 == 0. ) || ( a2 == 0. ) ) {                    /* Handle special case of 0 where accuracy can never be met. */
        x = 0.5 * ( x1 + x2 ); }
    else {
        if( ( a1 * a2 < 0. ) ) return( nfu_Okay );  /* Assume rounding error and no point needed as zero crossings are set in ptwXY_mul2_ptwXY. */
        a1 = std::sqrt( std::fabs( a1 ) );
        a2 = std::sqrt( std::fabs( a2 ) );
        x = ( a2 * x1 + a1 * x2 ) / ( a2 + a1 );
    }
    dx = x2 - x1;
    yp = ( u1 * v1 * ( x2 - x ) + u2 * v2 * ( x - x1 ) ) / dx;
    y = ( u1 * ( x2 - x ) + u2 * ( x - x1 ) ) * ( v1 * ( x2 - x ) + v2 * ( x - x1 ) ) / ( dx * dx );
    if( std::fabs( y - yp ) < std::fabs( y * n->accuracy ) ) return( nfu_Okay );
    if( ( status = ptwXY_setValueAtX( n, x, y ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_mul2_s_ptwXY( n, ptwXY1, ptwXY2, x, y, x2, y2, level ) ) != nfu_Okay ) return( status );
    status = ptwXY_mul2_s_ptwXY( n, ptwXY1, ptwXY2, x1, y1, x, y, level );
    return( status );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_div_ptwXY( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status, int safeDivide ) {

    int isNAN1, isNAN2;
    int64_t i, j, k, zeros = 0, length, iYs;
    double x1, x2, y1, y2, u1, u2, v1, v2, y, xz, nan = nfu_getNAN( ), s1, s2;
    ptwXYPoints *n = NULL;
    ptwXYPoint *p;

    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY1 ) ) != nfu_Okay ) return( NULL );
    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY2 ) ) != nfu_Okay ) return( NULL );
    *status = nfu_otherInterpolation;
    if( ptwXY1->interpolation == ptwXY_interpolationOther ) return( NULL );
    if( ptwXY2->interpolation == ptwXY_interpolationOther ) return( NULL );
    if( ( ptwXY1->interpolation == ptwXY_interpolationFlat ) || ( ptwXY1->interpolation == ptwXY_interpolationFlat ) )
        return( ptwXY_div_ptwXY_forFlats( ptwXY1, ptwXY2, status, safeDivide ) );
    
    if( ( *status = ptwXY_areDomainsMutual( ptwXY1, ptwXY2 ) ) != nfu_Okay ) return( NULL );
    if( ( n = ptwXY_union( ptwXY1, ptwXY2, status, ptwXY_union_fill | ptwXY_union_mergeClosePoints ) ) != NULL ) {
        for( i = 0, p = n->points; i < n->length; i++, p++ ) {
            if( ( *status = ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXY2, p->x, &y ) ) != nfu_Okay ) goto Err;
            if( y == 0. ) {
                if( p->y == 0. ) {
                    iYs = 0;
                    y1 = 0.;
                    y2 = 0.;
                    if( i > 0 ) {
                        if( ( *status = ptwXY_getSlopeAtX( ptwXY1, p->x, '-', &s1 ) ) != nfu_Okay ) {
                            if( *status != nfu_XOutsideDomain ) goto Err;
                            s1 = 0.;
                        }
                        if( ( *status = ptwXY_getSlopeAtX( ptwXY2, p->x, '-', &s2 ) ) != nfu_Okay ) goto Err;
                        if( s2 == 0. ) {
                            y1 = nan; }
                        else {
                            y1 = s1 / s2;
                        } 
                        iYs++;
                    }
                    if( i < ( n->length - 1 ) ) {
                        if( ( *status = ptwXY_getSlopeAtX( ptwXY1, p->x, '+', &s1 ) ) != nfu_Okay ) {
                            if( *status != nfu_XOutsideDomain ) goto Err;
                            s1 = 0.;
                        }
                        if( ( *status = ptwXY_getSlopeAtX( ptwXY2, p->x, '+', &s2 ) ) != nfu_Okay ) goto Err;
                        if( s2 == 0. ) {
                            y2 = nan; }
                        else {
                            y2 = s1 / s2;
                        }
                        iYs++;
                    }
                    p->y = ( y1 + y2 ) / iYs; 
                    if( nfu_isNAN( p->y ) ) zeros++; }
                else {
                    if( !safeDivide ) {
                        *status = nfu_divByZero;
                        goto Err;
                    }
                    zeros++;
                    p->y = nan;
                } }
            else {
                p->y /= y;
            }
        }
        length = n->length - 1;
        if( length > 0 ) {
            x2 = n->points[length].x;
            for( i = length - 1; i >= 0; i-- ) {             /* Find and add y zeros and NAN not currently in n's. */
                x1 = n->points[i].x;
                if( ( *status = ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXY1, x1, &u1 ) ) != nfu_Okay ) goto Err;
                if( ( *status = ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXY1, x2, &u2 ) ) != nfu_Okay ) goto Err;
                if( ( *status = ptwXY_getValueAtX( ptwXY2, x1, &v1 ) ) != nfu_Okay ) goto Err;
                if( ( *status = ptwXY_getValueAtX( ptwXY2, x2, &v2 ) ) != nfu_Okay ) goto Err;
                if( u1 * u2 < 0 ) {
                    xz = ( u1 * x2 - u2 * x1 ) / ( u1 - u2 );
                    if( ( *status = ptwXY_setValueAtX( n, xz, 0. ) ) != nfu_Okay ) goto Err;
                }
                if( v1 * v2 < 0 ) {
                    if( !safeDivide ) {
                        *status = nfu_divByZero;
                        goto Err;
                    }
                    zeros++;
                    xz = ( v1 * x2 - v2 * x1 ) / ( v1 - v2 );
                    if( ( *status = ptwXY_setValueAtX( n, xz, nan ) ) != nfu_Okay ) goto Err;
                }
                x2 = x1;
            }
            if( ( *status = ptwXY_simpleCoalescePoints( n ) ) != nfu_Okay ) goto Err;
            length = n->length;
            x2 = n->points[n->length-1].x;
            y2 = n->points[n->length-1].y;
            isNAN2 = nfu_isNAN( y2 );
            for( i = n->length - 2; i >= 0; i-- ) {             /* Make interpolation fit accuracy. Work backwards so new points will not mess up loop. */
                x1 = n->points[i].x;
                y1 = n->points[i].y;
                isNAN1 = nfu_isNAN( y1 );
                if( !isNAN1 || !isNAN2 ) {
                    if( ( *status = ptwXY_div_s_ptwXY( n, ptwXY1, ptwXY2, x1, y1, x2, y2, 0, isNAN1, isNAN2 ) ) != nfu_Okay ) goto Err;
                }
                x2 = x1;
                y2 = y1;
                isNAN2 = isNAN1;
            }
            ptwXY_update_biSectionMax( n, (double) length );
            if( zeros ) {
                if( ( *status = ptwXY_simpleCoalescePoints( n ) ) != nfu_Okay ) goto Err;
                for( i = 0; i < n->length; i++ ) if( !nfu_isNAN( n->points[i].y ) ) break;
                if( nfu_isNAN( n->points[0].y ) ) {                     /* Special case for first point. */
                    if( i == n->length ) {                              /* They are all nan's, what now? */
                        zeros = 0;
                        for( i = 0; i < n->length; i++ ) n->points[i].y = 0.; }
                    else {
                         n->points[0].y = 2. * n->points[i].y;
                        zeros--;
                    }
                }
                for( i = n->length - 1; i > 0; i-- ) if( !nfu_isNAN( n->points[i].y ) ) break;
                if( nfu_isNAN( n->points[n->length - 1].y ) ) {         /* Special case for last point. */
                    n->points[n->length - 1].y = 2. * n->points[i].y;
                    zeros--;
                }
                if( zeros ) {
                    for( i = 0; i < n->length; i++ ) if( nfu_isNAN( n->points[i].y ) ) break;
                    for( k = i + 1, j = i; k < n->length; k++ ) {
                    if( nfu_isNAN( n->points[k].y ) ) continue;
                        n->points[j] = n->points[k];
                        j++;
                    }
                    n->length = j;
                }
            }
        }
    }
    return( n );

Err:
    if( n ) ptwXY_free( n );
    return( NULL );
}
/*
************************************************************
*/
static nfu_status ptwXY_div_s_ptwXY( ptwXYPoints *n, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, double x1, double y1, double x2, double y2, 
        int level, int isNAN1, int isNAN2 ) {

    nfu_status status;
    double u1, u2, v1, v2, v, x, y, yp, dx, a1, a2;

    if( ( x2 - x1 ) < ClosestAllowXFactor * DBL_EPSILON * ( std::fabs( x1 ) + std::fabs( x2 ) ) ) return( nfu_Okay );
    if( level >= n->biSectionMax ) return( nfu_Okay );
    level++;
    if( ( status = ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXY1, x1, &u1 ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXY1, x2, &u2 ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_getValueAtX( ptwXY2, x1, &v1 ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_getValueAtX( ptwXY2, x2, &v2 ) ) != nfu_Okay ) return( status );
    if( isNAN1 ) {
        x = 0.5 * ( x1 + x2 );
        if( ( status = ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXY1, x, &u1 ) ) != nfu_Okay ) return( status );
        if( ( status = ptwXY_getValueAtX( ptwXY2, x, &v1 ) ) != nfu_Okay ) return( status );
        y = u1 / v1; }
    else if( isNAN2 ) {
        x = 0.5 * ( x1 + x2 );
        if( ( status = ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXY1, x, &u2 ) ) != nfu_Okay ) return( status );
        if( ( status = ptwXY_getValueAtX( ptwXY2, x, &v2 ) ) != nfu_Okay ) return( status );
        y = u2 / v2; }
    else {
        if( ( u1 == u2 ) || ( v1 == v2 ) ) return( nfu_Okay );
        if( ( y1 == 0. ) || ( y2 == 0. ) ) {                    /* Handle special case of 0 where accuracy can never be met. */
            x = 0.5 * ( x1 + x2 ); }
        else {
            if( ( u1 * u2 < 0. ) ) return( nfu_Okay );  /* Assume rounding error and no point needed. */
            a1 = std::sqrt( std::fabs( u1 ) );
            a2 = std::sqrt( std::fabs( u2 ) );
            x = ( a2 * x1 + a1 * x2 ) / ( a2 + a1 );
        }
        dx = x2 - x1;
        v = v1 * ( x2 - x ) + v2 * ( x - x1 );
        if( ( v1 == 0. ) || ( v2 == 0. ) || ( v == 0. ) ) return( nfu_Okay );     /* Probably not correct, but I had to do something. */
        yp = ( u1 / v1 * ( x2 - x ) + u2 / v2 * ( x - x1 ) ) / dx;
        y = ( u1 * ( x2 - x ) + u2 * ( x - x1 ) ) / v;
        if( std::fabs( y - yp ) < std::fabs( y * n->accuracy ) ) return( nfu_Okay );
    }
    if( ( status = ptwXY_setValueAtX( n, x, y ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_div_s_ptwXY( n, ptwXY1, ptwXY2, x, y, x2, y2, level, 0, isNAN2 ) ) != nfu_Okay ) return( status );
    status = ptwXY_div_s_ptwXY( n, ptwXY1, ptwXY2, x1, y1, x, y, level, isNAN1, 0 );
    return( status );
}
/*
************************************************************
*/
static ptwXYPoints *ptwXY_div_ptwXY_forFlats( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status, int safeDivide ) {

    int64_t i;
    ptwXYPoints *n = NULL;
    ptwXYPoint *p;
    double y;

    *status = nfu_invalidInterpolation;
    if( ptwXY1->interpolation != ptwXY_interpolationFlat ) return( NULL );
    if( ptwXY2->interpolation != ptwXY_interpolationFlat ) return( NULL );
    if( ( n = ptwXY_union( ptwXY1, ptwXY2, status, ptwXY_union_fill | ptwXY_union_mergeClosePoints ) ) != NULL ) {
        for( i = 0, p = n->points; i < n->length; i++, p++ ) {
            if( ( *status = ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXY2, p->x, &y ) ) != nfu_Okay ) goto Err;
            if( y == 0. ) {
                if( ( safeDivide ) && ( p->y == 0 ) ) {
                    *status = nfu_divByZero;
                    goto Err;
                } }
            else {
                p->y /= y;
            }
        }
    }
    return( n );

Err:
    if( n ) ptwXY_free( n );
    return( NULL );
}
/*
************************************************************
*/
static nfu_status ptwXY_getValueAtX_ignore_XOutsideDomainError( ptwXYPoints *ptwXY1, double x, double *y ) {

    nfu_status status = ptwXY_getValueAtX( ptwXY1, x, y );

    if( status == nfu_XOutsideDomain ) status = nfu_Okay;
    return( status );
}

#if defined __cplusplus
}
#endif
