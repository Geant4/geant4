/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

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

static nfu_status ptwXY_clip2( ptwXYPoints *ptwXY1, double y, double x1, double y1, double x2, double y2 );
static double ptwXY_thicken_linear_dx( int sectionSubdivideMax, double dxMax, double x1, double x2 );
static nfu_status ptwXY_thin2( ptwXYPoints *thinned, char *thin, double accuracy, int64_t i1, int64_t i2 );
/*
************************************************************
*/
nfu_status ptwXY_clip( ptwXYPoints *ptwXY1, double yMin, double yMax ) {
/*
    This function acts oddly for xy = [ [ 1, 0 ], [ 3, -2 ], [ 4, 1 ] ] and yMin = 0.2, why???????
    This function probably only works for linear, linear interpolation (mainly because of ptwXY_clip2).
*/
    int64_t i, j, n;
    double x2, y2;
    nfu_status status;
    ptwXYPoints *clipped;
    ptwXYPoint *points;

    if( ( status = ptwXY_simpleCoalescePoints( ptwXY1 ) ) != nfu_Okay ) return( status );
    if( ptwXY1->interpolation == ptwXY_interpolationOther ) return( nfu_otherInterpolation );
    n = ptwXY1->length;
    if( n > 0 ) {
        i = 0;
        if( ptwXY_getYMax( ptwXY1 ) < yMin ) i = 1;
        if( ptwXY_getYMin( ptwXY1 ) > yMax ) i = 1;
        if( i == 1 ) return( ptwXY_clear( ptwXY1 ) );
    }
    if( n == 1 ) {
        y2 = ptwXY1->points[0].y;
        if( y2 < yMin ) {
            ptwXY1->points[0].y = yMin; }
        else if( y2 > yMax ) {
            ptwXY1->points[0].y = yMax;
        } }
    else if( n > 1 ) {
        if( ( clipped = ptwXY_new( ptwXY1->interpolation, &(ptwXY1->interpolationOtherInfo), 
                ptwXY1->biSectionMax, ptwXY1->accuracy, n, 10, &status, ptwXY1->userFlag ) ) == NULL )
            return( ptwXY1->status = status );
        for( i = 0; i < n; i++ ) {
            x2 = ptwXY1->points[i].x;
            y2 = ptwXY1->points[i].y;
            if( y2 < yMin ) {
                if( i > 0 ) {
                    points = ptwXY_getPointAtIndex_Unsafely( clipped, clipped->length - 1 );
                    if( points->y > yMin ) {
                        if( ( status = ptwXY_clip2( clipped, yMin, points->x, points->y, x2, y2 ) ) != nfu_Okay ) goto Err;
                    }
                }
                if( ( status = ptwXY_setValueAtX( clipped, x2, yMin ) ) != nfu_Okay ) goto Err;
                j = i;
                for( i++; i < n; i++ ) if( !( ptwXY1->points[i].y < yMin ) ) break;
                if( i < n ) {
                    x2 = ptwXY1->points[i].x;
                    y2 = ptwXY1->points[i].y;
                    if( ( status = ptwXY_clip2( clipped, yMin, ptwXY1->points[i-1].x, ptwXY1->points[i-1].y, x2, y2 ) ) != nfu_Okay ) goto Err;
                    if( y2 > yMax ) {
                        if( ( status = ptwXY_clip2( clipped, yMax, ptwXY1->points[i-1].x, ptwXY1->points[i-1].y, x2, y2 ) ) != nfu_Okay ) goto Err;
                    } }
                else if( j != n - 1 ) {
                    if( ( status = ptwXY_setValueAtX( clipped, ptwXY1->points[n - 1].x, yMin ) ) != nfu_Okay ) goto Err;
                }
                i--; }
            else if( y2 > yMax ) {
                if( i > 0 ) {
                    points = ptwXY_getPointAtIndex_Unsafely( clipped, clipped->length - 1 );
                    if( points->y < yMax ) {
                        if( ( status = ptwXY_clip2( clipped, yMax, points->x, points->y, x2, y2 ) ) != nfu_Okay ) goto Err;
                    }
                }
                if( ( status = ptwXY_setValueAtX( clipped, x2, yMax ) ) != nfu_Okay ) goto Err;
                j = i;
                for( i++; i < n; i++ ) if( !( ptwXY1->points[i].y > yMax ) ) break;
                if( i < n ) {
                    x2 = ptwXY1->points[i].x;
                    y2 = ptwXY1->points[i].y;
                    if( ( status = ptwXY_clip2( clipped, yMax, ptwXY1->points[i-1].x, ptwXY1->points[i-1].y, x2, y2 ) ) != nfu_Okay ) goto Err;
                    if( y2 < yMin ) {
                        if( ( status = ptwXY_clip2( clipped, yMin, ptwXY1->points[i-1].x, ptwXY1->points[i-1].y, x2, y2 ) ) != nfu_Okay ) goto Err;
                    } }
                else if( j != n - 1 ) {
                    if( ( status = ptwXY_setValueAtX( clipped, ptwXY1->points[n - 1].x, yMax ) ) != nfu_Okay ) goto Err;
                }
                i--; }
            else {
                if( ( status = ptwXY_setValueAtX( clipped, x2, y2 ) ) != nfu_Okay ) goto Err;
            }
        }
        if( ( status = ptwXY_simpleCoalescePoints( clipped ) ) != nfu_Okay ) goto Err;
        ptwXY1->length = clipped->length;   /* The squeamish may want to skip the next few lines. */
        clipped->length = n;
        n = ptwXY1->allocatedSize;
        ptwXY1->allocatedSize = clipped->allocatedSize;
        clipped->allocatedSize = n;
        points = clipped->points;
        clipped->points = ptwXY1->points;
        ptwXY1->points = points;
        ptwXY_free( clipped );
    }

    return( ptwXY1->status );

Err:
    ptwXY_free( clipped );
    return( ptwXY1->status = status );
}
/*
************************************************************
*/
static nfu_status ptwXY_clip2( ptwXYPoints *clipped, double y, double x1, double y1, double x2, double y2 ) {

    double x;
    nfu_status status = nfu_Okay;

    x = ( y - y1 ) * ( x2 - x1 ) / ( y2 - y1 ) + x1;
    if( x <= x1 ) {
        x = x1; }
    else if( x >= x2 ) {
        x = x1; }
    else {
        status = ptwXY_setValueAtX( clipped, x, y );
    }
    return( status  );
}
/*
************************************************************
*/
nfu_status ptwXY_thicken( ptwXYPoints *ptwXY1, int sectionSubdivideMax, double dxMax, double fxMax ) {

    double x1, x2 = 0., y1, y2 = 0., fx = 1.1, x, dx, dxp, lfx, y;    /* fx initialized so compilers want complain. */
    int64_t i, notFirstPass = 0;
    int nfx, nDone, doLinear;
    nfu_status status;

    if( ptwXY1->interpolation == ptwXY_interpolationOther ) return( nfu_otherInterpolation );
    if( ( sectionSubdivideMax < 1 ) || ( dxMax < 0. ) || ( fxMax < 1. ) ) return( nfu_badInput );
    if( sectionSubdivideMax > ptwXY_sectionSubdivideMax ) sectionSubdivideMax = ptwXY_sectionSubdivideMax;
    if( ( status = ptwXY_simpleCoalescePoints( ptwXY1 ) ) != nfu_Okay ) return( status );
    for( i = ptwXY1->length - 1; i >= 0; i-- ) {
        x1 = ptwXY1->points[i].x;
        y1 = ptwXY1->points[i].y;
        if( notFirstPass ) {
            dx = ptwXY_thicken_linear_dx( sectionSubdivideMax, dxMax, x1, x2 );

            if( x1 == 0. ) {
                doLinear = 1; }
            else {
                fx = x2 / x1;
                if( fx > 0. ) {
                    lfx = G4Log( fx );
                    if( fxMax == 1. ) {
                        nfx = sectionSubdivideMax; }
                    else {
                        nfx = ( (int) ( lfx / G4Log( fxMax ) ) ) + 1;
                        if( nfx > sectionSubdivideMax ) nfx = sectionSubdivideMax;
                    }
                    if( nfx > 0 ) fx = G4Exp( lfx / nfx );
                    doLinear = 0;
                    if( dx < ( fx - 1 ) * x1 ) doLinear = 1; }
                else {
                    doLinear = 1;
                }
            }
            x = x1;
            dxp = dx;
            nDone = 0;
            while( 1 ) {
                if( doLinear ) {
                    x += dx; }
                else {
                    dx = ptwXY_thicken_linear_dx( sectionSubdivideMax - nDone, dxMax, x, x2 );
                    if( dx <= ( fx - 1 ) * x ) {
                        dxp = dx;
                        doLinear = 1;
                        continue;
                    }
                    dxp = ( fx - 1. ) * x;
                    x *= fx;
                }
                if( ( x2 - x ) < 0.05 * std::fabs( dxp ) ) break;
                if( ( status = ptwXY_interpolatePoint( ptwXY1->interpolation, x, &y, x1, y1, x2, y2 ) ) != nfu_Okay ) return( status );
                if( ( status = ptwXY_setValueAtX( ptwXY1, x, y ) ) != nfu_Okay ) return( status );
                nDone++;
            } // Loop checking, 11.06.2015, T. Koi
        } 
        notFirstPass = 1;
        x2 = x1;
        y2 = y1;
    }
    return( status );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_thin( ptwXYPoints *ptwXY1, double accuracy, nfu_status *status ) {

    int64_t i, j, length = ptwXY1->length;
    ptwXYPoints *thinned = NULL;
    double y1, y2, y3;
    char *thin = NULL;

    if( length < 3 ) return( ptwXY_clone( ptwXY1, status ) );   /* Logic below requires at least 2 points. */
    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY1 ) ) != nfu_Okay ) return( NULL );
    *status = nfu_otherInterpolation;
    if( ptwXY1->interpolation == ptwXY_interpolationOther ) return( NULL );
    if( accuracy < ptwXY1->accuracy ) accuracy = ptwXY1->accuracy;
    if( ( thinned = ptwXY_new( ptwXY1->interpolation, &(ptwXY1->interpolationOtherInfo), 
        ptwXY1->biSectionMax, accuracy, length, ptwXY1->overflowLength, status, ptwXY1->userFlag ) ) == NULL ) return( NULL );

    thinned->points[0] = ptwXY1->points[0];                     /* This sections removes middle point if surrounding points have the same y-value. */
    y1 = ptwXY1->points[0].y;
    y2 = ptwXY1->points[1].y;
    for( i = 2, j = 1; i < length; i++ ) {
        y3 = ptwXY1->points[i].y;
        if( ( y1 != y2 ) || ( y2 != y3 ) ) {
            thinned->points[j++] = ptwXY1->points[i - 1];
            y1 = y2;
            y2 = y3;
        }
    }
    thinned->points[j++] = ptwXY1->points[length - 1];

    if( ptwXY1->interpolation != ptwXY_interpolationFlat ) {    /* Now call ptwXY_thin2 for more thinning. */
        length = thinned->length = j;
        if( ( thin = (char *) nfu_calloc( 1, (size_t) length ) ) == NULL ) goto Err;
        if( ( *status = ptwXY_thin2( thinned, thin, accuracy, 0, length - 1 ) ) != nfu_Okay ) goto Err;
        for( j = 1; j < length; j++ ) if( thin[j] != 0 ) break;
        for( i = j + 1; i < length; i++ ) {
            if( thin[i] == 0 ) {
                thinned->points[j] = thinned->points[i];
                j++;
            }
        }
        nfu_free( thin );
    }
    thinned->length = j;

    return( thinned );

Err:
    ptwXY_free( thinned );
    if( thin != NULL ) nfu_free( thin );
    return( NULL );
}
/*
************************************************************
*/
static nfu_status ptwXY_thin2( ptwXYPoints *thinned, char *thin, double accuracy, int64_t i1, int64_t i2 ) {

    int64_t i, iMax = 0;
    double y, s, dY, dYMax = 0., dYR, dYRMax = 0;
    double x1 = thinned->points[i1].x, y1 = thinned->points[i1].y, x2 = thinned->points[i2].x, y2 = thinned->points[i2].y;
    nfu_status status = nfu_Okay;

    if( i1 + 1 >= i2 ) return( nfu_Okay );
    for( i = i1 + 1; i < i2; i++ ) {
        if( ( status = ptwXY_interpolatePoint( thinned->interpolation, thinned->points[i].x, &y, x1, y1, x2, y2 ) ) != nfu_Okay ) return( status );
        s = 0.5 * ( std::fabs( y ) + std::fabs( thinned->points[i].y ) );
        dY = std::fabs( y - thinned->points[i].y );
        dYR = 0;
        if( s != 0 ) dYR = dY / s;
        if( ( dYR > dYRMax ) || ( ( dYR >= 0.9999 * dYRMax ) && ( dY > dYMax ) ) ) {    /* The choice of 0.9999 is not exact science. */
            iMax = i;
            if( dY > dYMax ) dYMax = dY;
            if( dYR > dYRMax ) dYRMax = dYR;
        }
    }
    if( dYRMax < accuracy ) {
        for( i = i1 + 1; i < i2; i++ ) thin[i] = 1; }
    else {
        if( ( status = ptwXY_thin2( thinned, thin, accuracy, i1, iMax ) ) != nfu_Okay ) return( status );
        status = ptwXY_thin2( thinned, thin, accuracy, iMax, i2 );
    }
    return( status );
}
/*
************************************************************
*/
static double ptwXY_thicken_linear_dx( int sectionSubdivideMax, double dxMax, double x1, double x2 ) {

    int ndx;
    double dx = x2 - x1, dndx;

    if( dxMax == 0. ) {
        dx = ( x2 - x1 ) / sectionSubdivideMax; }
    else {
        dndx = dx / dxMax;
        ndx = (int) dndx;
        if( ( dndx  - ndx ) > 1e-6 ) ndx++;
        if( ndx > sectionSubdivideMax ) ndx = sectionSubdivideMax;
        if( ndx > 0 ) dx /= ndx;
    }

    return( dx );
}
/*
************************************************************
*/
nfu_status ptwXY_trim( ptwXYPoints *ptwXY ) {
/*
c   Remove extra zeros at beginning and end.
*/

    int64_t i, i1, i2;
    nfu_status status;

    if( ptwXY->status != nfu_Okay ) return( ptwXY->status );
    if( ( status = ptwXY_simpleCoalescePoints( ptwXY ) ) != nfu_Okay ) return( status );
    for( i1 = 0; i1 < ptwXY->length; i1++ ) {
        if( ptwXY->points[i1].y != 0 ) break;
    }
    if( i1 > 0 ) i1--;
    for( i2 = ptwXY->length - 1; i2 >= 0; i2-- ) {
        if( ptwXY->points[i2].y != 0 ) break;
    }
    i2++;
    if( i2 < ptwXY->length ) i2++;
    if( i2 > i1 ) {
        if( i1 > 0 ) {
            for( i = i1; i < i2; i++ ) ptwXY->points[i - i1] = ptwXY->points[i];
        }
        ptwXY->length = i2 - i1; }
    else if( i2 < i1 ) {                /* Remove all zeros between endpoints. */
        ptwXY->points[1] = ptwXY->points[ptwXY->length - 1];
        ptwXY->length = 2;
    }

    return( nfu_Okay );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_union( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status, int unionOptions ) {

    int64_t overflowSize, i, i1 = 0, i2 = 0, n1 = ptwXY1->length, n2 = ptwXY2->length, length;
    int fillWithFirst = unionOptions & ptwXY_union_fill, trim = unionOptions & ptwXY_union_trim;
    ptwXYPoints *n;
    double x1 = 0., x2 = 0., y1 = 0., y2 = 0., y, biSectionMax, accuracy;

    if( ( *status = ptwXY1->status ) != nfu_Okay ) return( NULL );
    if( ( *status = ptwXY2->status ) != nfu_Okay ) return( NULL );
    *status = nfu_otherInterpolation;
    if( ptwXY1->interpolation == ptwXY_interpolationOther ) return( NULL );
/*
*   Many other routines use the fact that ptwXY_union calls ptwXY_coalescePoints for ptwXY1 and ptwXY2 so do not change it.
*/
    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY1 ) ) != nfu_Okay ) return( NULL );
    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY2 ) ) != nfu_Okay ) return( NULL );

    if( ( n1 == 1 ) || ( n2 == 1 ) ) {
        *status = nfu_tooFewPoints;
        return( NULL );
    }
    if( trim ) {
        if( n1 > 0 ) {
            if( n2 > 0 ) {
                if( ptwXY1->points[0].x < ptwXY2->points[0].x ) {
                    while( i1 < n1 ) { // Loop checking, 11.05.2015, T. Koi
                        if( ptwXY1->points[i1].x >= ptwXY2->points[0].x ) break;
                        if( fillWithFirst ) {
                            if( i1 < ( ptwXY1->length - 1 ) ) {
                                x1 = ptwXY1->points[i1].x;
                                y1 = ptwXY1->points[i1].y;
                                x2 = ptwXY1->points[i1+1].x;
                                y2 = ptwXY1->points[i1+1].y; 
                            }
                        }
                        i1++;
                    } }
                else {
                    while( i2 < n2 ) { // Loop checking, 11.06.2015, T. Koi
                        if( ptwXY2->points[i2].x >= ptwXY1->points[0].x ) break;
                        i2++;
                    }
                }
                if( ptwXY1->points[n1-1].x > ptwXY2->points[n2-1].x ) {
                    while( i1 < n1 ) { // Loop checking, 11.06.2015, T. Koi
                        if( ptwXY1->points[n1-1].x <= ptwXY2->points[n2-1].x ) break;
                        n1--;
                    } }
                else {
                    while( i2 < n2 ) { // Loop checking, 11.06.2015, T. Koi
                        if( ptwXY2->points[n2-1].x <= ptwXY1->points[n1-1].x ) break;
                        n2--;
                    }
                } }
            else {
                n1 = 0;
            } }
        else {
            n2 = 0;
        }
    }
    overflowSize = ptwXY1->overflowAllocatedSize;
    if( overflowSize < ptwXY2->overflowAllocatedSize ) overflowSize = ptwXY2->overflowAllocatedSize;
    length = ( n1 - i1 ) + ( n2 - i2 );
    if( length == 0 ) length = ptwXY_minimumSize;
    biSectionMax = ptwXY1->biSectionMax;
    if( biSectionMax < ptwXY2->biSectionMax ) biSectionMax = ptwXY2->biSectionMax;
    accuracy = ptwXY1->accuracy;
    if( accuracy < ptwXY2->accuracy ) accuracy = ptwXY2->accuracy;
    n = ptwXY_new( ptwXY1->interpolation, NULL, biSectionMax, accuracy, length, overflowSize, status, ptwXY1->userFlag );
    if( n == NULL ) return( NULL );

    for( i = 0; ( i1 < n1 ) && ( i2 < n2 ); i++ ) {
        y = 0.;
        if( ptwXY1->points[i1].x <= ptwXY2->points[i2].x ) {
            n->points[i].x = ptwXY1->points[i1].x;
            if( fillWithFirst ) {
                y = ptwXY1->points[i1].y;
                if( i1 < ( ptwXY1->length - 1 ) ) {
                    x1 = ptwXY1->points[i1].x;
                    y1 = ptwXY1->points[i1].y;
                    x2 = ptwXY1->points[i1+1].x;
                    y2 = ptwXY1->points[i1+1].y; }
                else {
                    y1 = 0.;
                    y2 = 0.;
                }
            }
            if( ptwXY1->points[i1].x == ptwXY2->points[i2].x ) i2++;
            i1++; }
        else {
            n->points[i].x = ptwXY2->points[i2].x;
            if( fillWithFirst && ( ( y1 != 0. ) || ( y2 != 0. ) ) ) {
                if( ( *status = ptwXY_interpolatePoint( ptwXY1->interpolation, ptwXY2->points[i2].x, &y, x1, y1, x2, y2 ) ) != nfu_Okay ) {
                    ptwXY_free( n );
                    return( NULL );
                }
            }
            i2++;
        }
        n->points[i].y = y;
    }

    y = 0.;
    for( ; i1 < n1; i1++, i++ ) {
        n->points[i].x = ptwXY1->points[i1].x;
        if( fillWithFirst ) y = ptwXY1->points[i1].y;
        n->points[i].y = y;
    }
    for( ; i2 < n2; i2++, i++ ) {
        n->points[i].x = ptwXY2->points[i2].x;
        if( fillWithFirst && trim && ( n->points[i].x <= x2 ) ) {
            if( ( *status = ptwXY_interpolatePoint( ptwXY1->interpolation, n->points[i].x, &y, x1, y1, x2, y2 ) ) != nfu_Okay ) {
                ptwXY_free( n );
                return( NULL );
            }
        }
        n->points[i].y = y;
    }
    n->length = i;

    if( unionOptions & ptwXY_union_mergeClosePoints ) {
        if( ( *status = ptwXY_mergeClosePoints( n, 4 * DBL_EPSILON ) ) != nfu_Okay ) {
            ptwXY_free( n );
            return( NULL );
        }
    }
    return( n );
}
/*
************************************************************
*/
nfu_status ptwXY_scaleOffsetXAndY( ptwXYPoints *ptwXY, double xScale, double xOffset, double yScale, double yOffset ) {

    int64_t i1, length = ptwXY->length;
    ptwXYPoint *p1;
    nfu_status status;

    if( ptwXY->status != nfu_Okay ) return( ptwXY->status );
    if( xScale == 0 ) return( nfu_XNotAscending );

    if( ( status = ptwXY_simpleCoalescePoints( ptwXY ) ) != nfu_Okay ) return( status );

    for( i1 = 0, p1 = ptwXY->points; i1 < length; i1++, p1++ ) {
        p1->x = xScale * p1->x + xOffset;
        p1->y = yScale * p1->y + yOffset;
    }

    if( xScale < 0 ) {
        int64_t length_2 = length / 2;
        ptwXYPoint tmp, *p2;

        for( i1 = 0, p1 = ptwXY->points, p2 = &(ptwXY->points[length-1]); i1 < length_2; i1++ ) {
            tmp = *p1;
            *p1 = *p2;
            *p2 = tmp;
        }
    }

    return( ptwXY->status );
}

#if defined __cplusplus
}
#endif
