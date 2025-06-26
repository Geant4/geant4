/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <math.h>
#include <float.h>

#include "ptwXY.h"

typedef nfu_status (*interpolation_func)( statusMessageReporting *smr, ptwXYPoints *desc, double x1, double y1, 
        double x2, double y2, int depth );

static double ptwXY_flatInterpolationToLinear_eps( double px, double eps );
static nfu_status ptwXY_toOtherInterpolation2( statusMessageReporting *smr, ptwXYPoints *desc, ptwXYPoints *src, 
        interpolation_func func );
static nfu_status ptwXY_LogLogToLinLin( statusMessageReporting *smr, ptwXYPoints *desc, double x1, double y1, 
        double x2, double y2, int depth );
static nfu_status ptwXY_LogLinToLinLin( statusMessageReporting *smr, ptwXYPoints *desc, double x1, double y1, 
        double x2, double y2, int depth );
static nfu_status ptwXY_LinLogToLinLin( statusMessageReporting *smr, ptwXYPoints *desc, double x1, double y1, 
        double x2, double y2, int depth );
/*
************************************************************
*/
nfu_status ptwXY_interpolatePoint( statusMessageReporting *smr, ptwXY_interpolation interpolation, double x, double *y, 
        double x1, double y1, double x2, double y2 ) {

    nfu_status status = nfu_Okay;

    if( interpolation == ptwXY_interpolationOther ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_otherInterpolation, "Other interpolation not allow." );
        return( nfu_otherInterpolation );
    }
    if( ( x1 > x2 ) || ( x < x1 ) || ( x > x2 ) ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_invalidInterpolation, 
                "Interpolation point.x = %.17e not between function points x1 = %.17e and x2 = %.17e.", x, x1, x2 );
        return( nfu_invalidInterpolation );
    }
    if( y1 == y2 ) {
        *y = y1; }
    else if( x1 == x2 ) {
        *y = 0.5 * ( y1  + y2 ); }
    else if( x == x1 ) {
        *y = y1; }
    else if( x == x2 ) {
        *y = y2; }
    else {
        switch( interpolation ) {
        case ptwXY_interpolationLinLin :
            *y = ( y1 * ( x2 - x ) + y2 * ( x - x1 ) ) / ( x2 - x1 );
            break;
        case ptwXY_interpolationLinLog :
            if( ( x <= 0. ) || ( x1 <= 0. ) || ( x2 <= 0. ) ) {
                smr_setReportError2( smr, nfu_SMR_libraryID, nfu_invalidInterpolation, 
                        "For log(x), some x-values  less than equal to 0, x1 = %.17e, x = %.17e, x2 = %.17e.", x1, x, x2 );
                return( nfu_invalidInterpolation );
            }
            *y = ( y1 * log( x2 / x ) + y2 * log( x / x1 ) ) / log( x2 / x1 );
            break;
        case ptwXY_interpolationLogLin :
            if( ( y1 <= 0. ) || ( y2 <= 0. ) ) {
                smr_setReportError2( smr, nfu_SMR_libraryID, nfu_invalidInterpolation, 
                        "For log(y), some y-values less than equal to 0, y1 = %.17e, y2 = %.17e.", y1, y2 );
                return( nfu_invalidInterpolation );
            }
            *y = exp( ( log( y1 ) * ( x2 - x ) + log( y2 ) * ( x - x1 ) ) / ( x2 - x1 ) );
            break;
        case ptwXY_interpolationLogLog :
            if( ( x <= 0. ) || ( x1 <= 0. ) || ( x2 <= 0. ) ) {
                smr_setReportError2( smr, nfu_SMR_libraryID, nfu_invalidInterpolation, 
                        "For log(x), some x-values  less than equal to 0, x1 = %.17e, x = %.17e, x2 = %.17e.", x1, x, x2 );
                return( nfu_invalidInterpolation );
            }
            if( ( y1 <= 0. ) || ( y2 <= 0. ) ) {
                smr_setReportError2( smr, nfu_SMR_libraryID, nfu_invalidInterpolation, 
                        "For log(y), some y-values less than equal to 0, y1 = %.17e, y2 = %.17e.", y1, y2 );
                return( nfu_invalidInterpolation );
            }
            *y = exp( ( log( y1 ) * log( x2 / x ) + log( y2 ) * log( x / x1 ) ) / log( x2 / x1 ) );
            break;
        case ptwXY_interpolationFlat :
            *y = y1;
            break;
        default :
            smr_setReportError2( smr, nfu_SMR_libraryID, nfu_invalidInterpolation,
                    "Invalid interpolation token = %d.", interpolation );
            status = nfu_invalidInterpolation;
        }
    }
    return( status );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_flatInterpolationToLinear( statusMessageReporting *smr, ptwXYPoints *ptwXY, double lowerEps, double upperEps ) {

    int64_t i, length;
    double x;
    ptwXYPoints *n;
    ptwXYPoint *p1 = NULL, *p2 = NULL, *p3;

#define minEps 5e-16

    if( ptwXY_simpleCoalescePoints( smr, ptwXY ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }
        
    if( ptwXY->interpolation != ptwXY_interpolationFlat ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_invalidInterpolation, "Source interpolation not 'flat' but '%s'.",
                ptwXY->interpolationString );
        return( NULL );
    }

    if( ( lowerEps < 0 ) || ( upperEps < 0 ) || ( ( lowerEps == 0 ) && ( upperEps == 0 ) ) ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badInput, "Bad epsilons: lowerEps = %.17e and upperEps = %.17e.",
                lowerEps, upperEps );
        return( NULL );
    }
    if( ( lowerEps != 0 ) && ( lowerEps < minEps ) ) lowerEps = minEps;
    if( ( upperEps != 0 ) && ( upperEps < minEps ) ) upperEps = minEps;

    length = ptwXY->length * ( 1 + ( lowerEps == 0 ? 0 : 1 ) + ( lowerEps == 0 ? 0 : 1 ) );
    if( ( n = ptwXY_new( smr, ptwXY_interpolationLinLin, NULL, ptwXY->biSectionMax, ptwXY->accuracy, length, 
            ptwXY->overflowLength, ptwXY->userFlag ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }

    p3 = ptwXY->points;
    if( ptwXY->length > 0 ) {
        if( ptwXY_setValueAtX( smr, n, p3->x, p3->y ) != nfu_Okay ) goto Err;
    }
    for( i = 0; i < ptwXY->length; i++, p3++ ) {
        if( i > 1 ) {
            if( lowerEps > 0 ) {
                x = ptwXY_flatInterpolationToLinear_eps( p2->x, -lowerEps );
                if( x > p1->x ) {
                    if( ptwXY_setValueAtX( smr, n, x, p1->y ) != nfu_Okay ) goto Err;
                }
            }
            if( lowerEps == 0 ) if( ptwXY_setValueAtX( smr, n, p2->x, p1->y ) != nfu_Okay ) goto Err;
            if( upperEps == 0 ) if( ptwXY_setValueAtX( smr, n, p2->x, p2->y ) != nfu_Okay ) goto Err;
            if( upperEps > 0 ) {
                x = ptwXY_flatInterpolationToLinear_eps( p2->x, upperEps );
                if( x < p3->x ) {
                    if( ptwXY_setValueAtX( smr, n, x, p2->y ) != nfu_Okay ) goto Err;
                }
            }
        }
        p1 = p2;
        p2 = p3;
    }
    if( ptwXY->length > 1 ) {
        if( ( lowerEps != 0 ) && ( p1->y != p2->y ) ) {
            x = ptwXY_flatInterpolationToLinear_eps( p2->x, -lowerEps );
            if( x > p1->x ) {
                if( ptwXY_setValueAtX( smr, n, x, p1->y ) != nfu_Okay ) goto Err;
            }
        }
        if( ptwXY_setValueAtX( smr, n, p2->x, p2->y ) != nfu_Okay ) goto Err;
    }

    return( n );

Err:
    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    ptwXY_free( n );
    return( NULL );

#undef minEps
}
/*
************************************************************
*/
static double ptwXY_flatInterpolationToLinear_eps( double px, double eps ) {

    double x;

    if( px < 0 ) {
        x = ( 1 - eps ) * px; }
    else if( px > 0 ) {
        x = ( 1 + eps ) * px; }
    else {
        x = eps;
    }
    return( x );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_toOtherInterpolation( statusMessageReporting *smr, ptwXYPoints *ptwXY, 
        ptwXY_interpolation interpolationTo, double accuracy ) {
/*
*   This function only works when 'ptwXY->interpolation == interpolationTo' or when interpolationTo is ptwXY_interpolationLinLin.
*/
    int i1, logX = 0, logY = 0;
    ptwXYPoints *n1;
    interpolation_func func = NULL;

    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( NULL );
    }

    if( ptwXY->interpolation == interpolationTo ) {
        if( ( n1 = ptwXY_clone( smr, ptwXY ) ) == NULL ) smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( n1 ); }
    else {
        if( interpolationTo == ptwXY_interpolationLinLin ) {
            switch( ptwXY->interpolation ) {
            case ptwXY_interpolationLogLog :
                logX = logY = 1;
                func = ptwXY_LogLogToLinLin; break;
            case ptwXY_interpolationLogLin :
                logY = 1;
                func = ptwXY_LogLinToLinLin; break;
            case ptwXY_interpolationLinLog :
                logX = 1;
                func = ptwXY_LinLogToLinLin; break;
            case ptwXY_interpolationLinLin :        /* Stops compilers from complaining. */
            case ptwXY_interpolationFlat :
            case ptwXY_interpolationOther :
                break;
            }
        }
    }
    if( func == NULL ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_unsupportedInterpolationConversion,
                "Interpolation conversion from '%s' to %d not supported.", ptwXY->interpolationString, interpolationTo );
        return( NULL );
    }

    if( ( logX != 0 ) || ( logY != 0 ) ) {
        ptwXYPoint *point;

        if( ptwXY_simpleCoalescePoints( smr, ptwXY ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( NULL );
        }

        for( i1 = 0, point = ptwXY->points; i1 < ptwXY->length; ++i1, ++point ) {
            if( ( logX != 0 ) && ( point->x <= 0 ) ) {
                smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badLogValue, "At index %d x-value %.16e for log <= 0.", 
                        (int) i1, point->x );
                return( NULL );
            }
            if( ( logY != 0 ) && ( point->y <= 0 ) ) {
                smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badLogValue, "At index %d y-value %.16e for log <= 0.", 
                        (int) i1, point->y );
                return( NULL );
            }
        }
    }

    if( ( n1 = ptwXY_cloneToInterpolation( smr, ptwXY, interpolationTo ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }
    n1->accuracy = ptwXY_limitAccuracy( accuracy );

    if( ptwXY_toOtherInterpolation2( smr, n1, ptwXY, func ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        n1 = ptwXY_free( n1 ); }
    else {
        if( n1->accuracy < ptwXY->accuracy ) n1->accuracy = ptwXY->accuracy;
    }
    return( n1 );
}
/*
************************************************************
*/
static nfu_status ptwXY_toOtherInterpolation2( statusMessageReporting *smr, ptwXYPoints *desc, ptwXYPoints *src, 
        interpolation_func func ) {

    nfu_status status;
    int64_t i1;
    double x1, y1, x2, y2;

    if( ( status = ptwXY_simpleCoalescePoints( smr, src ) ) != nfu_Okay ) return( status );

    x1 = src->points[0].x;
    y1 = src->points[0].y;
    for( i1 = 1; i1 < src->length; i1++ ) {
        x2 = src->points[i1].x;
        y2 = src->points[i1].y;
        if( ( x1 != x2 ) && ( y1 != y2 ) ) {
            if( ( status = func( smr, desc, x1, y1, x2, y2, 0 ) ) != nfu_Okay ) break;
        }
        x1 = x2;
        y1 = y2;
    }
    return( status );
}
/*
************************************************************
*/
static nfu_status ptwXY_LogLogToLinLin( statusMessageReporting *smr, ptwXYPoints *desc, double x1, double y1, 
        double x2, double y2, int depth ) {

    nfu_status status = nfu_Okay;
    double x, y, u, u2 = x2 / x1, v2 = y2 / y1, logYXs, logXs = log( u2 ), logYs = log( v2 ), vLin, vLog, w;

    logYXs = logYs / logXs;

    if( depth > ptwXY_maxBiSectionMax ) return( nfu_Okay );
    if( fabs( logYXs  - 1 ) < 1e-5 ) {
        u = 0.5 * ( 1 + u2 );
        w = ( logYXs  - 1 ) * logXs;
        vLog = u * ( 1. + w * ( 1 + 0.5 * w ) ); }
    else {
        if( u2 > 10 ) {
            u = sqrt( u2 ); }
        else {
            u = logYXs * ( u2 - v2 ) / ( ( 1 - logYXs ) * ( v2 - 1 ) );
        }
        vLog = pow( u, logYXs );
    }
    vLin = ( u2 - u + v2 * ( u - 1 ) ) / ( u2 - 1 );
    if( fabs( vLog - vLin ) <= ( vLog * desc->accuracy ) ) return( status );
    x = x1 * u;
    y = y1 * vLog;
    if( ( status = ptwXY_setValueAtX( smr, desc, x, y ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_LogLogToLinLin( smr, desc, x1, y1, x, y, depth + 1 ) ) != nfu_Okay ) return( status );
    return( ptwXY_LogLogToLinLin( smr, desc, x, y, x2, y2, depth + 1 ) );
}
/*
************************************************************
*/
static nfu_status ptwXY_LogLinToLinLin( statusMessageReporting *smr, ptwXYPoints *desc, 
        double x1, double y1, double x2, double y2, int depth ) {

    nfu_status status = nfu_Okay;
    double x, y, logYs = log( y2 / y1 ), yLinLin;

    if( depth > ptwXY_maxBiSectionMax ) return( nfu_Okay );
    x = ( x2 - x1 ) / ( y2 - y1 ) * ( ( y2 - y1 ) / logYs - y1 ) + x1;
    y = y1 * exp( logYs / ( x2 - x1 ) * ( x - x1 ) );
    yLinLin = ( y1 * ( x2 - x ) + y2 * ( x - x1 ) ) / ( x2 - x1 );
    if( fabs( y - yLinLin ) <= ( y * desc->accuracy ) ) return( status );
    if( ( status = ptwXY_setValueAtX( smr, desc, x, y ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_LogLinToLinLin( smr, desc, x1, y1, x, y, depth + 1 ) ) != nfu_Okay ) return( status );
    return( ptwXY_LogLinToLinLin( smr, desc, x, y, x2, y2, depth + 1 ) );
}
/*
************************************************************
*/
static nfu_status ptwXY_LinLogToLinLin( statusMessageReporting *smr, ptwXYPoints *desc, double x1, double y1, 
        double x2, double y2, int depth ) {

    nfu_status status = nfu_Okay;
    double x = sqrt( x2 * x1 ), y, logXs = log( x2 / x1 ), yLinLin;

    if( depth > ptwXY_maxBiSectionMax ) return( nfu_Okay );
#if 0 /* The next line is very unstable at determineing x. Initial x must be chosen better. */
    x = ( y1 * x2 - y2 * x1 ) / ( y1 * logXs + ( y2 - y1 ) * ( log( x / x1 ) - 1 ) );
#endif
    y = ( y2 - y1 ) * log( x / x1 ) / logXs + y1;
    yLinLin = ( y1 * ( x2 - x ) + y2 * ( x - x1 ) ) / ( x2 - x1 );
    if( fabs( y - yLinLin ) <= fabs( y * desc->accuracy ) ) return( status );
    if( ( status = ptwXY_setValueAtX( smr, desc, x, y ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_LinLogToLinLin( smr, desc, x1, y1, x, y, depth + 1 ) ) != nfu_Okay ) return( status );
    return( ptwXY_LinLogToLinLin( smr, desc, x, y, x2, y2, depth + 1 ) );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_toUnitbase( statusMessageReporting *smr, ptwXYPoints *ptwXY, int scaleRange ) {

    int64_t i;
    ptwXYPoints *n;
    ptwXYPoint *p;
    double domainMin, domainMax, dx, inverseDx;

    if( ptwXY->length < 2 ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_tooFewPoints, "Too few points %d", (int) ptwXY->length );
        return( NULL );
    }
    if( ( n = ptwXY_clone( smr, ptwXY ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }

    domainMin = n->points[0].x;
    domainMax = n->points[n->length-1].x;
    dx = domainMax - domainMin;
    inverseDx = 1. / dx;
    for( i = 0, p = n->points; i < n->length; i++, p++ ) {
        p->x = ( p->x - domainMin ) * inverseDx;
        if( scaleRange ) p->y = p->y * dx;
    }
    n->points[n->length-1].x = 1.;                          /* Make sure last point is realy 1. */
    return( n );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_fromUnitbase( statusMessageReporting *smr, ptwXYPoints *ptwXY, double domainMin, double domainMax,
        int scaleRange ) {

    int64_t i, length;
    ptwXYPoints *n;
    ptwXYPoint *p, *p2;
    double dx, inverseDx, xLast = 0.;

    if( ptwXY->length < 2 ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_tooFewPoints, "Too few points %d", (int) ptwXY->length );
        return( NULL );
    }
    if( ( n = ptwXY_clone( smr, ptwXY ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }

    dx = domainMax - domainMin;
    inverseDx = 1. / dx;
    length = n->length;
    for( i = 0, p2 = p = n->points; i < length; ++i, ++p ) {
        p2->x = p->x * dx + domainMin;
        if( i > 0 ) {
            if( fabs( p2->x - xLast ) <= 10. * DBL_EPSILON * ( fabs( p2->x ) + fabs( xLast ) ) ) {
                --(n->length);
                continue;
            }
        }
        if( scaleRange ) p2->y = p->y * inverseDx;
        xLast = p2->x;
        ++p2;
    }
    n->points[n->length-1].x = domainMax;                          /* Make sure last point is realy domainMax. */
    return( n );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_unitbaseInterpolate( statusMessageReporting *smr, double w, double w1, ptwXYPoints *ptwXY1, 
        double w2, ptwXYPoints *ptwXY2, int scaleRange ) {
/* 
*   Should we not be checking the interpolation members???????
*/
    int64_t i;
    ptwXYPoints *n1 = NULL, *n2 = NULL, *a = NULL, *r = NULL;
    ptwXYPoint *p;
    double f, g, domainMin, domainMax;

    if( ( w < w1 ) || ( w > w2 ) ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_XOutsideDomain, "W value outside w-domain: (%.15e, %.15e)",
                w1, w2 );
        return( NULL );
    }
    if( w == w1 ) {
        if( ( n1 = ptwXY_clone( smr, ptwXY1 ) ) == NULL ) smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via (source1)." );
        return( n1 );
    }
    if( w == w2 ) {
        if( ( n1 = ptwXY_clone( smr, ptwXY2 ) ) == NULL ) smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via (source2)." );
        return( n1 );
    }
    if( ( n1 = ptwXY_toUnitbase( smr, ptwXY1, scaleRange ) ) == NULL ) goto Err;
    if( ( n2 = ptwXY_toUnitbase( smr, ptwXY2, scaleRange ) ) == NULL ) goto Err;
    f = ( w - w1 ) / ( w2 - w1 );
    g = 1. - f;
    for( i = 0, p = n1->points; i < n1->length; i++, p++ ) p->y *= g;
    for( i = 0, p = n2->points; i < n2->length; i++, p++ ) p->y *= f;
    if( ( a = ptwXY_add_ptwXY( smr, n1, n2 ) ) == NULL ) goto Err;

    domainMin = g * ptwXY1->points[0].x + f * ptwXY2->points[0].x;
    domainMax = g * ptwXY1->points[ptwXY1->length-1].x + f * ptwXY2->points[ptwXY2->length-1].x;
    if( ( r = ptwXY_fromUnitbase( smr, a, domainMin, domainMax, scaleRange ) ) == NULL ) goto Err;
    ptwXY_free( n1 );
    ptwXY_free( n2 );
    ptwXY_free( a );
    return( r );

Err:
    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    if( n1 != NULL ) ptwXY_free( n1 );
    if( n2 != NULL ) ptwXY_free( n2 );
    if( a  != NULL ) ptwXY_free( a );
    return( NULL );
}
