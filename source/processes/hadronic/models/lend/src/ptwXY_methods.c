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

static nfu_status ptwXY_clip2( statusMessageReporting *smr, ptwXYPoints *ptwXY1, double y, double x1, double y1, double x2, double y2 );
static double ptwXY_thicken_linear_dx( int sectionSubdivideMax, double dDomainMax, double x1, double x2 );
static nfu_status ptwXY_thin2( statusMessageReporting *smr, ptwXYPoints *thinned, char *thin, double accuracy, int64_t i1, int64_t i2 );
/*
************************************************************
*/
nfu_status ptwXY_clip( statusMessageReporting *smr, ptwXYPoints *ptwXY1, double rangeMin, double rangeMax ) {
/*
    This function acts oddly for xy = [ [ 1, 0 ], [ 3, -2 ], [ 4, 1 ] ] and rangeMin = 0.2, why???????
    This function probably only works for linear, linear interpolation (mainly because of ptwXY_clip2).
*/
    int64_t i, j, n;
    double x2, y2, _rangeMin, _rangeMax;
    ptwXYPoints *clipped;
    ptwXYPoint *points;

    if( ptwXY_simpleCoalescePoints( smr, ptwXY1 ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( ptwXY1->status );
    }

    if( ptwXY1->interpolation == ptwXY_interpolationOther ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_otherInterpolation, "Other interpolation not allowed." );
        return( nfu_otherInterpolation );
    }
    n = ptwXY1->length;
    if( n > 0 ) {
        i = 0;
        if( ptwXY_range( smr, ptwXY1, &_rangeMin, &_rangeMax ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( nfu_Error );
        }
        if( _rangeMax < rangeMin ) i = 1;
        if( _rangeMin > rangeMax ) i = 1;
        if( i == 1 ) {
            if( ptwXY_clear( smr, ptwXY1 ) != nfu_Okay ) smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( ptwXY1->status );
        }
    }
    if( n == 1 ) {
        y2 = ptwXY1->points[0].y;
        if( y2 < rangeMin ) {
            ptwXY1->points[0].y = rangeMin; }
        else if( y2 > rangeMax ) {
            ptwXY1->points[0].y = rangeMax;
        } }
    else if( n > 1 ) {
        if( ( clipped = ptwXY_new( smr, ptwXY1->interpolation, ptwXY1->interpolationString, 
                ptwXY1->biSectionMax, ptwXY1->accuracy, n, 10, ptwXY1->userFlag ) ) == NULL ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( ptwXY1->status );
        }
        for( i = 0; i < n; i++ ) {
            x2 = ptwXY1->points[i].x;
            y2 = ptwXY1->points[i].y;
            if( y2 < rangeMin ) {
                if( i > 0 ) {
                    points = ptwXY_getPointAtIndex_Unsafely( clipped, clipped->length - 1 );
                    if( points->y > rangeMin ) {
                        if( ptwXY_clip2( smr, clipped, rangeMin, points->x, points->y, x2, y2 ) != nfu_Okay ) goto Err;
                    }
                }
                if( ptwXY_setValueAtX( smr, clipped, x2, rangeMin ) != nfu_Okay ) goto Err;
                j = i;
                for( i++; i < n; i++ ) if( !( ptwXY1->points[i].y < rangeMin ) ) break;
                if( i < n ) {
                    x2 = ptwXY1->points[i].x;
                    y2 = ptwXY1->points[i].y;
                    if( ptwXY_clip2( smr, clipped, rangeMin, ptwXY1->points[i-1].x, ptwXY1->points[i-1].y, x2, y2 ) != nfu_Okay ) goto Err;
                    if( y2 > rangeMax ) {
                        if( ptwXY_clip2( smr, clipped, rangeMax, ptwXY1->points[i-1].x, ptwXY1->points[i-1].y, x2, y2 ) != nfu_Okay ) goto Err;
                    } }
                else if( j != n - 1 ) {
                    if( ptwXY_setValueAtX( smr, clipped, ptwXY1->points[n - 1].x, rangeMin ) != nfu_Okay ) goto Err;
                }
                i--; }
            else if( y2 > rangeMax ) {
                if( i > 0 ) {
                    points = ptwXY_getPointAtIndex_Unsafely( clipped, clipped->length - 1 );
                    if( points->y < rangeMax ) {
                        if( ptwXY_clip2( smr, clipped, rangeMax, points->x, points->y, x2, y2 ) != nfu_Okay ) goto Err;
                    }
                }
                if( ptwXY_setValueAtX( smr, clipped, x2, rangeMax ) != nfu_Okay ) goto Err;
                j = i;
                for( i++; i < n; i++ ) if( !( ptwXY1->points[i].y > rangeMax ) ) break;
                if( i < n ) {
                    x2 = ptwXY1->points[i].x;
                    y2 = ptwXY1->points[i].y;
                    if( ptwXY_clip2( smr, clipped, rangeMax, ptwXY1->points[i-1].x, ptwXY1->points[i-1].y, x2, y2 ) != nfu_Okay ) goto Err;
                    if( y2 < rangeMin ) {
                        if( ptwXY_clip2( smr, clipped, rangeMin, ptwXY1->points[i-1].x, ptwXY1->points[i-1].y, x2, y2 ) != nfu_Okay ) goto Err;
                    } }
                else if( j != n - 1 ) {
                    if( ptwXY_setValueAtX( smr, clipped, ptwXY1->points[n - 1].x, rangeMax ) != nfu_Okay ) goto Err;
                }
                i--; }
            else {
                if( ptwXY_setValueAtX( smr, clipped, x2, y2 ) != nfu_Okay ) goto Err;
            }
        }
        if( ptwXY_simpleCoalescePoints( smr, clipped ) != nfu_Okay ) goto Err;
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
    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    ptwXY_free( clipped );
    return( ptwXY1->status );
}
/*
************************************************************
*/
static nfu_status ptwXY_clip2( statusMessageReporting *smr, ptwXYPoints *clipped, double y, double x1, double y1, double x2, double y2 ) {

    double x;

    x = ( y - y1 ) * ( x2 - x1 ) / ( y2 - y1 ) + x1;
    if( x <= x1 ) {
        x = x1; }
    else if( x >= x2 ) {
        x = x1; }
    else {
        ptwXY_setValueAtX( smr, clipped, x, y );
    }
    return( clipped->status  );
}
/*
************************************************************
*/
nfu_status ptwXY_thicken( statusMessageReporting *smr, ptwXYPoints *ptwXY1, int sectionSubdivideMax, 
        double dDomainMax, double fDomainMax ) {

    double x1, x2 = 0., y1, y2 = 0., fx = 1.1, x, dx, dxp, lfx, y;    /* fx initialized so compilers want complain. */
    int64_t i, notFirstPass = 0;
    int nfx, nDone, doLinear;

    if( ptwXY1->interpolation == ptwXY_interpolationOther ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_otherInterpolation, "Other interpolation not allowed." );
        return( nfu_otherInterpolation );
    }
    if( ( sectionSubdivideMax < 1 ) || ( dDomainMax < 0. ) || ( fDomainMax < 1. ) ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_Error, 
                "ptwXY_thicken: One or more of the following not satisfied: sectionSubdivideMax = %d > 0, dDomainMax = %e >= 0.0, fDomainMax = %e >= 1.0",
                sectionSubdivideMax, dDomainMax, fDomainMax );
        return( nfu_badInput );
    }
    if( sectionSubdivideMax > ptwXY_sectionSubdivideMax ) sectionSubdivideMax = ptwXY_sectionSubdivideMax;
    if( ptwXY_simpleCoalescePoints( smr, ptwXY1 ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( ptwXY1->status );
    }
    for( i = ptwXY1->length - 1; i >= 0; i-- ) {
        x1 = ptwXY1->points[i].x;
        y1 = ptwXY1->points[i].y;
        if( notFirstPass ) {
            dx = ptwXY_thicken_linear_dx( sectionSubdivideMax, dDomainMax, x1, x2 );

            if( x1 == 0. ) {
                doLinear = 1; }
            else {
                fx = x2 / x1;
                if( fx > 0. ) {
                    lfx = log( fx );
                    if( fDomainMax == 1. ) {
                        nfx = sectionSubdivideMax; }
                    else {
                        nfx = ( (int) ( lfx / log( fDomainMax ) ) ) + 1;
                        if( nfx > sectionSubdivideMax ) nfx = sectionSubdivideMax;
                    }
                    if( nfx > 0 ) fx = exp( lfx / nfx );
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
                    dx = ptwXY_thicken_linear_dx( sectionSubdivideMax - nDone, dDomainMax, x, x2 );
                    if( dx <= ( fx - 1 ) * x ) {
                        dxp = dx;
                        doLinear = 1;
                        continue;
                    }
                    dxp = ( fx - 1. ) * x;
                    x *= fx;
                }
                if( ( x2 - x ) < 0.05 * fabs( dxp ) ) break;
                if( ( ptwXY1->status = ptwXY_interpolatePoint( smr, ptwXY1->interpolation, x, &y, x1, y1, x2, y2 ) ) != nfu_Okay ) {
                    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
                    return( ptwXY1->status );
                }
                if( ( ptwXY1->status = ptwXY_setValueAtX( smr, ptwXY1, x, y ) ) != nfu_Okay ) {
                    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
                    return( ptwXY1->status );
                }
                nDone++;
            }
        }
        notFirstPass = 1;
        x2 = x1;
        y2 = y1;
    }
    return( ptwXY1->status );
}
/*
************************************************************
*/
static double ptwXY_thicken_linear_dx( int sectionSubdivideMax, double dDomainMax, double x1, double x2 ) {

    int ndx;
    double dx = x2 - x1, dndx;

    if( dDomainMax == 0. ) {
        dx = ( x2 - x1 ) / sectionSubdivideMax; }
    else {
        dndx = dx / dDomainMax;
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
ptwXYPoints *ptwXY_thin( statusMessageReporting *smr, ptwXYPoints *ptwXY1, double accuracy ) {

    int64_t i, j, length = ptwXY1->length;
    ptwXYPoints *thinned = NULL;
    double y1, y2, y3, accuracyNew;
    char *thin = NULL;

    if( length < 3 ) {                          /* Logic below requires at least 2 points. */
        if( ( thinned = ptwXY_clone( smr, ptwXY1 ) ) == NULL )
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( thinned );
    }

    if( ptwXY_simpleCoalescePoints( smr, ptwXY1 ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }

    if( ptwXY1->interpolation == ptwXY_interpolationOther ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_otherInterpolation, "Other interpolation not allowed." );
        return( NULL );
    }

    accuracy = ptwXY_limitAccuracy( accuracy );
    accuracyNew = accuracy;
    if( accuracyNew < ptwXY1->accuracy ) accuracyNew = ptwXY1->accuracy;
    if( ( thinned = ptwXY_new( smr, ptwXY1->interpolation, ptwXY1->interpolationString, 
            ptwXY1->biSectionMax, accuracyNew, length, ptwXY1->overflowLength, ptwXY1->userFlag ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }

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
        if( ( thin = (char *) smr_malloc2( smr, (size_t) length, 1, "thin" ) ) == NULL ) goto Err2;
        if( ptwXY_thin2( smr, thinned, thin, accuracy, 0, length - 1 ) != nfu_Okay ) goto Err1;
        for( j = 1; j < length; j++ ) if( thin[j] != 0 ) break;
        for( i = j + 1; i < length; i++ ) {
            if( thin[i] == 0 ) {
                thinned->points[j] = thinned->points[i];
                j++;
            }
        }
        smr_freeMemory2( thin );
    }
    thinned->length = j;

    return( thinned );

Err1:
    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
Err2:
    ptwXY_free( thinned );
    if( thin != NULL ) smr_freeMemory2( thin );
    return( NULL );
}
/*
************************************************************
*/
static nfu_status ptwXY_thin2( statusMessageReporting *smr, ptwXYPoints *thinned, char *thin, double accuracy, int64_t i1, int64_t i2 ) {

    int64_t i, iMax = 0;
    double y, s, dRange, dRangeMax = 0., dRangeR, dRangeRMax = 0;
    double x1 = thinned->points[i1].x, y1 = thinned->points[i1].y, x2 = thinned->points[i2].x, y2 = thinned->points[i2].y;
    nfu_status status = nfu_Okay;

    if( i1 + 1 >= i2 ) return( nfu_Okay );
    for( i = i1 + 1; i < i2; i++ ) {
        if( ( thinned->status = ptwXY_interpolatePoint( smr, thinned->interpolation, thinned->points[i].x, &y, x1, y1, x2, y2 ) ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( thinned->status );
        }
        s = 0.5 * ( fabs( y ) + fabs( thinned->points[i].y ) );
        dRange = fabs( y - thinned->points[i].y );
        dRangeR = 0;
        if( s != 0 ) dRangeR = dRange / s;
        if( ( dRangeR > dRangeRMax ) || ( ( dRangeR >= 0.9999 * dRangeRMax ) && ( dRange > dRangeMax ) ) ) {
            iMax = i;							/* The choice of 0.9999 is not exact science. */
            if( dRange > dRangeMax ) dRangeMax = dRange;
            if( dRangeR > dRangeRMax ) dRangeRMax = dRangeR;
        }
    }
    if( dRangeRMax < accuracy ) {
        for( i = i1 + 1; i < i2; i++ ) thin[i] = 1; }
    else {
        if( ( status = ptwXY_thin2( smr, thinned, thin, accuracy, i1, iMax ) ) != nfu_Okay ) return( status );
        status = ptwXY_thin2( smr, thinned, thin, accuracy, iMax, i2 );
    }
    return( status );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_thinDomain( statusMessageReporting *smr, ptwXYPoints *ptwXY1, double epsilon ) {
/*
*   Thins domain points that are closer the '0.5 * (x[i+1] + x[i]) * epsilon'.
*/
    int64_t i1, i2, length = ptwXY1->length, lengthm1 = length - 1, thinnedLength = 0;
    ptwXYPoint *points;
    ptwXYPoints *thinned = NULL;
    double x1, x2, x3, dx, y2, half_epsilon = 0.5 * epsilon;

    if( ptwXY1->interpolation == ptwXY_interpolationFlat ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_invalidInterpolation, "Flat interpolation not allowed." );
        return( NULL );
    }

    if( ptwXY1->interpolation == ptwXY_interpolationOther ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_otherInterpolation, "Other interpolation not allowed." );
        return( NULL );
    }

    if( ptwXY_simpleCoalescePoints( smr, ptwXY1 ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }

    if( length > 1 ) {
        if( ( ptwXY1->points[length-1].x - ptwXY1->points[0].x ) < half_epsilon * ( fabs( ptwXY1->points[0].x ) + fabs( ptwXY1->points[0].x ) ) ) {
            smr_setReportError2( smr, nfu_SMR_libraryID, nfu_Error, "Domain (%.17e, %.17e) is less than epsilon = %.17e.", 
                    ptwXY1->points[0].x, ptwXY1->points[length-1].x, epsilon );
            return( NULL );
        }
    }

    if( ( length <= 2 ) || ( epsilon < 2 * DBL_EPSILON ) ) {
        if( ( thinned = ptwXY_clone( smr, ptwXY1 ) ) == NULL ) smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( thinned );
    }

    if( ( thinned = ptwXY_new( smr, ptwXY1->interpolation, ptwXY1->interpolationString, ptwXY1->biSectionMax, ptwXY1->accuracy, 
            length, ptwXY1->overflowAllocatedSize, ptwXY1->userFlag ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }

    points = thinned->points;
    *points = ptwXY1->points[0];
    ++points;
    ++thinnedLength;
    x1 = ptwXY1->points[0].x;
    x2 = x3 = x1;                                                       /* To stop some compilers from printing a warning. */
    for( i1 = 1; i1 < lengthm1; i1 = i2 ) {
        for( i2 = i1; i2 < length; ++i2 ) {                             /* Find next x3 that is  epsilon * ( x1 + x2 ) / 2 above x1. */
            x3 = ptwXY1->points[i2].x;
            if( ( x3 - x1 ) >= half_epsilon * ( fabs( x1 ) + fabs( x3 ) ) ) break;
            x2 = x3;
        }
        if( i1 == i2 ) {
            y2 = ptwXY1->points[i2].y;
            x2 = x3; }
        else {
            if( ( x3 - x1 ) > ( epsilon * ( fabs( x1 ) + fabs( x3 ) ) ) ) {
                dx = fabs( x2 * epsilon );
                x2 = x1 + dx;
                if( ptwXY_getValueAtX( smr, ptwXY1, x2, &y2 ) != nfu_Okay ) {
                    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
                    ptwXY_free( thinned );
                    return( NULL );
                }
                --i2; }
            else {
                if( i2 == length ) break;
                x2 = x3;
                y2 = ptwXY1->points[i2].y;
            }
        }
        points->x = x2;
        points->y = y2;
        ++points;
        ++thinnedLength;
        x1 = x2;
        ++i2;
    }

    x3 = ptwXY1->points[lengthm1].x;
    x2 = thinned->points[thinnedLength-1].x;
    if( ( x3 - x2 ) < ( half_epsilon * ( fabs( x2 ) + fabs( x3 ) ) ) ) {
        --points;
        --thinnedLength;
        x1 = thinned->points[thinnedLength-1].x;
        if( ( x3 - x1 ) > ( epsilon * ( fabs( x1 ) + fabs( x2 ) ) ) ) {
            dx = fabs( x3 * epsilon );
            x2 = x3 - dx;
            if( ptwXY_getValueAtX( smr, ptwXY1, x2, &y2 ) != nfu_Okay ) {
                smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
                ptwXY_free( thinned );
                return( NULL );
            }
            points->x = x2;
            points->y = y2;
            ++points;
            ++thinnedLength;
        }
    }
    points->x = x3;
    points->y = ptwXY1->points[lengthm1].y;
    ++thinnedLength;
    thinned->length = thinnedLength;

    return( thinned );
}
/*
************************************************************
*/
nfu_status ptwXY_trim( statusMessageReporting *smr, ptwXYPoints *ptwXY ) {
/*
c   Remove extra zeros at beginning and end.
*/
    int64_t i, i1, i2;

    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid destination." );
        return( ptwXY->status );
    }

    if( ptwXY_simpleCoalescePoints( smr, ptwXY ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( ptwXY->status );
    }

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
ptwXYPoints *ptwXY_union( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, int unionOptions ) {

    int64_t overflowSize, i, i1 = 0, i2 = 0, n1 = ptwXY1->length, n2 = ptwXY2->length, length;
    int fillWithFirst = unionOptions & ptwXY_union_fill, trim = unionOptions & ptwXY_union_trim;
    ptwXYPoints *n;
    double x1 = 0., x2 = 0., y1 = 0., y2 = 0., y, biSectionMax, accuracy;
/*
*   Many other routines use the fact that ptwXY_union calls ptwXY_coalescePoints for ptwXY1 and ptwXY2 so do not change it.
*/
    if( ptwXY_simpleCoalescePoints( smr, ptwXY1 ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }
    if( ptwXY_simpleCoalescePoints( smr, ptwXY2 ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }

    if( ptwXY1->interpolation == ptwXY_interpolationOther ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_otherInterpolation, "Other interpolation not allowed for source1." );
        return( NULL );
    }
    if( ptwXY2->interpolation == ptwXY_interpolationOther ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_otherInterpolation, "Other interpolation not allowed for source2." );
        return( NULL );
    }

    if( ( n1 == 1 ) || ( n2 == 1 ) ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_tooFewPoints, 
                "Too few point in one of the sources: len( source1 ) = %d, len( source2 ) = %d", (int) n1, (int) n2 );
        return( NULL );
    }
    if( trim ) {
        if( n1 > 0 ) {
            if( n2 > 0 ) {
                if( ptwXY1->points[0].x < ptwXY2->points[0].x ) {
                    while( i1 < n1 ) {
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
                    while( i2 < n2 ) {
                        if( ptwXY2->points[i2].x >= ptwXY1->points[0].x ) break;
                        i2++;
                    }
                }
                if( ptwXY1->points[n1-1].x > ptwXY2->points[n2-1].x ) {
                    while( i1 < n1 ) {
                        if( ptwXY1->points[n1-1].x <= ptwXY2->points[n2-1].x ) break;
                        n1--;
                    } }
                else {
                    while( i2 < n2 ) {
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
    n = ptwXY_new( smr, ptwXY1->interpolation, NULL, biSectionMax, accuracy, length, overflowSize, ptwXY1->userFlag );
    if( n == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }

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
                if( ptwXY_interpolatePoint( smr, ptwXY1->interpolation, ptwXY2->points[i2].x, &y, x1, y1, x2, y2 ) != nfu_Okay ) {
                    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
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
            if( ptwXY_interpolatePoint( smr, ptwXY1->interpolation, n->points[i].x, &y, x1, y1, x2, y2 ) != nfu_Okay ) {
                smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
                ptwXY_free( n );
                return( NULL );
            }
        }
        n->points[i].y = y;
    }
    n->length = i;

    if( unionOptions & ptwXY_union_mergeClosePoints ) {
        if( ptwXY_mergeClosePoints( smr, n, 4 * DBL_EPSILON ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            ptwXY_free( n );
            return( NULL );
        }
    }
    return( n );
}
/*
************************************************************
*/
nfu_status ptwXY_scaleOffsetXAndY( statusMessageReporting *smr, ptwXYPoints *ptwXY, double xScale, double xOffset, 
        double yScale, double yOffset ) {

    int64_t i1, length = ptwXY->length;
    ptwXYPoint *p1;
    nfu_status status;

    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid destination." );
        return( ptwXY->status );
    }

    if( xScale == 0 ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_XNotAscending, "xScale is 0 that will cause a non-ascending domain." );
        return( ptwXY->status = nfu_XNotAscending );
    }

    if( ( status = ptwXY_simpleCoalescePoints( smr, ptwXY ) ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( status );
    }

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
/*
************************************************************
*/
nfu_status ptwXY_scaleAndOffsetDomainWith_ptwXYs( statusMessageReporting *smr, ptwXYPoints *ptwXY, ptwXYPoints *offsetXY, ptwXYPoints *slopeXY, 
                int skipLastPoint ) {
/*
    This function only modifies ptwXY over the domain defined by offsetXY via the equation ptwXY = offsetXY + slopeXY * ptwXY. 
    slopeXY must have the same domain as offsetXY.  If skipLastPoint is non-zero, the last point modified is reverted by to its 
    origial value.  This allows one to call this function multiple times with abutting domains for the offset and slope without 
    getting weird jumps at the boundaries. For example, suppose ptwXY is initially defined to be 1 in the domain [0,10]. Let,
    ptwXY have the points [ [ 0, 1 ], [ 5, 1 ], [ 10, 1 ] ]. Calling this function with an offset and slope with the points
    [ [ 0, 0 ], [ 5, 0 ] ] and [ 0, 2 ], [ 5, 2 ], respectively would return ptwXY as [ [ 0, 2 ], [ 5, 2 ], [ 10, 1 ] ]
    if skipLastPoint is 0 and [ [ 0, 2 ], [ 5, 1 ], [ 10, 1 ] ] otherwise. Now calling the returned ptwXYPoints instances with an 
    offset and slope with the points [ [ 5, 0 ], [ 10, 0 ] ] and [ 5, 2 ], [ 10, 2 ], respectively would return ptwXY as 
    [ [ 0, 2 ], [ 5, 4 ], [ 10, 2 ] ] and [ [ 0, 2 ], [ 5, 2 ], [ 10, 2 ] ], respectively.
*/

    int64_t i1;
    ptwXYPoint *p1;
    nfu_status status1, status2;
    double offsetXYMin, offsetXYMax, slopeXYMin, slopeXYMax, domainMin, domainMax, domainMinXY, domainMaxXY;
    ptwXYPoints *ptwXY2 = NULL, *offsetXY2 = NULL, *slopeXY2 = NULL;
    ptwXYPoints *mulXY = NULL, *addXY = NULL;

    if( ptwXY_simpleCoalescePoints( smr, ptwXY ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( ptwXY->status );
    }

    status1 = ptwXY_domainMin( smr, offsetXY, &offsetXYMin );       /* Also verifies that offsetXY has no issues. */
    if( ( status1 != nfu_Okay ) && ( status1 != nfu_empty ) ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( status1 );
    }
    ptwXY_domainMax( smr, offsetXY, &offsetXYMax );                 /* If ptwXY_domainMin succeeded, this will succeed. */

    status2 = ptwXY_domainMin( smr, slopeXY, &slopeXYMin );
    if( ( status2 != nfu_Okay ) && ( status2 != nfu_empty ) ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( status2 );
    }
    ptwXY_domainMax( smr, slopeXY, &slopeXYMax );

    if( ( status1 == nfu_empty ) && ( status2 == nfu_empty ) ) return( nfu_Okay );

    if( ( offsetXYMin != slopeXYMin ) || ( offsetXYMax != slopeXYMax ) ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_domainsNotMutual, "Offset and slope do not have the same domain." );
        return( ptwXY->status = nfu_domainsNotMutual );
    }

    if( ptwXY->length == 0 ) return( nfu_Okay );

    ptwXY_domainMin( smr, ptwXY, &domainMinXY );
    ptwXY_domainMax( smr, ptwXY, &domainMaxXY );

    if( ( domainMinXY >= offsetXYMax ) || ( domainMaxXY <= offsetXYMin ) ) return( nfu_Okay );

    domainMin = ( domainMinXY > offsetXYMin ? domainMinXY : offsetXYMin );
    domainMax = ( domainMaxXY < offsetXYMax ? domainMaxXY : offsetXYMax );

    if( ( domainMinXY == offsetXYMin ) && ( domainMaxXY == offsetXYMax ) ) {
        ptwXY2 = ptwXY;
        offsetXY2 = offsetXY;
        slopeXY2 = slopeXY; }
    else {
        if( ( ptwXY2 = ptwXY_domainSlice( smr, ptwXY, domainMin, domainMax, 0, 1 ) ) == NULL ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            goto Err;
        }

        if( ( offsetXY2 = ptwXY_domainSlice( smr, offsetXY, domainMin, domainMax, 0, 1 ) ) == NULL ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            goto Err;
        }

        if( ( slopeXY2 = ptwXY_domainSlice( smr, slopeXY, domainMin, domainMax, 0, 1 ) ) == NULL ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            goto Err;
        }
    }

    if( ( mulXY = ptwXY_mul2_ptwXY( smr, ptwXY2, slopeXY2 ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        goto Err;
    }
    if( ( addXY = ptwXY_mul2_ptwXY( smr, mulXY, offsetXY2 ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        goto Err;
    }
    if( skipLastPoint != 0 ) addXY->points[addXY->length-1].y = ptwXY2->points[ptwXY2->length-1].y;

    if( domainMin > domainMinXY ) {
        for( i1 = 0, p1 = ptwXY2->points; i1 < ptwXY2->length; i1++, p1++ ) {
            if( p1->x >= domainMin ) break;
            if( ptwXY_setValueAtX( smr, addXY, p1->x, p1->y ) != nfu_Okay ) goto Err;
        }
    }

    if( domainMax < domainMaxXY ) {
        for( i1 = 0, p1 = ptwXY2->points; i1 < ptwXY2->length; i1++, p1++ ) if( p1->x > domainMax ) break;
        for( ; i1 < ptwXY2->length; i1++, p1++ ) {
            if( ptwXY_setValueAtX( smr, addXY, p1->x, p1->y ) != nfu_Okay ) goto Err;
        }
    }

    ptwXY_copy( smr, ptwXY, addXY );

TheEnd:
    if( offsetXY2 != offsetXY ) ptwXY_free( offsetXY2 );
    if( slopeXY2 != slopeXY ) ptwXY_free( slopeXY2 );
    if( ptwXY2 != ptwXY ) ptwXY_free( ptwXY2 );
    ptwXY_free( mulXY );
    ptwXY_free( addXY );

    return( ptwXY->status );

Err:
    ptwXY->status = nfu_Error;
    goto TheEnd;
}
