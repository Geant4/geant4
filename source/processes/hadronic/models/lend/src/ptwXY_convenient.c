/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "ptwXY.h"

static nfu_status ptwXY_createGaussianCenteredSigma1_2( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x1, double y1, 
        double x2, double y2, int addX1Point );
/*
************************************************************
*/
ptwXPoints *ptwXY_getXArray( statusMessageReporting *smr, ptwXYPoints *ptwXY ) {

    int64_t i, n;
    ptwXPoints *xArray;

    if( ptwXY_simpleCoalescePoints( smr, ptwXY ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }
    n = ptwXY->length;

    if( ( xArray = ptwX_new( smr, n ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }
    for( i = 0; i < n; i++ ) xArray->points[i] = ptwXY->points[i].x;
    xArray->length = n;

    return( xArray );
}
/*
************************************************************
*/
ptwXPoints *ptwXY_ysMappedToXs( statusMessageReporting *smr, ptwXYPoints *ptwXY, ptwXPoints *Xs, int64_t *offset ) {

    int64_t iXY, iX, nXY = ptwXY_length( NULL, ptwXY ), nX = ptwX_length( NULL, Xs );
    ptwXYPoint *point1, *point2;
    ptwXY_interpolation interpolation = ptwXY_getInterpolation( ptwXY );
    ptwXPoints *Ys = NULL;

    *offset = 0;

    if( ptwXY_simpleCoalescePoints( smr, ptwXY ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }

    if( Xs->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid destination." );
        return( NULL );
    }

    if( ( nXY == 1 ) || ( nX == 1 ) ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_tooFewPoints, "number of points less than 2: %lld  %lld", nXY, nX );
        return( NULL );
    }

    if( ( nXY == 0 ) || ( nX == 0 ) ) {
        if( ( Ys = ptwX_new( smr, 0 ) ) == NULL ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( NULL );
        }
        return( Ys );
    }

    point1 = &ptwXY->points[0];
    point2 = &ptwXY->points[nXY-1];

    for( iX = 0; iX < nX; ++iX ) {
        if( Xs->points[iX] >= point1->x ) break;
    }
    *offset = iX;

    for( iX = 0; iX < nX; ++iX ) {
        if( Xs->points[iX] > point2->x ) break;
    }
    nX = iX;
    iX = *offset;

    if( ( Ys = ptwX_new( smr, nX - iX ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }
    if( nX - iX < 2 ) return( Ys );

    for( iXY = 1; iXY < nXY; ++iXY ) {
        point2 = &ptwXY->points[iXY];
        if( point2->x >= Xs->points[iX] ) break;
        point1 = point2;
    }

    for( ; iXY < nXY; ++iXY ) {
        point2 = &ptwXY->points[iXY];

        while( iX < nX ) {
            double xValue = Xs->points[iX], yValue;

            if( xValue > point2->x ) break;

            if( ptwXY_interpolatePoint( smr, interpolation, xValue, &yValue, point1->x, point1->y, point2->x, point2->y ) != nfu_Okay ) {
                smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
                ptwX_free( Ys );
                return( NULL );
            }
            if( ptwX_setPointAtIndex( smr, Ys, ptwX_length( NULL, Ys ), yValue ) != nfu_Okay ) {
                smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
                ptwX_free( Ys );
                return( NULL );
            }
            ++iX;
        }
        point1 = point2;
    }

    return( Ys );
}

/*
************************************************************
*/
nfu_status ptwXY_mapToXsAndAdd( statusMessageReporting *a_smr, ptwXYPoints *a_ptwXY, int64_t a_offset, int64_t a_length, double const *a_Xs,
        double *a_results, double a_scaleFractor ) {

    int64_t offset, startIndex, length, length_m1;
    double x1, x2, y1, y2;
    ptwXYPoint *point;
    ptwXY_interpolation interpolation = ptwXY_getInterpolation( a_ptwXY );
    double xValue, yValue;

    if( a_offset < 0 ) a_offset = 0;                                    /* This and next line also ensure that a_length > 0. */
    if( a_offset >= a_length ) return( nfu_Okay );

    offset = a_offset;
    nfu_status status = ptwXY_startIndex( a_smr, a_ptwXY, a_Xs[a_offset], &startIndex, &length );
    if( status != nfu_Okay ) return( status );
    if( startIndex < 0 ) {
        if( startIndex == -2 ) {
            if( a_Xs[a_length-1] >= a_ptwXY->points[0].x ) startIndex = 0;      /* Case A. */
        }
        if( startIndex < 0 ) return( nfu_Okay );
    }

    length_m1 = length - 1;
    point = &a_ptwXY->points[startIndex];
    x1 = point->x;
    y1 = point->y;
    while( startIndex < length_m1 ) {
        ++startIndex;
        point = &a_ptwXY->points[startIndex];
        x2 = point->x;
        y2 = point->y;

        for( ; offset < a_length; ++offset ) {
            xValue = a_Xs[offset];

            if( xValue < x1 ) continue;                     /* Can happend per case A above. */
            if( xValue > x2 ) break;

            if( ( status = ptwXY_interpolatePoint( a_smr, interpolation, xValue, &yValue, x1, y1, x2, y2 ) ) != nfu_Okay ) {
                smr_setReportError2p( a_smr, nfu_SMR_libraryID, nfu_Error, "Via." );
                return( status );
            }
            a_results[offset] += a_scaleFractor * yValue;
        }
        x1 = x2;
        y1 = y2;
    }

    return( nfu_Okay );
}

/*
************************************************************
*/
nfu_status ptwXY_dullEdges( statusMessageReporting *smr, ptwXYPoints *ptwXY, double lowerEps, double upperEps, int positiveXOnly ) {

#define minEps 5e-16

    nfu_status status;
    double xm, xp, dx, y, x1, y1, x2, y2, sign;
    ptwXYPoint *p;

/* This routine can only be used for linear interpolation for the y-axes since for log interpolation, y cannot be 0. 
This needs to be fixed and documented. */
    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( ptwXY->status );
    }

    if( ptwXY->interpolation == ptwXY_interpolationOther ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_otherInterpolation, "Other interpolation not allowed." );
        return( nfu_otherInterpolation );
    }

    if( ptwXY->interpolation == ptwXY_interpolationFlat ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_invalidInterpolation, "Flat interpolation not allowed." );
        return( nfu_invalidInterpolation );
    }

    if( ptwXY->length < 2 ) return( nfu_Okay );

    if( lowerEps != 0. ) {
        if( fabs( lowerEps ) < minEps ) {
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
            dx = fabs( x1 * lowerEps );
            if( x1 == 0 ) dx = fabs( lowerEps );
            xm = x1 - dx;
            xp = x1 + dx;
            if( ( xp + dx ) < x2 ) {
                if( ( status = ptwXY_getValueAtX( smr, ptwXY, xp, &y ) ) != nfu_Okay ) return( ptwXY->status = status );
                if( ( status = ptwXY_setValueAtX( smr, ptwXY, xp,  y ) ) != nfu_Okay ) return( ptwXY->status = status ); }
            else {
                xp = x2;
                y = y2;
            }
            if( lowerEps > 0 ) {
                if( ( status = ptwXY_setValueAtX( smr, ptwXY, x1, 0. ) ) != nfu_Okay ) return( ptwXY->status = status ); }
            else {
                if( ( xm < 0. ) && ( x1 >= 0. ) && positiveXOnly ) {
                    if( ( status = ptwXY_setValueAtX( smr, ptwXY, x1, 0. ) ) != nfu_Okay ) return( ptwXY->status = status ); }
                else {
                    if( ( status = ptwXY_setValueAtX( smr, ptwXY, xm, 0. ) ) != nfu_Okay ) return( ptwXY->status = status );
                    if( ( status = ptwXY_interpolatePoint( smr, ptwXY->interpolation, x1, &y, xm, 0., xp, y )  ) != nfu_Okay )
                        return( ptwXY->status = status );
                    if( ( status = ptwXY_setValueAtX( smr, ptwXY, x1, y ) ) != nfu_Okay ) return( ptwXY->status = status );
                }
            }
        }
    }

    if( upperEps != 0. ) {
        if( fabs( upperEps ) < minEps ) {
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
            dx = fabs( x2 * upperEps );
            if( x2 == 0 ) dx = fabs( upperEps );
            xm = x2 - dx;
            xp = x2 + dx;
            if( ( xm - dx ) > x1 ) {
                if( ( status = ptwXY_getValueAtX( smr, ptwXY, xm, &y ) ) != nfu_Okay ) return( ptwXY->status = status );
                if( ( status = ptwXY_setValueAtX( smr, ptwXY, xm,  y ) ) != nfu_Okay ) return( ptwXY->status = status ); }
            else {
                xm = x1;
                y = y1;
            }
            if( upperEps < 0 ) {
                if( ( status = ptwXY_setValueAtX( smr, ptwXY, x2, 0. ) ) != nfu_Okay ) return( ptwXY->status = status ); }
            else {
                if( ( status = ptwXY_setValueAtX( smr, ptwXY, xp, 0. ) ) != nfu_Okay ) return( ptwXY->status = status );
                if( ( status = ptwXY_interpolatePoint( smr, ptwXY->interpolation, x2, &y, xm, y, xp, 0. )  ) != nfu_Okay )
                    return( ptwXY->status = status );
                if( ( status = ptwXY_setValueAtX( smr, ptwXY, x2, y ) ) != nfu_Okay ) return( ptwXY->status = status );
            }
        }
    }

    return( ptwXY->status );

#undef minEps
}
/*
************************************************************
*/
nfu_status ptwXY_mergeClosePoints( statusMessageReporting *smr, ptwXYPoints *ptwXY, double epsilon ) {

    int64_t i, i1, j, k, n = ptwXY->length;
    double x, y;
    ptwXYPoint *p1, *p2;

    if( n < 2 ) return( ptwXY->status );
    if( epsilon < 4 * DBL_EPSILON ) epsilon = 4 * DBL_EPSILON;
    if( ptwXY_simpleCoalescePoints( smr, ptwXY ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( ptwXY->status );
    }

    p2 = ptwXY->points;
    x = p2->x;
    for( i1 = 1, p2++; i1 < ( n - 1 ); i1++, p2++ ) {                 /* The first point shall remain the first point and all points close to it are deleted. */
        if( ( p2->x - x ) > 0.5 * epsilon * ( fabs( x ) + fabs( p2->x ) ) ) break;
    }
    if( i1 != 1 ) {
        for( i = i1, p1 = &(ptwXY->points[1]); i < n; i++, p1++, p2++ ) *p1 = *p2;
        n = ptwXY->length = ptwXY->length - i1 + 1;
    }

    p1 = &(ptwXY->points[n-1]);
    x = p1->x;
    for( i1 = n - 2, p1--; i1 > 0; i1--, p1-- ) {            /* The last point shall remain the last point and all points close to it are deleted. */
        if( x - p1->x > 0.5 * epsilon * ( fabs( x ) + fabs( p1->x ) ) ) break;
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
            if( ( p2->x - p1->x ) > 0.5 * epsilon * ( fabs( p2->x ) + fabs( p1->x ) ) ) break;
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
ptwXYPoints *ptwXY_intersectionWith_ptwX( statusMessageReporting *smr, ptwXYPoints *ptwXY, ptwXPoints *ptwX ) {

    int64_t i, i1, i2, lengthX = ptwX_length( smr, ptwX );
    double x, y, domainMin, domainMax;
    ptwXYPoints *n = NULL;

    if( ptwX->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid ptwXPoints." );
        return( NULL );
    }

    if( ptwXY_simpleCoalescePoints( smr, ptwXY ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }

    if( ptwXY->interpolation == ptwXY_interpolationOther ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_otherInterpolation, "Other interpolation not allowed." );
        return( NULL );
    }

    if( ( n = ptwXY_clone( smr, ptwXY ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }

    if( ptwXY->length == 0 ) return( n );
    domainMin = ptwXY->points[0].x;
    domainMax = ptwXY->points[ptwXY->length - 1].x;

    if( ( domainMin >= ptwX->points[lengthX-1] ) || ( domainMax <= ptwX->points[0] ) ) {  /* No overlap. */
        n->length = 0;
        return( n );
    }

    for( i = 0; i < lengthX; i++ ) {        /* Fill in ptwXY at x-points in ptwX. */
        x = ptwX->points[i];
        if( x <= domainMin ) continue;
        if( x >= domainMax ) break;
        if( ptwXY_getValueAtX( smr, ptwXY, x, &y ) != nfu_Okay ) goto Err;
        if( ptwXY_setValueAtX( smr, n, x, y ) != nfu_Okay ) goto Err;
    }
    if( ptwXY_simpleCoalescePoints( smr, n ) != nfu_Okay ) goto Err;

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
     smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
     ptwXY_free( n );
     return( NULL );
}
/*
************************************************************
*/
nfu_status ptwXY_areDomainsMutual( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 ) {

    nfu_status status = nfu_Okay;
    int64_t n1 = ptwXY1->length, n2 = ptwXY2->length;
    ptwXYPoint *xy1, *xy2;

    if( ptwXY1->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source1." );
        return( ptwXY1->status );
    }

    if( ptwXY2->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source2." );
        return( ptwXY2->status );
    }

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
nfu_status ptwXY_tweakDomainsToMutualify( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, 
        int epsilonFactor, double epsilon ) {

    nfu_status status = nfu_Okay;
    int64_t n1 = ptwXY1->length, n2 = ptwXY2->length;
    double sum, diff;
    ptwXYPoint *xy1, *xy2;

    epsilon = fabs( epsilon ) + fabs( epsilonFactor * DBL_EPSILON );

    if( ptwXY1->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source1." );
        return( ptwXY1->status );
    }
    if( ptwXY2->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source2." );
        return( ptwXY2->status );
    }

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
                sum = fabs( xy1->x ) + fabs( xy2->x );
                diff = fabs( xy2->x - xy1->x );
                if( diff > epsilon * sum ) {
                    status = nfu_domainsNotMutual; }
                else {
                    xy1->x = xy2->x;
                }
            } }
        else if( xy1->x > xy2->x ) {
            if( xy1->y != 0. ) {
                sum = fabs( xy1->x ) + fabs( xy2->x );
                diff = fabs( xy2->x - xy1->x );
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
                    sum = fabs( xy1->x ) + fabs( xy2->x );
                    diff = fabs( xy2->x - xy1->x );
                    if( diff > epsilon * sum ) {
                        status = nfu_domainsNotMutual; }
                    else {
                        xy2->x = xy1->x;
                    }
                } }
            else if( xy1->x > xy2->x ) {
                if( xy2->y != 0. ) {
                    sum = fabs( xy1->x ) + fabs( xy2->x );
                    diff = fabs( xy2->x - xy1->x );
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
nfu_status ptwXY_mutualifyDomains( statusMessageReporting *smr, ptwXYPoints *ptwXY1, double lowerEps1, double upperEps1, 
        int positiveXOnly1, ptwXYPoints *ptwXY2, double lowerEps2, double upperEps2, int positiveXOnly2 ) {

    nfu_status status;
    int64_t n1 = ptwXY1->length, n2 = ptwXY2->length;
    int code1, code2;
    ptwXYPoint *xy1, *xy2;

    switch( status = ptwXY_areDomainsMutual( smr, ptwXY1, ptwXY2 ) ) {
    case nfu_Okay :
    case nfu_empty :
        return( nfu_Okay );
    case nfu_domainsNotMutual :
        break;
    default :
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( status );
    }

    if( ptwXY1->interpolation == ptwXY_interpolationOther ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_otherInterpolation, "Other interpolation not allowed for source1." );
        return( nfu_otherInterpolation );
    }
    if( ptwXY2->interpolation == ptwXY_interpolationOther ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_otherInterpolation, "Other interpolation not allowed for source2." );
        return( nfu_otherInterpolation );
    }

    if( ptwXY1->interpolation == ptwXY_interpolationFlat ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_invalidInterpolation, "Flat interpolation not allowed for source1." );
        return( nfu_invalidInterpolation );
    }
    if( ptwXY2->interpolation == ptwXY_interpolationFlat ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_invalidInterpolation, "Flat interpolation not allowed for source2." );
        return( nfu_invalidInterpolation );
    }

    xy1 = ptwXY_getPointAtIndex_Unsafely( ptwXY1, 0 );
    xy2 = ptwXY_getPointAtIndex_Unsafely( ptwXY2, 0 );
    code1 = 0;
    if( xy1->x < xy2->x ) {
        lowerEps1 = 0.;
        if( xy2->y == 0. ) {
            lowerEps2 = 0.; }
        else {
            if( lowerEps2 == 0 ) {
                code1 = -1; }
            else {
                if( ( xy2->x - xy1->x ) < lowerEps2 * ( fabs( xy1->x ) + fabs( xy2->x ) ) ) {
                    lowerEps2 = 0;
                    xy1->x = xy2->x;
                    status = ptwXY_areDomainsMutual( smr, ptwXY1, ptwXY2 );
                }
            }
        } }
    else if( xy1->x > xy2->x ) {
        lowerEps2 = 0.;
        if( xy1->y == 0. ) {
            lowerEps1 = 0.; }
        else {
            if( lowerEps1 == 0 ) {
                code1 = 1; }
            else {
                if( ( xy1->x - xy2->x ) < lowerEps1 * ( fabs( xy1->x ) + fabs( xy2->x ) ) ) {
                    lowerEps1 = 0;
                    xy2->x = xy1->x;
                    status = ptwXY_areDomainsMutual( smr, ptwXY1, ptwXY2 );
                }
            }
        } }
    else {
        lowerEps1 = lowerEps2 = 0.;
    }

    xy1 = ptwXY_getPointAtIndex_Unsafely( ptwXY1, n1 - 1 );
    xy2 = ptwXY_getPointAtIndex_Unsafely( ptwXY2, n2 - 1 );
    code2 = 0;
    if( xy1->x < xy2->x ) {
        upperEps2 = 0.;
        if( xy1->y == 0. ) {
            upperEps1 = 0.; }
        else {
            if( upperEps1 == 0 ) {
                code2 = -1; }
            else {
                if( ( xy2->x - xy1->x ) < upperEps1 * ( fabs( xy1->x ) + fabs( xy2->x ) ) ) {
                    upperEps1 = 0;
                    xy2->x = xy1->x;
                    status = ptwXY_areDomainsMutual( smr, ptwXY1, ptwXY2 );
                }
            }
        } }
    else if( xy1->x > xy2->x ) {
        upperEps1 = 0.;
        if( xy2->y == 0. ) {
            upperEps2 = 0.; }
        else {
            if( upperEps2 == 0 ) {
                code2 = 1; }
            else {
                if( ( xy1->x - xy2->x ) < upperEps2 * ( fabs( xy1->x ) + fabs( xy2->x ) ) ) {
                    upperEps2 = 0;
                    xy1->x = xy2->x;
                    status = ptwXY_areDomainsMutual( smr, ptwXY1, ptwXY2 );
                }
            }
        } }
    else {
        upperEps1 = upperEps2 = 0.;
    }

    if( ( lowerEps1 != 0. ) || ( upperEps1 != 0. ) ) {
        if( ( status = ptwXY_dullEdges( smr, ptwXY1, lowerEps1, upperEps1, positiveXOnly1 ) ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( status );
        }
    }
    if( ( lowerEps2 != 0. ) || ( upperEps2 != 0. ) ) {
        if( ( status = ptwXY_dullEdges( smr, ptwXY2, lowerEps2, upperEps2, positiveXOnly2 ) ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( status );
        }
    }

    if( status == nfu_domainsNotMutual ) {
        char str[256] = "";

        if( code1 ==  1 ) strcat( str, " lowerEps1" );
        if( code1 == -1 ) strcat( str, " lowerEps2" );
        if( code2 ==  1 ) strcat( str, " upperEps2" );
        if( code2 == -1 ) strcat( str, " upperEps1" );
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badInput, 
                "The following inputs are 0 and must be a non 0 value: %s.", str );
        status = nfu_badInput;
    }

    return( status );
}
/*
************************************************************
*/
nfu_status ptwXY_copyToC_XY( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t index1, int64_t index2, 
        int64_t allocatedSize, int64_t *numberOfPoints, double *xys ) {

    int64_t i;
    double *d = xys;
    ptwXYPoint *pointFrom;

    if( ptwXY_simpleCoalescePoints( smr, ptwXY ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( ptwXY->status );
    }

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
nfu_status ptwXY_valuesToC_XsAndYs( statusMessageReporting *smr, ptwXYPoints *ptwXY, double **xs, double **ys ) {

    int64_t i1, length;
    double *xps, *yps;
    ptwXYPoint *pointFrom;

    if( ptwXY_simpleCoalescePoints( smr, ptwXY ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( ptwXY->status );
    }
    length = ptwXY_length( NULL, ptwXY );

    if( ( *xs = (double *) smr_malloc2( smr, (size_t) length * sizeof( double ), 0, "xs" ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( nfu_mallocError );
    }
    if( ( *ys = (double *) smr_malloc2( smr, (size_t) length * sizeof( double ), 0, "ys" ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        smr_freeMemory2( *xs );
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
ptwXYPoints *ptwXY_valueTo_ptwXY( statusMessageReporting *smr, double x1, double x2, double y ) {

    ptwXYPoints *n1;

    if( x1 >= x2 ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_XNotAscending,
                "X-values not ascend: x1 = %.17e, x2 = %.17e", x1, x2 );
        return( NULL );
    }
    if( ( n1 = ptwXY_new( smr, ptwXY_interpolationLinLin, NULL, ptwXY_maxBiSectionMax, ptwXY_minAccuracy, 2, 0, 0 ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }
    if( ptwXY_setValueAtX( smr, n1, x1, y ) != nfu_Okay ) goto Err;
    if( ptwXY_setValueAtX( smr, n1, x2, y ) != nfu_Okay ) goto Err;
    return( n1 );

Err:
    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    ptwXY_free( n1 );
    return( NULL );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_createGaussianCenteredSigma1( statusMessageReporting *smr, double accuracy ) {

    int64_t i, n;
    ptwXYPoint *pm, *pp;
    double x1, y1, x2, y2, accuracy2, rangeMin = 1e-10;
    ptwXYPoints *gaussian;

    if( accuracy < 1e-5 ) accuracy = 1e-5;
    if( accuracy > 1e-1 ) accuracy = 1e-1;
    if( ( gaussian = ptwXY_new( smr, ptwXY_interpolationLinLin, NULL, 1., accuracy, 200, 100, 0 ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }

    accuracy2 = accuracy = gaussian->accuracy;
    if( accuracy2 > 5e-3 ) accuracy2 = 5e-3;

    x1 = -sqrt( -2. * log( rangeMin ) );
    y1 = rangeMin;
    x2 = -5.2;
    y2 = exp( -0.5 * x2 * x2 );
    if( ptwXY_setValueAtX( smr, gaussian, x1, y1 ) != nfu_Okay ) goto Err;
    gaussian->accuracy = 20 * accuracy2;
    if( ptwXY_createGaussianCenteredSigma1_2( smr, gaussian, x1, y1, x2, y2, 1 ) != nfu_Okay ) goto Err;
    x1 = x2;
    y1 = y2;
    x2 = -4.;
    y2 = exp( -0.5 * x2 * x2 );
    gaussian->accuracy = 5 * accuracy2;
    if( ptwXY_createGaussianCenteredSigma1_2( smr, gaussian, x1, y1, x2, y2, 1 ) != nfu_Okay ) goto Err;
    x1 = x2;
    y1 = y2;
    x2 = -1;
    y2 = exp( -0.5 * x2 * x2 );
    gaussian->accuracy = accuracy;
    if( ptwXY_createGaussianCenteredSigma1_2( smr, gaussian, x1, y1, x2, y2, 1 ) != nfu_Okay ) goto Err;
    x1 = x2;
    y1 = y2;
    x2 =  0;
    y2 = exp( -0.5 * x2 * x2 );
    if( ptwXY_createGaussianCenteredSigma1_2( smr, gaussian, x1, y1, x2, y2, 1 ) != nfu_Okay ) goto Err;

    n = gaussian->length;
    if( ptwXY_coalescePoints( smr, gaussian, 2 * n + 1, NULL, 0 ) != nfu_Okay ) goto Err;
    if( ptwXY_setValueAtX( smr, gaussian, 0., 1. ) != nfu_Okay ) goto Err;
    pp = &(gaussian->points[gaussian->length]);
    for( i = 0, pm = pp - 2; i < n; i++, pp++, pm-- ) {
        *pp = *pm;
        pp->x *= -1;
    }
    gaussian->length = 2 * n + 1;

    return( gaussian );

Err:
    ptwXY_free( gaussian );
    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( NULL );
}
/*
************************************************************
*/
static nfu_status ptwXY_createGaussianCenteredSigma1_2( statusMessageReporting *smr, ptwXYPoints *ptwXY, 
        double x1, double y1, double x2, double y2, int addX1Point ) {

    nfu_status status = nfu_Okay;
    int morePoints = 0;
    double x = 0.5 * ( x1 + x2 );
    double y = exp( -0.5 * x * x ), rangeMin = ( y1 * ( x2 - x ) + y2 * ( x - x1 ) ) / ( x2 - x1 );

    if( fabs( y - rangeMin ) > y * ptwXY->accuracy ) morePoints = 1;
    if( morePoints && ( status = ptwXY_createGaussianCenteredSigma1_2( smr, ptwXY, x, y, x2, y2, 0 ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_setValueAtX( smr, ptwXY, x, y ) ) != nfu_Okay ) return( status );
    if( morePoints && ( status = ptwXY_createGaussianCenteredSigma1_2( smr, ptwXY, x1, y1, x, y, 0 ) ) != nfu_Okay ) return( status );
    if( addX1Point ) status = ptwXY_setValueAtX( smr, ptwXY, x1, y1 );
    return( status );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_createGaussian( statusMessageReporting *smr, double accuracy, double xCenter, double sigma, 
        double amplitude, double domainMin, double domainMax, double dullEps ) {

    int64_t i;
    ptwXYPoints *gaussian, *sliced;
    ptwXYPoint *point;

    if( ( gaussian = ptwXY_createGaussianCenteredSigma1( smr, accuracy ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }

    for( i = 0, point = gaussian->points; i < gaussian->length; i++, point++ ) {
        point->x = point->x * sigma + xCenter;
        point->y *= amplitude;
    }
    if( ( gaussian->points[0].x < domainMin ) || ( gaussian->points[gaussian->length - 1].x > domainMax ) ) {
        if( ( sliced = ptwXY_domainSlice( smr, gaussian, domainMin, domainMax, 10, 1 ) ) == NULL ) goto Err;
        ptwXY_free( gaussian );
        gaussian = sliced;
    }

    return( gaussian );

Err:
    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    ptwXY_free( gaussian );
    return( NULL );
}
