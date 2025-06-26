/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "ptwXY.h"

static nfu_status ptwXY_createFromFunctionBisect( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x1, double y1, 
        double x2, double y2, ptwXY_createFromFunction_callback func, void *argList, int level, int checkForRoots, double eps );
static nfu_status ptwXY_createFromFunctionZeroCrossing( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x1, double y1, 
        double x2, double y2, ptwXY_createFromFunction_callback func, void *argList, double eps );
static nfu_status ptwXY_applyFunction2( statusMessageReporting *smr, ptwXYPoints *ptwXY1, double y1, double y2, 
        ptwXYPoint *p1, ptwXYPoint *p2, ptwXY_applyFunction_callback func, void *argList, int level, int checkForRoots );
static nfu_status ptwXY_applyFunctionZeroCrossing( statusMessageReporting *smr, ptwXYPoints *ptwXY1, double y1, double y2, 
        ptwXYPoint *p1, ptwXYPoint *p2, ptwXY_applyFunction_callback func, void *argList );
/*
************************************************************
*/
double ptwXY_limitAccuracy( double accuracy ) {

    if( accuracy < ptwXY_minAccuracy ) accuracy = ptwXY_minAccuracy;
    if( accuracy > 1 ) accuracy = 1.;
    return( accuracy );
}
/*
************************************************************
*/
void ptwXY_update_biSectionMax( ptwXYPoints *ptwXY1, double oldLength ) {

    ptwXY1->biSectionMax = ptwXY1->biSectionMax - 1.442695 * log( ptwXY1->length / oldLength ); /* 1.442695 = 1 / log( 2. ) */
    if( ptwXY1->biSectionMax < 0 ) ptwXY1->biSectionMax = 0;
    if( ptwXY1->biSectionMax > ptwXY_maxBiSectionMax ) ptwXY1->biSectionMax = ptwXY_maxBiSectionMax;
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_createFromFunction( statusMessageReporting *smr, int n, double *xs, ptwXY_createFromFunction_callback func, 
        void *argList, double accuracy, int checkForRoots, int biSectionMax ) {

    int64_t i;
    double x1, y1, x2, y2, eps = ClosestAllowXFactor * DBL_EPSILON;
    ptwXYPoints *ptwXY;
    ptwXYPoint *p1, *p2;

    if( n < 2 ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_tooFewPoints, "Too few point = %d.", (int) n );
        return( NULL );
    }
    for( i = 1; i < n; i++ ) {
        if( xs[i-1] >= xs[i] ) {
            smr_setReportError2( smr, nfu_SMR_libraryID, nfu_XNotAscending, 
                    "Non-ascending domain values: xs[%d] = %.17e >= xs[%d] = %.17e.", 
                    (int) (i-1), xs[i-1], (int) i, xs[i] );
            return( NULL );
        }
    }

    x1 = xs[0];
    if( func( smr, x1, &y1, argList ) != nfu_Okay ) {
        return( NULL );
    }
    if( ( ptwXY = ptwXY_new( smr, ptwXY_interpolationLinLin, NULL, biSectionMax, accuracy, 500, 50, 0 ) ) == NULL ) goto Err;
    for( i = 1; i < n; i++ ) {
        if( ptwXY_setValueAtX_overrideIfClose( smr, ptwXY, x1, y1, eps, 0 ) != nfu_Okay ) goto Err;
        x2 = xs[i];
        if( func( smr, x2, &y2, argList ) != nfu_Okay ) goto Err;
        if( ptwXY_createFromFunctionBisect( smr, ptwXY, x1, y1, x2, y2, func, argList, 0, checkForRoots, eps ) != nfu_Okay ) goto Err;
        x1 = x2;
        y1 = y2;
    }
    if( ptwXY_setValueAtX_overrideIfClose( smr, ptwXY, x2, y2, eps, 1 ) != nfu_Okay ) goto Err;

    if( checkForRoots ) {
        if( ptwXY_simpleCoalescePoints( NULL, ptwXY ) != nfu_Okay ) goto Err;
        for( i = ptwXY->length - 1, p2 = NULL; i >= 0; i--, p2 = p1 ) { /* Work backward so lower points are still valid if a new point is added. */
            p1 = &(ptwXY->points[i]);
            if( p2 != NULL ) {
                if( ( p1->y * p2->y ) < 0. ) {
                    if( ptwXY_createFromFunctionZeroCrossing( smr, ptwXY, p1->x, p1->y, p2->x, p2->y, func, argList, eps ) != nfu_Okay ) goto Err;
                }
            }
        }
    }

    return( ptwXY );

Err:
    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    if( ptwXY != NULL ) ptwXY_free( ptwXY );
    return( NULL );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_createFromFunction2( statusMessageReporting *smr, ptwXPoints *xs, ptwXY_createFromFunction_callback func, 
        void *argList, double accuracy, int checkForRoots, int biSectionMax ) {

    ptwXYPoints *ptwXY = ptwXY_createFromFunction( smr, (int) xs->length, xs->points, func, argList, accuracy, 
            checkForRoots, biSectionMax );

    if( ptwXY == NULL ) smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( ptwXY );
}
/*
************************************************************
*/
static nfu_status ptwXY_createFromFunctionBisect( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x1, double y1, 
        double x2, double y2, ptwXY_createFromFunction_callback func, void *argList, int level, int checkForRoots, double eps ) {

    nfu_status status;
    double x, y, f;

    if( ( x2 - x1 ) < ClosestAllowXFactor * DBL_EPSILON * ( fabs( x1 ) + fabs( x2 ) ) ) return( nfu_Okay );
    if( level >= ptwXY->biSectionMax ) return( nfu_Okay );
    x = 0.5 * ( x1 + x2 );
    if( ( status = ptwXY_interpolatePoint( smr, ptwXY->interpolation, x, &y, x1, y1, x2, y2 ) ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( status );
    }
    if( ( status = func( smr, x, &f, argList ) ) != nfu_Okay ) return( status );
    if( fabs( f - y ) <= 0.8 * fabs( f * ptwXY->accuracy ) ) return( nfu_Okay );
    if( ptwXY_createFromFunctionBisect( smr, ptwXY, x1, y1, x, f, func, argList, level + 1, checkForRoots, eps ) ) return( status );
    if( ptwXY_setValueAtX_overrideIfClose( smr, ptwXY, x, f, eps, 0 ) != nfu_Okay ) return( status );
    return( ptwXY_createFromFunctionBisect( smr, ptwXY, x, f, x2, y2, func, argList, level + 1, checkForRoots, eps ) );
}
/*
************************************************************
*/
static nfu_status ptwXY_createFromFunctionZeroCrossing( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x1, double y1, 
        double x2, double y2, ptwXY_createFromFunction_callback func, void *argList, double eps ) {

    int i;
    double x = 0, y;        /* Initialize x so some compilers do not complain. */
    nfu_status status;

    for( i = 0; i < 16; i++ ) {
        if( y2 == y1 ) break;
        x = ( y2 * x1 - y1 * x2 ) / ( y2 - y1 );
        if( x <= x1 ) x = x1 + 0.1 * ( x2 - x1 );
        if( x >= x2 ) x =  x2 - 0.1 * ( x2 - x1 );
        if( ( status = func( smr, x, &y, argList ) ) != nfu_Okay ) return( status );
        if( y == 0 ) break;
        if( y1 * y < 0 ) {
            x2 = x;
            y2 = y; }
        else {
            x1 = x;
            y1 = y;
        }
    }
    return( ptwXY_setValueAtX_overrideIfClose( smr, ptwXY, x, 0., eps, 1 ) );
}
/*
************************************************************
*/
nfu_status ptwXY_applyFunction( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXY_applyFunction_callback func, 
        void *argList, int checkForRoots ) {

    int64_t i, originalLength = ptwXY1->length, notFirstPass = 0;
    double y1, y2 = 0;
    nfu_status status;
    ptwXYPoint p1, p2;

    checkForRoots = checkForRoots && ptwXY1->biSectionMax;

    if( ptwXY1->interpolation == ptwXY_interpolationOther ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_otherInterpolation, "Other interpolation not allowed." );
        return( ptwXY1->status = nfu_otherInterpolation );
    }
    if( ptwXY1->interpolation == ptwXY_interpolationFlat ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_flatInterpolation, "Flat interpolation not allowed." );
        return( ptwXY1->status = nfu_flatInterpolation );
    }

    if( ptwXY_simpleCoalescePoints( smr, ptwXY1 ) != nfu_Okay ) goto Err;
    for( i = originalLength - 1; i >= 0; i-- ) {
        y1 = ptwXY1->points[i].y;
        if( ( status = func( smr, &(ptwXY1->points[i]), argList ) ) != nfu_Okay ) {
            if( ptwXY1->status == nfu_Okay ) ptwXY1->status = status;
            return( status );
        }
        p1 = ptwXY1->points[i];
        if( notFirstPass ) {
            if( ptwXY_applyFunction2( smr, ptwXY1, y1, y2, &p1, &p2, func, argList, 0, checkForRoots ) != nfu_Okay ) goto Err;
        }
        notFirstPass = 1;
        p2 = p1;
        y2 = y1;
    }
    ptwXY_update_biSectionMax( ptwXY1, (double) originalLength );
    return( status );

Err:
    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    if( ptwXY1->status == nfu_Okay ) ptwXY1->status = nfu_Error;
    return( nfu_Error );
}
/*
************************************************************
*/
static nfu_status ptwXY_applyFunction2( statusMessageReporting *smr, ptwXYPoints *ptwXY1, double y1, double y2, 
        ptwXYPoint *p1, ptwXYPoint *p2, ptwXY_applyFunction_callback func, void *argList, int level, int checkForRoots ) {

    double y;
    ptwXYPoint p;
    nfu_status status;

    if( ( p2->x - p1->x ) < ClosestAllowXFactor * DBL_EPSILON * ( fabs( p1->x ) + fabs( p2->x ) ) ) return( nfu_Okay );
    if( level >= ptwXY1->biSectionMax ) goto checkForZeroCrossing;
    p.x = 0.5 * ( p1->x + p2->x );
    if( ( status = ptwXY_interpolatePoint( smr, ptwXY1->interpolation, p.x, &y, p1->x, y1, p2->x, y2 ) ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( status );
    }
    p.y = y;
    if( ( status = func( smr, &p, argList ) ) != nfu_Okay ) return( status );
    if( fabs( ( p.x - p1->x ) * ( p2->y - p1->y ) + ( p2->x - p1->x ) * ( p1->y - p.y ) ) <= 0.8 * fabs( ( p2->x - p1->x ) * p.y * ptwXY1->accuracy ) ) 
        goto checkForZeroCrossing;
    if( ( status = ptwXY_setValueAtX( smr, ptwXY1, p.x, p.y ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_applyFunction2( smr, ptwXY1, y1, y, p1, &p, func, argList, level + 1, checkForRoots ) ) ) return( status );
    return( ptwXY_applyFunction2( smr, ptwXY1, y, y2, &p, p2, func, argList, level + 1, checkForRoots ) );

checkForZeroCrossing:
    if( checkForRoots && ( ( p1->y * p2->y ) < 0. ) )
        return( ptwXY_applyFunctionZeroCrossing( smr, ptwXY1, y1, y2, p1, p2, func, argList ) );
    return( nfu_Okay );
}
/*
************************************************************
*/
static nfu_status ptwXY_applyFunctionZeroCrossing( statusMessageReporting *smr, ptwXYPoints *ptwXY1, double y1, double y2, 
        ptwXYPoint *p1, ptwXYPoint *p2, ptwXY_applyFunction_callback func, void *argList ) {

    int i;
    double y, x1 = p1->x, x2 = p2->x, nY1 = p1->y, nY2 = p2->y, refY = 0.5 * ( fabs( p1->y ) + fabs( p2->y ) );
    ptwXYPoint p = { 0.5 * ( p1->x + p2->x ), 0.0 };
    nfu_status status;

    for( i = 0; i < 6; i++ ) {
        if( nY2 == nY1 ) break;
        p.x = ( nY2 * x1 - nY1 * x2 ) / ( nY2 - nY1 );
        if( p.x <= x1 ) p.x = 0.5 * ( x1 + x2 );
        if( p.x >= x2 ) p.x = 0.5 * ( x1 + x2 );
        if( ( status = ptwXY_interpolatePoint( smr, ptwXY1->interpolation, p.x, &y, p1->x, y1, p2->x, y2 ) ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( status );
        }
        p.y = y;
        if( ( status = func( smr, &p, argList ) ) != nfu_Okay ) return( status );
        if( p.y == 0 ) break;
        if( 0.5 * refY < fabs( p.y ) ) break;
        refY = fabs( p.y );
        if( p1->y * p.y < 0 ) {
            x2 = p.x;
            nY2 = p.y; }
        else {
            x1 = p.x;
            nY1 = p.y;
        }
    }
    return( ptwXY_setValueAtX( smr, ptwXY1, p.x, 0. ) );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_fromString( statusMessageReporting *smr, char const *str, char sep, ptwXY_interpolation interpolation, 
        char const *interpolationString, double biSectionMax, double accuracy, char **endCharacter, int useSystem_strtod ) {

    int64_t numberConverted;
    double  *doublePtr;
    ptwXYPoints *ptwXY = NULL;

    if( ( doublePtr = nfu_stringToListOfDoubles( smr, str, sep, &numberConverted, endCharacter, useSystem_strtod ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }
    if( ( numberConverted % 2 ) == 0 ) {
        ptwXY = ptwXY_create( smr, interpolation, interpolationString, biSectionMax, accuracy, numberConverted / 2, 10, numberConverted / 2, doublePtr, 0 ); }
    else {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_oddNumberOfValues, "Odd number = %d of float for ptwXY.", (int) numberConverted );
    }
    smr_freeMemory2( doublePtr );
    return( ptwXY );
}
/*
************************************************************
*/
void ptwXY_showInteralStructure( ptwXYPoints *ptwXY, FILE *f, int printPointersAsNull ) {

    int64_t i, n1;
    ptwXYPoint *point = ptwXY->points;
    ptwXYOverflowPoint *overflowPoint;

    n1 = ptwXY_getNonOverflowLength( NULL, ptwXY );

    fprintf( f, "status = %d  interpolation = %d  length = %d  allocatedSize = %d\n", 
        (int) ptwXY->status, (int) ptwXY->interpolation, (int) ptwXY->length, (int) ptwXY->allocatedSize );
    fprintf( f, "userFlag = %d  biSectionMax = %.8e  accuracy = %.2e  minFractional_dx = %.6e\n", 
        ptwXY->userFlag, ptwXY->biSectionMax, ptwXY->accuracy, ptwXY->minFractional_dx );
    fprintf( f, "interpolationString = %s\n", ptwXY->interpolationString );
    fprintf( f, "  overflowLength = %d  overflowAllocatedSize = %d  mallocFailedSize = %d\n", 
        (int) ptwXY->overflowLength, (int) ptwXY->overflowAllocatedSize, (int) ptwXY->mallocFailedSize );
    fprintf( f, "  Points data, points = %20p\n", ( printPointersAsNull ? NULL : ptwXY->points ) );
    for( i = 0; i < n1; i++,  point++ ) fprintf( f, "    %14.7e %14.7e\n", point->x, point->y );
    fprintf( f, "  Overflow points data; %20p\n", ( printPointersAsNull ? NULL : &(ptwXY->overflowHeader) ) );
    for( overflowPoint = ptwXY->overflowHeader.next; overflowPoint != &(ptwXY->overflowHeader); overflowPoint = overflowPoint->next ) {
        fprintf( f, "    %14.7e %14.7e %8d %20p %20p %20p\n", overflowPoint->point.x, overflowPoint->point.y, (int) overflowPoint->index, 
            ( printPointersAsNull ? NULL : overflowPoint ), ( printPointersAsNull ? NULL : overflowPoint->prior ), 
            ( printPointersAsNull ? NULL : overflowPoint->next ) );
    }
    fprintf( f, "  Points in order\n" );
    for( i = 0; i < ptwXY->length; i++ ) {
        point = ptwXY_getPointAtIndex_Unsafely( ptwXY, i );
        fprintf( f, "    %14.7e %14.7e\n", point->x, point->y );
    }
}
/*
************************************************************
*/
void ptwXY_simpleWrite( ptwXYPoints *ptwXY, FILE *f, char const *format ) {

    int64_t i;
    ptwXYPoint *point;

    for( i = 0; i < ptwXY->length; i++ ) {
        point = ptwXY_getPointAtIndex_Unsafely( ptwXY, i );
        fprintf( f, format, point->x, point->y );
    }
}
/*
************************************************************
*/
void ptwXY_simplePrint( ptwXYPoints *ptwXY, char const *format ) {

    ptwXY_simpleWrite( ptwXY, stdout, format );
}
