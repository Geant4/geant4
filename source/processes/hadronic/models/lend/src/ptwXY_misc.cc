/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <float.h>

#include "ptwXY.h"

#if defined __cplusplus
#include <cmath>
#include "G4Log.hh"
namespace GIDI {
using namespace GIDI;
#endif

static nfu_status ptwXY_createFromFunctionBisect( ptwXYPoints *ptwXY, double x1, double y1, double x2, double y2, ptwXY_createFromFunction_callback func,
        void *argList, int level, int checkForRoots, double eps );
static nfu_status ptwXY_createFromFunctionZeroCrossing( ptwXYPoints *ptwXY, double x1, double y1, double x2, double y2, 
        ptwXY_createFromFunction_callback func, void *argList, double eps );
static nfu_status ptwXY_applyFunction2( ptwXYPoints *ptwXY1, double y1, double y2, ptwXYPoint *p1, ptwXYPoint *p2, ptwXY_applyFunction_callback func, 
    void *argList, int level, int checkForRoots );
static nfu_status ptwXY_applyFunctionZeroCrossing( ptwXYPoints *ptwXY1, double y1, double y2, ptwXYPoint *p1, ptwXYPoint *p2, 
    ptwXY_applyFunction_callback func, void *argList );
/*
************************************************************
*/
void ptwXY_update_biSectionMax( ptwXYPoints *ptwXY1, double oldLength ) {

    ptwXY1->biSectionMax = ptwXY1->biSectionMax - 1.442695 * G4Log( ptwXY1->length / oldLength ); /* 1.442695 = 1 / std::log( 2. ) */
    if( ptwXY1->biSectionMax < 0 ) ptwXY1->biSectionMax = 0;
    if( ptwXY1->biSectionMax > ptwXY_maxBiSectionMax ) ptwXY1->biSectionMax = ptwXY_maxBiSectionMax;
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_createFromFunction( int n, double *xs, ptwXY_createFromFunction_callback func, void *argList, double accuracy, int checkForRoots, 
    int biSectionMax, nfu_status *status ) {

    int64_t i;
    double x1, y1, x2, y2, eps = ClosestAllowXFactor * DBL_EPSILON;
    ptwXYPoints *ptwXY;
    ptwXYPoint *p1, *p2;

    *status = nfu_Okay;
    if( n < 2 ) { *status = nfu_tooFewPoints; return( NULL ); }
    for( i = 1; i < n; i++ ) {
        if( xs[i-1] >= xs[i] ) *status = nfu_XNotAscending;
    }
    if( *status == nfu_XNotAscending ) return( NULL );

    x1 = xs[0];
    if( ( *status = func( x1, &y1, argList ) ) != nfu_Okay ) return( NULL );
    if( ( ptwXY = ptwXY_new( ptwXY_interpolationLinLin, NULL, biSectionMax, accuracy, 500, 50, status, 0 ) ) == NULL ) return( NULL );
    for( i = 1; i < n; i++ ) {
        if( ( *status = ptwXY_setValueAtX_overrideIfClose( ptwXY, x1, y1, eps, 0 ) ) != nfu_Okay ) goto err;
        x2 = xs[i];
        if( ( *status = func( x2, &y2, argList ) ) != nfu_Okay ) goto err;
        if( ( *status = ptwXY_createFromFunctionBisect( ptwXY, x1, y1, x2, y2, func, argList, 0, checkForRoots, eps ) ) != nfu_Okay ) goto err;
        x1 = x2;
        y1 = y2;
    }
    if( ( *status = ptwXY_setValueAtX_overrideIfClose( ptwXY, x2, y2, eps, 1 ) ) != nfu_Okay ) goto err;

    if( checkForRoots ) {
        if( ( *status = ptwXY_simpleCoalescePoints( ptwXY ) ) != nfu_Okay ) goto err;
        for( i = ptwXY->length - 1, p2 = NULL; i >= 0; i--, p2 = p1 ) { /* Work backward so lower points are still valid if a new point is added. */
            p1 = &(ptwXY->points[i]);
            if( p2 != NULL ) {
                if( ( p1->y * p2->y ) < 0. ) {
                    if( ( *status = ptwXY_createFromFunctionZeroCrossing( ptwXY, p1->x, p1->y, p2->x, p2->y, func, argList, eps ) ) != nfu_Okay ) goto err;
                }
            }
        }
    }

    return( ptwXY );

err:
    ptwXY_free( ptwXY );
    return( NULL );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_createFromFunction2( ptwXPoints *xs, ptwXY_createFromFunction_callback func, void *argList, double accuracy, int checkForRoots, 
    int biSectionMax, nfu_status *status ) {

    return( ptwXY_createFromFunction( (int) xs->length, xs->points, func, argList, accuracy, checkForRoots, biSectionMax, status ) );
}
/*
************************************************************
*/
static nfu_status ptwXY_createFromFunctionBisect( ptwXYPoints *ptwXY, double x1, double y1, double x2, double y2, ptwXY_createFromFunction_callback func,
        void *argList, int level, int checkForRoots, double eps ) {

    nfu_status status;
    double x, y, f;

    if( ( x2 - x1 ) < ClosestAllowXFactor * DBL_EPSILON * ( std::fabs( x1 ) + std::fabs( x2 ) ) ) return( nfu_Okay );
    if( level >= ptwXY->biSectionMax ) return( nfu_Okay );
    x = 0.5 * ( x1 + x2 );
    if( ( status = ptwXY_interpolatePoint( ptwXY->interpolation, x, &y, x1, y1, x2, y2 ) ) != nfu_Okay ) return( status );
    if( ( status = func( x, &f, argList ) ) != nfu_Okay ) return( status );
    if( std::fabs( f - y ) <= 0.8 * std::fabs( f * ptwXY->accuracy ) ) return( nfu_Okay );
    if( ( status = ptwXY_createFromFunctionBisect( ptwXY, x1, y1, x, f, func, argList, level + 1, checkForRoots, eps ) ) ) return( status );
    if( ( status = ptwXY_setValueAtX_overrideIfClose( ptwXY, x, f, eps, 0 ) ) != nfu_Okay ) return( status );
    return( ptwXY_createFromFunctionBisect( ptwXY, x, f, x2, y2, func, argList, level + 1, checkForRoots, eps ) );
}
/*
************************************************************
*/
static nfu_status ptwXY_createFromFunctionZeroCrossing( ptwXYPoints *ptwXY, double x1, double y1, double x2, double y2, 
        ptwXY_createFromFunction_callback func, void *argList, double eps ) {

    //For coverity #63077
    if ( y2 == y1 ) return ( nfu_badInput );

    int i;
    double x, y;
    nfu_status status;

    for( i = 0; i < 16; i++ ) {
        if( y2 == y1 ) break;
        x = ( y2 * x1 - y1 * x2 ) / ( y2 - y1 );
        if( x <= x1 ) x = x1 + 0.1 * ( x2 - x1 );
        if( x >= x2 ) x =  x2 - 0.1 * ( x2 - x1 );
        if( ( status = func( x, &y, argList ) ) != nfu_Okay ) return( status );
        if( y == 0 ) break;
        if( y1 * y < 0 ) {
            x2 = x;
            y2 = y; }
        else {
            x1 = x;
            y1 = y;
        }
    }
    return( ptwXY_setValueAtX_overrideIfClose( ptwXY, x, 0., eps, 1 ) );
}
/*
************************************************************
*/
nfu_status ptwXY_applyFunction( ptwXYPoints *ptwXY1, ptwXY_applyFunction_callback func, void *argList, int checkForRoots ) {

    int64_t i, originalLength = ptwXY1->length, notFirstPass = 0;
    double y1, y2 = 0;
    nfu_status status;
    ptwXYPoint p1, p2;

    checkForRoots = checkForRoots && ptwXY1->biSectionMax;
    if( ptwXY1->status != nfu_Okay ) return( ptwXY1->status );
    if( ptwXY1->interpolation == ptwXY_interpolationOther ) return( nfu_otherInterpolation );
    if( ptwXY1->interpolation == ptwXY_interpolationFlat ) return( nfu_invalidInterpolation );
    if( ( status = ptwXY_simpleCoalescePoints( ptwXY1 ) ) != nfu_Okay ) return( status );
    for( i = originalLength - 1; i >= 0; i-- ) {
        y1 = ptwXY1->points[i].y;
        if( ( status = func( &(ptwXY1->points[i]), argList ) ) != nfu_Okay ) return( status );
        p1 = ptwXY1->points[i];
        if( notFirstPass ) {
            if( ( status = ptwXY_applyFunction2( ptwXY1, y1, y2, &p1, &p2, func, argList, 0, checkForRoots ) ) != nfu_Okay ) return( status );
        }
        notFirstPass = 1;
        p2 = p1;
        y2 = y1;
    }
    ptwXY_update_biSectionMax( ptwXY1, (double) originalLength );
    return( status );
}
/*
************************************************************
*/
static nfu_status ptwXY_applyFunction2( ptwXYPoints *ptwXY1, double y1, double y2, ptwXYPoint *p1, ptwXYPoint *p2, ptwXY_applyFunction_callback func, 
        void *argList, int level, int checkForRoots ) {

    double y;
    ptwXYPoint p;
    nfu_status status;

    if( ( p2->x - p1->x ) < ClosestAllowXFactor * DBL_EPSILON * ( std::fabs( p1->x ) + std::fabs( p2->x ) ) ) return( nfu_Okay );
    if( level >= ptwXY1->biSectionMax ) goto checkForZeroCrossing;
    p.x = 0.5 * ( p1->x + p2->x );
    if( ( status = ptwXY_interpolatePoint( ptwXY1->interpolation, p.x, &y, p1->x, y1, p2->x, y2 ) ) != nfu_Okay ) return( status );
    p.y = y;
    if( ( status = func( &p, argList ) ) != nfu_Okay ) return( status );
    if( std::fabs( ( p.x - p1->x ) * ( p2->y - p1->y ) + ( p2->x - p1->x ) * ( p1->y - p.y ) ) <= 0.8 * std::fabs( ( p2->x - p1->x ) * p.y * ptwXY1->accuracy ) ) 
        goto checkForZeroCrossing;
    if( ( status = ptwXY_setValueAtX( ptwXY1, p.x, p.y ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_applyFunction2( ptwXY1, y1, y, p1, &p, func, argList, level + 1, checkForRoots ) ) ) return( status );
    return( ptwXY_applyFunction2( ptwXY1, y, y2, &p, p2, func, argList, level + 1, checkForRoots ) );

checkForZeroCrossing:
    if( checkForRoots && ( ( p1->y * p2->y ) < 0. ) ) return( ptwXY_applyFunctionZeroCrossing( ptwXY1, y1, y2, p1, p2, func, argList ) );
    return( nfu_Okay );
}
/*
************************************************************
*/
static nfu_status ptwXY_applyFunctionZeroCrossing( ptwXYPoints *ptwXY1, double y1, double y2, ptwXYPoint *p1, ptwXYPoint *p2, 
        ptwXY_applyFunction_callback func, void *argList ) {

    int i;
    double y, x1 = p1->x, x2 = p2->x, nY1 = p1->y, nY2 = p2->y, refY = 0.5 * ( std::fabs( p1->y ) + std::fabs( p2->y ) );
    ptwXYPoint p;
    nfu_status status;

    //For coverity #63074
    if ( nY2 == nY1 ) return ( nfu_badInput );

    for( i = 0; i < 6; i++ ) {
        if( nY2 == nY1 ) break;
        p.x = ( nY2 * x1 - nY1 * x2 ) / ( nY2 - nY1 );
        if( p.x <= x1 ) p.x = 0.5 * ( x1 + x2 );
        if( p.x >= x2 ) p.x = 0.5 * ( x1 + x2 );
        if( ( status = ptwXY_interpolatePoint( ptwXY1->interpolation, p.x, &y, p1->x, y1, p2->x, y2 ) ) != nfu_Okay ) return( status );
        p.y = y;
        if( ( status = func( &p, argList ) ) != nfu_Okay ) return( status );
        if( p.y == 0 ) break;
        if( 0.5 * refY < std::fabs( p.y ) ) break;
        refY = std::fabs( p.y );
        if( p1->y * p.y < 0 ) {
            x2 = p.x;
            nY2 = p.y; }
        else {
            x1 = p.x;
            nY1 = p.y;
        }
    }
    return( ptwXY_setValueAtX( ptwXY1, p.x, 0. ) );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_fromString( char const *str, ptwXY_interpolation interpolation, ptwXY_interpolationOtherInfo const *interpolationOtherInfo,
    double biSectionMax, double accuracy, char **endCharacter, nfu_status *status ) {

    int64_t numberConverted;
    double  *doublePtr;
    ptwXYPoints *ptwXY = NULL;

    if( ( *status = nfu_stringToListOfDoubles( str, &numberConverted, &doublePtr, endCharacter ) ) != nfu_Okay ) return( NULL );
    *status = nfu_oddNumberOfValues;
    if( ( numberConverted % 2 ) == 0 )
        ptwXY = ptwXY_create( interpolation, interpolationOtherInfo, biSectionMax, accuracy, numberConverted, 10, numberConverted / 2, doublePtr, status, 0 );
    nfu_free( doublePtr );
    return( ptwXY );
}
/*
************************************************************
*/
void ptwXY_showInteralStructure( ptwXYPoints *ptwXY, FILE *f, int printPointersAsNull ) {

    int64_t i, n = ptwXY_getNonOverflowLength( ptwXY );
    ptwXYPoint *point = ptwXY->points;
    ptwXYOverflowPoint *overflowPoint;

    fprintf( f, "status = %d  interpolation = %d  length = %d  allocatedSize = %d\n", 
        (int) ptwXY->status, (int) ptwXY->interpolation, (int) ptwXY->length, (int) ptwXY->allocatedSize );
    fprintf( f, "userFlag = %d  biSectionMax = %.8e  accuracy = %.2e  minFractional_dx = %.6e\n", 
        ptwXY->userFlag, ptwXY->biSectionMax, ptwXY->accuracy, ptwXY->minFractional_dx );
    fprintf( f, "interpolationString = %s\n", ptwXY->interpolationOtherInfo.interpolationString );
    fprintf( f, "getValueFunc is NULL = %d. argList is NULL = %d.\n", 
        ptwXY->interpolationOtherInfo.getValueFunc == NULL, ptwXY->interpolationOtherInfo.argList == NULL );
    fprintf( f, "  overflowLength = %d  overflowAllocatedSize = %d  mallocFailedSize = %d\n", 
        (int) ptwXY->overflowLength, (int) ptwXY->overflowAllocatedSize, (int) ptwXY->mallocFailedSize );
    fprintf( f, "  Points data, points = %20p\n", ( printPointersAsNull ? NULL : (void*)ptwXY->points ) );
    for( i = 0; i < n; i++,  point++ ) fprintf( f, "    %14.7e %14.7e\n", point->x, point->y );
    fprintf( f, "  Overflow points data; %20p\n", ( printPointersAsNull ? NULL : (void*)&(ptwXY->overflowHeader) ) );
    for( overflowPoint = ptwXY->overflowHeader.next; overflowPoint != &(ptwXY->overflowHeader); overflowPoint = overflowPoint->next ) {
        fprintf( f, "    %14.7e %14.7e %8d %20p %20p %20p\n", overflowPoint->point.x, overflowPoint->point.y, (int) overflowPoint->index, 
    (void*) ( printPointersAsNull ? NULL : overflowPoint ), (void*) ( printPointersAsNull ? NULL : overflowPoint->prior ), 
    (void*) ( printPointersAsNull ? NULL : overflowPoint->next ) );
    }
    fprintf( f, "  Points in order\n" );
    for( i = 0; i < ptwXY->length; i++ ) {
        point = ptwXY_getPointAtIndex( ptwXY, i );
        fprintf( f, "    %14.7e %14.7e\n", point->x, point->y );
    }
}
/*
************************************************************
*/
void ptwXY_simpleWrite( ptwXYPoints *ptwXY, FILE *f, char *format ) {

    int64_t i;
    ptwXYPoint *point;

    for( i = 0; i < ptwXY->length; i++ ) {
        point = ptwXY_getPointAtIndex( ptwXY, i );
        fprintf( f, format, point->x, point->y );
    }
}
/*
************************************************************
*/
void ptwXY_simplePrint( ptwXYPoints *ptwXY, char *format ) {

    ptwXY_simpleWrite( ptwXY, stdout, format );
}

#if defined __cplusplus
}
#endif
