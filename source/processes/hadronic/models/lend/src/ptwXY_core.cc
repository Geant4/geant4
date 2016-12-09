/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "ptwXY.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

static char const linLinInterpolationString[] = "linear,linear";
static char const linLogInterpolationString[] = "linear,log";
static char const logLinInterpolationString[] = "log,linear";
static char const logLogInterpolationString[] = "log,log";
static char const flatInterpolationString[] = "flat";

static void ptwXY_initialOverflowPoint( ptwXYOverflowPoint *overflowPoint, ptwXYOverflowPoint *prior, ptwXYOverflowPoint *next );
static nfu_status ptwXY_mergeFrom( ptwXYPoints *ptwXY, int incY, int length, double *xs, double *ys );
static int ptwXY_mergeCompareFunction( void const *x1p, void const *x2p );
/*
************************************************************
*/
ptwXYPoints *ptwXY_new( ptwXY_interpolation interpolation, ptwXY_interpolationOtherInfo const *interpolationOtherInfo, double biSectionMax, 
    double accuracy, int64_t primarySize, int64_t secondarySize, nfu_status *status, int userFlag ) {

    ptwXYPoints *ptwXY = (ptwXYPoints *) nfu_calloc( sizeof( ptwXYPoints ), 1 );

    *status = nfu_mallocError;
    if( ptwXY == NULL ) return( NULL );
    ptwXY_setup( ptwXY, interpolation, interpolationOtherInfo, biSectionMax, accuracy, primarySize, 
        secondarySize, userFlag );
    if( ( *status = ptwXY->status ) != nfu_Okay ) {
        ptwXY = (ptwXYPoints *) nfu_free( ptwXY );
    }
    return( ptwXY );
}
/*
************************************************************
*/
nfu_status ptwXY_setup( ptwXYPoints *ptwXY, ptwXY_interpolation interpolation, ptwXY_interpolationOtherInfo const *interpolationOtherInfo,
    double biSectionMax, double accuracy, int64_t primarySize, int64_t secondarySize, int userFlag ) {

    ptwXY->status = nfu_Okay;
    ptwXY->typeX =  ptwXY_sigma_none;
    ptwXY->typeY =  ptwXY_sigma_none;
    ptwXY->interpolation = interpolation;
    ptwXY->interpolationOtherInfo.interpolationString = NULL;
    ptwXY->interpolationOtherInfo.getValueFunc = NULL;
    ptwXY->interpolationOtherInfo.argList = NULL;
    switch( interpolation ) {
    case ptwXY_interpolationLinLin :
        ptwXY->interpolationOtherInfo.interpolationString = linLinInterpolationString; break;
    case ptwXY_interpolationLinLog :
        ptwXY->interpolationOtherInfo.interpolationString = linLogInterpolationString; break;
    case ptwXY_interpolationLogLin :
        ptwXY->interpolationOtherInfo.interpolationString = logLinInterpolationString; break;
    case ptwXY_interpolationLogLog :
        ptwXY->interpolationOtherInfo.interpolationString = logLogInterpolationString; break;
    case ptwXY_interpolationFlat :
        ptwXY->interpolationOtherInfo.interpolationString = flatInterpolationString; break;
    case ptwXY_interpolationOther :         /* For ptwXY_interpolationOther, interpolationOtherInfo and interpolationString must be defined. */
        if( interpolationOtherInfo == NULL ) {              
                ptwXY->status = nfu_otherInterpolation; }
        else {
            if( interpolationOtherInfo->interpolationString == NULL ) {
                ptwXY->status = nfu_otherInterpolation; }
            else {
                if( ( ptwXY->interpolationOtherInfo.interpolationString = strdup( interpolationOtherInfo->interpolationString ) ) == NULL ) {
                    ptwXY->status = nfu_mallocError;
                }
            }
            ptwXY->interpolationOtherInfo.getValueFunc = interpolationOtherInfo->getValueFunc;
            ptwXY->interpolationOtherInfo.argList = interpolationOtherInfo->argList;
        }
    }
    ptwXY->userFlag = 0;
    ptwXY_setUserFlag( ptwXY, userFlag );
    ptwXY->biSectionMax = ptwXY_maxBiSectionMax;
    ptwXY_setBiSectionMax( ptwXY, biSectionMax );
    ptwXY->accuracy = ptwXY_minAccuracy;
    ptwXY_setAccuracy( ptwXY, accuracy );

    ptwXY->length = 0;
    ptwXY->allocatedSize = 0;
    ptwXY->overflowLength = 0;
    ptwXY->overflowAllocatedSize = 0;
    ptwXY->mallocFailedSize = 0;

    ptwXY_initialOverflowPoint( &(ptwXY->overflowHeader), &(ptwXY->overflowHeader), &(ptwXY->overflowHeader) );

    ptwXY->points = NULL;
    ptwXY->overflowPoints = NULL;

    ptwXY_reallocatePoints( ptwXY, primarySize, 0 );
    ptwXY_reallocateOverflowPoints( ptwXY, secondarySize );
    if( ptwXY->status != nfu_Okay ) ptwXY_release( ptwXY );
    return( ptwXY->status );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_create( ptwXY_interpolation interpolation, ptwXY_interpolationOtherInfo const *interpolationOtherInfo,
    double biSectionMax, double accuracy, int64_t primarySize, int64_t secondarySize, int64_t length, double const *xy, 
    nfu_status *status, int userFlag ) {

    ptwXYPoints *ptwXY;

    if( primarySize < length ) primarySize = length;
    if( ( ptwXY = ptwXY_new( interpolation, interpolationOtherInfo, biSectionMax, accuracy, primarySize, 
            secondarySize, status, userFlag ) ) != NULL ) {
        if( ( *status = ptwXY_setXYData( ptwXY, length, xy ) ) != nfu_Okay ) {
            ptwXY = ptwXY_free( ptwXY );
        }
    }
    return( ptwXY );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_createFrom_Xs_Ys( ptwXY_interpolation interpolation, ptwXY_interpolationOtherInfo const *interpolationOtherInfo,
    double biSectionMax, double accuracy, int64_t primarySize, int64_t secondarySize, int64_t length, double const *Xs, 
    double const *Ys, nfu_status *status, int userFlag ) {

    int i;
    ptwXYPoints *ptwXY;

    if( primarySize < length ) primarySize = length;
    if( ( ptwXY = ptwXY_new( interpolation, interpolationOtherInfo, biSectionMax, accuracy, primarySize, 
            secondarySize, status, userFlag ) ) != NULL ) {
        for( i = 0; i < length; i++ ) {
            ptwXY->points[i].x = Xs[i];
            ptwXY->points[i].y = Ys[i];
        }
        ptwXY->length = length;
    }

    return( ptwXY );
}
/*
************************************************************
*/
nfu_status ptwXY_copy( ptwXYPoints *dest, ptwXYPoints *src ) {

    int64_t i, nonOverflowLength = ptwXY_getNonOverflowLength( src );
    ptwXYPoint *pointFrom, *pointTo;
    ptwXYOverflowPoint *o, *overflowHeader = &(src->overflowHeader);

    if( dest->status != nfu_Okay ) return( dest->status );
    if( src->status != nfu_Okay ) return( src->status );

    ptwXY_clear( dest );
    if( dest->interpolation == ptwXY_interpolationOther ) {
        if( dest->interpolationOtherInfo.interpolationString != NULL ) {
            dest->interpolationOtherInfo.interpolationString = (char const *) nfu_free( (void *) dest->interpolationOtherInfo.interpolationString );
        }
    }
    dest->interpolation = ptwXY_interpolationLinLin; /* This and prior lines are in case interpolation is 'other' and ptwXY_reallocatePoints fails. */
    if( dest->allocatedSize < src->length ) ptwXY_reallocatePoints( dest, src->length, 0 );
    if( dest->status != nfu_Okay ) return( dest->status );
    dest->interpolation = src->interpolation;
    if( dest->interpolation == ptwXY_interpolationOther ) {
        if( src->interpolationOtherInfo.interpolationString != NULL ) {
            if( ( dest->interpolationOtherInfo.interpolationString = strdup( src->interpolationOtherInfo.interpolationString ) ) == NULL ) 
                return( dest->status = nfu_mallocError );
        } }
    else {
        dest->interpolationOtherInfo.interpolationString = src->interpolationOtherInfo.interpolationString;
    }
    dest->interpolationOtherInfo.getValueFunc = src->interpolationOtherInfo.getValueFunc;
    dest->interpolationOtherInfo.argList = src->interpolationOtherInfo.argList;
    dest->userFlag = src->userFlag;
    dest->biSectionMax = src->biSectionMax;
    dest->accuracy = src->accuracy;
    dest->minFractional_dx = src->minFractional_dx;
    pointFrom = src->points;
    o = src->overflowHeader.next;
    pointTo = dest->points;
    i = 0;
    while( o != overflowHeader ) {
        if( i < nonOverflowLength ) {
            if( pointFrom->x < o->point.x ) {
                *pointTo = *pointFrom;
                i++;
                pointFrom++; }
            else {
                *pointTo = o->point;
                o = o->next;
            } }
        else {
            *pointTo = o->point;
            o = o->next;
        }
        pointTo++;
    } // Loop checking, 11.06.2015, T. Koi
    for( ; i < nonOverflowLength; i++, pointFrom++, pointTo++ ) *pointTo = *pointFrom;
    dest->length = src->length;
    return( dest->status );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_clone( ptwXYPoints *ptwXY, nfu_status *status ) {

    return( ptwXY_slice( ptwXY, 0, ptwXY->length, ptwXY->overflowAllocatedSize, status ) );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_cloneToInterpolation( ptwXYPoints *ptwXY, ptwXY_interpolation interpolationTo, nfu_status *status ) {

    ptwXYPoints *n1;

    if( interpolationTo == ptwXY_interpolationOther ) {
        *status = nfu_otherInterpolation;
        return( NULL );
    }
    if( ( n1 = ptwXY_clone( ptwXY, status ) ) != NULL ) {
        if( n1->interpolation == ptwXY_interpolationOther ) nfu_free( (void *) n1->interpolationOtherInfo.interpolationString );
        n1->interpolation = interpolationTo;
        switch( interpolationTo ) {
            case ptwXY_interpolationLinLin :
                n1->interpolationOtherInfo.interpolationString = linLinInterpolationString; break;
            case ptwXY_interpolationLinLog :
                n1->interpolationOtherInfo.interpolationString = linLogInterpolationString; break;
            case ptwXY_interpolationLogLin :
                n1->interpolationOtherInfo.interpolationString = logLinInterpolationString; break;
            case ptwXY_interpolationLogLog :
                n1->interpolationOtherInfo.interpolationString = logLogInterpolationString; break;
            case ptwXY_interpolationFlat :
                n1->interpolationOtherInfo.interpolationString = flatInterpolationString; break;
            case ptwXY_interpolationOther :     /* Does not happen, but needed to stop compilers from complaining. */
                break;
        }
        n1->interpolationOtherInfo.getValueFunc = NULL;
        n1->interpolationOtherInfo.argList = NULL;
    }
    return( n1 );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_slice( ptwXYPoints *ptwXY, int64_t index1, int64_t index2, int64_t secondarySize, nfu_status *status ) {

    int64_t i, length;
    ptwXYPoints *n;

    *status = nfu_badSelf;
    if( ptwXY->status != nfu_Okay ) return( NULL );

    *status = nfu_badIndex;
    if( index2 < index1 ) return( NULL );
    if( index1 < 0 ) index1 = 0;
    if( index2 > ptwXY->length ) index2 = ptwXY->length;

    length = index2 - index1;
    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY ) ) != nfu_Okay ) return( NULL );
    if( ( n = ptwXY_new( ptwXY->interpolation, &(ptwXY->interpolationOtherInfo), ptwXY->biSectionMax, 
        ptwXY->accuracy, length, secondarySize, status, ptwXY->userFlag ) ) == NULL ) return( NULL );

    *status = n->status = ptwXY->status;
    for( i = index1; i < index2; i++ ) n->points[i - index1] = ptwXY->points[i];
    n->length = length;
    return( n );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_xSlice( ptwXYPoints *ptwXY, double xMin, double xMax, int64_t secondarySize, int fill, nfu_status *status ) {

    int64_t i, i1, i2;
    double y;
    ptwXYPoints *n = NULL;

    if( ( *status = ptwXY->status ) != nfu_Okay ) return( NULL );

    if( ( ptwXY->length == 0 ) || ( ptwXY_getXMin( ptwXY ) >= xMax ) || ( ptwXY_getXMax( ptwXY ) <= xMin ) ) {
        n = ptwXY_new( ptwXY->interpolation, &(ptwXY->interpolationOtherInfo), ptwXY->biSectionMax, 
            ptwXY->accuracy, 0, secondarySize, status, ptwXY->userFlag ); }
    else {
        if( ( n = ptwXY_clone( ptwXY, status ) ) == NULL ) return( NULL );
        if( ( n->points[0].x < xMin ) || ( n->points[n->length - 1].x > xMax ) ) {
            if( fill && ( n->points[n->length - 1].x > xMax ) ) {
                if( ( *status = ptwXY_getValueAtX( n, xMax, &y ) ) != nfu_Okay ) goto Err;
                if( ( *status = ptwXY_setValueAtX( n, xMax,  y ) ) != nfu_Okay ) goto Err;
            }
            if( fill && ( n->points[0].x < xMin ) ) {
                if( ( *status = ptwXY_getValueAtX( n, xMin, &y ) ) != nfu_Okay ) goto Err;
                if( ( *status = ptwXY_setValueAtX( n, xMin,  y ) ) != nfu_Okay ) goto Err;
            }
            ptwXY_coalescePoints( n, n->length + n->overflowAllocatedSize, NULL, 0 );
            for( i1 = 0; i1 < n->length; i1++ ) if( n->points[i1].x >= xMin ) break;
            for( i2 = n->length - 1; i2 > 0; i2-- ) if( n->points[i2].x <= xMax ) break;
            i2++;
            if( i1 > 0 ) {
                for( i = i1; i < i2; i++ ) n->points[i- i1] = n->points[i];
            }
            n->length = i2 - i1;
        }
    }
    return( n );

Err:
    if( n != NULL ) ptwXY_free( n );
    return( NULL );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_xMinSlice( ptwXYPoints *ptwXY, double xMin, int64_t secondarySize, int fill, nfu_status *status ) {

    double xMax = 1.1 * xMin + 1;

    if( xMin < 0 ) xMax = 0.9 * xMin + 1;
    if( ptwXY->length > 0 ) xMax = ptwXY_getXMax( ptwXY );
    return( ptwXY_xSlice( ptwXY, xMin, xMax, secondarySize, fill, status ) );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_xMaxSlice( ptwXYPoints *ptwXY, double xMax, int64_t secondarySize, int fill, nfu_status *status ) {

    double xMin = 0.9 * xMax - 1;

    if( xMax < 0 ) xMin = 1.1 * xMax - 1;
    if( ptwXY->length > 0 ) xMin = ptwXY_getXMin( ptwXY );
    return( ptwXY_xSlice( ptwXY, xMin, xMax, secondarySize, fill, status ) );
}
/*
************************************************************
*/
ptwXY_interpolation ptwXY_getInterpolation( ptwXYPoints *ptwXY ) {

    return( ptwXY->interpolation );
}
/*
************************************************************
*/
char const *ptwXY_getInterpolationString( ptwXYPoints *ptwXY ) {

    return( ptwXY->interpolationOtherInfo.interpolationString );
}
/*
************************************************************
*/
nfu_status ptwXY_getStatus( ptwXYPoints *ptwXY ) {

    return( ptwXY->status );
}
/*
************************************************************
*/
int ptwXY_getUserFlag( ptwXYPoints *ptwXY ) {

    return( ptwXY->userFlag );
}
/*
************************************************************
*/
void ptwXY_setUserFlag( ptwXYPoints *ptwXY, int userFlag ) {

    ptwXY->userFlag = userFlag;
}
/*
************************************************************
*/
double ptwXY_getAccuracy( ptwXYPoints *ptwXY ) {

    return( ptwXY->accuracy );
}
/*
************************************************************
*/
double ptwXY_setAccuracy( ptwXYPoints *ptwXY, double accuracy ) {

    if( accuracy < ptwXY_minAccuracy ) accuracy = ptwXY_minAccuracy;
    if( accuracy < ptwXY->accuracy ) accuracy = ptwXY->accuracy;
    if( accuracy > 1 ) accuracy = 1.;
    ptwXY->accuracy = accuracy;
    return( ptwXY->accuracy );
}
/*
************************************************************
*/
double ptwXY_getBiSectionMax( ptwXYPoints *ptwXY ) {

    return( ptwXY->biSectionMax );
}
/*
************************************************************
*/
double ptwXY_setBiSectionMax( ptwXYPoints *ptwXY, double biSectionMax ) {

    if( biSectionMax < 0 ) {
        biSectionMax = 0; }
    else if( biSectionMax > ptwXY_maxBiSectionMax ) {
        biSectionMax = ptwXY_maxBiSectionMax;
    }
    ptwXY->biSectionMax = biSectionMax;
    return( ptwXY->biSectionMax );
}
/*
************************************************************
*/
nfu_status ptwXY_reallocatePoints( ptwXYPoints *ptwXY, int64_t size, int forceSmallerResize ) {
/*
*   This is for allocating/reallocating the primary data memory.
*/
    if( ptwXY->status != nfu_Okay ) return( ptwXY->status );

    if( size < ptwXY_minimumSize ) size = ptwXY_minimumSize;                      /* ptwXY_minimumSize must be > 0. */
    if( size < ptwXY->length ) size = ptwXY->length;
    if( size != ptwXY->allocatedSize ) {
        if( size > ptwXY->allocatedSize ) {                                        /* Increase size of allocated points. */
            ptwXY->points = (ptwXYPoint *) nfu_realloc( (size_t) size * sizeof( ptwXYPoint ), ptwXY->points ); }
        else if( ( ptwXY->allocatedSize > 2 * size ) || forceSmallerResize ) {     /* Decrease size, if at least 1/2 size reduction or if forced to. */
            ptwXY->points = (ptwXYPoint *) nfu_realloc( (size_t) size * sizeof( ptwXYPoint ), ptwXY->points ); }
        else {
            size = ptwXY->allocatedSize;                                           /* Size is < ptwXY->allocatedSize, but realloc not called. */
        }
        if( ptwXY->points == NULL ) {
            ptwXY->length = 0;
            ptwXY->mallocFailedSize = size;
            size = 0;
            ptwXY->status = nfu_mallocError;
        }
        ptwXY->allocatedSize = size;
    }
    return( ptwXY->status );
}
/*
************************************************************
*/
nfu_status ptwXY_reallocateOverflowPoints( ptwXYPoints *ptwXY, int64_t size ) {
/*
*   This is for allocating/reallocating the secondary data memory.
*/
    nfu_status status = nfu_Okay;

    if( ptwXY->status != nfu_Okay ) return( ptwXY->status );

    if( size < ptwXY_minimumOverflowSize ) size = ptwXY_minimumOverflowSize;      /* ptwXY_minimumOverflowSize must be > 0. */
    if( size < ptwXY->overflowLength ) status = ptwXY_coalescePoints( ptwXY, ptwXY->length + ptwXY->overflowAllocatedSize, NULL, 0 );
    if( status == nfu_Okay ) {
        if( size != ptwXY->overflowAllocatedSize ) {
            ptwXY->overflowPoints = (ptwXYOverflowPoint *) nfu_realloc( (size_t) size * sizeof( ptwXYOverflowPoint ), ptwXY->overflowPoints );
            if( ptwXY->overflowPoints == NULL ) {
                ptwXY->length = 0;
                ptwXY->overflowLength = 0;
                ptwXY->mallocFailedSize = size;
                size = 0;
                ptwXY->status = nfu_mallocError;
            }
        }
        ptwXY->overflowAllocatedSize = size; }
    else {
        ptwXY->status = status;
    }
    return( ptwXY->status );
}
/*
************************************************************
*/
nfu_status ptwXY_coalescePoints( ptwXYPoints *ptwXY, int64_t size, ptwXYPoint *newPoint, int forceSmallerResize ) {

    int addNewPoint;
    int64_t length = ptwXY->length + ( ( newPoint != NULL ) ? 1 : 0 );
    ptwXYOverflowPoint *last = ptwXY->overflowHeader.prior;
    ptwXYPoint *pointsFrom, *pointsTo;

    if( ptwXY->status != nfu_Okay ) return( ptwXY->status );
    if( ptwXY->overflowLength == 0 ) return( nfu_Okay );

    if( size < length ) size = length;
    if( size > ptwXY->allocatedSize ) {
        if( ptwXY_reallocatePoints( ptwXY, size, forceSmallerResize ) != nfu_Okay ) return( ptwXY->status );
    }
    pointsFrom = &(ptwXY->points[ptwXY_getNonOverflowLength( ptwXY ) - 1]);
    pointsTo = &(ptwXY->points[length - 1]);
    while( last != &(ptwXY->overflowHeader) ) {
        addNewPoint = 0;
        if( newPoint != NULL ) {
            if( ( pointsFrom >= ptwXY->points ) && ( pointsFrom->x > last->point.x ) ) {
                if( newPoint->x > pointsFrom->x ) addNewPoint = 1; }
            else {
                if( newPoint->x > last->point.x ) addNewPoint = 1;
            }
            if( addNewPoint == 1 ) {
                *pointsTo = *newPoint;
                newPoint = NULL;
            }
        }
        if( addNewPoint == 0 ) {
            if( ( pointsFrom >= ptwXY->points ) && ( pointsFrom->x > last->point.x ) ) {
                *pointsTo = *pointsFrom;
                pointsFrom--; }
            else {
                *pointsTo = last->point;
                last = last->prior;
            }
        }
        pointsTo--;
    }  // Loop checking, 11.06.2015, T. Koi
    while( ( newPoint != NULL ) && ( pointsFrom >= ptwXY->points ) ) { 
        if( newPoint->x > pointsFrom->x ) {
            *pointsTo = *newPoint;
             newPoint = NULL; }
         else {
            *pointsTo = *pointsFrom;
            pointsFrom--;
         }
         pointsTo--;
    } // Loop checking, 11.06.2015, T. Koi
    if( newPoint != NULL ) *pointsTo = *newPoint;
    ptwXY->overflowHeader.prior = &(ptwXY->overflowHeader);
    ptwXY->overflowHeader.next = &(ptwXY->overflowHeader);
    ptwXY->length = length;
    ptwXY->overflowLength = 0;
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwXY_simpleCoalescePoints( ptwXYPoints *ptwXY ) {

    return( ptwXY_coalescePoints( ptwXY, ptwXY->length, NULL, 0 ) );
}
/*
************************************************************
*/
nfu_status ptwXY_clear( ptwXYPoints *ptwXY ) {

    if( ptwXY->status != nfu_Okay ) return( ptwXY->status );

    ptwXY->length = 0;
    ptwXY->overflowLength = 0;
    ptwXY->overflowHeader.prior = &(ptwXY->overflowHeader);
    ptwXY->overflowHeader.next = &(ptwXY->overflowHeader);
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwXY_release( ptwXYPoints *ptwXY ) {
/*
*   Note, this routine does not free ptwXY (i.e., it does not undo all of ptwXY_new).
*/

    if( ptwXY->interpolation == ptwXY_interpolationOther ) {
        if( ptwXY->interpolationOtherInfo.interpolationString != NULL ) 
            ptwXY->interpolationOtherInfo.interpolationString = (char const *) nfu_free( (void *) ptwXY->interpolationOtherInfo.interpolationString );
    }
    ptwXY->interpolation = ptwXY_interpolationLinLin;
    ptwXY->interpolationOtherInfo.getValueFunc = NULL;
    ptwXY->interpolationOtherInfo.argList = NULL;
    ptwXY->length = 0;
    ptwXY->allocatedSize = 0;
    ptwXY->points = (ptwXYPoint *) nfu_free( ptwXY->points );

    ptwXY->overflowLength = 0;
    ptwXY->overflowAllocatedSize = 0;
    ptwXY->overflowPoints = (ptwXYOverflowPoint *) nfu_free( ptwXY->overflowPoints );

    return( nfu_Okay );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_free( ptwXYPoints *ptwXY ) {

    if( ptwXY != NULL ) ptwXY_release( ptwXY );
    nfu_free( ptwXY );
    return( (ptwXYPoints *) NULL );
}
/*
************************************************************
*/
int64_t ptwXY_length( ptwXYPoints *ptwXY ) {

    return( ptwXY->length );
}
/*
************************************************************
*/
int64_t ptwXY_getNonOverflowLength( ptwXYPoints const *ptwXY ) {

    return( ptwXY->length - ptwXY->overflowLength );
}
/*
************************************************************
*/
nfu_status ptwXY_setXYData( ptwXYPoints *ptwXY, int64_t length, double const *xy ) {

    nfu_status status = nfu_Okay;
    int64_t i;
    ptwXYPoint *p;
    double const *d = xy;
    double xOld = 0.;

    if( length > ptwXY->allocatedSize ) {
        status = ptwXY_reallocatePoints( ptwXY, length, 0 );
        if( status != nfu_Okay ) return( status );
    }
    for( i = 0, p = ptwXY->points; i < length; i++, p++ ) {
        if( i != 0 ) {
            if( *d <= xOld ) {
                status = nfu_XNotAscending;
                ptwXY->status = nfu_XNotAscending;
                length = 0;
                break;
            }
        }
        xOld = *d;
        p->x = *(d++);
        p->y = *(d++);
    }
    ptwXY->overflowHeader.next = &(ptwXY->overflowHeader);
    ptwXY->overflowHeader.prior = &(ptwXY->overflowHeader);
    ptwXY->overflowLength = 0;
    ptwXY->length = length;
    return( ptwXY->status = status );
}
/*
************************************************************
*/  
nfu_status ptwXY_setXYDataFromXsAndYs( ptwXYPoints *ptwXY, int64_t length, double const *x, double const *y ) {

    nfu_status status;
    int64_t i;
    ptwXYPoint *p;
    double xOld = 0.;

    if( ( status = ptwXY_clear( ptwXY ) ) != nfu_Okay ) return( status );
    if( length > ptwXY->allocatedSize ) {
        if( ( status = ptwXY_reallocatePoints( ptwXY, length, 0 ) ) != nfu_Okay ) return( status );
    }
    for( i = 0, p = ptwXY->points; i < length; i++, p++, x++, y++ ) {
        if( i != 0 ) {
            if( *x <= xOld ) {
                status = ptwXY->status = nfu_XNotAscending;
                length = 0;
                break;
            }
        }
        xOld = *x;
        p->x = *x;
        p->y = *y;
    }
    ptwXY->length = length;
    return( status );
}
/*
************************************************************
*/  
nfu_status ptwXY_deletePoints( ptwXYPoints *ptwXY, int64_t i1, int64_t i2 ) {

    int64_t n = ptwXY->length - ( i2 - i1 );

    if( ( ptwXY->status = ptwXY_simpleCoalescePoints( ptwXY ) ) != nfu_Okay ) return( ptwXY->status );
    if( ( i1 < 0 ) || ( i1 > i2 ) || ( i2 > ptwXY->length ) ) return( nfu_badIndex );
    if( i1 != i2 ) {
        for( ; i2 < ptwXY->length; i1++, i2++ ) ptwXY->points[i1] = ptwXY->points[i2];
        ptwXY->length = n;
    }
    return( ptwXY->status );
}
/*
************************************************************
*/
ptwXYPoint *ptwXY_getPointAtIndex( ptwXYPoints *ptwXY, int64_t index ) {

    if( ptwXY->status != nfu_Okay ) return( NULL );
    if( ( index < 0 ) || ( index >= ptwXY->length ) ) return( NULL );
    return( ptwXY_getPointAtIndex_Unsafely( ptwXY, index ) );
}
/*
************************************************************
*/
ptwXYPoint *ptwXY_getPointAtIndex_Unsafely( ptwXYPoints *ptwXY, int64_t index ) {

    int64_t i;
    ptwXYOverflowPoint *overflowPoint;

    for( overflowPoint = ptwXY->overflowHeader.next, i = 0; overflowPoint != &(ptwXY->overflowHeader); overflowPoint = overflowPoint->next, i++ ) {
        if( overflowPoint->index == index ) return( &(overflowPoint->point) );
        if( overflowPoint->index > index ) break;
    }
    return( &(ptwXY->points[index - i]) );
}
/*
************************************************************
*/
nfu_status ptwXY_getXYPairAtIndex( ptwXYPoints *ptwXY, int64_t index, double *x, double *y ) {

    ptwXYPoint *p = ptwXY_getPointAtIndex( ptwXY, index );

    if( p == NULL ) return( nfu_badIndex );
    *x = p->x;
    *y = p->y;
    return( nfu_Okay );
}
/*
************************************************************
*/
ptwXY_lessEqualGreaterX ptwXY_getPointsAroundX( ptwXYPoints *ptwXY, double x, ptwXYOverflowPoint *lessThanEqualXPoint, ptwXYOverflowPoint *greaterThanXPoint ) {

    int closeIsEqual;
    ptwXYPoint *closePoint;

    return( ptwXY_getPointsAroundX_closeIsEqual( ptwXY, x, lessThanEqualXPoint, greaterThanXPoint, 0, &closeIsEqual, &closePoint ) );
}
/*
************************************************************
*/
ptwXY_lessEqualGreaterX ptwXY_getPointsAroundX_closeIsEqual( ptwXYPoints *ptwXY, double x, ptwXYOverflowPoint *lessThanEqualXPoint, 
        ptwXYOverflowPoint *greaterThanXPoint, double eps, int *closeIsEqual, ptwXYPoint **closePoint ) {

    int64_t overflowIndex, nonOverflowLength = ptwXY_getNonOverflowLength( ptwXY );
    int64_t indexMin, indexMid, indexMax;
    ptwXY_dataFrom xMinFrom, xMaxFrom;
    double xMin = ptwXY_getXMinAndFrom( ptwXY, &xMinFrom ), xMax = ptwXY_getXMaxAndFrom( ptwXY, &xMaxFrom );
    ptwXYOverflowPoint *overflowPoint, *overflowHeader = &(ptwXY->overflowHeader);
    ptwXY_lessEqualGreaterX status = ptwXY_lessEqualGreaterX_empty;
    ptwXYPoint *lowerPoint = NULL, *upperPoint = NULL;

    ptwXY_initialOverflowPoint( lessThanEqualXPoint, overflowHeader, NULL );
    ptwXY_initialOverflowPoint( greaterThanXPoint, overflowHeader, NULL );
    if( ptwXY->length != 0 ) {
        if( x < xMin ) {
            status = ptwXY_lessEqualGreaterX_lessThan;
            if( xMinFrom == ptwXY_dataFrom_Points ) {
                greaterThanXPoint->prior = overflowHeader;
                greaterThanXPoint->index = 0;
                greaterThanXPoint->point = ptwXY->points[0];
                *closePoint = &(ptwXY->points[0]); }
            else {
                *greaterThanXPoint = *(overflowHeader->next);
                *closePoint = &(overflowHeader->next->point);
            } }
        else if( x > xMax ) {
            status = ptwXY_lessEqualGreaterX_greater;
            if( xMaxFrom == ptwXY_dataFrom_Points ) {
                lessThanEqualXPoint->prior = overflowHeader->prior;
                lessThanEqualXPoint->index = nonOverflowLength - 1;
                lessThanEqualXPoint->point = ptwXY->points[lessThanEqualXPoint->index];
                *closePoint = &(ptwXY->points[lessThanEqualXPoint->index]); }
            else {
                *lessThanEqualXPoint = *(overflowHeader->prior);
                *closePoint = &(overflowHeader->prior->point);
            } }
        else {                                                  /* xMin <= x <= xMax */
            status = ptwXY_lessEqualGreaterX_between;           /* Default for this condition, can only be between or equal. */
            for( overflowPoint = overflowHeader->next, overflowIndex = 0; overflowPoint != overflowHeader; 
                overflowPoint = overflowPoint->next, overflowIndex++ ) if( overflowPoint->point.x > x ) break;
            overflowPoint = overflowPoint->prior;
            if( ( overflowPoint != overflowHeader ) && ( overflowPoint->point.x == x ) ) {
                status = ptwXY_lessEqualGreaterX_equal;
                *lessThanEqualXPoint = *overflowPoint; }
            else if( ptwXY->length == 1 ) {                    /* If here and length = 1, then ptwXY->points[0].x == x. */
                status = ptwXY_lessEqualGreaterX_equal;
                lessThanEqualXPoint->index = 0;
                lessThanEqualXPoint->point = ptwXY->points[0]; }
            else {                                              /* ptwXY->length > 1 */
                indexMin = 0;
                indexMax = nonOverflowLength - 1;
                indexMid = ( indexMin + indexMax ) >> 1;
                while( ( indexMin != indexMid ) && ( indexMid != indexMax ) ) {
                    if( ptwXY->points[indexMid].x > x ) {
                        indexMax = indexMid; }
                    else {
                        indexMin = indexMid;
                    }
                    indexMid = ( indexMin + indexMax ) >> 1;
                } // Loop checking, 11.06.2015, T. Koi
                if( ptwXY->points[indexMin].x == x ) {
                    status = ptwXY_lessEqualGreaterX_equal;
                    lessThanEqualXPoint->index = indexMin;
                    lessThanEqualXPoint->point = ptwXY->points[indexMin]; }
                else if( ptwXY->points[indexMax].x == x ) {
                    status = ptwXY_lessEqualGreaterX_equal;
                    lessThanEqualXPoint->index = indexMax;
                    lessThanEqualXPoint->point = ptwXY->points[indexMax]; }
                else {
                    if( ptwXY->points[indexMin].x > x ) indexMax = 0;
                    if( ptwXY->points[indexMax].x < x ) indexMin = indexMax;
                    if( ( overflowPoint == overflowHeader ) ||     /* x < xMin of overflow points. */
                            ( ( ptwXY->points[indexMin].x > overflowPoint->point.x ) && ( ptwXY->points[indexMin].x < x ) ) ) {
                        if( overflowPoint != overflowHeader ) lessThanEqualXPoint->prior = overflowPoint;
                        lowerPoint = &(ptwXY->points[indexMin]);
                        lessThanEqualXPoint->index = indexMin;
                        lessThanEqualXPoint->point = ptwXY->points[indexMin]; }
                    else {
                        lowerPoint = &(overflowPoint->point);
                        *lessThanEqualXPoint = *overflowPoint;
                    }
                    if( ( overflowPoint->next == overflowHeader ) ||     /* x > xMax of overflow points. */
                            ( ( ptwXY->points[indexMax].x < overflowPoint->next->point.x ) && ( ptwXY->points[indexMax].x > x ) ) ) {
                        upperPoint = &(ptwXY->points[indexMax]);
                        greaterThanXPoint->index = indexMax;
                        greaterThanXPoint->point = ptwXY->points[indexMax]; }
                    else {
                        upperPoint = &(overflowPoint->next->point);
                        *greaterThanXPoint = *(overflowPoint->next);
                    }
                }
            }
        }
    }

    *closeIsEqual = 0;
    if( eps > 0 ) {
        double absX = std::fabs( x );

        if( status == ptwXY_lessEqualGreaterX_lessThan ) {
            if( absX < std::fabs( greaterThanXPoint->point.x ) ) absX = std::fabs( greaterThanXPoint->point.x );
            if( ( greaterThanXPoint->point.x - x ) < eps * absX ) *closeIsEqual = 1; }
        else if( status == ptwXY_lessEqualGreaterX_greater ) {
            if( absX < std::fabs( lessThanEqualXPoint->point.x ) ) absX = std::fabs( lessThanEqualXPoint->point.x );
            if( ( x - lessThanEqualXPoint->point.x ) < eps * absX ) *closeIsEqual = -1; }
        else if( status == ptwXY_lessEqualGreaterX_between ) {
            if( ( x - lessThanEqualXPoint->point.x ) < ( greaterThanXPoint->point.x - x ) ) {   /* x is closer to lower point. */
                *closePoint = lowerPoint;
                if( absX < std::fabs( lessThanEqualXPoint->point.x ) ) absX = std::fabs( lessThanEqualXPoint->point.x );
                if( ( x - lessThanEqualXPoint->point.x ) < eps * absX ) *closeIsEqual = -1; }
            else {                                                                              /* x is closer to upper point. */
                *closePoint = upperPoint;
                if( absX < std::fabs( greaterThanXPoint->point.x ) ) absX = std::fabs( greaterThanXPoint->point.x );
                if( ( greaterThanXPoint->point.x - x ) < eps * absX ) *closeIsEqual = 1;
            } }
        else if( status == ptwXY_lessEqualGreaterX_equal ) {
            *closeIsEqual = 1;
        }
    }
    return( status );
}
/*
************************************************************
*/
nfu_status ptwXY_getValueAtX( ptwXYPoints *ptwXY, double x, double *y ) {

    nfu_status status = nfu_XOutsideDomain;
    ptwXYOverflowPoint lessThanEqualXPoint, greaterThanXPoint;
    ptwXY_lessEqualGreaterX legx = ptwXY_getPointsAroundX( ptwXY, x, &lessThanEqualXPoint, &greaterThanXPoint );

    *y = 0.;
    if( ptwXY->status != nfu_Okay ) return( ptwXY->status );
    switch( legx ) {
    case ptwXY_lessEqualGreaterX_empty :
    case ptwXY_lessEqualGreaterX_lessThan :
    case ptwXY_lessEqualGreaterX_greater :
        break;
    case ptwXY_lessEqualGreaterX_equal :
        status = nfu_Okay;
        *y = lessThanEqualXPoint.point.y;
        break;
    case ptwXY_lessEqualGreaterX_between :
        if( ptwXY->interpolationOtherInfo.getValueFunc != NULL ) {
            status = ptwXY->interpolationOtherInfo.getValueFunc( ptwXY->interpolationOtherInfo.argList, x, y, 
                lessThanEqualXPoint.point.x, lessThanEqualXPoint.point.y, greaterThanXPoint.point.x, greaterThanXPoint.point.y ); }
        else {
            status = ptwXY_interpolatePoint( ptwXY->interpolation, x, y, lessThanEqualXPoint.point.x, lessThanEqualXPoint.point.y, 
                greaterThanXPoint.point.x, greaterThanXPoint.point.y );
        }
        break;
    }
    return( status );
}
/*
************************************************************
*/
nfu_status ptwXY_setValueAtX( ptwXYPoints *ptwXY, double x, double y ) {

    return( ptwXY_setValueAtX_overrideIfClose( ptwXY, x, y, 0., 0 ) );
}
/*
************************************************************
*/
nfu_status ptwXY_setValueAtX_overrideIfClose( ptwXYPoints *ptwXY, double x, double y, double eps, int override ) {

    int closeIsEqual;
    int64_t nonOverflowLength = ptwXY_getNonOverflowLength( ptwXY ), i;
    nfu_status status = nfu_Okay;
    ptwXY_lessEqualGreaterX legx;
    ptwXYPoint *point = NULL, newPoint = { x, y };
    ptwXYOverflowPoint *overflowPoint, *p, *overflowHeader = &(ptwXY->overflowHeader);
    ptwXYOverflowPoint lessThanEqualXPoint, greaterThanXPoint;
    ptwXYPoint *closePoint = NULL;

    if( ptwXY->status != nfu_Okay ) return( ptwXY->status );

    legx = ptwXY_getPointsAroundX_closeIsEqual( ptwXY, x, &lessThanEqualXPoint, &greaterThanXPoint, eps, &closeIsEqual, &closePoint );
    switch( legx ) {
    case ptwXY_lessEqualGreaterX_lessThan :
    case ptwXY_lessEqualGreaterX_greater :
    case ptwXY_lessEqualGreaterX_between :
        if( closeIsEqual ) {
            if( !override ) return( status );
            point = closePoint;
            legx = ptwXY_lessEqualGreaterX_equal;
            x = point->x; }
        else {
            if( ( legx == ptwXY_lessEqualGreaterX_greater ) && ( nonOverflowLength < ptwXY->allocatedSize ) ) {
                point = &(ptwXY->points[nonOverflowLength]); }
            else {
                if( ptwXY->overflowLength == ptwXY->overflowAllocatedSize ) 
                    return( ptwXY_coalescePoints( ptwXY, ptwXY->length + ptwXY->overflowAllocatedSize, &newPoint, 0 ) );
                overflowPoint = &(ptwXY->overflowPoints[ptwXY->overflowLength]);
                if( legx == ptwXY_lessEqualGreaterX_lessThan ) {
                    overflowPoint->prior = greaterThanXPoint.prior;
                    overflowPoint->index = 0; }
                else {                                      /* Between or greater and must go in overflow area. */
                    if( legx == ptwXY_lessEqualGreaterX_greater ) {
                        overflowPoint->prior = overflowHeader->prior;
                        overflowPoint->index = ptwXY->length; }
                    else {
                        overflowPoint->prior = lessThanEqualXPoint.prior;
                        if( lessThanEqualXPoint.next != NULL ) {
                            if( lessThanEqualXPoint.point.x < x ) 
                                overflowPoint->prior = lessThanEqualXPoint.prior->next;
                            i = 1; }
                        else {
                            for( p = overflowHeader->next, i = 1; p != overflowHeader; p = p->next, i++ ) 
                                if( p->point.x > x ) break;
                        }
                        overflowPoint->index = lessThanEqualXPoint.index + i;
                    }
                }
                overflowPoint->next = overflowPoint->prior->next;
                overflowPoint->prior->next = overflowPoint;
                overflowPoint->next->prior = overflowPoint;
                point = &(overflowPoint->point);
                for( overflowPoint = overflowPoint->next; overflowPoint != overflowHeader; overflowPoint = overflowPoint->next ) {
                    overflowPoint->index++;
                }
                ptwXY->overflowLength++;
            }
        }
        break;
    case ptwXY_lessEqualGreaterX_empty :
        point = ptwXY->points;                 /* ptwXY_minimumSize must be > 0 so there is always space here. */
        break;
    case ptwXY_lessEqualGreaterX_equal :
        if( closeIsEqual && !override ) return( status );
        if( lessThanEqualXPoint.next == NULL ) {
            point = &(ptwXY->points[lessThanEqualXPoint.index]); }
        else {
            point = &(lessThanEqualXPoint.prior->next->point);
        }
        break;
    }
    if( status == nfu_Okay ) {
        point->x = x;
        point->y = y;
        if( legx != ptwXY_lessEqualGreaterX_equal ) ptwXY->length++;
    }
    return( status );
}
/*
************************************************************
*/
nfu_status ptwXY_mergeFromXsAndYs( ptwXYPoints *ptwXY, int length, double *xs, double *ys ) {

    return( ptwXY_mergeFrom( ptwXY, 1, length, xs, ys ) );
}
/*
************************************************************
*/
nfu_status ptwXY_mergeFromXYs( ptwXYPoints *ptwXY, int length, double *xys ) {

    int i;
    double *xs, *p1, *p2;
    nfu_status status;

    if( length < 0 ) return( nfu_badInput );
    if( length == 0 ) return( nfu_Okay );
    if( ( xs = (double *) nfu_malloc( length * sizeof( double ) ) ) == NULL ) return( nfu_mallocError );
    for( i = 0, p1 = xs, p2 = xys; i < length; i++, p1++, p2 += 2 ) *p1 = *p2;
    status = ptwXY_mergeFrom( ptwXY, 2, length, xs, xys );
    nfu_free( xs );

    return( status );
}
/*
************************************************************
*/
static nfu_status ptwXY_mergeFrom( ptwXYPoints *ptwXY, int /*incY*/, int length, double *xs, double *ys ) {

    int i, j,  n;
    double *sortedXs, *p1, *p2;
    nfu_status status;
    ptwXYPoint *point1, *point2;

    if( length < 0 ) return( nfu_badInput );
    if( length == 0 ) return( nfu_Okay );
    if( ( status = ptwXY_simpleCoalescePoints( ptwXY ) ) != nfu_Okay ) return( status );

    //if( ( sortedXs = (double *) nfu_malloc( length * sizeof( double * ) ) ) == NULL ) return( nfu_mallocError );
    //TK fixed for Coverity 63081 
    if( ( sortedXs = (double *) nfu_malloc( length * sizeof( double ) ) ) == NULL ) return( nfu_mallocError );

    for( i = 0, p1 = sortedXs, p2 = xs; i < length; i++, p1++, p2++ ) *p1 = *p2;
    //qsort( sortedXs, length, sizeof( double * ), ptwXY_mergeCompareFunction );
    //TK fixed for Coverity 63079
    qsort( sortedXs, length, sizeof( double ), ptwXY_mergeCompareFunction );

    for( i = 0, p1 = sortedXs, j = 0, n = 0; i < length; i++, p1++, n++ ) {
        for( ; j < ptwXY->length; j++, n++ ) {
            if( *p1 <= ptwXY->points[j].x ) break;
        }
        if( j == ptwXY->length ) break;             /* Completed all ptwXY points. */
    }
    n += (int) ( ( length - i ) + ( ptwXY->length - j ) );

    if( ( status = ptwXY_reallocatePoints( ptwXY, n, 0 ) ) == nfu_Okay ) {
        point1 = &(ptwXY->points[n-1]);
        point2 = &(ptwXY->points[length-1]);
        for( i = 0, j = 0, p1 = &(sortedXs[length-1]); ( i < length ) && ( j < length ) && ( n > 0 ); n--, point1-- ) {
            if( *p1 >= point2->x ) {
                point1->x = *p1;
                point1->y = ys[(int)(p1 - xs)];
                if( *p1 >= point2->x ) {
                    point2++;
                    j++;
                }
                p1--;
                i--; }
            else {
                *point1 = *point2;
                point2++;
                j++;
            }
        }
        for( ; i < length; i++, p1--, point1-- ) {
            point1->x = *p1;
            point1->y = ys[(int)(p1 - xs)];
        }
        for( ; j < length; j++, point1--, point2-- ) *point1 = *point2;
    }
    nfu_free( sortedXs );

    return( status );
}
/*
************************************************************
*/
static int ptwXY_mergeCompareFunction( void const *x1p, void const *x2p ) {

    double d1 = *((double *) x1p), d2 = *((double *) x2p);

    if( d1  < d2 ) return( -1 );
    if( d1 == d2 ) return(  0 );
    return( 1 );
}
/*
************************************************************
*/
nfu_status ptwXY_appendXY( ptwXYPoints *ptwXY, double x, double y ) {

    int64_t nonOverflowLength = ptwXY_getNonOverflowLength( ptwXY );
    ptwXY_dataFrom dataFrom;

    if( ptwXY->length != 0 ) {
        double xMax = ptwXY_getXMaxAndFrom( ptwXY, &dataFrom );
        if( xMax >= x ) return( nfu_XNotAscending );
    }

    if( nonOverflowLength < ptwXY->allocatedSize ) {      /* Room at end of points. Also handles the case when length = 0. */
        ptwXY->points[nonOverflowLength].x = x;
        ptwXY->points[nonOverflowLength].y = y; }
    else {
        if( ptwXY->overflowLength == ptwXY->overflowAllocatedSize ) {
            ptwXYPoint newPoint = { x, y };
            return( ptwXY_coalescePoints( ptwXY, ptwXY->length + ptwXY->overflowAllocatedSize, &newPoint, 0 ) ); }
        else {                                              /* Add to end of overflow. */
            ptwXYOverflowPoint *overflowPoint = &(ptwXY->overflowPoints[ptwXY->overflowLength]);

            overflowPoint->prior = ptwXY->overflowHeader.prior;
            overflowPoint->next = overflowPoint->prior->next;
            overflowPoint->index = ptwXY->length;
            overflowPoint->prior->next = overflowPoint;
            overflowPoint->next->prior = overflowPoint;
            overflowPoint->point.x = x;
            overflowPoint->point.y = y;
            ptwXY->overflowLength++;
        }
    }
    ptwXY->length++;
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwXY_setXYPairAtIndex( ptwXYPoints *ptwXY, int64_t index, double x, double y ) {

    int64_t i, ip1;
    ptwXYOverflowPoint *overflowPoint, *pm1, *pp1;

    if( ptwXY->status != nfu_Okay ) return( ptwXY->status );

    if( ( index < 0 ) || ( index >= ptwXY->length ) ) return( nfu_badIndex );
    for( overflowPoint = ptwXY->overflowHeader.next, i = 0; overflowPoint != &(ptwXY->overflowHeader); overflowPoint = overflowPoint->next, i++ ) {
        if( overflowPoint->index >= index ) break;
    }
    ip1 = i;
    pm1 = pp1 = overflowPoint;
    if( overflowPoint->index == index ) {                                           /* Note, if overflowPoint is header, then its index = -1. */
        pp1 = overflowPoint->next;
        ip1++;
    }
    if( ( pp1 != &(ptwXY->overflowHeader) ) && ( pp1->index == ( index + 1 ) ) ) {     /* This if and else check that x < element[index+1]'s x values. */
        if( pp1->point.x <= x ) return( nfu_badIndexForX ); }
    else {
        if( ( ( index + 1 ) < ptwXY->length ) && ( ptwXY->points[index + 1 - ip1].x <= x ) ) return( nfu_badIndexForX );
    }
    if( overflowPoint != &(ptwXY->overflowHeader) ) pm1 = overflowPoint->prior;
    if( ( pm1 != &(ptwXY->overflowHeader) ) && ( pm1->index == ( index - 1 ) ) ) {     /* This if and else check that x > element[index-1]'s x values. */
        if( pm1->point.x >= x ) return( nfu_badIndexForX ); }
    else {
        if( ( ( index - 1 ) >= 0 ) && ( ptwXY->points[index - 1 - i].x >= x ) ) return( nfu_badIndexForX );
    }
    if( ( overflowPoint != &(ptwXY->overflowHeader) ) && ( overflowPoint->index == index ) ) {
        overflowPoint->point.x = x;
        overflowPoint->point.y = y; }
    else {
        index -= i;
        ptwXY->points[index].x = x;
        ptwXY->points[index].y = y;
    }
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwXY_getSlopeAtX( ptwXYPoints *ptwXY, double x, const char side, double *slope ) {

    nfu_status status  = nfu_Okay;
    ptwXYOverflowPoint lessThanEqualXPoint, greaterThanXPoint;
    ptwXY_lessEqualGreaterX legx = ptwXY_getPointsAroundX( ptwXY, x, &lessThanEqualXPoint, &greaterThanXPoint );
    ptwXYPoint *point;

    *slope = 0.;
    if( ( side != '-' ) && ( side != '+' ) ) return( nfu_badInput );

    switch( legx ) {
    case ptwXY_lessEqualGreaterX_empty :
    case ptwXY_lessEqualGreaterX_lessThan :
    case ptwXY_lessEqualGreaterX_greater :
        status = nfu_XOutsideDomain;
        break;
    case ptwXY_lessEqualGreaterX_between :
        *slope = ( greaterThanXPoint.point.y - lessThanEqualXPoint.point.y ) / 
            ( greaterThanXPoint.point.x - lessThanEqualXPoint.point.x );
        break;
    case ptwXY_lessEqualGreaterX_equal :
        if( side == '-' ) {
            if( lessThanEqualXPoint.index == 0 ) {
                status = nfu_XOutsideDomain; }
            else {
                point = ptwXY_getPointAtIndex_Unsafely( ptwXY, lessThanEqualXPoint.index - 1 );
                *slope = ( lessThanEqualXPoint.point.y - point->y ) / ( lessThanEqualXPoint.point.x - point->x );
            } }
        else {
            if( lessThanEqualXPoint.index == ( ptwXY->length - 1 ) ) {
                status = nfu_XOutsideDomain; }
            else {
                point = ptwXY_getPointAtIndex_Unsafely( ptwXY, lessThanEqualXPoint.index + 1 );
                *slope = ( point->y - lessThanEqualXPoint.point.y ) / ( point->x - lessThanEqualXPoint.point.x );
            }
        }
    }

    return( status );
}
/*
************************************************************
*/
double ptwXY_getXMinAndFrom( ptwXYPoints *ptwXY, ptwXY_dataFrom *dataFrom ) {

    int64_t nonOverflowLength = ptwXY_getNonOverflowLength( ptwXY );
    double xMin = nfu_getNAN( );

    *dataFrom = ptwXY_dataFrom_Unknown;
    if( ptwXY->overflowLength > 0 ) {
        *dataFrom = ptwXY_dataFrom_Overflow;
        xMin = ptwXY->overflowHeader.next->point.x;
        if( nonOverflowLength >= 0 ) {
            if( xMin > ptwXY->points[0].x ) {
                *dataFrom = ptwXY_dataFrom_Points;
                xMin = ptwXY->points[0].x;
            }
        } }
    else if( nonOverflowLength > 0 ) {
        *dataFrom = ptwXY_dataFrom_Points;
        xMin = ptwXY->points[0].x;
    }
    return( xMin );
}
/*
************************************************************
*/
double ptwXY_getXMin( ptwXYPoints *ptwXY ) {

    ptwXY_dataFrom dataFrom;

    return( ptwXY_getXMinAndFrom( ptwXY, &dataFrom ) );
}
/*
************************************************************
*/
double ptwXY_getXMaxAndFrom( ptwXYPoints *ptwXY, ptwXY_dataFrom *dataFrom ) {

    int64_t nonOverflowLength = ptwXY_getNonOverflowLength( ptwXY );
    double xMax = nfu_getNAN( );

    *dataFrom = ptwXY_dataFrom_Unknown;
    if( ptwXY->overflowLength > 0 ) {
        *dataFrom = ptwXY_dataFrom_Overflow;
        xMax = ptwXY->overflowHeader.prior->point.x;
        if( ( nonOverflowLength > 0 ) ) {
            if( xMax < ptwXY->points[nonOverflowLength-1].x ) {
                *dataFrom = ptwXY_dataFrom_Points;
                xMax = ptwXY->points[nonOverflowLength-1].x;
            }
        } }
    else if( ptwXY->length > 0 ) {
        *dataFrom = ptwXY_dataFrom_Points;
        xMax = ptwXY->points[nonOverflowLength-1].x;
    }
    return( xMax );
}
/*
************************************************************
*/
double ptwXY_getXMax( ptwXYPoints *ptwXY ) {

    ptwXY_dataFrom dataFrom;

    return( ptwXY_getXMaxAndFrom( ptwXY, &dataFrom ) );
}
/*
************************************************************
*/
double ptwXY_getYMin( ptwXYPoints *ptwXY ) {

    int64_t i, n = ptwXY_getNonOverflowLength( ptwXY  );
    ptwXYPoint *p = ptwXY->points;
    ptwXYOverflowPoint *overflowPoint = ptwXY->overflowHeader.next;
    double yMin;

    if( ptwXY->length == 0 ) return( 0. );
    if( n > 0 ) {
        yMin = p->y;
        for( i = 1, p++; i < n; i++, p++ ) yMin = ( ( yMin < p->y ) ? yMin : p->y ); }
    else {
        yMin = overflowPoint->point.y;
    }
    for( ; overflowPoint != &(ptwXY->overflowHeader); overflowPoint = overflowPoint->next ) 
        yMin = ( ( yMin < overflowPoint->point.y ) ? yMin : overflowPoint->point.y );
    return( yMin );
}
/*
************************************************************
*/
double ptwXY_getYMax( ptwXYPoints *ptwXY ) {

    int64_t i, n = ptwXY_getNonOverflowLength( ptwXY  );
    ptwXYPoint *p = ptwXY->points;
    ptwXYOverflowPoint *overflowPoint = ptwXY->overflowHeader.next;
    double yMax;

    if( ptwXY->length == 0 ) return( 0. );
    if( n > 0 ) {
        yMax = p->y;
        for( i = 1, p++; i < n; i++, p++ ) yMax = ( ( yMax > p->y ) ? yMax : p->y ); }
    else {
        yMax = overflowPoint->point.y;
    }
    for( ; overflowPoint != &(ptwXY->overflowHeader); overflowPoint = overflowPoint->next )
        yMax = ( ( yMax > overflowPoint->point.y ) ? yMax : overflowPoint->point.y );
    return( yMax );
}
/*
************************************************************
*/
static void ptwXY_initialOverflowPoint( ptwXYOverflowPoint *overflowPoint, ptwXYOverflowPoint *prior, ptwXYOverflowPoint *next ) {

    overflowPoint->prior = prior;
    overflowPoint->next = next;
    overflowPoint->index = -1;
    overflowPoint->point.x = 0.;
    overflowPoint->point.y = 0.;
}

#if defined __cplusplus
}
#endif
