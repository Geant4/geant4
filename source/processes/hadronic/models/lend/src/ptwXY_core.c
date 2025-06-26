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

#include "ptwXY.h"

/* Need to change these when conversion to Fudge3 is complete. */
static char const linLinInterpolationString[] = "lin-lin";
static char const logLinInterpolationString[] = "log-lin";
static char const linLogInterpolationString[] = "lin-log";
static char const logLogInterpolationString[] = "log-log";
static char const flatInterpolationString[] = "flat";

static nfu_status ptwXY_mergeFrom( statusMessageReporting *smr, ptwXYPoints *ptwXY, int incY, int length, double *xs, double *ys );
static void ptwXY_initialOverflowPoint( ptwXYOverflowPoint *overflowPoint, ptwXYOverflowPoint *prior, ptwXYOverflowPoint *next );
/*
************************************************************
*/
ptwXYPoints *ptwXY_new( statusMessageReporting *smr, ptwXY_interpolation interpolation, char const *interpolationString, 
        double biSectionMax, double accuracy, int64_t primarySize, int64_t secondarySize, int userFlag ) {

    ptwXYPoints *ptwXY = (ptwXYPoints *) smr_malloc2( smr, sizeof( ptwXYPoints ), 1, "ptwXY" );

    if( ptwXY == NULL ) return( NULL );
    if( ptwXY_initialize( smr, ptwXY, interpolation, interpolationString, biSectionMax, accuracy, primarySize, secondarySize, userFlag ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        smr_freeMemory2( ptwXY );
    }
    return( ptwXY );
}

/*
************************************************************
*/
ptwXYPoints *ptwXY_new2( statusMessageReporting *smr, ptwXY_interpolation interpolation, int64_t primarySize, int64_t secondarySize ) {

    char const *interpolationString = ptwXY_interpolationToString( interpolation );

    return( ptwXY_new( smr, interpolation, interpolationString, 12, 1e-3, primarySize, secondarySize, 0 ) );
}
/*
************************************************************
*/
nfu_status ptwXY_initialize( statusMessageReporting *smr, ptwXYPoints *ptwXY, ptwXY_interpolation interpolation, 
        char const *interpolationString, double biSectionMax, double accuracy, int64_t primarySize, int64_t secondarySize, 
        int userFlag ) {

    ptwXY->status = nfu_Okay;
    ptwXY->interpolation = interpolation;
    ptwXY->interpolationString = NULL;
    switch( interpolation ) {
    case ptwXY_interpolationLinLin :
        ptwXY->interpolationString = linLinInterpolationString; break;
    case ptwXY_interpolationLogLin :
        ptwXY->interpolationString = logLinInterpolationString; break;
    case ptwXY_interpolationLinLog :
        ptwXY->interpolationString = linLogInterpolationString; break;
    case ptwXY_interpolationLogLog :
        ptwXY->interpolationString = logLogInterpolationString; break;
    case ptwXY_interpolationFlat :
        ptwXY->interpolationString = flatInterpolationString; break;
    case ptwXY_interpolationOther :         /* For ptwXY_interpolationOther, interpolationString must be defined. */
        if( interpolationString == NULL ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_otherInterpolation,
                    "Invalid other interplation. interpolationString is NULL, it must be a defined string" );
            ptwXY->status = nfu_Error; }
        else {
            if( ( ptwXY->interpolationString = smr_allocateCopyString2( smr, interpolationString, "interpolationString" ) ) == NULL ) 
                    ptwXY->status = nfu_Error;
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

    if( ptwXY_reallocatePoints( smr, ptwXY, primarySize, 0 ) == nfu_Okay )
        ptwXY_reallocateOverflowPoints( smr, ptwXY, secondarySize );
    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        ptwXY_release( smr, ptwXY );
    }
    return( ptwXY->status );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_create( statusMessageReporting *smr, ptwXY_interpolation interpolation, char const *interpolationString,
    double biSectionMax, double accuracy, int64_t primarySize, int64_t secondarySize, int64_t length, double const *xy, 
    int userFlag ) {

    ptwXYPoints *ptwXY;

    if( primarySize < length ) primarySize = length;
    if( ( ptwXY = ptwXY_new( smr, interpolation, interpolationString, biSectionMax, accuracy, primarySize, 
            secondarySize, userFlag ) ) != NULL ) {
        if( ptwXY_setXYData( smr, ptwXY, length, xy ) != nfu_Okay ) ptwXY = ptwXY_free( ptwXY );
    }

    if( ptwXY == NULL ) smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( ptwXY );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_create2( statusMessageReporting *smr, ptwXY_interpolation interpolation, 
        int64_t primarySize, int64_t secondarySize, int64_t length, double const *xy, int userFlag ) {

    char const *interpolationString = ptwXY_interpolationToString( interpolation );

    return( ptwXY_create( smr, interpolation, interpolationString, 12, 1e-3, primarySize, secondarySize, length, xy, userFlag ) );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_createFrom_Xs_Ys( statusMessageReporting *smr, ptwXY_interpolation interpolation, char const *interpolationString,
    double biSectionMax, double accuracy, int64_t primarySize, int64_t secondarySize, int64_t length, double const *Xs, 
    double const *Ys, int userFlag ) {

    int i;
    ptwXYPoints *ptwXY;

    if( primarySize < length ) primarySize = length;
    if( ( ptwXY = ptwXY_new( smr, interpolation, interpolationString, biSectionMax, accuracy, primarySize, 
            secondarySize, userFlag ) ) != NULL ) {
        for( i = 0; i < length; i++ ) {
            ptwXY->points[i].x = Xs[i];
            ptwXY->points[i].y = Ys[i];
        }
        ptwXY->length = length;
    }

    if( ptwXY == NULL ) smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( ptwXY );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_createFrom_Xs_Ys2( statusMessageReporting *smr, ptwXY_interpolation interpolation, int64_t primarySize, 
    int64_t secondarySize, int64_t length, double const *Xs, double const *Ys, int userFlag ) {

    char const *interpolationString = ptwXY_interpolationToString( interpolation );

    return( ptwXY_createFrom_Xs_Ys( smr, interpolation, interpolationString, 12, 1e-3, primarySize, secondarySize, length, Xs, Ys, userFlag ) );
}
/*
************************************************************
*/
nfu_status ptwXY_copy( statusMessageReporting *smr, ptwXYPoints *dest, ptwXYPoints *src ) {

    int64_t i, nonOverflowLength;
    ptwXYPoint *pointFrom, *pointTo;
    ptwXYOverflowPoint *o, *overflowHeader = &(src->overflowHeader);

    if( dest->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid destination." );
        return( dest->status );
    }
    if( src->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( src->status );
    }

    nonOverflowLength = ptwXY_getNonOverflowLength( smr, src );     /* No need to check return value. */

    if( ptwXY_clear( smr, dest ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( dest->status );
    }
    if( dest->interpolation == ptwXY_interpolationOther ) {
        if( dest->interpolationString != NULL ) {
            dest->interpolationString = (char const *) smr_freeMemory2( dest->interpolationString );
        }
    }
    dest->interpolation = ptwXY_interpolationLinLin; /* This and prior lines are in case interpolation is 'other' and ptwXY_reallocatePoints fails. */
    if( dest->allocatedSize < src->length ) {
        if( ptwXY_reallocatePoints( smr, dest, src->length, 0 ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( dest->status );
        }
    }
    if( dest->status != nfu_Okay ) return( dest->status );
    dest->interpolation = src->interpolation;
    if( dest->interpolation == ptwXY_interpolationOther ) {
        if( src->interpolationString != NULL ) {
            if( ( dest->interpolationString = smr_allocateCopyString2( smr, src->interpolationString, "interpolationString" ) ) == NULL ) {
                smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
                return( dest->status = nfu_Error );
            }
        } }
    else {
        dest->interpolationString = src->interpolationString;
    }
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
    }
    for( ; i < nonOverflowLength; i++, pointFrom++, pointTo++ ) *pointTo = *pointFrom;
    dest->length = src->length;
    return( dest->status );
}
/*
************************************************************
*/
nfu_status ptwXY_copyPointsOnly( statusMessageReporting *smr, ptwXYPoints *dest, ptwXYPoints *src ) {

    int64_t i, nonOverflowLength;
    ptwXYPoint *pointFrom, *pointTo;
    ptwXYOverflowPoint *o, *overflowHeader = &(src->overflowHeader);

    if( dest->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid destination." );
        return( dest->status );
    }
    if( src->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( src->status );
    }

    nonOverflowLength = ptwXY_getNonOverflowLength( smr, src );     /* No need to check return value. */

    if( ptwXY_clear( smr, dest ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( dest->status );
    }

    if( dest->allocatedSize < src->length ) {
        if( ptwXY_reallocatePoints( smr, dest, src->length, 0 ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( dest->status );
        }
    }
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
    }
    for( ; i < nonOverflowLength; i++, pointFrom++, pointTo++ ) *pointTo = *pointFrom;
    dest->length = src->length;
    return( dest->status );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_clone( statusMessageReporting *smr, ptwXYPoints *ptwXY ) {

    ptwXYPoints *ptwXY2 = ptwXY_slice( smr, ptwXY, 0, ptwXY->length, ptwXY->overflowAllocatedSize );

    if( ptwXY2 == NULL ) smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( ptwXY2 );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_clone2( statusMessageReporting *smr, ptwXYPoints const *ptwXY ) {

    int64_t length = ptwXY->length;
    ptwXYPoints *ptwXY2 = NULL;
    ptwXYPoint *pointsFrom, *pointsTo;
    ptwXYOverflowPoint *last = ptwXY->overflowHeader.prior;

    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( NULL );
    }

    ptwXY2 = ptwXY_new( smr, ptwXY->interpolation, ptwXY->interpolationString,
            ptwXY->biSectionMax, ptwXY->accuracy, length, ptwXY->overflowAllocatedSize, ptwXY->userFlag );
    if( ptwXY2 == NULL ) smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );

    pointsFrom = &(ptwXY->points[ptwXY_getNonOverflowLength( smr, ptwXY ) - 1]);
    pointsTo = &(ptwXY2->points[length - 1]);
    while( last != &(ptwXY->overflowHeader) ) {
        if( ( pointsFrom >= ptwXY->points ) && ( pointsFrom->x > last->point.x ) ) {
            *pointsTo = *pointsFrom;
            --pointsFrom; }
        else {
            *pointsTo = last->point;
            last = last->prior;
        }
        --pointsTo;
    }

    for( ; pointsFrom >= ptwXY->points; --pointsFrom, --pointsTo ) *pointsTo = *pointsFrom;
    ptwXY2->length = length;

    return( ptwXY2 );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_cloneToInterpolation( statusMessageReporting *smr, ptwXYPoints *ptwXY, ptwXY_interpolation interpolationTo ) {

    ptwXYPoints *n1;

    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( NULL );
    }

/* Other interpolation should probably be allowed. */
    if( interpolationTo == ptwXY_interpolationOther ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_otherInterpolation, "Other interpolation not allowed." );
        return( NULL );
    }
    if( ( n1 = ptwXY_clone( smr, ptwXY ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." ); }
    else {
        if( n1->interpolation == ptwXY_interpolationOther ) smr_freeMemory2( n1->interpolationString );
        n1->interpolation = interpolationTo;
        switch( interpolationTo ) {
            case ptwXY_interpolationLinLin :
                n1->interpolationString = linLinInterpolationString; break;
            case ptwXY_interpolationLogLin :
                n1->interpolationString = logLinInterpolationString; break;
            case ptwXY_interpolationLinLog :
                n1->interpolationString = linLogInterpolationString; break;
            case ptwXY_interpolationLogLog :
                n1->interpolationString = logLogInterpolationString; break;
            case ptwXY_interpolationFlat :
                n1->interpolationString = flatInterpolationString; break;
            case ptwXY_interpolationOther :     /* Does not happen, but needed to stop compilers from complaining. */
                break;
        }
    }
    return( n1 );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_slice( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t index1, int64_t index2, int64_t secondarySize ) {

    int64_t i, length;
    ptwXYPoints *n;

    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( NULL );
    }

    if( ( index1 < 0 ) || ( index2 < index1 ) || ( index2 > ptwXY->length ) ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badIndex, "Indices = %d, %d out of bounds: length = %d",
                (int) index1, (int) index2, (int) ptwXY->length );
        return( NULL );
    }

    length = index2 - index1;
    if( ptwXY_simpleCoalescePoints( smr, ptwXY ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }
    if( ( n = ptwXY_new( smr, ptwXY->interpolation, ptwXY->interpolationString, ptwXY->biSectionMax, 
            ptwXY->accuracy, length, secondarySize, ptwXY->userFlag ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }

    for( i = index1; i < index2; i++ ) n->points[i - index1] = ptwXY->points[i];
    n->length = length;
    return( n );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_domainSlice( statusMessageReporting *smr, ptwXYPoints *ptwXY, double domainMin, double domainMax, 
        int64_t secondarySize, int fill ) {

    int64_t i, i1, i2;
    double y, _domainMin, _domainMax;
    ptwXYPoints *n = NULL;

    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( NULL );
    }

    if( ptwXY_domainMin( smr, ptwXY, &_domainMin ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }
    if( ptwXY_domainMax( smr, ptwXY, &_domainMax ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }

    if( ( ptwXY->length == 0 ) || ( _domainMin >= domainMax ) || ( _domainMax <= domainMin ) ) {
        if( ( n = ptwXY_new( smr, ptwXY->interpolation, ptwXY->interpolationString, ptwXY->biSectionMax, 
                ptwXY->accuracy, 0, secondarySize, ptwXY->userFlag ) ) == NULL ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        } }
    else {
        if( ( n = ptwXY_clone( smr, ptwXY ) ) == NULL ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( NULL );
        }
        if( ( n->points[0].x < domainMin ) || ( n->points[n->length - 1].x > domainMax ) ) {
            if( fill && ( n->points[n->length - 1].x > domainMax ) ) {
                if( ptwXY_getValueAtX( smr, n, domainMax, &y ) != nfu_Okay ) goto Err;
                if( ptwXY_setValueAtX( smr, n, domainMax,  y ) != nfu_Okay ) goto Err;
            }
            if( fill && ( n->points[0].x < domainMin ) ) {
                if( ptwXY_getValueAtX( smr, n, domainMin, &y ) != nfu_Okay ) goto Err;
                if( ptwXY_setValueAtX( smr, n, domainMin,  y ) != nfu_Okay ) goto Err;
            }
            if( ptwXY_coalescePoints( smr, n, n->length + n->overflowAllocatedSize, NULL, 0 ) != nfu_Okay ) goto Err;
            for( i1 = 0; i1 < n->length; i1++ ) if( n->points[i1].x >= domainMin ) break;
            for( i2 = n->length - 1; i2 > 0; i2-- ) if( n->points[i2].x <= domainMax ) break;
            i2++;
            if( i1 > 0 ) {
                for( i = i1; i < i2; i++ ) n->points[i- i1] = n->points[i];
            }
            n->length = i2 - i1;
        }
    }
    return( n );

Err:
    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    if( n != NULL ) ptwXY_free( n );
    return( NULL );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_domainMinSlice( statusMessageReporting *smr, ptwXYPoints *ptwXY, double domainMin, int64_t secondarySize, int fill ) {

    double domainMax = 1.1 * domainMin + 1;
    ptwXYPoints *ptwXY2;

    if( domainMin < 0 ) domainMax = 0.9 * domainMin + 1;
    if( ptwXY->length > 0 ) {
        if( ptwXY_domainMax( smr, ptwXY, &domainMax ) != nfu_Okay )
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    }
    if( ( ptwXY2 = ptwXY_domainSlice( smr, ptwXY, domainMin, domainMax, secondarySize, fill ) ) == NULL )
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( ptwXY2 );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_domainMaxSlice( statusMessageReporting *smr, ptwXYPoints *ptwXY, double domainMax, int64_t secondarySize, int fill ) {

    double domainMin = 0.9 * domainMax - 1;
    ptwXYPoints *ptwXY2;

    if( domainMax < 0 ) domainMin = 1.1 * domainMax - 1;
    if( ptwXY->length > 0 ) {
        if( ptwXY_domainMin( smr, ptwXY, &domainMin ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( NULL );
        }
    }
    if( ( ptwXY2 = ptwXY_domainSlice( smr, ptwXY, domainMin, domainMax, secondarySize, fill ) ) == NULL )
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( ptwXY2 );
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

    return( ptwXY->interpolationString );
}
/*
************************************************************
*/
nfu_status ptwXY_setInterpolationString( ptwXYPoints *ptwXY, char const *interpolationString ) {

    ptwXY_interpolation interpolation = ptwXY_stringToInterpolation( interpolationString );

    if( interpolation == ptwXY_interpolationOther ) return( nfu_invalidInterpolation );

    ptwXY->interpolation = interpolation;
    ptwXY->interpolationString = ptwXY_interpolationToString( interpolation );
    return( nfu_Okay );
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

    accuracy = ptwXY_limitAccuracy( accuracy );
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
nfu_status ptwXY_reallocatePoints( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t size, int forceSmallerResize ) {
/*
*   This is for allocating/reallocating the primary data memory.
*/
    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( ptwXY->status );
    }

    if( size < ptwXY_minimumSize ) size = ptwXY_minimumSize;                      /* ptwXY_minimumSize must be > 0. */
    if( size < ptwXY->length ) size = ptwXY->length;
    if( size != ptwXY->allocatedSize ) {
        if( size > ptwXY->allocatedSize ) {                                        /* Increase size of allocated points. */
            ptwXY->points = (ptwXYPoint *) smr_realloc2( smr, ptwXY->points, (size_t) size * sizeof( ptwXYPoint ), "ptwXY->points" ); }
        else if( ( ptwXY->allocatedSize > 2 * size ) || forceSmallerResize ) {     /* Decrease size, if at least 1/2 size reduction or if forced to. */
            ptwXY->points = (ptwXYPoint *) smr_realloc2( smr, ptwXY->points, (size_t) size * sizeof( ptwXYPoint ), "ptwXY->points" ); }
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
nfu_status ptwXY_reallocateOverflowPoints( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t size ) {
/*
*   This is for allocating/reallocating the secondary data memory.
*/
    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( ptwXY->status );
    }

    if( size < ptwXY_minimumOverflowSize ) size = ptwXY_minimumOverflowSize;      /* ptwXY_minimumOverflowSize must be > 0. */
    if( size < ptwXY->overflowLength ) {
        if( ptwXY_coalescePoints( smr, ptwXY, ptwXY->length + ptwXY->overflowAllocatedSize, NULL, 0 ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( ptwXY->status );
        }
    }
    if( size != ptwXY->overflowAllocatedSize ) {
        ptwXY->overflowPoints = (ptwXYOverflowPoint *) smr_realloc2( smr, ptwXY->overflowPoints, 
                (size_t) size * sizeof( ptwXYOverflowPoint ), "ptwXY->overflowPoints" );
        if( ptwXY->overflowPoints == NULL ) {
            ptwXY->length = 0;
            ptwXY->overflowLength = 0;
            ptwXY->mallocFailedSize = size;
            size = 0;
            ptwXY->status = nfu_mallocError;
        }
    }
    ptwXY->overflowAllocatedSize = size;
    return( ptwXY->status );
}
/*
************************************************************
*/
nfu_status ptwXY_coalescePoints( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t size, 
        ptwXYPoint *newPoint, int forceSmallerResize ) {

    int addNewPoint;
    int64_t length = ptwXY->length + ( ( newPoint != NULL ) ? 1 : 0 );
    ptwXYOverflowPoint *last = ptwXY->overflowHeader.prior;
    ptwXYPoint *pointsFrom, *pointsTo;

    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( ptwXY->status );
    }
    if( ptwXY->overflowLength == 0 ) return( nfu_Okay );

    if( size < length ) size = length;
    if( size > ptwXY->allocatedSize ) {
        if( ptwXY_reallocatePoints( smr, ptwXY, size, forceSmallerResize ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( ptwXY->status );
        }
    }
    pointsFrom = &(ptwXY->points[ptwXY_getNonOverflowLength( smr, ptwXY ) - 1]);
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
    }
    while( ( newPoint != NULL ) && ( pointsFrom >= ptwXY->points ) ) {
        if( newPoint->x > pointsFrom->x ) {
            *pointsTo = *newPoint;
             newPoint = NULL; }
         else {
            *pointsTo = *pointsFrom;
            pointsFrom--;
         }
         pointsTo--;
    }
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
nfu_status ptwXY_simpleCoalescePoints( statusMessageReporting *smr, ptwXYPoints *ptwXY ) {

    if( ptwXY_coalescePoints( smr, ptwXY, ptwXY->length, NULL, 0 ) != nfu_Okay )
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( ptwXY->status );
}
/*
************************************************************
*/
nfu_status ptwXY_clear( statusMessageReporting *smr, ptwXYPoints *ptwXY ) {

    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( ptwXY->status );
    }

    ptwXY->length = 0;
    ptwXY->overflowLength = 0;
    ptwXY->overflowHeader.prior = &(ptwXY->overflowHeader);
    ptwXY->overflowHeader.next = &(ptwXY->overflowHeader);
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwXY_release( statusMessageReporting *smr, ptwXYPoints *ptwXY ) {
/*
*   Note, this routine does not free ptwXY (i.e., it does not undo all of ptwXY_new).
*/

    if( ptwXY->interpolation == ptwXY_interpolationOther ) {
        if( ptwXY->interpolationString != NULL ) 
            ptwXY->interpolationString = (char const *) smr_freeMemory2( ptwXY->interpolationString );
    }
    ptwXY->interpolation = ptwXY_interpolationLinLin;
    ptwXY->length = 0;
    ptwXY->allocatedSize = 0;
    ptwXY->points = (ptwXYPoint *) smr_freeMemory2( ptwXY->points );

    ptwXY->overflowLength = 0;
    ptwXY->overflowAllocatedSize = 0;
    ptwXY->overflowPoints = (ptwXYOverflowPoint *) smr_freeMemory2( ptwXY->overflowPoints );

    return( nfu_Okay );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_free( ptwXYPoints *ptwXY ) {

    if( ptwXY != NULL ) {
        ptwXY_release( NULL, ptwXY );
        smr_freeMemory2( ptwXY );
    }
    return( (ptwXYPoints *) NULL );
}
/*
************************************************************
*/
int64_t ptwXY_length( statusMessageReporting *smr, ptwXYPoints *ptwXY ) {

    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( -ptwXY->status );
    }

    return( ptwXY->length );
}
/*
************************************************************
*/
int64_t ptwXY_getNonOverflowLength( statusMessageReporting *smr, ptwXYPoints const *ptwXY ) {

    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( -ptwXY->status );
    }

    return( ptwXY->length - ptwXY->overflowLength );
}

/*
************************************************************
*/
nfu_status ptwXY_startIndex( statusMessageReporting *a_smr, ptwXYPoints *a_ptwXY, double a_x, int64_t *a_startIndex, int64_t *a_length ) {
/*
    Sets *a_startIndex to -2 if a_x < domainMin, -1 if a_x > domainMax, otherwise to the lowest index in a_ptwXY->points 
where a_x >= a_ptwXY->points[*a_startIndex]. The logic below guarantees that *a_startIndex < (*a_length - 1 ). For example, if
a_x == domainMax, then *a_startIndex = *a_length - 2, or the next to the last point.
*/

    int64_t lower = 0, mid, upper;
    ptwXYPoint *point;
    *a_length = ptwXY_length( NULL, a_ptwXY );

    if( ptwXY_simpleCoalescePoints( a_smr, a_ptwXY ) != nfu_Okay ) { 
        smr_setReportError2p( a_smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( a_ptwXY->status );
    }

    if( *a_length < 2 ) {
        smr_setReportError2( a_smr, nfu_SMR_libraryID, nfu_tooFewPoints, "number of points = %lld < 2", *a_length );
        return( nfu_tooFewPoints );
    }

    *a_startIndex = -2;
    point = &a_ptwXY->points[lower];
    if( a_x < point->x ) return( nfu_Okay );

    upper = *a_length - 1;
    *a_startIndex = -1;
    point = &a_ptwXY->points[upper];
    if( a_x > point->x ) return( nfu_Okay );

    while( 1 ) {
        mid = ( lower + upper ) >> 1;
        if( mid == lower ) break;
        point = &a_ptwXY->points[mid];
        if( a_x < point->x ) {
            upper = mid; }
        else {
            lower = mid;
        }
    }

    *a_startIndex = mid;
    return( nfu_Okay );
}

/*
************************************************************
*/
nfu_status ptwXY_setXYData( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t length, double const *xy ) {

    int64_t index;
    ptwXYPoint *p;
    double const *d = xy;
    double priorX = 0.;

    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( ptwXY->status );
    }

    if( length > ptwXY->allocatedSize ) {
        if( ptwXY_reallocatePoints( smr, ptwXY, length, 0 ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( ptwXY->status );
        }
    }
    for( index = 0, p = ptwXY->points; index < length; index++, p++ ) {
        if( index != 0 ) {
            if( *d <= priorX ) {
                smr_setReportError2( smr, nfu_SMR_libraryID, nfu_XNotAscending, 
                        "X value at index = %d of %.17e is <= prior value of %.17e", (int) index, *d, priorX );
                ptwXY->status = nfu_XNotAscending;
                length = 0;
                break;
            }
        }
        priorX = *d;
        p->x = *(d++);
        p->y = *(d++);
    }
    ptwXY->overflowHeader.next = &(ptwXY->overflowHeader);
    ptwXY->overflowHeader.prior = &(ptwXY->overflowHeader);
    ptwXY->overflowLength = 0;
    ptwXY->length = length;
    return( ptwXY->status );
}
/*
************************************************************
*/  
nfu_status ptwXY_setXYDataFromXsAndYs( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t length, 
        double const *x, double const *y ) {

    int64_t i;
    ptwXYPoint *p;
    double xOld = 0.;

    if( ptwXY_clear( smr, ptwXY ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( ptwXY->status );
    }

    if( length > ptwXY->allocatedSize ) {
        if( ptwXY_reallocatePoints( smr, ptwXY, length, 0 ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( ptwXY->status );
        }
    }
    for( i = 0, p = ptwXY->points; i < length; i++, p++, x++, y++ ) {
        if( i != 0 ) {
            if( *x <= xOld ) {
                ptwXY->status = nfu_XNotAscending;
                length = 0;
                break;
            }
        }
        xOld = *x;
        p->x = *x;
        p->y = *y;
    }
    ptwXY->length = length;
    return( ptwXY->status );
}
/*
************************************************************
*/  
nfu_status ptwXY_deletePoints( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t i1, int64_t i2 ) {

    int64_t n = ptwXY->length - ( i2 - i1 );

    if( ptwXY_simpleCoalescePoints( smr, ptwXY ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( ptwXY->status );
    }

    if( ( i1 < 0 ) || ( i2 < i1 ) || ( i2 > ptwXY->length ) ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badIndex, "Indices = %d, %d out of bounds: length = %d",
                (int) i1, (int) i2, (int) ptwXY->length );
        return( ptwXY->status = nfu_badIndex );
    }

    if( i1 != i2 ) {
        for( ; i2 < ptwXY->length; i1++, i2++ ) ptwXY->points[i1] = ptwXY->points[i2];
        ptwXY->length = n;
    }
    return( ptwXY->status );
}
/*
************************************************************
*/
nfu_status ptwXY_getLowerIndexBoundingX( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x, int64_t *index ) {

    int64_t i1, length = ptwXY->length;

    *index = -1;

    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( ptwXY->status );
    }

    if( ptwXY_simpleCoalescePoints( smr, ptwXY ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( ptwXY->status );
    }

    if( x < ptwXY->points[0].x ) return( nfu_Okay );
    if( x > ptwXY->points[length-1].x ) return( nfu_Okay );
    for( i1 = 1; i1 < length; ++i1 ) {
        if( x < ptwXY->points[i1].x ) break;
    }
    *index = i1 - 1;
    return( ptwXY->status );
}
/*
************************************************************
*/
ptwXYPoint *ptwXY_getPointAtIndex( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t index ) {

    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( NULL );
    }

    if( ( index < 0 ) || ( index >= ptwXY->length ) ) return( NULL );
    return( ptwXY_getPointAtIndex_Unsafely( ptwXY, index ) );
}
/*
************************************************************
*/
ptwXYPoint *ptwXY_getPointAtIndex_Unsafely( ptwXYPoints const *ptwXY, int64_t index ) {

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
nfu_status ptwXY_getXYPairAtIndex( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t index, double *x, double *y ) {

    ptwXYPoint *p = ptwXY_getPointAtIndex( smr, ptwXY, index );

    if( p == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( ptwXY->status );
    }
    *x = p->x;
    *y = p->y;
    return( nfu_Okay );
}
/*
************************************************************
*/
ptwXY_lessEqualGreaterX ptwXY_getPointsAroundX( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x,  
        ptwXYOverflowPoint *lessThanEqualXPoint, ptwXYOverflowPoint *greaterThanXPoint ) {

    int closeIsEqual;
    ptwXYPoint *closePoint;
    ptwXY_lessEqualGreaterX lessEqualGreaterX;

    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( ptwXY_lessEqualGreaterX_Error );
    }

    lessEqualGreaterX = ptwXY_getPointsAroundX_closeIsEqual( smr, ptwXY, x, lessThanEqualXPoint, greaterThanXPoint, 
            0, &closeIsEqual, &closePoint );
    if( lessEqualGreaterX == ptwXY_lessEqualGreaterX_Error ) smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( lessEqualGreaterX );
}
/*
************************************************************
*/
ptwXY_lessEqualGreaterX ptwXY_getPointsAroundX_closeIsEqual( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x, 
        ptwXYOverflowPoint *lessThanEqualXPoint, ptwXYOverflowPoint *greaterThanXPoint, double eps, int *closeIsEqual, 
        ptwXYPoint **closePoint ) {

    int64_t overflowIndex, nonOverflowLength = ptwXY_getNonOverflowLength( smr, ptwXY );
    int64_t indexMin, indexMid, indexMax;
    ptwXY_dataFrom domainMinFrom, domainMaxFrom;
    double domainMin, domainMax;
    ptwXYOverflowPoint *overflowPoint, *overflowHeader = &(ptwXY->overflowHeader);
    ptwXY_lessEqualGreaterX status = ptwXY_lessEqualGreaterX_empty;
    ptwXYPoint *lowerPoint = NULL, *upperPoint = NULL;

    if( nonOverflowLength < 0 ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( ptwXY_lessEqualGreaterX_Error );
    }

    *closeIsEqual = 0;
    if( ptwXY->length == 0 ) return( status );

    if( ptwXY_domainMinAndFrom( smr, ptwXY, &domainMinFrom, &domainMin ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( ptwXY_lessEqualGreaterX_Error );
    }
        
    if( ptwXY_domainMaxAndFrom( smr, ptwXY, &domainMaxFrom, &domainMax ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( ptwXY_lessEqualGreaterX_Error );
    }

    ptwXY_initialOverflowPoint( lessThanEqualXPoint, overflowHeader, NULL );
    ptwXY_initialOverflowPoint( greaterThanXPoint, overflowHeader, NULL );
    if( x < domainMin ) {
        status = ptwXY_lessEqualGreaterX_lessThan;
        if( domainMinFrom == ptwXY_dataFrom_Points ) {
            greaterThanXPoint->prior = overflowHeader;
            greaterThanXPoint->index = 0;
            greaterThanXPoint->point = ptwXY->points[0];
            *closePoint = &(ptwXY->points[0]); }
        else {
            *greaterThanXPoint = *(overflowHeader->next);
            *closePoint = &(overflowHeader->next->point);
        } }
    else if( x > domainMax ) {
        status = ptwXY_lessEqualGreaterX_greater;
        if( domainMaxFrom == ptwXY_dataFrom_Points ) {
            lessThanEqualXPoint->prior = overflowHeader->prior;
            lessThanEqualXPoint->index = nonOverflowLength - 1;
            lessThanEqualXPoint->point = ptwXY->points[lessThanEqualXPoint->index];
            *closePoint = &(ptwXY->points[lessThanEqualXPoint->index]); }
        else {
            *lessThanEqualXPoint = *(overflowHeader->prior);
            *closePoint = &(overflowHeader->prior->point);
        } }
    else {                                                  /* domainMin <= x <= domainMax */
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
            }
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
                if( ( overflowPoint == overflowHeader ) ||     /* x < domainMin of overflow points. */
                        ( ( ptwXY->points[indexMin].x > overflowPoint->point.x ) && ( ptwXY->points[indexMin].x < x ) ) ) {
                    if( overflowPoint != overflowHeader ) lessThanEqualXPoint->prior = overflowPoint;
                    lowerPoint = &(ptwXY->points[indexMin]);
                    lessThanEqualXPoint->index = indexMin;
                    lessThanEqualXPoint->point = ptwXY->points[indexMin]; }
                else {
                    lowerPoint = &(overflowPoint->point);
                    *lessThanEqualXPoint = *overflowPoint;
                }
                if( ( overflowPoint->next == overflowHeader ) ||     /* x > domainMax of overflow points. */
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

    if( eps > 0 ) {
        double absX = fabs( x );

        if( status == ptwXY_lessEqualGreaterX_lessThan ) {
            if( absX < fabs( greaterThanXPoint->point.x ) ) absX = fabs( greaterThanXPoint->point.x );
            if( ( greaterThanXPoint->point.x - x ) < eps * absX ) *closeIsEqual = 1; }
        else if( status == ptwXY_lessEqualGreaterX_greater ) {
            if( absX < fabs( lessThanEqualXPoint->point.x ) ) absX = fabs( lessThanEqualXPoint->point.x );
            if( ( x - lessThanEqualXPoint->point.x ) < eps * absX ) *closeIsEqual = -1; }
        else if( status == ptwXY_lessEqualGreaterX_between ) {
            if( ( x - lessThanEqualXPoint->point.x ) < ( greaterThanXPoint->point.x - x ) ) {   /* x is closer to lower point. */
                *closePoint = lowerPoint;
                if( absX < fabs( lessThanEqualXPoint->point.x ) ) absX = fabs( lessThanEqualXPoint->point.x );
                if( ( x - lessThanEqualXPoint->point.x ) < eps * absX ) *closeIsEqual = -1; }
            else {                                                                              /* x is closer to upper point. */
                *closePoint = upperPoint;
                if( absX < fabs( greaterThanXPoint->point.x ) ) absX = fabs( greaterThanXPoint->point.x );
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
nfu_status ptwXY_getValueAtX( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x, double *y ) {

    nfu_status status = nfu_XOutsideDomain;
    ptwXYOverflowPoint lessThanEqualXPoint, greaterThanXPoint;
    ptwXY_lessEqualGreaterX legx = ptwXY_getPointsAroundX( smr, ptwXY, x, &lessThanEqualXPoint, &greaterThanXPoint );

    *y = 0.;
    switch( legx ) {
    case ptwXY_lessEqualGreaterX_Error :
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( ptwXY->status );
    case ptwXY_lessEqualGreaterX_empty :
    case ptwXY_lessEqualGreaterX_lessThan :
    case ptwXY_lessEqualGreaterX_greater :
        break;
    case ptwXY_lessEqualGreaterX_equal :
        status = nfu_Okay;
        *y = lessThanEqualXPoint.point.y;
        break;
    case ptwXY_lessEqualGreaterX_between :
        status = ptwXY_interpolatePoint( smr, ptwXY->interpolation, x, y, lessThanEqualXPoint.point.x, lessThanEqualXPoint.point.y, 
                greaterThanXPoint.point.x, greaterThanXPoint.point.y );
        break;
    }
    return( status );
}
/*
************************************************************
*/
nfu_status ptwXY_setValueAtX( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x, double y ) {

    if( ptwXY_setValueAtX_overrideIfClose( smr, ptwXY, x, y, 0., 0 ) != nfu_Okay )
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( ptwXY->status );
}
/*
************************************************************
*/
nfu_status ptwXY_setValueAtX_overrideIfClose( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x, double y, 
        double eps, int override ) {

    int closeIsEqual;
    int64_t nonOverflowLength = ptwXY_getNonOverflowLength( smr, ptwXY ), i;
    ptwXY_lessEqualGreaterX legx;
    ptwXYPoint *point = NULL, newPoint = { x, y };
    ptwXYOverflowPoint *overflowPoint, *p, *overflowHeader = &(ptwXY->overflowHeader);
    ptwXYOverflowPoint lessThanEqualXPoint, greaterThanXPoint;
    ptwXYPoint *closePoint;

    if( nonOverflowLength < 0 ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( nfu_Error );
    }

    legx = ptwXY_getPointsAroundX_closeIsEqual( smr, ptwXY, x, &lessThanEqualXPoint, &greaterThanXPoint, eps, &closeIsEqual, &closePoint );
    switch( legx ) {
    case ptwXY_lessEqualGreaterX_Error :
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( nfu_Error );
    case ptwXY_lessEqualGreaterX_lessThan :
    case ptwXY_lessEqualGreaterX_greater :
    case ptwXY_lessEqualGreaterX_between :
        if( closeIsEqual ) {
            if( !override ) return( nfu_Okay );
            point = closePoint;
            legx = ptwXY_lessEqualGreaterX_equal;
            x = point->x; }
        else {
            if( ( legx == ptwXY_lessEqualGreaterX_greater ) && ( nonOverflowLength < ptwXY->allocatedSize ) ) {
                point = &(ptwXY->points[nonOverflowLength]); }
            else {
                if( ptwXY->overflowLength == ptwXY->overflowAllocatedSize ) {
                    if( ptwXY_coalescePoints( smr, ptwXY, ptwXY->length + ptwXY->overflowAllocatedSize, &newPoint, 0 ) != nfu_Okay )
                        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
                    return( ptwXY->status );
                }
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
        if( closeIsEqual && !override ) return( nfu_Okay );
        if( lessThanEqualXPoint.next == NULL ) {
            point = &(ptwXY->points[lessThanEqualXPoint.index]); }
        else {
            point = &(lessThanEqualXPoint.prior->next->point);
        }
        break;
    }

    point->x = x;
    point->y = y;
    if( legx != ptwXY_lessEqualGreaterX_equal ) ptwXY->length++;
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwXY_mergeFromXsAndYs( statusMessageReporting *smr, ptwXYPoints *ptwXY, int length, double *xs, double *ys ) {

    if( ptwXY_mergeFrom( smr, ptwXY, 1, length, xs, ys ) != nfu_Okay )
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( ptwXY->status );
}
/*
************************************************************
*/
nfu_status ptwXY_mergeFromXYs( statusMessageReporting *smr, ptwXYPoints *ptwXY, int length, double *xys ) {

    int i;
    double *xs, *p1, *p2;

    if( length == 0 ) return( nfu_Okay );
    if( length < 0 ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badInput, "Negative length = %d.", length );
        return( nfu_badInput );
    }

    if( ( xs = (double *) smr_malloc2( smr, length * sizeof( double ), 0, "xs" ) ) == NULL ) return( nfu_mallocError );
    for( i = 0, p1 = xs, p2 = xys; i < length; i++, p1++, p2 += 2 ) *p1 = *p2;
    if( ptwXY_mergeFrom( smr, ptwXY, 2, length, xs, xys ) != nfu_Okay )
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    smr_freeMemory2( xs );

    return( ptwXY->status );
}
/*
************************************************************
*/
static nfu_status ptwXY_mergeFrom( statusMessageReporting *smr, ptwXYPoints *ptwXY, int incY, int length, double *xs, double *ys ) {

    int i1, j1,  n1 = 0;
    double *p1, priorX;
    ptwXYPoint *point1, *point2;

    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( ptwXY->status );
    }

    if( length == 0 ) return( nfu_Okay );
    if( length < 0 ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badInput, "Negative length = %d.", length );
        return( nfu_badInput );
    }

    if( ptwXY_simpleCoalescePoints( smr, ptwXY ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( ptwXY->status );
    }

    if( xs[0] < 0 ) {
        priorX = 1.1 * xs[0]; }
    else {
        priorX = 0.9 * xs[0] - 1;
    }
    for( i1 = 0, p1 = xs; i1 < length; ++i1, ++p1 ) {
        if( *p1 <= priorX ) return( nfu_XNotAscending );
        priorX = *p1;
    }

    for( i1 = 0, p1 = xs, j1 = 0; i1 < length; ++i1, ++p1 ) { /* Count the number of x-values same in xs and ptwXY. */
        for( ; j1 < ptwXY->length; ++j1 ) {
            if( *p1 <= ptwXY->points[j1].x ) break;
        }
        if( j1 == ptwXY->length ) break;             /* Completed all ptwXY points. */
        if( *p1 == ptwXY->points[j1].x ) ++n1;
    }
    n1 = length + (int) ptwXY->length - n1;

    if( ptwXY_reallocatePoints( smr, ptwXY, n1, 0 ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        n1 = 0; }
    else {
        point1 = &(ptwXY->points[n1-1]);            /* Go backwards through arrays. */
        point2 = &(ptwXY->points[ptwXY->length-1]);
        p1 = &(xs[length-1]);
        for( i1 = length - 1, j1 = (int) ptwXY->length - 1; ( i1 >= 0 ) && ( j1 >= 0 ); --point1 ) {
            if( *p1 >= point2->x ) {
                point1->x = *p1;
                point1->y = ys[i1];
                if( *p1 == point2->x ) {
                    --point2;
                    --j1;
                }
                --p1;
                --i1; }
            else {
                *point1 = *point2;
                --point2;
                --j1;
            }
        }
        for( ; i1 >= 0; --i1, --p1, --point1 ) {
            point1->x = *p1;
            point1->y = ys[i1];
        }
    }
    ptwXY->length = n1;

    return( ptwXY->status );
}
/*
************************************************************
*/
nfu_status ptwXY_appendXY( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x, double y ) {

    int64_t nonOverflowLength = ptwXY_getNonOverflowLength( smr, ptwXY );
    ptwXY_dataFrom dataFrom;

    if( nonOverflowLength < 0 ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( nfu_Error );
    }

    if( ptwXY->length != 0 ) {
        double domainMax;
        nfu_status status;

        if( ( status = ptwXY_domainMaxAndFrom( smr, ptwXY, &dataFrom, &domainMax ) ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( status );
        }
        if( domainMax >= x ) {
            smr_setReportError2( smr, nfu_SMR_libraryID, nfu_XNotAscending, "domainMax = %17.e >= x = %.17e", domainMax, x );
            return( ptwXY->status = nfu_XNotAscending );
        }
    }

    if( nonOverflowLength < ptwXY->allocatedSize ) {      /* Room at end of points. Also handles the case when length = 0. */
        ptwXY->points[nonOverflowLength].x = x;
        ptwXY->points[nonOverflowLength].y = y; }
    else {
        if( ptwXY->overflowLength == ptwXY->overflowAllocatedSize ) {
            ptwXYPoint newPoint = { x, y };
            if( ptwXY_coalescePoints( smr, ptwXY, ptwXY->length + ptwXY->overflowAllocatedSize, &newPoint, 0 ) != nfu_Okay ) {
                smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
                return( ptwXY->status ); } }
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
nfu_status ptwXY_setXYPairAtIndex( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t index, double x, double y ) {

    int64_t i, ip1;
    ptwXYOverflowPoint *overflowPoint, *pm1, *pp1;

    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( ptwXY->status );
    }

    if( ( index < 0 ) || ( index >= ptwXY->length ) ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badIndex, "Index = %d, out of bounds: length = %d",
                (int) index, (int) ptwXY->length );
        return( ptwXY->status = nfu_badIndex );
    }

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
nfu_status ptwXY_getSlopeAtX( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x, const char side, double *slope ) {

    nfu_status status  = nfu_Okay;
    ptwXYOverflowPoint lessThanEqualXPoint = { NULL, NULL, 0, {0.0, 0.0}}, greaterThanXPoint;
    ptwXY_lessEqualGreaterX legx = ptwXY_getPointsAroundX( smr, ptwXY, x, &lessThanEqualXPoint, &greaterThanXPoint );
    ptwXYPoint *point;
    greaterThanXPoint = lessThanEqualXPoint;        /* Done to stop a compiler from complaining. */

    *slope = 0.;

    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( ptwXY->status );
    }
    if( ( side != '-' ) && ( side != '+' ) ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid side = '%c'.", side );
        return( nfu_badInput );
    }

    switch( legx ) {
    case ptwXY_lessEqualGreaterX_Error :
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( nfu_Error );
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
nfu_status ptwXY_domainMinAndFrom( statusMessageReporting *smr, ptwXYPoints *ptwXY, ptwXY_dataFrom *dataFrom, double *domainMin ) {

    int64_t nonOverflowLength = ptwXY_getNonOverflowLength( smr, ptwXY );

    *domainMin  = 0;
    *dataFrom = ptwXY_dataFrom_Unknown;

    if( nonOverflowLength < 0 ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( nfu_Error );
    }

    if( ptwXY->overflowLength > 0 ) {
        *dataFrom = ptwXY_dataFrom_Overflow;
        *domainMin = ptwXY->overflowHeader.next->point.x;
        if( nonOverflowLength >= 0 ) {
            if( *domainMin > ptwXY->points[0].x ) {
                *dataFrom = ptwXY_dataFrom_Points;
                *domainMin = ptwXY->points[0].x;
            }
        } }
    else if( nonOverflowLength > 0 ) {
        *dataFrom = ptwXY_dataFrom_Points;
        *domainMin = ptwXY->points[0].x; }
    else {
        return( nfu_empty );
    }
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwXY_domainMin( statusMessageReporting *smr, ptwXYPoints *ptwXY, double *domainMin ) {

    ptwXY_dataFrom dataFrom;
    nfu_status status;
    
    if( ( status = ptwXY_domainMinAndFrom( smr, ptwXY, &dataFrom, domainMin ) ) != nfu_Okay ) {
        if( status == nfu_empty ) return( status );
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( nfu_Error );
    }

    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwXY_domainMaxAndFrom( statusMessageReporting *smr, ptwXYPoints *ptwXY, ptwXY_dataFrom *dataFrom, double *domainMax ) {

    int64_t nonOverflowLength = ptwXY_getNonOverflowLength( smr, ptwXY );

    *domainMax  = 0;
    *dataFrom = ptwXY_dataFrom_Unknown;

    if( nonOverflowLength < 0 ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( nfu_Error );
    }

    if( ptwXY->overflowLength > 0 ) {
        *dataFrom = ptwXY_dataFrom_Overflow;
        *domainMax = ptwXY->overflowHeader.prior->point.x;
        if( ( nonOverflowLength > 0 ) ) {
            if( *domainMax < ptwXY->points[nonOverflowLength-1].x ) {
                *dataFrom = ptwXY_dataFrom_Points;
                *domainMax = ptwXY->points[nonOverflowLength-1].x;
            }
        } }
    else if( ptwXY->length > 0 ) {
        *dataFrom = ptwXY_dataFrom_Points;
        *domainMax = ptwXY->points[nonOverflowLength-1].x; }
    else {
        return( nfu_empty );
    }
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwXY_domainMax( statusMessageReporting *smr, ptwXYPoints *ptwXY, double *domainMax ) {

    ptwXY_dataFrom dataFrom;
    nfu_status status;

    if( ( status = ptwXY_domainMaxAndFrom( smr, ptwXY, &dataFrom, domainMax ) ) != nfu_Okay ) {
        if( status == nfu_empty ) return( status );
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( nfu_Error );
    }
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwXY_range( statusMessageReporting *smr, ptwXYPoints *ptwXY, double *rangeMin, double *rangeMax ) {

    int64_t i, nonOverflowLength = ptwXY_getNonOverflowLength( smr, ptwXY  );
    ptwXYPoint *p = ptwXY->points;
    ptwXYOverflowPoint *overflowPoint = ptwXY->overflowHeader.next;

    *rangeMin = *rangeMax = 0.;

    if( nonOverflowLength < 0 ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( nfu_Error );
    }

    if( ptwXY->length == 0 ) return( nfu_empty );
    if( nonOverflowLength > 0 ) {
        *rangeMin = *rangeMax = p->y;
        for( i = 1, p++; i < nonOverflowLength; i++, p++ ) {
            *rangeMin = ( ( *rangeMin < p->y ) ? *rangeMin : p->y );
            *rangeMax = ( ( *rangeMax > p->y ) ? *rangeMax : p->y );
        } }
    else {
        *rangeMin = *rangeMax = overflowPoint->point.y;
    }
    for( ; overflowPoint != &(ptwXY->overflowHeader); overflowPoint = overflowPoint->next ) {
        *rangeMin = ( ( *rangeMin < overflowPoint->point.y ) ? *rangeMin : overflowPoint->point.y );
        *rangeMax = ( ( *rangeMax < overflowPoint->point.y ) ? *rangeMax : overflowPoint->point.y );
    }
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwXY_rangeMin( statusMessageReporting *smr, ptwXYPoints *ptwXY, double *rangeMin ) {

    int64_t i, nonOverflowLength = ptwXY_getNonOverflowLength( smr, ptwXY  );
    ptwXYPoint *p = ptwXY->points;
    ptwXYOverflowPoint *overflowPoint = ptwXY->overflowHeader.next;

    *rangeMin = 0.;

    if( nonOverflowLength < 0 ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( nfu_Error );
    }

    if( ptwXY->length == 0 ) return( nfu_empty );
    if( nonOverflowLength > 0 ) {
        *rangeMin = p->y;
        for( i = 1, p++; i < nonOverflowLength; i++, p++ ) *rangeMin = ( ( *rangeMin < p->y ) ? *rangeMin : p->y ); }
    else {
        *rangeMin = overflowPoint->point.y;
    }
    for( ; overflowPoint != &(ptwXY->overflowHeader); overflowPoint = overflowPoint->next ) 
        *rangeMin = ( ( *rangeMin < overflowPoint->point.y ) ? *rangeMin : overflowPoint->point.y );
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwXY_rangeMax( statusMessageReporting *smr, ptwXYPoints *ptwXY, double *rangeMax ) {

    int64_t i, nonOverflowLength = ptwXY_getNonOverflowLength( smr, ptwXY  );
    ptwXYPoint *p = ptwXY->points;
    ptwXYOverflowPoint *overflowPoint = ptwXY->overflowHeader.next;

    *rangeMax = 0.;

    if( nonOverflowLength < 0 ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( nfu_Error );
    }

    if( ptwXY->length == 0 ) return( nfu_empty );
    if( nonOverflowLength > 0 ) {
        *rangeMax = p->y;
        for( i = 1, p++; i < nonOverflowLength; i++, p++ ) *rangeMax = ( ( *rangeMax > p->y ) ? *rangeMax : p->y ); }
    else {
        *rangeMax = overflowPoint->point.y;
    }
    for( ; overflowPoint != &(ptwXY->overflowHeader); overflowPoint = overflowPoint->next )
        *rangeMax = ( ( *rangeMax > overflowPoint->point.y ) ? *rangeMax : overflowPoint->point.y );
    return( nfu_Okay );
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
/*
************************************************************
*/
char const *ptwXY_interpolationToString( ptwXY_interpolation interpolation ) {

    switch( interpolation ) {
    case ptwXY_interpolationLinLin : return( linLinInterpolationString );
    case ptwXY_interpolationLogLin : return( logLinInterpolationString );
    case ptwXY_interpolationLinLog : return( linLogInterpolationString );
    case ptwXY_interpolationLogLog : return( logLogInterpolationString );
    case ptwXY_interpolationFlat :   return( flatInterpolationString );
    default :
        break;
    }
    return( NULL );
}
/*
************************************************************
*/
ptwXY_interpolation ptwXY_stringToInterpolation( char const *interpolationString ) {

    if( strcmp( interpolationString, "" ) == 0 ) return( ptwXY_interpolationLinLin );
    if( strcmp( interpolationString, linLinInterpolationString ) == 0 ) return( ptwXY_interpolationLinLin );
    if( strcmp( interpolationString, logLinInterpolationString ) == 0 ) return( ptwXY_interpolationLogLin );
    if( strcmp( interpolationString, linLogInterpolationString ) == 0 ) return( ptwXY_interpolationLinLog );
    if( strcmp( interpolationString, logLogInterpolationString ) == 0 ) return( ptwXY_interpolationLogLog );
    if( strcmp( interpolationString, flatInterpolationString ) == 0 ) return( ptwXY_interpolationFlat );
    return( ptwXY_interpolationOther );
}
