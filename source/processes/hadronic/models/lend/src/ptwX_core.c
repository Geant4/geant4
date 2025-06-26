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

#include "ptwX.h"

static int ptwX_sort_descending( void const *p1, void const *p2 );
static int ptwX_sort_ascending( void const *p1, void const *p2 );
/*
************************************************************
*/
ptwXPoints *ptwX_new( statusMessageReporting *smr, int64_t size ) {

    ptwXPoints *ptwX = (ptwXPoints *) smr_malloc2( smr, sizeof( ptwXPoints ), 1, "ptwX" );

    if( ptwX == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }

    if( ptwX_initialize( smr, ptwX, size ) != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        smr_freeMemory2( ptwX );
    }
    return( ptwX );
}
/*
************************************************************
*/
nfu_status ptwX_initialize( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t size ) {

    ptwX->status = nfu_Okay;
    ptwX->length = 0;
    ptwX->allocatedSize = 0;
    ptwX->mallocFailedSize = 0;
    ptwX->points = NULL;
    if( ptwX_reallocatePoints( smr, ptwX, size, 0 ) != nfu_Okay )
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( ptwX->status );
}
/*
************************************************************
*/
ptwXPoints *ptwX_create( statusMessageReporting *smr, int64_t size, int64_t length, double const *xs ) {

    ptwXPoints *ptwX = ptwX_new( smr, size );

    if( ptwX == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." ); }
    else {
        if( ptwX_setData( smr, ptwX, length, xs ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            ptwX = ptwX_free( ptwX );
        }
    }
    return( ptwX );
}
/*
************************************************************
*/
ptwXPoints *ptwX_createLine( statusMessageReporting *smr, int64_t size, int64_t length, double slope, double offset ) {

    int64_t i1;
    double *p1;
    ptwXPoints *ptwX;

    if( size < length ) size = length;
    if( ( ptwX = ptwX_new( smr, size ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." ); }
    else {
        for( i1 = 0, p1 = ptwX->points; i1 < length; i1++, p1++ ) *p1 = slope * i1 + offset;
        ptwX->length = length;
    }
    return( ptwX );
}
/*
************************************************************
*/
nfu_status ptwX_copy( statusMessageReporting *smr, ptwXPoints *dest, ptwXPoints *src ) {

    if( dest->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid destination." );
        return( nfu_badSelf );
    }
    if( src->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    if( ptwX_clear( smr, dest ) ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( src->status );
    }
    if( ptwX_setData( smr, dest, src->length, src->points ) ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( src->status );
    }
    return( nfu_Okay );
}
/*
************************************************************
*/
ptwXPoints *ptwX_clone( statusMessageReporting *smr, ptwXPoints *ptwX ) {

    ptwXPoints *clone = ptwX_slice( smr, ptwX, 0, ptwX->length );

    if( clone == NULL ) smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( clone );
}
/*
************************************************************
*/
ptwXPoints *ptwX_slice( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t index1, int64_t index2 ) {

    int64_t i1, i2, length;
    ptwXPoints *n1;

    if( ptwX->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( NULL );
    }

    if( index1 < 0 ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badIndex, "negative index1 = %d.", (int) index1 );
        return( NULL );
    }
    if( index2 < index1 ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badIndex, "index1 = %d greater than index2 = %d", 
                (int) index1, (int) index2 );
        return( NULL );
    }
    if( index2 > ptwX->length ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badIndex, "index2 = %d greater than length = %d.", 
                (int) index2, (int) ptwX->length );
        return( NULL );
    }

    length = ( index2 - index1 );
    if( ( n1 = ptwX_new( smr, length ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }
    for( i2 = 0, i1 = index1; i1 < index2; i1++, i2++ ) n1->points[i2] = ptwX->points[i1];
    n1->length = length;

    return( n1 );
}
/*
************************************************************
*/
nfu_status ptwX_reallocatePoints( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t size, int forceSmallerResize ) {

    if( ptwX->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    if( size < ptwX_minimumSize ) size = ptwX_minimumSize;                      /* ptwX_minimumSize must be > 0 for other routines to work properly. */
    if( size < ptwX->length ) size = ptwX->length;
    if( size != ptwX->allocatedSize ) {
        if( size > ptwX->allocatedSize ) {                                      /* Increase size of allocated points. */
            ptwX->points = (double *) smr_realloc2( smr, ptwX->points, (size_t) size * sizeof( double ), "ptwX->points" );
            if( ptwX->points == NULL ) smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." ); }
        else if( ( ptwX->allocatedSize > 2 * size ) || forceSmallerResize ) {   /* Decrease size, if at least 1/2 size reduction or if forced to. */
            ptwX->points = (double *) smr_realloc2( smr, ptwX->points, (size_t) size * sizeof( double ), "ptwX->points" );
            if( ptwX->points == NULL ) smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        }
        if( ptwX->points == NULL ) {
            ptwX->mallocFailedSize = size;
            size = 0;
            ptwX->status = nfu_mallocError;
        }
        ptwX->allocatedSize = size;
    }

    return( ptwX->status );
}
/*
************************************************************
*/
nfu_status ptwX_clear( statusMessageReporting *smr, ptwXPoints *ptwX ) {

    ptwX->length = 0;
    ptwX->status = nfu_Okay;

    return( ptwX->status );
}
/*
************************************************************
*/
nfu_status ptwX_release( statusMessageReporting *smr, ptwXPoints *ptwX ) {

    ptwX->status = nfu_Okay;
    ptwX->length = 0;
    ptwX->allocatedSize = 0;
    smr_freeMemory2( ptwX->points );

    return( nfu_Okay );
}
/*
************************************************************
*/
ptwXPoints *ptwX_free( ptwXPoints *ptwX ) {

    if( ptwX != NULL ) ptwX_release( NULL, ptwX );
    smr_freeMemory2( ptwX );
    return( ptwX );
}
/*
************************************************************
*/
int64_t ptwX_length( statusMessageReporting *smr, ptwXPoints *ptwX ) {

    if( ptwX->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( -nfu_badSelf );
    }

    return( ptwX->length );
}
/*
************************************************************
*/
nfu_status ptwX_setData( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t length, double const *xs ) {

    int64_t  i;

    if( ptwX->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    if( length > ptwX->allocatedSize ) {
        if( ptwX_reallocatePoints( smr, ptwX, length, 0 ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( ptwX->status );
        }
    }
    for( i = 0; i < length; i++ ) ptwX->points[i] = xs[i];
    ptwX->length = length;

    return( ptwX->status );
}
/*
************************************************************
*/ 
nfu_status ptwX_deletePoints( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t i1, int64_t i2 ) {

    if( ptwX->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    if( ( i1 < 0 ) || ( i1 > i2 ) || ( i2 > ptwX->length ) ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badIndex, "index1 = %d, index2 = %d and length = %d", 
                (int) i1, (int) i2, (int) ptwX->length );
        return( nfu_badIndex );
    }
    if( i1 != i2 ) {
        int64_t n1 = ptwX->length - ( i2 - i1 );

        for( ; i2 < ptwX->length; i1++, i2++ ) ptwX->points[i1] = ptwX->points[i2];
        ptwX->length = n1;
    }
    return( ptwX->status );
}
/*
************************************************************
*/
double *ptwX_getPointAtIndex( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t index ) {

    if( ptwX->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( NULL );
    }

    if( ( index < 0 ) || ( index >= ptwX->length ) ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badIndex, "Index = %d out of bounds: length = %d", 
                (int) index, (int) ptwX->length );
        return( NULL );
    }
    return( &(ptwX->points[index]) );
}
/*
************************************************************
*/
double ptwX_getPointAtIndex_Unsafely( ptwXPoints *ptwX, int64_t index ) {

    return( ptwX->points[index] );
}
/*
************************************************************
*/
nfu_status ptwX_setPointAtIndex( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t index, double x ) {

    if( ptwX->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    if( ( index < 0 ) || ( index > ptwX->length ) ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badIndex, "Index = %d out of bounds: length = %d", 
                (int) index, (int) ptwX->length );
        return( nfu_badIndex );
    }

    if( index == ptwX->allocatedSize ) {
        if( ptwX_reallocatePoints( smr, ptwX, ptwX->allocatedSize + 10, 0 ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( ptwX->status );
        }
    }
    ptwX->points[index] = x;
    if( index == ptwX->length ) ptwX->length++;
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwX_insertPointsAtIndex( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t index, int64_t n1, double const *xs ) {

    int64_t i1, i2, n1p, size = n1 + ptwX->length;

    if( ptwX->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    if( n1 < 1 ) return( nfu_Okay );        /* No points to insert. */

    if( ( index < 0 ) || ( index > ptwX->length ) ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badIndex, "Index = %d out of bounds: length = %d", 
                (int) index, (int) ptwX->length );
        return( ptwX->status = nfu_Error );
    }

    if( size > ptwX->allocatedSize ) {
        if( ptwX_reallocatePoints( smr, ptwX, size, 0 ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( ptwX->status );
        }
    }

    for( i1 = ptwX->length - 1, i2 = size - 1, n1p = ptwX->length - index; n1p > 0; i1--, i2--, n1p-- )
        ptwX->points[i2] = ptwX->points[i1];
    for( i1 = 0, i2 = index; i1 < n1; i1++, i2++ ) ptwX->points[i2] = xs[i1];
    ptwX->length += n1;

    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwX_ascendingOrder( statusMessageReporting *smr, ptwXPoints *ptwX, int *order ) {
/*
*    Returns -1 list is descending, 1 if ascending and 0 otherwise (i.e., mixed).
*/
    int64_t i1;
    double x1, x2;

    if( ptwX->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    *order = 1;
    if( ptwX->length < 2 ) return( nfu_Okay );

    if( ( x1 = ptwX->points[0] ) < ( x2 = ptwX->points[1] ) ) {     /* Check for ascending order. */
        for( i1 = 2; i1 < ptwX->length; i1++ ) {
            x1 = x2;
            x2 = ptwX->points[i1];
            if( x2 <= x1 ) {
                *order = 0;
                return( nfu_Okay );
            }
        } }
    else {
        *order = -1;                                                 /* Check for descending order. */
        for( i1 = 1; i1 < ptwX->length; i1++ ) {
            x2 = ptwX->points[i1];
            if( x1 <= x2 ) {
                *order = 0;
                return( nfu_Okay );
            }
            x1 = x2;
        }
    }
    return( nfu_Okay );
}
/*
************************************************************
*/
ptwXPoints *ptwX_fromString( statusMessageReporting *smr, char const *str, char sep, char **endCharacter ) {

    int64_t numberConverted;
    double  *doublePtr;
    ptwXPoints *ptwX = NULL;

    if( ( doublePtr = nfu_stringToListOfDoubles( smr, str, sep, &numberConverted, endCharacter, 1 ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }
    if( ( ptwX = ptwX_create( smr, numberConverted, numberConverted, doublePtr ) ) == NULL )
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    smr_freeMemory2( doublePtr );
    return( ptwX );
}
/*
************************************************************
*/
int ptwX_countOccurrences( statusMessageReporting *smr, ptwXPoints *ptwX, double value ) {

    int count;
    int64_t i1;

    if( ptwX->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( -nfu_badSelf );
    }

    count = 0;
    for( i1 = 0; i1 < ptwX->length; i1++ ) {
        if( ptwX->points[i1] == value ) ++count;
    }
    return( count );
}
/*
************************************************************
*/
nfu_status ptwX_reverse( statusMessageReporting *smr, ptwXPoints *ptwX ) {

    int64_t i1, i2 = ptwX->length - 1, n1 = ptwX->length / 2;
    double tmp;

    if( ptwX->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    for( i1 = 0; i1 < n1; i1++, i2-- ) {
        tmp = ptwX->points[i1];
        ptwX->points[i1] = ptwX->points[i2];
        ptwX->points[i2] = tmp;
    }
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwX_sort( statusMessageReporting *smr, ptwXPoints *ptwX, enum ptwX_sort_order order ) {

    int (*cmp)( void const *, void const * ) = ptwX_sort_descending;

    if( ptwX->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    if( order == ptwX_sort_order_ascending ) cmp = ptwX_sort_ascending;
    qsort( ptwX->points, (size_t) ptwX->length, sizeof( ptwX->points[0] ), cmp );
    return( nfu_Okay );
}
/*
************************************************************
*/
static int ptwX_sort_descending( void const *p1, void const *p2 ) { return( -ptwX_sort_ascending( p1, p2 ) ); }
static int ptwX_sort_ascending( void const *p1, void const *p2 ) {

    double *d1 = (double *) p1, *d2 = (double *) p2;

    if( *d1 < *d2 ) return( -1 );
    if( *d1 == *d2 ) return( 0 );
    return( 1 );
}
/*
************************************************************
*/
nfu_status ptwX_closesDifference( statusMessageReporting *smr, ptwXPoints *ptwX, double value, int64_t *index, double *difference ) {

    return( ptwX_closesDifferenceInRange( smr, ptwX, 0, ptwX->length, value, index, difference ) );
}
/*
************************************************************
*/
nfu_status ptwX_closesDifferenceInRange( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t i1, int64_t i2, 
        double value, int64_t *index, double *difference ) {
/*
*   Finds the closes datum to value. If *difference is zero, datum is same as value.
*/
    double d1;

    if( ptwX->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    if( i1 < 0 ) i1 = 0;
    if( i2 > ptwX->length ) i2 = ptwX->length;
    if( i1 >= i2 ) return( nfu_Okay );
    *index = i1;
    *difference = value - ptwX->points[i1];
    for( i1++; i1 < i2; i1++ ) {
        d1 = value - ptwX->points[i1];
        if( fabs( *difference ) > fabs( d1 ) ) {
            *index = i1;
            *difference = d1;
        }
    }
    return( nfu_Okay );
}
/*
************************************************************
*/
ptwXPoints *ptwX_unique( statusMessageReporting *smr, ptwXPoints *ptwX, int order ) {
/*
*   Returns a new ptwXPoints instance that is a unique list of the values in ptwX.
*   If order < 0 order is descending, if order > 0 order is ascending, otherwise, order is the same as ptwX.
*/
    int64_t i1, i2, n1 = 0;
    double x1, *p2;
    ptwXPoints *ptwX2 = NULL;

    if( ptwX->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( NULL );
    }

    if( order == 0 ) {
        if( ( ptwX2 = ptwX_new( smr, ptwX->length ) ) == NULL ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( NULL );
        }
        for( i1 = 0; i1 < ptwX->length; i1++ ) {
            x1 = ptwX->points[i1];
            for( i2 = 0, p2 = ptwX2->points; i2 < ptwX2->length; i2++, p2++ ) {
                if( *p2 == x1 ) break;
            }
            if( i2 == ptwX2->length ) {
                ptwX2->points[ptwX2->length] = x1;
                ptwX2->length++;
            }
        } }
    else {
        enum ptwX_sort_order sort_order = ( order > 0 ) ? ptwX_sort_order_ascending : ptwX_sort_order_descending;

        if( ( ptwX2 = ptwX_clone( smr, ptwX ) ) == NULL ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( NULL );
        }
        if( ptwX_sort( smr, ptwX2, sort_order ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            goto err;
        }

        if( ptwX2->length > 1 ) {
            x1 = ptwX2->points[n1];                 /* n1 is initially 0. */
            n1++;
            for( i1 = 1; i1 < ptwX2->length; i1++ ) {
                if( x1 != ptwX2->points[i1] ) {
                    x1 = ptwX2->points[i1];
                    ptwX2->points[n1] = x1;
                    n1++;
                }
            }
            ptwX2->length = n1;
        }
    }
    return( ptwX2 );

err:
    if( ptwX2 != NULL ) ptwX_free( ptwX2 );
    return( NULL );
}
/*
************************************************************
*/
nfu_status ptwX_abs( statusMessageReporting *smr, ptwXPoints *ptwX ) {

    int64_t i1;
    double *p1;

    if( ptwX->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    for( i1 = 0, p1 = ptwX->points; i1 < ptwX->length; i1++, p1++ ) *p1 = fabs( *p1 );
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwX_neg( statusMessageReporting *smr, ptwXPoints *ptwX ) {

    if( ptwX_slopeOffset( smr, ptwX, -1, 0 ) != nfu_Okay )
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( ptwX->status );
}
/*
************************************************************
*/
nfu_status ptwX_add_double( statusMessageReporting *smr, ptwXPoints *ptwX, double value ) {

    if( ptwX_slopeOffset( smr, ptwX, 1, value ) != nfu_Okay )
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( ptwX->status );
}
/*
************************************************************
*/
nfu_status ptwX_mul_double( statusMessageReporting *smr, ptwXPoints *ptwX, double value ) {

    if( ptwX_slopeOffset( smr, ptwX, value, 0 ) != nfu_Okay )
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( ptwX->status );
}
/*
************************************************************
*/
nfu_status ptwX_slopeOffset( statusMessageReporting *smr, ptwXPoints *ptwX, double slope, double offset ) {

    int64_t i1;
    double *p1;

    if( ptwX->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    for( i1 = 0, p1 = ptwX->points; i1 < ptwX->length; i1++, p1++ ) *p1 = slope * *p1 + offset;
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwX_add_ptwX( statusMessageReporting *smr, ptwXPoints *ptwX1, ptwXPoints *ptwX2 ) {

    int64_t i1;
    double *p1 = ptwX1->points, *p2 = ptwX2->points;

    if( ptwX1->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source1." );
        return( nfu_badSelf );
    }
    if( ptwX2->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source2." );
        return( nfu_badSelf );
    }

    if( ptwX1->length != ptwX2->length ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_domainsNotMutual, 
                "length of source1 = %d not the same as length of source2 = %d.", (int) ptwX1->length, (int) ptwX2->length );
        return( nfu_domainsNotMutual );
    }

    for( i1 = 0; i1 < ptwX1->length; i1++, p1++, p2++ ) *p1 += *p2;
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwX_sub_ptwX( statusMessageReporting *smr, ptwXPoints *ptwX1, ptwXPoints *ptwX2 ) {

    int64_t i1;
    double *p1 = ptwX1->points, *p2 = ptwX2->points;

    if( ptwX1->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source1" );
        return( nfu_badSelf );
    }
    if( ptwX2->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source2." );
        return( nfu_badSelf );
    }

    if( ptwX1->length != ptwX2->length ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_domainsNotMutual, 
                "length of source1 = %d not the same as length of source2 = %d.", (int) ptwX1->length, (int) ptwX2->length );
        return( nfu_domainsNotMutual );
    }

    for( i1 = 0; i1 < ptwX1->length; i1++, p1++, p2++ ) *p1 -= *p2;
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwX_range( statusMessageReporting *smr, ptwXPoints *ptwX, double *rangeMin, double *rangeMax ) {

    int64_t i1, n1 = ptwX->length;
    *rangeMin = *rangeMax = 0;
    double *p1 = ptwX->points;

    if( ptwX->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    if( n1 > 0 ) {
        *rangeMin = *rangeMax = *(p1++);
        for( i1 = 1; i1 < n1; ++i1, ++p1 ) {
            if( *p1 < *rangeMin ) *rangeMin = *p1;
            if( *p1 > *rangeMax ) *rangeMax = *p1;
        }
    }
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwX_compare( statusMessageReporting *smr, ptwXPoints *ptwX1, ptwXPoints *ptwX2, int *comparison ) {

    int64_t i1, n1 = ptwX1->length, n2 = ptwX2->length, nn = n1;
    double *p1 = ptwX1->points, *p2 = ptwX2->points;

    *comparison = 0;
    if( ptwX1->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source1." );
        return( nfu_badSelf );
    }
    if( ptwX2->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source2." );
        return( nfu_badSelf );
    }

    if( nn > n2 ) nn = n2;
    for( i1 = 0; i1 < nn; i1++, p1++, p2++ ) {
        if( *p1 == *p2 ) continue;
        *comparison = 1;
        if( *p1 < *p2 ) *comparison = -1;
        return( nfu_Okay );
    }
    if( n1 < n2 ) {
        *comparison = -1; }
    else if( n1 > n2 ) {
        *comparison = 1;
    }
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwX_close( statusMessageReporting *smr, ptwXPoints *ptwX1, ptwXPoints *ptwX2, int epsilonFactor, double epsilon,
        int *index ) {
/*
*   Returns the index where ptwX1 and ptwX2 differ significantly as determined by epsilonFactor and epsilon.
*/

    int64_t i1, n1 = ptwX1->length;
    double larger;
    double *p1 = ptwX1->points, *p2 = ptwX2->points;

    *index = -1;
    if( ptwX1->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source1." );
        return( nfu_badSelf );
    }
    if( ptwX2->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source2." );
        return( nfu_badSelf );
    }

    if( ptwX1->length != ptwX2->length ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_domainsNotMutual, 
                "length of source1 = %d not the same as length of source2 = %d.", (int) ptwX1->length, (int) ptwX2->length );
        return( nfu_domainsNotMutual );
    }

    epsilon = fabs( epsilon ) + abs( epsilonFactor ) * DBL_EPSILON;

    for( i1 = 0; i1 < n1; i1++, p1++, p2++ ) {
        larger = fabs( *p1 );
        if( fabs( *p2 ) > larger ) larger = fabs( *p2 );
        if( fabs( *p2 - *p1 ) > epsilon * larger ) break;
    }
    *index = (int) i1;
    return( nfu_Okay );
}
