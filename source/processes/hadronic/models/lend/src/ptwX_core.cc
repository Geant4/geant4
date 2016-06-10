/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <cstdlib>   //for std::abs(int)
#include <float.h>

#include "ptwX.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

static int ptwX_sort_descending( void const *p1, void const *p2 );
static int ptwX_sort_ascending( void const *p1, void const *p2 );
/*
************************************************************
*/
ptwXPoints *ptwX_new( int64_t size, nfu_status *status ) {

    ptwXPoints *ptwX = (ptwXPoints *) nfu_calloc( sizeof( ptwXPoints ), 1 );

    *status = nfu_mallocError;
    if( ptwX == NULL ) return( NULL );
    ptwX_setup( ptwX, size );
    if( ( *status = ptwX->status ) != nfu_Okay ) ptwX = (ptwXPoints *) nfu_free( ptwX );
    return( ptwX );
}
/*
************************************************************
*/
nfu_status ptwX_setup( ptwXPoints *ptwX, int64_t size ) {

    ptwX->status = nfu_Okay;
    ptwX->length = 0;
    ptwX->allocatedSize = 0;
    ptwX->mallocFailedSize = 0;
    ptwX->points = NULL;
    ptwX_reallocatePoints( ptwX, size, 0 );
    return( ptwX->status );
}
/*
************************************************************
*/
ptwXPoints *ptwX_create( int64_t size, int64_t length, double const *xs, nfu_status *status ) {

    ptwXPoints *ptwX = ptwX_new( size, status );

    if( ptwX != NULL ) {
        if( ( *status = ptwX_setData( ptwX, length, xs ) ) != nfu_Okay ) ptwX = ptwX_free( ptwX );
    }
    return( ptwX );
}
/*
************************************************************
*/
ptwXPoints *ptwX_createLine( int64_t size, int64_t length, double slope, double offset, nfu_status *status ) {

    int64_t i1;
    double *p1;
    ptwXPoints *ptwX;

    if( size < length ) size = length;
    if( ( ptwX = ptwX_new( size, status ) ) != NULL ) {
        for( i1 = 0, p1 = ptwX->points; i1 < length; i1++, p1++ ) *p1 = slope * i1 + offset;
        ptwX->length = length;
    }
    return( ptwX );
}
/*
************************************************************
*/
nfu_status ptwX_copy( ptwXPoints *dest, ptwXPoints *src ) {

    if( dest->status == nfu_Okay ) return( dest->status );
    if( src->status == nfu_Okay ) return( src->status );
    ptwX_clear( dest );
    return( ptwX_setData( dest, src->length, src->points ) );
}
/*
************************************************************
*/
ptwXPoints *ptwX_clone( ptwXPoints *ptwX, nfu_status *status ) {

    return( ptwX_slice( ptwX, 0, ptwX->length, status ) );
}
/*
************************************************************
*/
ptwXPoints *ptwX_slice( ptwXPoints *ptwX, int64_t index1, int64_t index2, nfu_status *status ) {

    int64_t i, j, length;
    ptwXPoints *n;

    *status = nfu_badSelf;
    if( ptwX->status != nfu_Okay ) return( NULL );
    *status = nfu_badIndex;
    if( index1 < 0 ) return( NULL );
    if( index2 < index1 ) return( NULL );
    if( index2 > ptwX->length ) return( NULL );
    length = ( index2 - index1 );
    if( ( n = ptwX_new( length, status ) ) == NULL ) return( n );
    *status = n->status;
    for( j = 0, i = index1; i < index2; i++, j++ ) n->points[j] = ptwX->points[i];
    n->length = length;
    return( n );
}
/*
************************************************************
*/
nfu_status ptwX_reallocatePoints( ptwXPoints *ptwX, int64_t size, int forceSmallerResize ) {

    if( size < ptwX_minimumSize ) size = ptwX_minimumSize;                        /* ptwX_minimumSize must be > 0 for other routines to work properly. */
    if( size < ptwX->length ) size = ptwX->length;
    if( size != ptwX->allocatedSize ) {
        if( size > ptwX->allocatedSize ) {                                         /* Increase size of allocated points. */
             ptwX->points = (double *) nfu_realloc( (size_t) size * sizeof( double ), ptwX->points ); }
        else if( ( ptwX->allocatedSize > 2 * size ) || forceSmallerResize ) {      /* Decrease size, if at least 1/2 size reduction or if forced to. */
            ptwX->points = (double *) nfu_realloc( (size_t) size * sizeof( double ), ptwX->points );
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
nfu_status ptwX_clear( ptwXPoints *ptwX ) {

    ptwX->length = 0;
    return( ptwX->status );
}
/*
************************************************************
*/
nfu_status ptwX_release( ptwXPoints *ptwX ) {

    ptwX->length = 0;
    ptwX->allocatedSize = 0;
    ptwX->points = (double *) nfu_free( ptwX->points );

    return( nfu_Okay );
}
/*
************************************************************
*/
ptwXPoints *ptwX_free( ptwXPoints *ptwX ) {

    if( ptwX != NULL ) ptwX_release( ptwX );
    return( (ptwXPoints *) nfu_free( ptwX ) );
}
/*
************************************************************
*/
int64_t ptwX_length( ptwXPoints *ptwX ) {

    return( ptwX->length );
}
/*
************************************************************
*/
nfu_status ptwX_setData( ptwXPoints *ptwX, int64_t length, double const *xs ) {

    int64_t  i;

    if( ptwX->status != nfu_Okay ) return( ptwX->status );

    if( length > ptwX->allocatedSize ) {
        ptwX_reallocatePoints( ptwX, length, 0 );
        if( ptwX->status != nfu_Okay ) return( ptwX->status );
    }
    for( i = 0; i < length; i++ ) ptwX->points[i] = xs[i];
    ptwX->length = length;

    return( ptwX->status );
}
/*
************************************************************
*/ 
nfu_status ptwX_deletePoints( ptwXPoints *ptwX, int64_t i1, int64_t i2 ) {

    int64_t n = ptwX->length - ( i2 - i1 );

    if( ptwX->status != nfu_Okay ) return( ptwX->status );
    if( ( i1 < 0 ) || ( i1 > i2 ) || ( i2 > ptwX->length ) ) return( nfu_badIndex );
    if( i1 != i2 ) {
        for( ; i2 < ptwX->length; i1++, i2++ ) ptwX->points[i1] = ptwX->points[i2];
        ptwX->length = n;
    }
    return( ptwX->status );
}
/*
************************************************************
*/
double *ptwX_getPointAtIndex( ptwXPoints *ptwX, int64_t index ) {

    if( ptwX->status != nfu_Okay ) return( NULL );
    if( ( index < 0 ) || ( index >= ptwX->length ) ) return( NULL );
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
nfu_status ptwX_setPointAtIndex( ptwXPoints *ptwX, int64_t index, double x ) {

    nfu_status status;

    if( ptwX->status != nfu_Okay ) return( ptwX->status );
    if( ( index < 0 ) || ( index > ptwX->length ) ) return( nfu_badIndex );
    if( index == ptwX->allocatedSize ) {
        if( ( status = ptwX_reallocatePoints( ptwX, ptwX->allocatedSize + 10, 0 ) ) != nfu_Okay ) return( status );
    }
    ptwX->points[index] = x;
    if( index == ptwX->length ) ptwX->length++;
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwX_insertPointsAtIndex( ptwXPoints *ptwX, int64_t index, int64_t n1, double const *xs ) {

    nfu_status status;
    int64_t i1, i2, n1p, size = n1 + ptwX->length;

    if( ptwX->status != nfu_Okay ) return( ptwX->status );
    if( n1 < 1 ) return( nfu_Okay );
    if( ( index < 0 ) || ( index > ptwX->length ) ) return( nfu_badIndex );
    if( size > ptwX->allocatedSize ) {
        if( ( status = ptwX_reallocatePoints( ptwX, size, 0 ) ) != nfu_Okay ) return( status );
    }
    for( i1 = ptwX->length - 1, i2 = size - 1, n1p = ptwX->length - index + 1; n1p > 0; i1--, i2--, n1p-- ) ptwX->points[i2] = ptwX->points[i1];
    for( i1 = 0, i2 = index; i1 < n1; i1++, i2++ ) ptwX->points[i2] = xs[i1];
    ptwX->length += n1;
    return( nfu_Okay );
}
/*
************************************************************
*/
int ptwX_ascendingOrder( ptwXPoints *ptwX ) {
/*
*    Returns -1 list is descending, 1 if ascending and 0 otherwise (i.e., mixed).
*/
    int order = 1;
    int64_t i;
    double x1, x2;

    if( ptwX->length < 2 ) return( 0 );

    if( ( x1 = ptwX->points[0] ) < ( x2 = ptwX->points[1] ) ) {     /* Check for ascending order. */
        for( i = 2; i < ptwX->length; i++ ) {
            x1 = x2;
            x2 = ptwX->points[i];
            if( x2 <= x1 ) return( 0 );
        } }
    else {
        if( x1 == x2 ) return( 0 );
        order = -1;                                                 /* Check for descending order. */
        for( i = 2; i < ptwX->length; i++ ) {
            x1 = x2;
            x2 = ptwX->points[i];
            if( x1 <= x2 ) return( 0 );
        }
    }
    return( order );
}
/*
************************************************************
*/
ptwXPoints *ptwX_fromString( char const *str, char **endCharacter, nfu_status *status ) {

    int64_t numberConverted;
    double  *doublePtr;
    ptwXPoints *ptwX = NULL;

    if( ( *status = nfu_stringToListOfDoubles( str, &numberConverted, &doublePtr, endCharacter ) ) != nfu_Okay ) return( NULL );
    ptwX = ptwX_create( numberConverted, numberConverted, doublePtr, status );
    nfu_free( doublePtr );
    return( ptwX );
}
/*
************************************************************
*/
nfu_status ptwX_countOccurrences( ptwXPoints *ptwX, double value, int *count ) {

    int64_t i1;

    *count = 0;
    for( i1 = 0; i1 < ptwX->length; i1++ ) {
        if( ptwX->points[i1] == value ) (*count)++;
    }
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwX_reverse( ptwXPoints *ptwX ) {

    int64_t i1, i2 = ptwX->length - 1, n1 = ptwX->length / 2;
    double tmp;

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
nfu_status ptwX_sort( ptwXPoints *ptwX, enum ptwX_sort_order order ) {

    int (*cmp)( void const *, void const * ) = ptwX_sort_descending;

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
nfu_status ptwX_closesDifference( ptwXPoints *ptwX, double value, int64_t *index, double *difference ) {

    return( ptwX_closesDifferenceInRange( ptwX, 0, ptwX->length, value, index, difference ) );
}
/*
************************************************************
*/
nfu_status ptwX_closesDifferenceInRange( ptwXPoints *ptwX, int64_t i1, int64_t i2, double value, int64_t *index, double *difference ) {
/*
*   Finds the closes datum to value. If *difference is zero, datum is same as value.
*/
    double d1;

    *index = -1;
    *difference = -1;
    if( ptwX->status != nfu_Okay ) return( ptwX->status );
    if( i1 < 0 ) i1 = 0;
    if( i2 > ptwX->length ) i2 = ptwX->length;
    if( i1 >= i2 ) return( nfu_Okay );
    *index = i1;
    *difference = value - ptwX->points[i1];
    for( i1++; i1 < i2; i1++ ) {
        d1 = value - ptwX->points[i1];
        if( std::fabs( *difference ) > std::fabs( d1 ) ) {
            *index = i1;
            *difference = d1;
        }
    }
    return( nfu_Okay );
}
/*
************************************************************
*/
ptwXPoints *ptwX_unique( ptwXPoints *ptwX, int order, nfu_status *status ) {
/*
*   If order < 0 order is descending, if order > 0 order is ascending, otherwise, order is the same as ptwX.
*/
    int64_t i1, i2, n1 = 0;
    double x1, *p2;
    ptwXPoints *ptwX2 = NULL;

    if( order == 0 ) {
        if( ( ptwX2 = ptwX_new( ptwX->length, status ) ) == NULL ) return( NULL );
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
        if( ( ptwX2 = ptwX_clone( ptwX, status ) ) == NULL ) return( NULL );
        if( ( *status = ptwX_sort( ptwX2, ptwX_sort_order_ascending ) ) != nfu_Okay ) goto err;

        if( ptwX2->length > 1 ) {
            x1 = ptwX2->points[n1];
            n1++;
            for( i1 = 1; i1 < ptwX2->length; i1++ ) {
                if( x1 != ptwX2->points[i1] ) {
                    x1 = ptwX2->points[i1];
                    ptwX2->points[n1] = x1;
                    n1++;
                }
            }
            ptwX2->length = n1;
            if( order < 0 ) {
                if( ( *status = ptwX_sort( ptwX2, ptwX_sort_order_descending ) ) != nfu_Okay ) goto err;
            }
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
nfu_status ptwX_abs( ptwXPoints *ptwX ) {

    int64_t i1;
    double *p1;

    if( ptwX->status != nfu_Okay ) return( ptwX->status );
    for( i1 = 0, p1 = ptwX->points; i1 < ptwX->length; i1++, p1++ ) *p1 = std::fabs( *p1 );
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwX_neg( ptwXPoints *ptwX ) {

    return( ptwX_slopeOffset( ptwX, -1, 0 ) );
}
/*
************************************************************
*/
nfu_status ptwX_add_double( ptwXPoints *ptwX, double value ) {

    return( ptwX_slopeOffset( ptwX, 1, value ) );
}
/*
************************************************************
*/
nfu_status ptwX_mul_double( ptwXPoints *ptwX, double value ) {

    return( ptwX_slopeOffset( ptwX, value, 0 ) );
}
/*
************************************************************
*/
nfu_status ptwX_slopeOffset( ptwXPoints *ptwX, double slope, double offset ) {

    int64_t i1;
    double *p1;

    if( ptwX->status != nfu_Okay ) return( ptwX->status );
    for( i1 = 0, p1 = ptwX->points; i1 < ptwX->length; i1++, p1++ ) *p1 = slope * *p1 + offset;
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwX_add_ptwX( ptwXPoints *ptwX1, ptwXPoints *ptwX2 ) {

    int64_t i1;
    double *p1 = ptwX1->points, *p2 = ptwX2->points;

    if( ptwX1->status != nfu_Okay ) return( ptwX1->status );
    if( ptwX2->status != nfu_Okay ) return( ptwX2->status );
    if( ptwX1->length != ptwX2->length ) return( nfu_domainsNotMutual );

    for( i1 = 0; i1 < ptwX1->length; i1++, p1++, p2++ ) *p1 += *p2;
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwX_sub_ptwX( ptwXPoints *ptwX1, ptwXPoints *ptwX2 ) {

    int64_t i1;
    double *p1 = ptwX1->points, *p2 = ptwX2->points;

    if( ptwX1->status != nfu_Okay ) return( ptwX1->status );
    if( ptwX2->status != nfu_Okay ) return( ptwX2->status );
    if( ptwX1->length != ptwX2->length ) return( nfu_domainsNotMutual );

    for( i1 = 0; i1 < ptwX1->length; i1++, p1++, p2++ ) *p1 -= *p2;
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwX_xMinMax( ptwXPoints *ptwX, double *xMin, double *xMax ) {

    int64_t i1, n1 = ptwX->length;
    *xMin = *xMax = 0;
    double *p1 = ptwX->points;

    if( ptwX->status != nfu_Okay ) return( ptwX->status );
    if( n1 > 0 ) {
        *xMin = *xMax = *(p1++);
        for( i1 = 1; i1 < n1; ++i1, ++p1 ) {
            if( *p1 < *xMin ) *xMin = *p1;
            if( *p1 > *xMax ) *xMax = *p1;
        }
    }
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwX_compare( ptwXPoints *ptwX1, ptwXPoints *ptwX2, int *comparison ) {

    int64_t i1, n1 = ptwX1->length, n2 = ptwX2->length, nn = n1;
    double *p1 = ptwX1->points, *p2 = ptwX2->points;

    *comparison = 0;
    if( ptwX1->status != nfu_Okay ) return( ptwX1->status );
    if( ptwX2->status != nfu_Okay ) return( ptwX2->status );
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
int ptwX_close( ptwXPoints *ptwX1, ptwXPoints *ptwX2, int epsilonFactor, double epsilon, nfu_status *status ) {

    int64_t i1, n1 = ptwX1->length;
    double larger;
    double *p1 = ptwX1->points, *p2 = ptwX2->points;

    epsilon = std::fabs( epsilon ) + std::abs( epsilonFactor ) * DBL_EPSILON;

    *status = ptwX1->status;
    if( ptwX1->status != nfu_Okay ) return( -1 );
    *status = ptwX2->status;
    if( ptwX2->status != nfu_Okay ) return( -1 );
    *status = nfu_domainsNotMutual;
    if( n1 != ptwX2->length ) return( -1 );

    *status = nfu_Okay;
    for( i1 = 0; i1 < n1; i1++, p1++, p2++ ) {
        larger = std::fabs( *p1 );
        if( std::fabs( *p2 ) > larger ) larger = std::fabs( *p2 );
        if( std::fabs( *p2 - *p1 ) > epsilon * larger ) return( (int) ( i1 + 1 ) );
    }
    return( 0 );
}

#if defined __cplusplus
}
#endif
