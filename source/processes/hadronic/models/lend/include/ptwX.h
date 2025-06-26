/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef ptwX_h_included
#define ptwX_h_included

#include <stdio.h>
#include <stdint.h>

#include <nf_utilities.h>

#if defined __cplusplus
    extern "C" {
#endif

#define ptwX_minimumSize 10

enum ptwX_sort_order { ptwX_sort_order_descending, ptwX_sort_order_ascending };

typedef
    struct ptwXPoints_s {
        nfu_status status;
        int64_t length;
        int64_t allocatedSize;
        int64_t mallocFailedSize;
        double *points;
    } ptwXPoints;

/*
* Routines in ptwX_core.c
*/
ptwXPoints *ptwX_new( statusMessageReporting *smr, int64_t size );
nfu_status ptwX_initialize( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t size );
ptwXPoints *ptwX_create( statusMessageReporting *smr, int64_t size, int64_t length, double const *xs );
ptwXPoints *ptwX_createLine( statusMessageReporting *smr, int64_t size, int64_t length, double slope, double offset );
nfu_status ptwX_copy( statusMessageReporting *smr, ptwXPoints *dest, ptwXPoints *src );
ptwXPoints *ptwX_clone( statusMessageReporting *smr, ptwXPoints *ptwX );
ptwXPoints *ptwX_slice( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t index1, int64_t index2 );
nfu_status ptwX_reallocatePoints( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t size, int forceSmallerResize );
nfu_status ptwX_clear( statusMessageReporting *smr, ptwXPoints *ptwX );
nfu_status ptwX_release( statusMessageReporting *smr, ptwXPoints *ptwX );
ptwXPoints *ptwX_free( ptwXPoints *ptwX );

int64_t ptwX_length( statusMessageReporting *smr, ptwXPoints *ptwX );
nfu_status ptwX_setData( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t length, double const *xs );
nfu_status ptwX_deletePoints( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t i1, int64_t i2 );
double *ptwX_getPointAtIndex( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t index );
double ptwX_getPointAtIndex_Unsafely( ptwXPoints *ptwX, int64_t index );
nfu_status ptwX_setPointAtIndex( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t index, double x );
nfu_status ptwX_insertPointsAtIndex( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t index, int64_t n1, double const *xs );
nfu_status ptwX_ascendingOrder( statusMessageReporting *smr, ptwXPoints *ptwX, int *order );
ptwXPoints *ptwX_fromString( statusMessageReporting *smr, char const *str, char sep, char **endCharacter );
int ptwX_countOccurrences( statusMessageReporting *smr, ptwXPoints *ptwX, double value );
nfu_status ptwX_reverse( statusMessageReporting *smr, ptwXPoints *ptwX );
nfu_status ptwX_sort( statusMessageReporting *smr, ptwXPoints *ptwX, enum ptwX_sort_order order );
nfu_status ptwX_closesDifference( statusMessageReporting *smr, ptwXPoints *ptwX, double value, int64_t *index, double *difference );
nfu_status ptwX_closesDifferenceInRange( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t i1, int64_t i2, 
        double value, int64_t *index, double *difference );
ptwXPoints *ptwX_unique( statusMessageReporting *smr, ptwXPoints *ptwX, int order );

nfu_status ptwX_abs( statusMessageReporting *smr, ptwXPoints *ptwX );
nfu_status ptwX_neg( statusMessageReporting *smr, ptwXPoints *ptwX );
nfu_status ptwX_add_double( statusMessageReporting *smr, ptwXPoints *ptwX, double value );
nfu_status ptwX_mul_double( statusMessageReporting *smr, ptwXPoints *ptwX, double value );
nfu_status ptwX_slopeOffset( statusMessageReporting *smr, ptwXPoints *ptwX, double slope, double offset );
nfu_status ptwX_add_ptwX( statusMessageReporting *smr, ptwXPoints *ptwX1, ptwXPoints *ptwX2 );
nfu_status ptwX_sub_ptwX( statusMessageReporting *smr, ptwXPoints *ptwX1, ptwXPoints *ptwX2 );

nfu_status ptwX_range( statusMessageReporting *smr, ptwXPoints *ptwX, double *rangeMin, double *rangeMax );

nfu_status ptwX_compare( statusMessageReporting *smr, ptwXPoints *ptwX1, ptwXPoints *ptwX2, int *comparison );
nfu_status ptwX_close( statusMessageReporting *smr, ptwXPoints *ptwX1, ptwXPoints *ptwX2, int epsilonFactor, double epsilon,
        int *index );

/*
* Routines in ptwX_misc.c
*/
nfu_status ptwX_simpleWrite( statusMessageReporting *smr, ptwXPoints const *ptwX, FILE *f, char const *format );
nfu_status ptwX_simplePrint( statusMessageReporting *smr, ptwXPoints const *ptwX, char const *format );

#if defined __cplusplus
    }
#endif

#endif          /* End of ptwX_h_included. */
