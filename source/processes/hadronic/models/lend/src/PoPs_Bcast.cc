#ifdef PoPs_MPI
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "PoPs.h"
#include "PoPs_private.h"
#include "PoPs_Bcast_private.h"

#define NumberOfBcastArrays 3

enum PoPs_Bcast_mode { PoPs_Bcast_mode_count, PoPs_Bcast_mode_pack, PoPs_Bcast_mode_unpack };

typedef struct PoPs_Bcast_info {
    enum PoPs_Bcast_mode mode;
    int int_count, char_count, double_count;
    int *int_array;
    char *char_array;
    double *double_array;
} PoPs_Bcast_info;

static int PoPs_Bcast3( statusMessageReporting *smr, MPI_Comm comm, PoPs_Bcast_info *info, unitsDB *unitsRoot, PoPs *popsRoot );
static int PoPs_Bcast_PoPs( statusMessageReporting *smr, PoPs_Bcast_info *info, int index, PoPs *popsRoot );
static int PoPs_Bcast_PoPs2( statusMessageReporting *smr, PoPs_Bcast_info *info, PoP *pop );
static int PoPs_Bcast_int( statusMessageReporting *smr, PoPs_Bcast_info *info, int *value );
static int PoPs_Bcast_charAllocate( statusMessageReporting *smr, PoPs_Bcast_info *info, char **value );
static int PoPs_Bcast_double( statusMessageReporting *smr, PoPs_Bcast_info *info, double *value );
/*
========================================================================
*/
int PoPs_Bcast2( statusMessageReporting *smr, MPI_Comm comm, int bossRank, unitsDB *unitsRoot, PoPs *popsRoot ) {

    int myRank, status;
    int description[NumberOfBcastArrays]; 
    PoPs_Bcast_info info = { PoPs_Bcast_mode_count, 0, 0, 0, NULL, NULL, NULL };

    if( ( status = MPI_Errhandler_set( comm, MPI_ERRORS_RETURN ) ) != 0 ) return( status );
/*          New way but not on all systems yet.
    if( ( status = MPI_Comm_set_errhandler( comm, MPI_ERRORS_RETURN ) ) != 0 ) return( status );
*/
    if( ( status = MPI_Comm_rank( comm, &myRank ) ) != 0 ) return( status );

    if( myRank == bossRank ) {
        info.mode = PoPs_Bcast_mode_count;
        if( ( status = PoPs_Bcast3( smr, comm, &info, unitsRoot, popsRoot ) ) != 0 ) return( status );
        description[0] = info.int_count;
        description[1] = info.char_count;
        description[2] = info.double_count;
        if( ( info.int_array = (int *) smr_malloc2( smr, info.int_count * sizeof( int ), 1, "info.int_array" ) ) == NULL ) goto err;
        if( ( info.char_array = (char *) smr_malloc2( smr, info.char_count * sizeof( char ), 1, "info.char_array" ) ) == NULL ) goto err;
        if( ( info.double_array = (double *) smr_malloc2( smr, info.double_count * sizeof( double ), 1, "info.double_array" ) ) == NULL ) goto err;

        info.mode = PoPs_Bcast_mode_pack;
        info.int_count = 0;
        info.char_count = 0;
        info.double_count = 0;
        if( ( status = PoPs_Bcast3( smr, comm, &info, unitsRoot, popsRoot ) ) != 0 ) return( status );
        if( info.int_count != description[0] ) {
            smr_setReportError2( smr, PoPs_smr_ID, 1, "int counting count = %d != packing count = %d", info.int_count, description[0] );
            goto err;
        }
        if( info.char_count != description[1] ) {
            smr_setReportError2( smr, PoPs_smr_ID, 1, "char counting count = %d != packing count = %d", info.char_count, description[1] );
            goto err;
        }
        if( info.double_count != description[2] ) {
            smr_setReportError2( smr, PoPs_smr_ID, 1, "double counting count = %d != packing count = %d", info.double_count, description[2] );
            goto err;
        }
    }

    if( ( status = MPI_Bcast( description, NumberOfBcastArrays, MPI_INT, bossRank, comm ) ) != 0 ) goto err;

    if( myRank != bossRank ) {
        if( ( info.int_array = (int *) smr_malloc2( smr, description[0] * sizeof( int ), 1, "info.int_array (2)" ) ) == NULL ) goto err;
        if( ( info.char_array = (char *) smr_malloc2( smr, description[1] * sizeof( char ), 1, "info.char_array (2)" ) ) == NULL ) goto err;
        if( ( info.double_array = (double *) smr_malloc2( smr, description[2] * sizeof( double ), 1, "info.double_array (2)" ) ) == NULL ) goto err;
    }
    if( ( status = MPI_Bcast( info.int_array, description[0], MPI_INT, bossRank, comm ) ) != 0 ) goto err;
    if( ( status = MPI_Bcast( info.char_array, description[1], MPI_CHAR, bossRank, comm ) ) != 0 ) goto err;
    if( ( status = MPI_Bcast( info.double_array, description[2], MPI_DOUBLE, bossRank, comm ) ) != 0 ) goto err;

    if( myRank != bossRank ) {
        info.mode = PoPs_Bcast_mode_unpack;
        if( ( status = PoPs_Bcast3( smr, comm, &info, unitsRoot, popsRoot ) ) != 0 ) goto err;
    }

    if( info.int_array != NULL ) smr_freeMemory( (void **) &(info.int_array) );
    if( info.char_array != NULL ) smr_freeMemory( (void **) &(info.char_array) );
    if( info.double_array != NULL ) smr_freeMemory( (void **) &(info.double_array) );

    return( 0 );

err:
    if( info.int_array != NULL ) smr_freeMemory( (void **) &(info.int_array) );
    if( info.char_array != NULL ) smr_freeMemory( (void **) &(info.char_array) );
    if( info.double_array != NULL ) smr_freeMemory( (void **) &(info.double_array) );
    if( unitsRoot->unsorted != NULL ) smr_freeMemory( (void **) &(unitsRoot->unsorted) );
    if( popsRoot->pops != NULL ) smr_freeMemory( (void **) &(popsRoot->pops) );
    if( popsRoot->sorted != NULL ) smr_freeMemory( (void **) &(popsRoot->sorted) );
    return( -1 );
}
/*
========================================================================
*/
static int PoPs_Bcast3( statusMessageReporting *smr, MPI_Comm comm, PoPs_Bcast_info *info, unitsDB *unitsRoot, PoPs *popsRoot ) {

    int i, status, numberOfUnits, numberOfParticles;

    if( info->mode == PoPs_Bcast_mode_unpack ) PoPs_releasePrivate( smr );
    if( ( status = PoPs_Bcast_int( smr, info, &(unitsRoot->numberOfUnits) ) ) != 0 ) return( status );
    numberOfUnits = unitsRoot->numberOfUnits;
    if( info->mode == PoPs_Bcast_mode_unpack ) {
        unitsRoot->allocated = unitsRoot->numberOfUnits;
        unitsRoot->numberOfUnits = 0;
        if( ( unitsRoot->unsorted = (char const **) smr_malloc2( smr, unitsRoot->allocated * sizeof( char const ** ), 1, "unitsRoot->unsorted" ) ) == NULL ) return( -1 );
    }
    for( i = 0; i < numberOfUnits; i++ ) {
        if( ( status = PoPs_Bcast_charAllocate( smr, info, (char **) &(unitsRoot->unsorted[i]) ) ) != 0 ) return( status );
        if( info->mode == PoPs_Bcast_mode_unpack ) unitsRoot->numberOfUnits++;
    }

    if( ( status = PoPs_Bcast_int( smr, info, &(popsRoot->numberOfParticles) ) ) != 0 ) return( status );
    numberOfParticles = popsRoot->numberOfParticles;
    if( info->mode == PoPs_Bcast_mode_unpack ) {
        popsRoot->allocated = popsRoot->numberOfParticles;
        popsRoot->numberOfParticles = 0;
        if( ( popsRoot->pops = (PoP **) smr_malloc2( smr, popsRoot->allocated * sizeof( PoP * ), 1, "popsRoot->pops" ) ) == NULL ) return( -1 );
        if( ( popsRoot->sorted = (PoP **) smr_malloc2( smr, popsRoot->allocated * sizeof( PoP * ), 1, "popsRoot->unsorted" ) ) == NULL ) return( -1 );
    }
    for( i = 0; i < numberOfParticles; i++ ) {
        if( ( status = PoPs_Bcast_PoPs( smr, info, i, popsRoot ) ) != 0 ) return( status );
    }
    return( 0 );
}
/*
========================================================================
*/
static int PoPs_Bcast_PoPs( statusMessageReporting *smr, PoPs_Bcast_info *info, int index, PoPs *popsRoot ) {

    int status;
    PoP pop;

    if( info->mode != PoPs_Bcast_mode_unpack ) return( PoPs_Bcast_PoPs2( smr, info, popsRoot->pops[index] ) );
    if( ( status = PoPs_Bcast_PoPs2( smr, info, &pop ) ) != 0 ) return( status );
    return( 0 );
}
/*
========================================================================
*/
static int PoPs_Bcast_PoPs2( statusMessageReporting *smr, PoPs_Bcast_info *info, PoP *pop ) {

    int status, n = 0;

    if( ( status = PoPs_Bcast_int( smr, info, &(pop->index) ) ) != 0 ) return( status );
    if( ( status = PoPs_Bcast_int( smr, info, &(pop->properIndex) ) ) != 0 ) return( status );
    if( ( status = PoPs_Bcast_int( smr, info, &(pop->aliasIndex) ) ) != 0 ) return( status );   /* Not needed, see below. */
    if( ( status = PoPs_Bcast_int( smr, info, (int *) &(pop->genre) ) ) != 0 ) return( status );

    if( ( status = PoPs_Bcast_int( smr, info, &(pop->Z) ) ) != 0 ) return( status );
    if( ( status = PoPs_Bcast_int( smr, info, &(pop->A) ) ) != 0 ) return( status );
    if( ( status = PoPs_Bcast_int( smr, info, &(pop->l) ) ) != 0 ) return( status );
    if( ( status = PoPs_Bcast_double( smr, info, &(pop->mass) ) ) != 0 ) return( status );

    if( info->mode == PoPs_Bcast_mode_pack ) {
        n = -1;
        if( pop->massUnit != NULL ) {
            if( ( n = unitsDB_index( smr, pop->massUnit ) ) < 0 ) return( n );
        }
    }
    if( ( status = PoPs_Bcast_int( smr, info, &n ) ) != 0 ) return( status );
    if( ( status = PoPs_Bcast_charAllocate( smr, info, (char **) &(pop->name) ) ) != 0 ) return( status );

    if( info->mode == PoPs_Bcast_mode_unpack ) {
        pop->aliasIndex = -1;       /* Reset here as it will be set in PoPs_addParticleIfNeeded via PoPs_copyAddParticleIfNeeded. */

        if( n < 0 ) {
            pop->massUnit = NULL; }
        else {
            if( ( pop->massUnit = unitsDB_stringFromIndex( smr, n ) ) == NULL ) goto err;
        }
        if( PoPs_copyAddParticleIfNeeded( smr, pop ) == NULL ) goto err;

        if( pop->name != NULL ) smr_freeMemory( (void **) &(pop->name) );
    }

    return( 0 );

err:
    if( info->mode == PoPs_Bcast_mode_unpack ) {
        if( pop->name != NULL ) smr_freeMemory( (void **) &(pop->name) );
    }
    return( -1 );
}
/*
========================================================================
*/
static int PoPs_Bcast_int( statusMessageReporting *smr, PoPs_Bcast_info *info, int *value ) {

    if( info->mode == PoPs_Bcast_mode_pack ) {
        info->int_array[info->int_count] = *value; }
    else if( info->mode == PoPs_Bcast_mode_unpack ) {
        *value = info->int_array[info->int_count];
    }
    info->int_count++;
    return( 0 );
}
/*
========================================================================
*/
static int PoPs_Bcast_charAllocate( statusMessageReporting *smr, PoPs_Bcast_info *info, char **value ) {

    int i, n = 0, status;

    if( info->mode != PoPs_Bcast_mode_unpack ) {
        n = (int) strlen( *value ) + 1;
        if( ( status = PoPs_Bcast_int( smr, info, &n ) ) != 0 ) return( status );
        if( info->mode == PoPs_Bcast_mode_pack ) {
            for( i = 0; i < n; i++ ) info->char_array[info->char_count + i] = (*value)[i];
        } }
    else {
        if( ( status = PoPs_Bcast_int( smr, info, &n ) ) != 0 ) return( status  );
        if( ( *value = (char *) smr_malloc2( smr, n * sizeof( char ), 0, "*value" ) ) == NULL ) return( -1 );
        for( i = 0; i < n; i++ ) (*value)[i] = info->char_array[info->char_count + i];
    }
    info->char_count += n;
 
    return( 0 );
}
/*
========================================================================
*/
static int PoPs_Bcast_double( statusMessageReporting *smr, PoPs_Bcast_info *info, double *value ) {

    if( info->mode == PoPs_Bcast_mode_pack ) {
        info->double_array[info->double_count] = *value; }
    else if( info->mode == PoPs_Bcast_mode_unpack ) {
        *value = info->double_array[info->double_count];
    }
    info->double_count++;
    return( 0 );
}
#endif      /* End of #ifdef PoPs_MPI */
